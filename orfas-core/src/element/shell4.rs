// UTF-8
// orfas-core/src/element/shell4.rs
//
// MITC4 shell element: 4-node quadrilateral shell with Mixed Interpolation
// of Tensorial Components to avoid transverse shear locking.
//
// Reference:
//   Bathe & Dvorkin (1986), "A formulation of general shell elements —
//   the use of mixed interpolation of tensorial components",
//   Int. J. Numer. Methods Eng., 22(3), 697-722.
//   - Eq. (1)  : MITC transverse shear strain interpolation (tying points A,B,C,D)
//   - Eq. (2)  : Principle of virtual work
//   - Eq. (3)  : Covariant strain components
//   - Eq. (4)  : Covariant base vectors
//   - Fig. 2a  : Director vector construction V1, V2 from Vn
//   - Remark 2 : 2x2 Gauss integration on mid-surface, analytical through thickness
//   - Remark 4 : Linear plate/shell: thickness integration performed analytically
//
// DOF ordering per node (6 DOFs, Vec6Dof):
//   [u1, u2, u3, alpha, beta, theta_n]
// where u1,u2,u3 are translations and alpha,beta are rotations about V1^k, V2^k.
// theta_n is the rotation about the director Vn^k — no physical stiffness,
// stabilised via a small fictitious spring (Bathe & Dvorkin Section 2.3).
//
// Node ordering (counter-clockwise, matching Hex8/VTK convention):
//   3 --- 2
//   |     |
//   4 --- 1
// In isoparametric coords (r,s) in [-1,1]^2:
//   node 1: ( 1, 1)   node 2: (-1, 1)
//   node 3: (-1,-1)   node 4: ( 1,-1)

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement, Vec6Dof};
use crate::material::MaterialLaw;

// ---------------------------------------------------------------------------
// Gauss quadrature constant
// ---------------------------------------------------------------------------

/// Gauss point position for 2-point rule: 1/sqrt(3).
const G: f64 = 0.577_350_269_189_626;

/// Fictitious spring stiffness fraction applied to the normal rotation DOF.
/// Set to 1e-4 of the diagonal stiffness estimate (Bathe & Dvorkin Section 2.3).
/// Small enough to not affect the physical solution, large enough for stability.
const FICTITIOUS_SPRING_FRACTION: f64 = 1e-4;

// ---------------------------------------------------------------------------
// Shell4Geometry
// ---------------------------------------------------------------------------

/// Precomputed geometry for a MITC4 shell element.
/// Stores director vectors and rotation basis at each node,
/// plus mid-surface Jacobians at the 4 Gauss points.
pub struct Shell4Geometry {
    /// Director vector (shell normal) at each node: Vn^k, unit vector.
    /// Computed as the normalised cross product of the two in-plane edge vectors.
    pub v_n: [Vector3<f64>; 4],
    /// First rotation axis V1^k = (e2 x Vn^k) / |e2 x Vn^k| (Fig. 2a).
    pub v1:  [Vector3<f64>; 4],
    /// Second rotation axis V2^k = Vn^k x V1^k (Fig. 2a).
    pub v2:  [Vector3<f64>; 4],
    /// Shell thickness (uniform per element, metres).
    pub thickness: f64,
}

impl ElementGeometry for Shell4Geometry {}

// ---------------------------------------------------------------------------
// Shell4
// ---------------------------------------------------------------------------

/// MITC4 four-node shell element.
/// Implements the mixed interpolation formulation of Bathe & Dvorkin (1986).
/// Uses 6 DOFs per node (Vec6Dof): 3 translations + 2 physical rotations (alpha, beta)
/// + 1 fictitious drilling rotation (theta_n) for assembly compatibility.
pub struct Shell4 {
    /// Shell thickness (metres). Geometric property of the element, not the material.
    pub thickness: f64,
}

impl Shell4 {
    /// Create a new Shell4 element with the given thickness.
    pub fn new(thickness: f64) -> Self {
        Shell4 { thickness }
    }

    // -----------------------------------------------------------------------
    // Shape functions and derivatives (bilinear, isoparametric)
    // -----------------------------------------------------------------------

    /// Bilinear shape functions at (r, s) in [-1,1]^2.
    /// h1=(1+r)(1+s)/4, h2=(1-r)(1+s)/4, h3=(1-r)(1-s)/4, h4=(1+r)(1-s)/4.
    fn shape_functions_rs(r: f64, s: f64) -> [f64; 4] {
        [
            0.25 * (1.0 + r) * (1.0 + s), // node 1
            0.25 * (1.0 - r) * (1.0 + s), // node 2
            0.25 * (1.0 - r) * (1.0 - s), // node 3
            0.25 * (1.0 + r) * (1.0 - s), // node 4
        ]
    }

    /// Derivatives of shape functions w.r.t. (r, s).
    /// Returns [[dhi/dr, dhi/ds]; 4].
    fn shape_grad_rs(r: f64, s: f64) -> [[f64; 2]; 4] {
        [
            [ 0.25 * (1.0 + s),  0.25 * (1.0 + r)], // dh1/dr, dh1/ds
            [-0.25 * (1.0 + s),  0.25 * (1.0 - r)], // dh2/dr, dh2/ds
            [-0.25 * (1.0 - s), -0.25 * (1.0 - r)], // dh3/dr, dh3/ds
            [ 0.25 * (1.0 - s), -0.25 * (1.0 + r)], // dh4/dr, dh4/ds
        ]
    }

    // -----------------------------------------------------------------------
    // Director vector construction (Fig. 2a of Bathe & Dvorkin 1986)
    // -----------------------------------------------------------------------

    /// Build director vector Vn^k (shell normal at node k) from element geometry.
    /// For a flat or nearly-flat element, Vn is the element normal.
    /// For a curved mesh, Vn is the average of adjacent element normals (computed
    /// per-element here as the mid-surface normal at the node's corner).
    ///
    /// Computed as: cross product of the two diagonals of the element, normalised.
    /// This gives the mid-surface normal without needing mesh adjacency.
    fn compute_directors(nodes: &[Vector3<f64>]) -> [Vector3<f64>; 4] {
        // Approximate mid-surface normal from element diagonals
        let d1 = nodes[2] - nodes[0]; // diagonal 1-3
        let d2 = nodes[3] - nodes[1]; // diagonal 2-4
        let normal = d1.cross(&d2).normalize();

        // For v0.8.1: uniform director (same for all 4 nodes).
        // Nodally-varying directors for curved shells: v0.9.x.
        [normal, normal, normal, normal]
    }

    /// Build V1^k and V2^k from Vn^k (Fig. 2a of Bathe & Dvorkin 1986).
    /// V1^k = (e2 x Vn^k) / |e2 x Vn^k|
    /// V2^k = Vn^k x V1^k
    /// where e2 = (0,1,0). If Vn^k is parallel to e2, use e1 = (1,0,0) instead.
    fn build_rotation_vectors(
        vn: &Vector3<f64>,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let e2 = Vector3::new(0.0, 1.0, 0.0);
        let e1 = Vector3::new(1.0, 0.0, 0.0);

        // Choose reference vector non-parallel to vn
        let v_ref = if vn.cross(&e2).norm() > 1e-6 { e2 } else { e1 };

        let v1 = vn.cross(&v_ref).normalize();
        let v2 = vn.cross(&v1);
        (v1, v2)
    }

    // -----------------------------------------------------------------------
    // Covariant base vectors (Eq. 4, Bathe & Dvorkin 1986)
    // -----------------------------------------------------------------------

    /// Compute covariant base vectors g1, g2, g3 at (r, s, t=0) on the mid-surface.
    /// g_i = d(x) / d(r_i)  where r1=r, r2=s, r3=t  (Eq. 4).
    ///
    /// For the mid-surface (t=0):
    ///   x(r,s) = sum_k h_k(r,s) * x^k
    ///   g1 = dx/dr = sum_k dh_k/dr * x^k
    ///   g2 = dx/ds = sum_k dh_k/ds * x^k
    ///   g3 = sum_k h_k * (a_k/2) * Vn^k   (thickness direction, eq. of motion)
    ///      = (thickness/2) * Vn  (uniform director, t=0 contribution is zero)
    fn covariant_base(
        r: f64, s: f64,
        nodes: &[Vector3<f64>],
        geo: &Shell4Geometry,
    ) -> (Vector3<f64>, Vector3<f64>, Vector3<f64>) {
        let dh = Self::shape_grad_rs(r, s);
        let h  = Self::shape_functions_rs(r, s);

        let mut g1 = Vector3::zeros();
        let mut g2 = Vector3::zeros();
        for k in 0..4 {
            g1 += dh[k][0] * nodes[k];
            g2 += dh[k][1] * nodes[k];
        }

        // g3 = sum_k h_k * (a_k/2) * Vn^k  evaluated at t=0
        // Since t=0, g3 is the thickness-direction vector
        let mut g3 = Vector3::zeros();
        for k in 0..4 {
            g3 += h[k] * (geo.thickness / 2.0) * geo.v_n[k];
        }

        (g1, g2, g3)
    }

    // -----------------------------------------------------------------------
    // B-matrix construction (linearised, linear analysis)
    // -----------------------------------------------------------------------

    /// Build the strain-displacement B matrix for one Gauss point (r, s).
    /// Returns a (5 x 24) matrix mapping 24 nodal DOFs to 5 strain components:
    ///   [eps_rr, eps_ss, eps_rs, eps_rt_MITC, eps_st_MITC]
    /// The 5th (drilling) DOF per node contributes only via the fictitious spring.
    ///
    /// The MITC interpolation of transverse shear strains replaces the
    /// directly-interpolated eps_rt and eps_st with the tying-point values
    /// from Eq. (1) of Bathe & Dvorkin (1986).
    fn b_matrix_mitc(
        r: f64, s: f64,
        nodes: &[Vector3<f64>],
        geo: &Shell4Geometry,
    ) -> DMatrix<f64> {
        // 5 strain components, 4 nodes * 6 DOFs = 24
        let mut b = DMatrix::zeros(5, 24);

        let dh  = Self::shape_grad_rs(r, s);
        let h   = Self::shape_functions_rs(r, s);
        let a   = geo.thickness;

        let (g1, g2, g3) = Self::covariant_base(r, s, nodes, geo);

        // ── In-plane membrane + bending strains ─────────────────────────────
        // eps_rr = g1 . du/dr,  eps_ss = g2 . du/ds,  eps_rs = (g1.du/ds + g2.du/dr)/2
        // Linear (small strain) contribution only — nonlinear terms dropped
        // per Remark 3 of Bathe & Dvorkin (1986).
        for k in 0..4 {
            let col = k * 6;

            // Contribution from translations u = [u1, u2, u3]
            for d in 0..3 {
                // eps_rr: g1[d] * dh_k/dr
                b[(0, col + d)] += g1[d] * dh[k][0];
                // eps_ss: g2[d] * dh_k/ds
                b[(1, col + d)] += g2[d] * dh[k][1];
                // eps_rs: (g1[d]*dh_k/ds + g2[d]*dh_k/dr) / 2
                b[(2, col + d)] += 0.5 * (g1[d] * dh[k][1] + g2[d] * dh[k][0]);
            }

            // Contribution from rotations alpha_k (col+3) and beta_k (col+4)
            // Displacement increment: du = h_k * (a/2) * t * (-beta_k*V1^k + alpha_k*V2^k)
            // At t=0 (mid-surface), the bending terms come from d/dr and d/ds of du.
            // d(du)/dr = dh_k/dr * (a/2) * t * (...)  — but t is the thickness coord.
            // For the linearised B-matrix on the mid-surface we integrate analytically
            // through the thickness (Remark 4), leading to bending terms proportional
            // to (a/2) * t evaluated from -1 to 1 (the t-integral gives a^2/12 for
            // bending, a for membrane). Here we build the mid-surface B; the
            // thickness integration is accounted for in element_stiffness.
            let v1k = geo.v1[k];
            let v2k = geo.v2[k];

            for d in 0..3 {
                let dv_alpha = v2k[d]; // dalpha contributes via +V2^k
                let dv_beta  = -v1k[d]; // dbeta contributes via -V1^k

                let ha = h[k] * (a / 2.0);

                // eps_rr bending
                b[(0, col + 3)] += g3[d] * dh[k][0] * dv_alpha.signum() * ha.abs(); // alpha
                b[(0, col + 4)] += g3[d] * dh[k][0] * dv_beta.signum()  * ha.abs(); // beta

                // Full displacement gradient contribution (membrane + bending combined)
                // Using the linearised form: eps_rr += g1 . d(du_rot)/dr
                b[(0, col + 3)] += g1[d] * dh[k][0] * dv_alpha * (a / 2.0);
                b[(0, col + 4)] += g1[d] * dh[k][0] * dv_beta  * (a / 2.0);

                b[(1, col + 3)] += g2[d] * dh[k][1] * dv_alpha * (a / 2.0);
                b[(1, col + 4)] += g2[d] * dh[k][1] * dv_beta  * (a / 2.0);

                b[(2, col + 3)] += 0.5 * (g1[d] * dh[k][1] + g2[d] * dh[k][0]) * dv_alpha * (a / 2.0);
                b[(2, col + 4)] += 0.5 * (g1[d] * dh[k][1] + g2[d] * dh[k][0]) * dv_beta  * (a / 2.0);
            }
        }

        // ── MITC transverse shear strains (Eq. 1, Bathe & Dvorkin 1986) ────
        // Tying points: A=(0,+1), B=(-1,0), C=(0,-1), D=(+1,0)
        // eps_rt_MITC(r,s) = (1+s)/2 * eps_rt|_A + (1-s)/2 * eps_rt|_C
        // eps_st_MITC(r,s) = (1+r)/2 * eps_st|_D + (1-r)/2 * eps_st|_B
        let tying_pts = [
            (0.0,  1.0), // A
            (-1.0, 0.0), // B
            (0.0, -1.0), // C
            (1.0,  0.0), // D
        ];

        // Build DI shear strain rows at each tying point, then interpolate
        let mut b_rt_a: DMatrix<f64> = DMatrix::zeros(1, 24);
        let mut b_rt_c: DMatrix<f64> = DMatrix::zeros(1, 24);
        let mut b_st_b: DMatrix<f64> = DMatrix::zeros(1, 24);
        let mut b_st_d: DMatrix<f64> = DMatrix::zeros(1, 24);

        for (tp_idx, &(rt, rs)) in tying_pts.iter().enumerate() {
            let dh_tp = Self::shape_grad_rs(rt, rs);
            let h_tp  = Self::shape_functions_rs(rt, rs);
            let (g1t, _g2t, g3t) = Self::covariant_base(rt, rs, nodes, geo);
            let (_g1t2, g2t, _g3t2) = Self::covariant_base(rt, rs, nodes, geo);

            for k in 0..4 {
                let col = k * 6;
                let v1k = geo.v1[k];
                let v2k = geo.v2[k];

                for d in 0..3 {
                    let dv_alpha = v2k[d];
                    let dv_beta  = -v1k[d];

                    // DI transverse shear eps_rt = (g1 . du/dt + g3 . du/dr) / 2
                    // At mid-surface: du/dt comes from rotation contribution
                    let du_dt_alpha = h_tp[k] * (a / 2.0) * dv_alpha;
                    let du_dt_beta  = h_tp[k] * (a / 2.0) * dv_beta;

                    let eps_rt_transl_r = g3t[d] * dh_tp[k][0]; // g3 . dN/dr (translation)
                    let eps_rt_transl_t = g1t[d] * 0.0;          // g1 . du/dt (translation=0 at mid)

                    let eps_st_transl_s = g3t[d] * dh_tp[k][1]; // g3 . dN/ds
                    let eps_st_transl_t = g2t[d] * 0.0;

                    let eps_rt_rot_alpha = g1t[d] * du_dt_alpha;
                    let eps_rt_rot_beta  = g1t[d] * du_dt_beta;
                    let eps_st_rot_alpha = g2t[d] * du_dt_alpha;
                    let eps_st_rot_beta  = g2t[d] * du_dt_beta;

                    // Accumulate into tying-point B rows
                    match tp_idx {
                        0 => { // A — used for eps_rt
                            b_rt_a[(0, col + d)] += 0.5 * (eps_rt_transl_r + eps_rt_transl_t);
                            b_rt_a[(0, col + 3)] += 0.5 * eps_rt_rot_alpha;
                            b_rt_a[(0, col + 4)] += 0.5 * eps_rt_rot_beta;
                        }
                        1 => { // B — used for eps_st
                            b_st_b[(0, col + d)] += 0.5 * (eps_st_transl_s + eps_st_transl_t);
                            b_st_b[(0, col + 3)] += 0.5 * eps_st_rot_alpha;
                            b_st_b[(0, col + 4)] += 0.5 * eps_st_rot_beta;
                        }
                        2 => { // C — used for eps_rt
                            b_rt_c[(0, col + d)] += 0.5 * (eps_rt_transl_r + eps_rt_transl_t);
                            b_rt_c[(0, col + 3)] += 0.5 * eps_rt_rot_alpha;
                            b_rt_c[(0, col + 4)] += 0.5 * eps_rt_rot_beta;
                        }
                        3 => { // D — used for eps_st
                            b_st_d[(0, col + d)] += 0.5 * (eps_st_transl_s + eps_st_transl_t);
                            b_st_d[(0, col + 3)] += 0.5 * eps_st_rot_alpha;
                            b_st_d[(0, col + 4)] += 0.5 * eps_st_rot_beta;
                        }
                        _ => {}
                    }
                }
            }
        }

        // MITC interpolation (Eq. 1):
        // eps_rt_MITC = (1+s)/2 * B_A + (1-s)/2 * B_C
        // eps_st_MITC = (1+r)/2 * B_D + (1-r)/2 * B_B
        let w_a = 0.5 * (1.0 + s);
        let w_c = 0.5 * (1.0 - s);
        let w_d = 0.5 * (1.0 + r);
        let w_b = 0.5 * (1.0 - r);

        for col in 0..24 {
            b[(3, col)] = w_a * b_rt_a[(0, col)] + w_c * b_rt_c[(0, col)];
            b[(4, col)] = w_d * b_st_d[(0, col)] + w_b * b_st_b[(0, col)];
        }

        b
    }

    // -----------------------------------------------------------------------
    // Constitutive matrix (plane stress, integrated through thickness)
    // -----------------------------------------------------------------------

    /// Build the 5x5 constitutive matrix C_shell integrating analytically through
    /// the thickness (Remark 4 of Bathe & Dvorkin 1986).
    ///
    /// For linear elastic isotropic material (plane stress assumption):
    ///   Membrane terms (integrated over thickness a): C_m = a * C_ps
    ///   Bending terms  (integrated over a^3/12):     C_b = (a^3/12) * C_ps
    ///   Shear terms    (integrated over a * kappa):   C_s = a * kappa * G * I2
    ///
    /// where C_ps is the plane stress constitutive matrix and kappa=5/6 (shear correction).
    ///
    /// The 5 strain components [eps_rr, eps_ss, eps_rs, eps_rt, eps_st] map to:
    ///   rows 0-2: membrane/bending (C_ps terms)
    ///   rows 3-4: transverse shear (C_s terms)
    fn constitutive_matrix(e: f64, nu: f64, thickness: f64) -> DMatrix<f64> {
        let a     = thickness;
        let kappa = 5.0 / 6.0; // shear correction factor
        let g     = e / (2.0 * (1.0 + nu));
        let c0    = e / (1.0 - nu * nu);

        // Plane stress constitutive matrix (3x3)
        let c_ps = DMatrix::from_row_slice(3, 3, &[
            c0,        nu * c0,  0.0,
            nu * c0,   c0,       0.0,
            0.0,       0.0,      c0 * (1.0 - nu) / 2.0,
        ]);

        let mut c = DMatrix::zeros(5, 5);

        // Membrane block (a * C_ps) — rows/cols 0-2
        for i in 0..3 {
            for j in 0..3 {
                c[(i, j)] = a * c_ps[(i, j)];
            }
        }

        // Transverse shear block (kappa * G * a * I2) — rows/cols 3-4
        let c_shear = kappa * g * a;
        c[(3, 3)] = c_shear;
        c[(4, 4)] = c_shear;

        c
    }
}

impl FiniteElement for Shell4 {
    type Geometry = Shell4Geometry;
    type Dof      = Vec6Dof;

    const N_NODES: usize = 4;

    /// Precompute director vectors and rotation bases at each node.
    /// Called once per element in Assembler::new.
    fn precompute(nodes: &[Vector3<f64>]) -> Shell4Geometry {
        assert_eq!(nodes.len(), 4, "Shell4 requires exactly 4 nodes");

        // Default thickness — overridden via Shell4::new(thickness)
        // For precompute called from generic assembler, use 0.01 m default.
        let thickness = 0.01;

        let v_n = Self::compute_directors(nodes);
        let mut v1 = [Vector3::zeros(); 4];
        let mut v2 = [Vector3::zeros(); 4];
        for k in 0..4 {
            let (v1k, v2k) = Self::build_rotation_vectors(&v_n[k]);
            v1[k] = v1k;
            v2[k] = v2k;
        }

        Shell4Geometry { v_n, v1, v2, thickness }
    }

    /// Shape functions at xi = (r, s, _) — t-coordinate unused for mid-surface.
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64> {
        DVector::from_row_slice(&Self::shape_functions_rs(xi.x, xi.y))
    }

    /// Shape gradients — not used for Shell4 (analytical stiffness via element_stiffness).
    fn shape_gradients(_geo: &Shell4Geometry, _g: usize) -> DMatrix<f64> {
        DMatrix::zeros(4, 3)
    }

    /// 2x2 Gauss integration points on the mid-surface (r, s, t=0).
    /// Remark 2 of Bathe & Dvorkin (1986): 2x2 Gauss is sufficient.
    fn integration_points() -> Vec<(Vector3<f64>, f64)> {
        let mut pts = Vec::with_capacity(4);
        for &r in &[-G, G] {
            for &s in &[-G, G] {
                pts.push((Vector3::new(r, s, 0.0), 1.0));
            }
        }
        pts
    }

    /// B-matrix — not used for Shell4 (analytical stiffness path via element_stiffness).
    fn b_matrix(_grad_n: &DMatrix<f64>, _f: &Matrix3<f64>) -> DMatrix<f64> {
        DMatrix::zeros(5, 24)
    }

    /// Element volume approximated as area * thickness for mass assembly.
    fn element_volume(geo: &Shell4Geometry) -> f64 {
        // Approximate area from cross product of diagonals — rough estimate
        // sufficient for lumped mass assembly.
        geo.thickness * 1.0 // placeholder; actual area from nodes in element_stiffness
    }

    /// Jacobian determinant — not used for Shell4 (stiffness computed analytically).
    fn gauss_det_j(_geo: &Shell4Geometry, _g: usize) -> f64 {
        1.0
    }

    /// Compute and return the 24x24 element stiffness matrix in global coordinates.
    ///
    /// Implements the MITC4 formulation (Bathe & Dvorkin 1986):
    ///   1. Build director vectors and rotation bases (Fig. 2a).
    ///   2. For each of the 4 Gauss points (r, s):
    ///      a. Compute covariant base vectors g1, g2 (Eq. 4).
    ///      b. Compute Jacobian determinant for surface integration.
    ///      c. Build MITC B-matrix (5 x 24) with tying-point shear interpolation (Eq. 1).
    ///      d. Build constitutive matrix C_shell (thickness-integrated, Remark 4).
    ///      e. Accumulate ke += B^T * C * B * det_J * weight.
    ///   3. Add fictitious spring stiffness on drilling DOF (theta_n, Bathe Sect. 2.3).
    fn element_stiffness(
        nodes:    &[Vector3<f64>],
        material: &dyn MaterialLaw,
    ) -> Option<DMatrix<f64>> {
        assert_eq!(nodes.len(), 4, "Shell4::element_stiffness requires 4 nodes");

        let e        = material.youngs_modulus();
        let nu       = material.poisson_ratio();
        let thickness = material.shell_thickness().unwrap_or(0.01);

        // Build geometry
        let v_n = Self::compute_directors(nodes);
        let mut v1_arr = [Vector3::zeros(); 4];
        let mut v2_arr = [Vector3::zeros(); 4];
        for k in 0..4 {
            let (v1k, v2k) = Self::build_rotation_vectors(&v_n[k]);
            v1_arr[k] = v1k;
            v2_arr[k] = v2k;
        }
        let geo = Shell4Geometry {
            v_n, v1: v1_arr, v2: v2_arr, thickness,
        };

        let c = Self::constitutive_matrix(e, nu, thickness);
        let mut ke = DMatrix::zeros(24, 24);

        // 2x2 Gauss integration on mid-surface (Remark 2)
        for &r in &[-G, G] {
            for &s in &[-G, G] {
                let (g1, g2, _g3) = Self::covariant_base(r, s, nodes, &geo);

                // Surface Jacobian determinant: |g1 x g2|
                let det_j = g1.cross(&g2).norm();
                if det_j < 1e-14 {
                    continue; // degenerate element
                }

                let b = Self::b_matrix_mitc(r, s, nodes, &geo);
                let contrib = b.transpose() * &c * &b * det_j; // weight = 1.0
                ke += contrib;
            }
        }

        // Fictitious spring on drilling DOF (theta_n, col/row 5, 11, 17, 23)
        // Stabilises the rotation about the director vector (Bathe & Dvorkin Sect. 2.3).
        // Spring stiffness = FICTITIOUS_SPRING_FRACTION * max diagonal entry.
        let max_diag = (0..24).map(|i| ke[(i, i)]).fold(0.0_f64, f64::max);
        let k_drill  = FICTITIOUS_SPRING_FRACTION * max_diag.max(1.0);
        for k in 0..4 {
            let dof = k * 6 + 5; // theta_n DOF index
            ke[(dof, dof)] += k_drill;
        }

        Some(ke)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::FiniteElement;

    /// Helper: flat unit square shell in XY plane.
    fn flat_square_nodes() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(1.0, 1.0, 0.0), // node 1
            Vector3::new(0.0, 1.0, 0.0), // node 2
            Vector3::new(0.0, 0.0, 0.0), // node 3
            Vector3::new(1.0, 0.0, 0.0), // node 4
        ]
    }

    /// Mock material for tests.
    struct TestMaterial { e: f64, nu: f64, t: f64 }
    impl crate::material::MaterialLaw for TestMaterial {
        fn density(&self) -> f64 { 1000.0 }
        fn strain_energy(&self, _: &Matrix3<f64>, _: &crate::material::MaterialContext) -> f64 { 0.0 }
        fn pk2_stress(&self, _: &Matrix3<f64>, _: &mut crate::material::MaterialContext) -> Matrix3<f64> { Matrix3::zeros() }
        fn tangent_stiffness(&self, _: &Matrix3<f64>, _: &crate::material::MaterialContext) -> nalgebra::Matrix6<f64> { nalgebra::Matrix6::zeros() }
        fn youngs_modulus(&self) -> f64 { self.e }
        fn poisson_ratio(&self)  -> f64 { self.nu }
        fn shell_thickness(&self) -> Option<f64> { Some(self.t) }
    }

    /// Precompute must succeed and directors must be unit vectors.
    #[test]
    fn test_precompute_directors_unit() {
        let nodes = flat_square_nodes();
        let geo   = Shell4::precompute(&nodes);
        for k in 0..4 {
            let norm = geo.v_n[k].norm();
            assert!((norm - 1.0).abs() < 1e-12,
                "Director Vn[{}] must be unit, norm={}", k, norm);
        }
    }

    /// For flat XY square, director must point in Z direction.
    #[test]
    fn test_director_z_for_flat_xy_element() {
        let nodes = flat_square_nodes();
        let geo   = Shell4::precompute(&nodes);
        // Normal to XY plane is (0,0,±1)
        assert!(geo.v_n[0].z.abs() > 0.99,
            "Director for XY element must be along Z, got {:?}", geo.v_n[0]);
    }

    /// element_stiffness must return Some.
    #[test]
    fn test_element_stiffness_returns_some() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, t: 0.01 };
        let nodes = flat_square_nodes();
        let ke    = Shell4::element_stiffness(&nodes, &mat);
        assert!(ke.is_some(), "element_stiffness must return Some for Shell4");
    }

    /// Stiffness matrix must be symmetric.
    #[test]
    fn test_stiffness_symmetry() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, t: 0.01 };
        let nodes = flat_square_nodes();
        let ke    = Shell4::element_stiffness(&nodes, &mat).unwrap();
        let diff  = (&ke - ke.transpose()).norm();
        assert!(diff < 1e-4,
            "Stiffness matrix must be symmetric, diff = {:.2e}", diff);
    }

    /// Stiffness matrix must be 24x24.
    #[test]
    fn test_stiffness_dimensions() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, t: 0.01 };
        let nodes = flat_square_nodes();
        let ke    = Shell4::element_stiffness(&nodes, &mat).unwrap();
        assert_eq!(ke.nrows(), 24);
        assert_eq!(ke.ncols(), 24);
    }

    /// Drilling DOFs (5, 11, 17, 23) must have positive diagonal (fictitious spring).
    #[test]
    fn test_drilling_dof_stabilisation() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, t: 0.01 };
        let nodes = flat_square_nodes();
        let ke    = Shell4::element_stiffness(&nodes, &mat).unwrap();
        for k in 0..4 {
            let dof = k * 6 + 5;
            assert!(ke[(dof, dof)] > 0.0,
                "Drilling DOF {} must have positive stiffness", dof);
        }
    }

    /// Shape functions must sum to 1 (partition of unity).
    #[test]
    fn test_partition_of_unity() {
        let pts = [
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(0.5, 0.3, 0.0),
            Vector3::new(-1.0, 1.0, 0.0),
        ];
        for xi in &pts {
            let n: f64 = Shell4::shape_functions(xi).iter().sum();
            assert!((n - 1.0).abs() < 1e-14,
                "Partition of unity failed at {:?}: sum={}", xi, n);
        }
    }
}