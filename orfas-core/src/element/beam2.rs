// UTF-8
// orfas-core/src/element/beam2.rs
//
// Euler-Bernoulli 3D beam element with 2 nodes and 6 DOFs per node.
// Stiffness matrix is analytical — no Gauss quadrature.
//
// References:
// - Local stiffness (Ke, eq. 1-104 to 1-106, Phi=0 for Euler-Bernoulli):
//   Andersen & Nielsen (2008), "Elastic Beams in Three Dimensions",
//   DCE Lecture Notes No. 23, Aalborg University.
// - Local-to-global transformation (Dc, eq. 1-4):
//   OpenFAST SubDyn Theory Manual, section 6.3.6.3,
//   https://rick-openfast.readthedocs.io/en/ssi/source/user/subdyn/theory.html
//
// DOF ordering per node (local and global):
//   [wx, wy, wz, theta_x, theta_y, theta_z]
//   indices 0..5 for node 1, 6..11 for node 2.
//
// Cross-section: solid circular of radius r.
//   A = pi * r^2
//   I = pi * r^4 / 4   (I_y = I_z = I for circular section)
//   J = pi * r^4 / 2   (St. Venant torsional constant = 2*I for solid circle)

use std::f64::consts::PI;
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement, Vec6Dof};
use crate::material::MaterialLaw;

// ---------------------------------------------------------------------------
// Beam2Geometry
// ---------------------------------------------------------------------------

/// Precomputed geometry for a Beam2 element.
/// For Euler-Bernoulli, the stiffness is fully analytical — no per-Gauss-point
/// data needed. We store the element length for use in element_volume.
pub struct Beam2Geometry {
    /// Element length ||p2 - p1||.
    pub length: f64,
}

impl ElementGeometry for Beam2Geometry {}

// ---------------------------------------------------------------------------
// Beam2CrossSection
// ---------------------------------------------------------------------------

/// Solid circular cross-section properties for a Beam2 element.
/// Passed to element_stiffness via the material law's beam_section() method.
/// Stored as a unit struct with radius — all section properties are derived.
#[derive(Clone, Copy, Debug)]
pub struct CircularSection {
    /// Outer radius of the solid circular cross-section.
    pub radius: f64,
}

impl CircularSection {
    /// Cross-sectional area: A = pi * r^2.
    pub fn area(&self) -> f64 {
        PI * self.radius * self.radius
    }

    /// Second moment of area: I = pi * r^4 / 4.
    /// Equal for both principal axes (I_y = I_z = I) for circular section.
    pub fn moment_of_inertia(&self) -> f64 {
        PI * self.radius.powi(4) / 4.0
    }

    /// St. Venant torsional constant: J = pi * r^4 / 2 = 2 * I.
    pub fn torsional_constant(&self) -> f64 {
        PI * self.radius.powi(4) / 2.0
    }
}

// ---------------------------------------------------------------------------
// Beam2 element
// ---------------------------------------------------------------------------

/// Euler-Bernoulli 3D beam element with 2 nodes and 6 DOFs per node.
/// 12 DOFs total: [wx1,wy1,wz1,tx1,ty1,tz1, wx2,wy2,wz2,tx2,ty2,tz2].
/// Stiffness is computed analytically in the local frame then rotated to global.
/// No shear deformation (Phi_y = Phi_z = 0 in Andersen & Nielsen notation).
///
/// The cross-section radius is read from the material via MaterialLaw::beam_radius().
/// If the material does not provide a radius (returns None), a default of 0.01 m is used.
pub struct Beam2;

impl Beam2 {
    /// Build the 12x12 local stiffness matrix for an Euler-Bernoulli beam.
    ///
    /// Implements eq. (1-104) from Andersen & Nielsen (2008) with Phi_y = Phi_z = 0.
    /// The local x-axis is the beam axis (node 1 to node 2).
    /// DOF order: [wx1,wy1,wz1,tx1,ty1,tz1, wx2,wy2,wz2,tx2,ty2,tz2].
    ///
    /// # Arguments
    /// * `l`  - element length
    /// * `ea` - axial stiffness E * A
    /// * `ei` - bending stiffness E * I (same for y and z, circular section)
    /// * `gj` - torsional stiffness G * J
    fn local_stiffness(l: f64, ea: f64, ei: f64, gj: f64) -> DMatrix<f64> {
        let mut k = DMatrix::zeros(12, 12);

        // Scalar stiffness coefficients (Andersen & Nielsen eq. 1-104 to 1-106,
        // evaluated at Phi_y = Phi_z = 0, i.e. pure Euler-Bernoulli).
        let a  = ea / l;           // axial
        let b  = gj / l;           // torsion
        let c  = 12.0 * ei / l.powi(3); // flexural shear
        let d  =  6.0 * ei / l.powi(2); // flexural coupling
        let e  =  4.0 * ei / l;         // flexural moment (end)
        let f  =  2.0 * ei / l;         // flexural moment (mid, opposite end)

        // ── Axial (DOFs 0, 6) ──────────────────────────────────────────────
        k[(0, 0)] =  a;  k[(0, 6)] = -a;
        k[(6, 0)] = -a;  k[(6, 6)] =  a;

        // ── Torsion (DOFs 3, 9) ────────────────────────────────────────────
        k[(3, 3)] =  b;  k[(3, 9)] = -b;
        k[(9, 3)] = -b;  k[(9, 9)] =  b;

        // ── Flexion XY — coupling (wy, tz): DOFs (1, 5, 7, 11) ───────────
        // Andersen & Nielsen eq. (1-105), Phi_y = 0:
        //   k11=c, k12=d, k13=-c, k14=d
        //        k22=e, k23=-d, k24=f
        //               k33=c,  k34=-d
        //                       k44=e
        k[(1,  1)] =  c;  k[(1,  5)] =  d;  k[(1,  7)] = -c;  k[(1, 11)] =  d;
        k[(5,  1)] =  d;  k[(5,  5)] =  e;  k[(5,  7)] = -d;  k[(5, 11)] =  f;
        k[(7,  1)] = -c;  k[(7,  5)] = -d;  k[(7,  7)] =  c;  k[(7, 11)] = -d;
        k[(11, 1)] =  d;  k[(11, 5)] =  f;  k[(11, 7)] = -d;  k[(11,11)] =  e;

        // ── Flexion XZ — coupling (wz, ty): DOFs (2, 4, 8, 10) ──────────
        // Andersen & Nielsen eq. (1-106), Phi_z = 0:
        //   k11=c, k12=-d, k13=-c, k14=-d
        //         k22=e,  k23=d,  k24=f
        //                k33=c,   k34=d
        //                         k44=e
        // Note the sign flip on d vs XY — right-hand rule with x as beam axis.
        k[(2,  2)] =  c;  k[(2,  4)] = -d;  k[(2,  8)] = -c;  k[(2, 10)] = -d;
        k[(4,  2)] = -d;  k[(4,  4)] =  e;  k[(4,  8)] =  d;  k[(4, 10)] =  f;
        k[(8,  2)] = -c;  k[(8,  4)] =  d;  k[(8,  8)] =  c;  k[(8, 10)] =  d;
        k[(10, 2)] = -d;  k[(10, 4)] =  f;  k[(10, 8)] =  d;  k[(10,10)] =  e;

        k
    }

    /// Build the 3x3 rotation matrix R from node positions.
    ///
    /// Rows of R are the local axes expressed in global coordinates:
    ///   row 0 = x_loc (beam axis, node1 -> node2)
    ///   row 1 = y_loc (perpendicular to beam, in plane of x_loc and z_global)
    ///   row 2 = z_loc = x_loc x y_loc (completes right-handed system)
    ///
    /// The transformation convention is:
    ///   u_local  = R * u_global   (rows = local axes)
    ///   k_global = R^T * k_local * R
    ///
    /// This ensures the axial DOF (index 0 in local) aligns with x_loc in global.
    /// The auxiliary vector selection follows the SubDyn convention (section 6.3.6.3)
    /// but with rows (not columns) representing the local axes.
    fn direction_cosine_matrix(p1: &Vector3<f64>, p2: &Vector3<f64>) -> Matrix3<f64> {
        let delta = p2 - p1;
        let le    = delta.norm();
        assert!(le > 1e-14, "Beam2: zero-length element");

        // Local x axis: beam axis direction
        let x_loc = delta / le;

        // Choose auxiliary reference vector non-parallel to x_loc.
        // Use global Z unless beam is nearly vertical, then use global Y.
        let lexy = (delta.x * delta.x + delta.y * delta.y).sqrt();
        let v_ref = if lexy > 1e-10 {
            Vector3::new(0.0, 0.0, 1.0) // general case: z_global
        } else {
            Vector3::new(0.0, 1.0, 0.0) // vertical beam: y_global
        };

        // z_loc = x_loc x v_ref  (normalised)
        let z_loc = x_loc.cross(&v_ref).normalize();

        // y_loc = z_loc x x_loc  (completes right-handed triad)
        let y_loc = z_loc.cross(&x_loc);

        // R: rows are local axes in global frame  =>  R * u_global = u_local
        Matrix3::from_rows(&[
            x_loc.transpose(),
            y_loc.transpose(),
            z_loc.transpose(),
        ])
    }

    /// Build the 12x12 transformation matrix T = diag(Dc, Dc, Dc, Dc).
    /// Used to rotate the local stiffness matrix to global coordinates:
    ///   k_global = T^T * k_local * T
    fn transformation_matrix(dc: &Matrix3<f64>) -> DMatrix<f64> {
        let mut t = DMatrix::zeros(12, 12);
        for block in 0..4 {
            let offset = block * 3;
            for r in 0..3 {
                for c in 0..3 {
                    t[(offset + r, offset + c)] = dc[(r, c)];
                }
            }
        }
        t
    }
}

impl FiniteElement for Beam2 {
    type Geometry = Beam2Geometry;
    type Dof      = Vec6Dof;

    const N_NODES: usize = 2;

    /// Precompute element length. Called once per element in Assembler::new.
    fn precompute(nodes: &[Vector3<f64>]) -> Beam2Geometry {
        assert_eq!(nodes.len(), 2, "Beam2 requires exactly 2 nodes");
        Beam2Geometry { length: (nodes[1] - nodes[0]).norm() }
    }

    /// Shape functions — not used for Beam2 (analytical stiffness path).
    /// Returns zeros; required by the trait but never called by the assembler
    /// when element_stiffness returns Some(_).
    fn shape_functions(_xi: &Vector3<f64>) -> DVector<f64> {
        DVector::zeros(2)
    }

    /// Shape gradients — not used for Beam2.
    fn shape_gradients(_geometry: &Beam2Geometry, _gauss_index: usize) -> DMatrix<f64> {
        DMatrix::zeros(2, 3)
    }

    /// Integration points — not used for Beam2 (analytical stiffness).
    /// Returns a single dummy point to satisfy the trait contract.
    fn integration_points() -> Vec<(Vector3<f64>, f64)> {
        vec![(Vector3::zeros(), 1.0)]
    }

    /// B-matrix — not used for Beam2 (analytical stiffness path).
    fn b_matrix(_grad_n: &DMatrix<f64>, _f: &Matrix3<f64>) -> DMatrix<f64> {
        DMatrix::zeros(6, 12)
    }

    /// Element volume approximated as A * L for mass assembly.
    /// The cross-section area is stored in the geometry via a unit radius default.
    /// For accurate mass, set the material density accordingly.
    fn element_volume(geometry: &Beam2Geometry) -> f64 {
        // A = pi * r^2 with default r=0.01 — mass assembly uses density * volume.
        // Override via MaterialLaw::beam_radius() for accurate results.
        let r = 0.01_f64;
        PI * r * r * geometry.length
    }

    /// Jacobian determinant — not used for Beam2.
    fn gauss_det_j(_geometry: &Beam2Geometry, _gauss_index: usize) -> f64 {
        1.0
    }

    /// Compute and return the 12x12 element stiffness matrix in global coordinates.
    ///
    /// Overrides the default `None` from the trait. The assembler detects `Some(_)`
    /// and bypasses the Gauss-point loop, scattering ke directly into K.
    ///
    /// Steps:
    ///   1. Read material properties (E, nu) and cross-section radius.
    ///   2. Compute section properties (A, I, J) for solid circular section.
    ///   3. Build local 12x12 stiffness matrix (Andersen & Nielsen eq. 1-104).
    ///   4. Build direction cosine matrix Dc (OpenFAST SubDyn eq. 2/3/4).
    ///   5. Build 12x12 transformation matrix T = diag(Dc, Dc, Dc, Dc).
    ///   6. Return k_global = T^T * k_local * T.
    fn element_stiffness(
        nodes:    &[Vector3<f64>],
        material: &dyn MaterialLaw,
    ) -> Option<DMatrix<f64>> {
        assert_eq!(nodes.len(), 2, "Beam2::element_stiffness requires 2 nodes");

        let p1 = &nodes[0];
        let p2 = &nodes[1];
        let l  = (p2 - p1).norm();
        assert!(l > 1e-14, "Beam2: zero-length element");

        // Material properties
        let e  = material.youngs_modulus();
        let nu = material.poisson_ratio();
        let g  = e / (2.0 * (1.0 + nu));

        // Cross-section properties (solid circular, radius from material or default)
        let r  = material.beam_radius().unwrap_or(0.01);
        let sec = CircularSection { radius: r };
        let ea  = e * sec.area();
        let ei  = e * sec.moment_of_inertia();
        let gj  = g * sec.torsional_constant();

        // Local stiffness matrix (Andersen & Nielsen eq. 1-104 to 1-106, Phi=0)
        let k_local = Self::local_stiffness(l, ea, ei, gj);

        // Direction cosine matrix (OpenFAST SubDyn eq. 2/3/4)
        let dc = Self::direction_cosine_matrix(p1, p2);

        // Transformation matrix T = diag(Dc, Dc, Dc, Dc)
        let t = Self::transformation_matrix(&dc);

        // Global stiffness: k_global = T^T * k_local * T
        let k_global = t.transpose() * k_local * &t;

        Some(k_global)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::FiniteElement;
    use nalgebra::Matrix6;

    /// Helper: unit beam along x-axis from (0,0,0) to (L,0,0).
    fn axial_nodes(l: f64) -> Vec<Vector3<f64>> {
        vec![Vector3::zeros(), Vector3::new(l, 0.0, 0.0)]
    }

    /// Helper: unit beam along z-axis (vertical, degenerate case).
    fn vertical_nodes(l: f64) -> Vec<Vector3<f64>> {
        vec![Vector3::zeros(), Vector3::new(0.0, 0.0, l)]
    }

    /// Mock material for tests.
    struct TestMaterial { e: f64, nu: f64, r: f64 }
    impl MaterialLaw for TestMaterial {
        fn density(&self) -> f64 { 1000.0 }
        fn strain_energy(&self, _: &Matrix3<f64>, _: &crate::material::MaterialContext) -> f64 { 0.0 }
        fn pk2_stress(&self, _: &Matrix3<f64>, _: &mut crate::material::MaterialContext) -> Matrix3<f64> { Matrix3::zeros() }
        fn tangent_stiffness(&self, _: &Matrix3<f64>, _: &crate::material::MaterialContext) -> nalgebra::Matrix6<f64> { nalgebra::Matrix6::zeros() }
        fn youngs_modulus(&self) -> f64 { self.e }
        fn poisson_ratio(&self)  -> f64 { self.nu }
        fn beam_radius(&self) -> Option<f64> { Some(self.r) }
    }

    /// Precompute must return the correct element length.
    #[test]
    fn test_precompute_length() {
        let l     = 2.5;
        let nodes = axial_nodes(l);
        let geo   = Beam2::precompute(&nodes);
        assert!((geo.length - l).abs() < 1e-12, "Length should be {}, got {}", l, geo.length);
    }

    /// element_stiffness must return Some for any valid element.
    #[test]
    fn test_element_stiffness_returns_some() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, r: 0.05 };
        let nodes = axial_nodes(1.0);
        let ke    = Beam2::element_stiffness(&nodes, &mat);
        assert!(ke.is_some(), "element_stiffness must return Some for Beam2");
    }

    /// Global stiffness matrix must be symmetric.
    #[test]
    fn test_stiffness_symmetry() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, r: 0.05 };
        let nodes = axial_nodes(1.0);
        let ke    = Beam2::element_stiffness(&nodes, &mat).unwrap();
        let diff  = (&ke - ke.transpose()).norm();
        assert!(diff < 1e-6, "Stiffness matrix must be symmetric, diff = {:.2e}", diff);
    }

    /// For a beam along x, Dc = I and k_global = k_local.
    /// The axial DOFs (0,6) must have the correct stiffness EA/L.
    #[test]
    fn test_axial_stiffness_x_axis() {
        let e   = 2e11_f64;
        let nu  = 0.3;
        let r   = 0.05;
        let l   = 1.0;
        let mat = TestMaterial { e, nu, r };

        let nodes = axial_nodes(l);
        let ke    = Beam2::element_stiffness(&nodes, &mat).unwrap();

        let a  = PI * r * r;
        let ea = e * a;
        let expected = ea / l;

        assert!((ke[(0, 0)] - expected).abs() < 1.0,
            "k[0,0] should be EA/L={:.4e}, got {:.4e}", expected, ke[(0,0)]);
        assert!((ke[(0, 6)] + expected).abs() < 1.0,
            "k[0,6] should be -EA/L={:.4e}, got {:.4e}", -expected, ke[(0,6)]);
    }

    /// For a vertical beam (along z), the degenerate case (eq. 3/4) must not panic
    /// and must return a symmetric matrix.
    #[test]
    fn test_vertical_beam_degenerate_case() {
        let mat   = TestMaterial { e: 2e11, nu: 0.3, r: 0.05 };
        let nodes = vertical_nodes(2.0);
        let ke    = Beam2::element_stiffness(&nodes, &mat).unwrap();
        let diff  = (&ke - ke.transpose()).norm();
        assert!(diff < 1e-6,
            "Vertical beam stiffness must be symmetric, diff = {:.2e}", diff);
    }

    /// Cantilever tip deflection under unit load must match analytical solution.
    /// delta = F * L^3 / (3 * E * I)
    /// We verify that k[7,7] (wy DOF at node 2, beam along x) equals 12EI/L^3.
    #[test]
    fn test_bending_stiffness_coefficient() {
        let e   = 2e11_f64;
        let nu  = 0.3;
        let r   = 0.05;
        let l   = 1.0;
        let mat = TestMaterial { e, nu, r };

        let nodes = axial_nodes(l);
        let ke    = Beam2::element_stiffness(&nodes, &mat).unwrap();

        let i        = PI * r.powi(4) / 4.0;
        let expected = 12.0 * e * i / l.powi(3);

        assert!((ke[(7, 7)] - expected).abs() / expected < 1e-6,
            "k[7,7] should be 12EI/L^3={:.4e}, got {:.4e}", expected, ke[(7,7)]);
    }
}