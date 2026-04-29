// UTF-8
// orfas-core/src/element/hex8.rs
//
// Trilinear hexahedral element with 8 nodes (1 per corner).
// Reference domain: (r, s, t) in [-1, 1]^3.
// Integration: 2x2x2 Gauss quadrature (8 points, weight=1 each).
//
// Node ordering (VTK convention):
//   0: (-1,-1,-1)   4: (-1,-1,+1)
//   1: (+1,-1,-1)   5: (+1,-1,+1)
//   2: (+1,+1,-1)   6: (+1,+1,+1)
//   3: (-1,+1,-1)   7: (-1,+1,+1)
//
// Shape functions:
//   Ni(r,s,t) = 1/8 * (1 + ri*r)(1 + si*s)(1 + ti*t)
// where (ri, si, ti) are the reference coordinates of node i.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement, Vec3Dof};

// ---------------------------------------------------------------------------
// Gauss quadrature constants
// ---------------------------------------------------------------------------

/// Gauss point position along one axis for 2-point rule.
/// g = 1/sqrt(3) — exact for polynomials up to degree 3.
const G: f64 = 0.577_350_269_189_626;

/// Node reference coordinates (ri, si, ti) in [-1,+1]^3.
/// VTK ordering: bottom face (t=-1) CCW then top face (t=+1) CCW.
const NODE_RST: [[f64; 3]; 8] = [
    [-1.0, -1.0, -1.0], // 0
    [ 1.0, -1.0, -1.0], // 1
    [ 1.0,  1.0, -1.0], // 2
    [-1.0,  1.0, -1.0], // 3
    [-1.0, -1.0,  1.0], // 4
    [ 1.0, -1.0,  1.0], // 5
    [ 1.0,  1.0,  1.0], // 6
    [-1.0,  1.0,  1.0], // 7
];

// ---------------------------------------------------------------------------
// Hex8Geometry
// ---------------------------------------------------------------------------

/// Precomputed geometry for a Hex8 element.
/// Unlike Tet4, the Jacobian varies across the element — it is evaluated
/// and stored at each of the 8 Gauss points during precompute.
/// gauss_grad[g] is the (8 x 3) physical-space gradient matrix at Gauss point g.
/// det_j[g] is the Jacobian determinant at Gauss point g.
pub struct Hex8Geometry {
    /// Physical-space shape function gradients at each Gauss point.
    /// gauss_grad[g][(i, k)] = dNi/dxk evaluated at Gauss point g.
    pub gauss_grad: Vec<DMatrix<f64>>,
    /// Jacobian determinant at each Gauss point.
    /// Used to scale the integrand: integral ~ sum_g f(g) * det_J(g) * w(g).
    pub det_j: Vec<f64>,
}

impl ElementGeometry for Hex8Geometry {}

// ---------------------------------------------------------------------------
// Hex8 element
// ---------------------------------------------------------------------------

/// Trilinear hexahedron with 8 nodes.
/// More accurate than Tet4 for structured meshes — no shear locking
/// for pure translational DOFs with 2x2x2 quadrature.
/// Preferred over Tet4 when the geometry permits regular hexahedral meshing.
pub struct Hex8;

impl Hex8 {
    /// Evaluate the 8 shape functions at reference point (r, s, t).
    /// Ni = 1/8 * (1 + ri*r)(1 + si*s)(1 + ti*t).
    fn shape_functions_rst(r: f64, s: f64, t: f64) -> [f64; 8] {
        let mut n = [0.0f64; 8];
        for i in 0..8 {
            let [ri, si, ti] = NODE_RST[i];
            n[i] = 0.125 * (1.0 + ri * r) * (1.0 + si * s) * (1.0 + ti * t);
        }
        n
    }

    /// Evaluate shape function gradients w.r.t. reference coordinates (r, s, t).
    /// Returns an (8 x 3) matrix: row i = [dNi/dr, dNi/ds, dNi/dt].
    /// dNi/dr = 1/8 * ri*(1 + si*s)(1 + ti*t), and cyclically for s and t.
    fn shape_grad_rst(r: f64, s: f64, t: f64) -> DMatrix<f64> {
        let mut grad = DMatrix::zeros(8, 3);
        for i in 0..8 {
            let [ri, si, ti] = NODE_RST[i];
            grad[(i, 0)] = 0.125 * ri * (1.0 + si * s) * (1.0 + ti * t); // dNi/dr
            grad[(i, 1)] = 0.125 * si * (1.0 + ri * r) * (1.0 + ti * t); // dNi/ds
            grad[(i, 2)] = 0.125 * ti * (1.0 + ri * r) * (1.0 + si * s); // dNi/dt
        }
        grad
    }

    /// Compute physical-space gradients and Jacobian determinant at one point.
    /// J[dim, k] = sum_i dNi/drk * x_i[dim]  (3x3 Jacobian matrix).
    /// grad_phys_i = J^-T * [dNi/dr, dNi/ds, dNi/dt].
    fn compute_gauss_grad(
        nodes: &[Vector3<f64>],
        r: f64, s: f64, t: f64,
    ) -> (DMatrix<f64>, f64) {
        let grad_ref = Self::shape_grad_rst(r, s, t);

        // Build the 3x3 Jacobian matrix: J[dim, k] = sum_i dNi/drk * x_i[dim]
        let mut j = Matrix3::<f64>::zeros();
        for i in 0..8 {
            let p = nodes[i];
            for dim in 0..3 {
                for k in 0..3 {
                    j[(dim, k)] += grad_ref[(i, k)] * p[dim];
                }
            }
        }

        let det_j = j.determinant();
        assert!(det_j > 0.0, "Hex8: negative or zero Jacobian — check node ordering");

        let j_inv = j.try_inverse()
            .expect("Hex8: singular Jacobian — degenerate element");

        // Physical gradients: grad_phys_i = J^-T * grad_ref_i
        let mut grad_phys = DMatrix::zeros(8, 3);
        for i in 0..8 {
            let g_ref = Vector3::new(grad_ref[(i, 0)], grad_ref[(i, 1)], grad_ref[(i, 2)]);
            let g_phys = j_inv.transpose() * g_ref;
            grad_phys[(i, 0)] = g_phys.x;
            grad_phys[(i, 1)] = g_phys.y;
            grad_phys[(i, 2)] = g_phys.z;
        }

        (grad_phys, det_j)
    }
}

impl FiniteElement for Hex8 {
    type Geometry = Hex8Geometry;
    type Dof      = Vec3Dof;

    const N_NODES: usize = 8;

    /// Precompute physical-space gradients and Jacobian determinants
    /// at each of the 8 Gauss points. Called once per element in Assembler::new.
    fn precompute(nodes: &[Vector3<f64>]) -> Hex8Geometry {
        assert_eq!(nodes.len(), 8, "Hex8 requires exactly 8 nodes");

        let gauss_pts      = Self::integration_points();
        let mut gauss_grad = Vec::with_capacity(8);
        let mut det_j      = Vec::with_capacity(8);

        for (xi, _w) in &gauss_pts {
            let (grad, dj) = Self::compute_gauss_grad(nodes, xi.x, xi.y, xi.z);
            gauss_grad.push(grad);
            det_j.push(dj);
        }

        Hex8Geometry { gauss_grad, det_j }
    }

    /// Trilinear shape functions at xi = (r, s, t) in [-1,+1]^3.
    /// Returns a vector of length 8.
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64> {
        DVector::from_row_slice(&Self::shape_functions_rst(xi.x, xi.y, xi.z))
    }

    /// Physical-space gradients at Gauss point gauss_index.
    /// Returns an (8 x 3) matrix: row i = [dNi/dx, dNi/dy, dNi/dz].
    fn shape_gradients(geometry: &Hex8Geometry, gauss_index: usize) -> DMatrix<f64> {
        geometry.gauss_grad[gauss_index].clone()
    }

    /// 2x2x2 Gauss quadrature rule.
    /// 8 points at (+-G, +-G, +-G) with weight 1.0 each.
    /// Exact for trilinear functions on undistorted elements.
    /// xi is stored as Vector3(r, s, t) for compatibility with the FiniteElement trait.
    fn integration_points() -> Vec<(Vector3<f64>, f64)> {
        let mut pts = Vec::with_capacity(8);
        for &r in &[-G, G] {
            for &s in &[-G, G] {
                for &t in &[-G, G] {
                    pts.push((Vector3::new(r, s, t), 1.0));
                }
            }
        }
        pts
    }

    /// Build the 6x24 B matrix from physical-space shape function gradients.
    /// Same Voigt convention as Tet4: [eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz].
    /// F is unused — Hex8 uses the same linear B-matrix as Tet4 (updated Lagrangian
    /// formulation is deferred; nonlinear effects enter via material tangent only).
    fn b_matrix(grad_n: &DMatrix<f64>, _f: &Matrix3<f64>) -> DMatrix<f64> {
        let n_nodes = grad_n.nrows(); // 8
        let n_dofs  = n_nodes * 3;   // 24
        let mut b   = DMatrix::zeros(6, n_dofs);

        for i in 0..n_nodes {
            let bi  = grad_n[(i, 0)]; // dNi/dx
            let ci  = grad_n[(i, 1)]; // dNi/dy
            let di  = grad_n[(i, 2)]; // dNi/dz
            let col = i * 3;

            b[(0, col)]     = bi; // eps_xx = dux/dx
            b[(1, col + 1)] = ci; // eps_yy = duy/dy
            b[(2, col + 2)] = di; // eps_zz = duz/dz

            b[(3, col)]     = ci; // eps_xy: dux/dy
            b[(3, col + 1)] = bi; // eps_xy: duy/dx
            b[(4, col + 1)] = di; // eps_yz: duy/dz
            b[(4, col + 2)] = ci; // eps_yz: duz/dy
            b[(5, col)]     = di; // eps_xz: dux/dz
            b[(5, col + 2)] = bi; // eps_xz: duz/dx
        }
        b
    }

    /// Total element volume = sum_g det_J(g) * w(g).
    /// For 2x2x2 Gauss with w=1, this is just the sum of det_J values.
    fn element_volume(geometry: &Hex8Geometry) -> f64 {
        geometry.det_j.iter().sum()
    }

    /// Jacobian determinant at Gauss point gauss_index.
    /// Already scaled correctly — no extra factor needed (unlike Tet conventions).
    fn gauss_det_j(geometry: &Hex8Geometry, gauss_index: usize) -> f64 {
        geometry.det_j[gauss_index]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Unit cube node positions in VTK order.
    fn unit_hex8_nodes() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0), // 0
            Vector3::new(1.0, 0.0, 0.0), // 1
            Vector3::new(1.0, 1.0, 0.0), // 2
            Vector3::new(0.0, 1.0, 0.0), // 3
            Vector3::new(0.0, 0.0, 1.0), // 4
            Vector3::new(1.0, 0.0, 1.0), // 5
            Vector3::new(1.0, 1.0, 1.0), // 6
            Vector3::new(0.0, 1.0, 1.0), // 7
        ]
    }

    /// Shape functions must sum to 1 at any reference point (partition of unity).
    #[test]
    fn test_partition_of_unity() {
        let test_pts = vec![
            Vector3::new( 0.0,  0.0,  0.0),
            Vector3::new( 0.5,  0.3, -0.2),
            Vector3::new(-1.0, -1.0, -1.0),
            Vector3::new( 1.0,  1.0,  1.0),
        ];
        for xi in test_pts {
            let n: f64 = Hex8::shape_functions(&xi).iter().sum();
            assert!((n - 1.0).abs() < 1e-14, "Partition of unity failed at {:?}: sum={}", xi, n);
        }
    }

    /// Shape function Ni must equal 1 at node i and 0 at all other nodes.
    #[test]
    fn test_kronecker_delta() {
        for i in 0..8 {
            let [ri, si, ti] = NODE_RST[i];
            let xi = Vector3::new(ri, si, ti);
            let n  = Hex8::shape_functions(&xi);
            assert!((n[i] - 1.0).abs() < 1e-14, "N{}({:?}) should be 1, got {}", i, xi, n[i]);
            for j in 0..8 {
                if j != i {
                    assert!(n[j].abs() < 1e-14, "N{}({:?}) should be 0, got {}", j, xi, n[j]);
                }
            }
        }
    }

    /// Volume of a unit cube must be 1.0.
    #[test]
    fn test_unit_cube_volume() {
        let nodes = unit_hex8_nodes();
        let geo   = Hex8::precompute(&nodes);
        let vol   = Hex8::element_volume(&geo);
        assert!((vol - 1.0).abs() < 1e-12, "Unit cube volume should be 1.0, got {}", vol);
    }

    /// Sum of shape function gradients over all nodes must be zero (completeness).
    /// sum_i dNi/dxk = 0 for k=0,1,2 at any Gauss point.
    #[test]
    fn test_gradient_sum_zero() {
        let nodes = unit_hex8_nodes();
        let geo   = Hex8::precompute(&nodes);
        for g in 0..8 {
            let grad = Hex8::shape_gradients(&geo, g);
            for col in 0..3 {
                let sum: f64 = (0..8).map(|i| grad[(i, col)]).sum();
                assert!(
                    sum.abs() < 1e-12,
                    "Gradient sum col {} at Gauss {} should be 0, got {}",
                    col, g, sum
                );
            }
        }
    }

    /// For a rectangular element (lx, ly, lz), volume must equal lx*ly*lz.
    #[test]
    fn test_scaled_element_volume() {
        let lx = 2.0; let ly = 3.0; let lz = 0.5;
        let nodes = vec![
            Vector3::new( 0.0,  0.0,  0.0),
            Vector3::new(  lx,  0.0,  0.0),
            Vector3::new(  lx,   ly,  0.0),
            Vector3::new( 0.0,   ly,  0.0),
            Vector3::new( 0.0,  0.0,   lz),
            Vector3::new(  lx,  0.0,   lz),
            Vector3::new(  lx,   ly,   lz),
            Vector3::new( 0.0,   ly,   lz),
        ];
        let geo = Hex8::precompute(&nodes);
        let vol = Hex8::element_volume(&geo);
        assert!(
            (vol - lx * ly * lz).abs() < 1e-10,
            "Volume should be {}, got {}", lx * ly * lz, vol
        );
    }
}