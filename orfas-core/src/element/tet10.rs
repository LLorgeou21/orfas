// UTF-8
// orfas-core/src/element/tet10.rs
// Quadratic tetrahedron with 10 nodes (4 corner + 6 midside).
// Uses 4-point Keast quadrature (exact for polynomials up to degree 2).
//
// Node ordering follows basix convention:
//   0,1,2,3 — corner nodes
//   4 — midside of edge (2,3)
//   5 — midside of edge (1,3)
//   6 — midside of edge (1,2)
//   7 — midside of edge (0,3)
//   8 — midside of edge (0,2)
//   9 — midside of edge (0,1)
//
// Barycentric coordinates: xi = (L2, L3, L4), L1 = 1 - L2 - L3 - L4.
// Jacobian convention: J[dim, k] = sum_i dNi/drk * xi[dim]
// where r = (L2, L3, L4) are independent variables and
// dNi/drk = dNi/dLk - dNi/dL1 (chain rule with L1 = 1 - r - s - t).
// This convention matches basix and gives det_J = 1 for the unit tetrahedron.

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement};

// ---------------------------------------------------------------------------
// Keast 4-point quadrature rule for tetrahedron
// ---------------------------------------------------------------------------
// Keast (1986), CMAME. Points are permutations of (b,a,a,a) in barycentric.
// Weights w = 1/4 sum to 1 — must be divided by 6 for volume scaling.
// a = (5 - sqrt(5)) / 20, b = (5 + 3*sqrt(5)) / 20.

const KEAST4_A: f64 = 0.138_196_601_125_011;
const KEAST4_B: f64 = 0.585_410_196_624_969;
const KEAST4_W: f64 = 0.25;

// ---------------------------------------------------------------------------
// Tet10Geometry
// ---------------------------------------------------------------------------

/// Precomputed geometry for a quadratic tetrahedron (Tet10).
/// Gradients vary inside the element — evaluated at each Gauss point.
/// gauss_grad[g] is the (10 x 3) physical-space gradient matrix at point g.
/// det_j[g] is the jacobian determinant at point g.
pub struct Tet10Geometry {
    pub gauss_grad: Vec<DMatrix<f64>>,
    pub det_j:      Vec<f64>,
}

impl ElementGeometry for Tet10Geometry {}

// ---------------------------------------------------------------------------
// Tet10 element
// ---------------------------------------------------------------------------

/// Quadratic tetrahedron with 10 nodes.
/// Eliminates shear locking present in Tet4 for bending-dominated problems.
/// Validated against basix (FEniCS): gradients match to machine precision.
pub struct Tet10;

impl Tet10 {
    /// Shape functions in barycentric coordinates (L1,L2,L3,L4).
    /// Corner nodes: Ni = Li*(2*Li-1). Midside nodes: Nij = 4*Li*Lj.
    /// Node-to-barycentric mapping follows basix ordering.
    fn shape_functions_bary(l1: f64, l2: f64, l3: f64, l4: f64) -> [f64; 10] {
        [
            l1 * (2.0 * l1 - 1.0), // N0 — corner L1
            l2 * (2.0 * l2 - 1.0), // N1 — corner L2
            l3 * (2.0 * l3 - 1.0), // N2 — corner L3
            l4 * (2.0 * l4 - 1.0), // N3 — corner L4
            4.0 * l3 * l4,          // N4 — mid(2,3)
            4.0 * l2 * l4,          // N5 — mid(1,3)
            4.0 * l2 * l3,          // N6 — mid(1,2)
            4.0 * l1 * l4,          // N7 — mid(0,3)
            4.0 * l1 * l3,          // N8 — mid(0,2)
            4.0 * l1 * l2,          // N9 — mid(0,1)
        ]
    }

    /// Gradients of shape functions w.r.t. barycentric coords (L1,L2,L3,L4).
    /// Row i = [dNi/dL1, dNi/dL2, dNi/dL3, dNi/dL4].
    fn shape_grad_bary(l1: f64, l2: f64, l3: f64, l4: f64) -> [[f64; 4]; 10] {
        [
            [4.0*l1 - 1.0,  0.0,           0.0,           0.0          ], // N0
            [0.0,           4.0*l2 - 1.0,  0.0,           0.0          ], // N1
            [0.0,           0.0,           4.0*l3 - 1.0,  0.0          ], // N2
            [0.0,           0.0,           0.0,           4.0*l4 - 1.0 ], // N3
            [0.0,           0.0,           4.0*l4,        4.0*l3       ], // N4 mid(2,3)
            [0.0,           4.0*l4,        0.0,           4.0*l2       ], // N5 mid(1,3)
            [0.0,           4.0*l3,        4.0*l2,        0.0          ], // N6 mid(1,2)
            [4.0*l4,        0.0,           0.0,           4.0*l1       ], // N7 mid(0,3)
            [4.0*l3,        0.0,           4.0*l1,        0.0          ], // N8 mid(0,2)
            [4.0*l2,        4.0*l1,        0.0,           0.0          ], // N9 mid(0,1)
        ]
    }

    /// Compute physical-space gradients and jacobian determinant at one point.
    /// Uses r=(L2,L3,L4) as independent variables with L1 = 1-r-s-t.
    /// dNi/dr = dNi/dL2 - dNi/dL1 (chain rule).
    /// J[dim, k] = sum_i dNi/drk * xi[dim].
    /// grad_phys_i = J^-T * [dNi/dr, dNi/ds, dNi/dt].
    fn compute_gauss_grad(
        nodes: &[Vector3<f64>],
        l1: f64, l2: f64, l3: f64, l4: f64,
    ) -> (DMatrix<f64>, f64) {
        let gb = Self::shape_grad_bary(l1, l2, l3, l4);

        // Build jacobian J[dim, k] using the chain rule dNi/drk = dNi/dLk - dNi/dL1
        let mut j = Matrix3::<f64>::zeros();
        for i in 0..10 {
            let p     = nodes[i];
            let dn_dr = gb[i][1] - gb[i][0]; // dNi/dL2 - dNi/dL1
            let dn_ds = gb[i][2] - gb[i][0]; // dNi/dL3 - dNi/dL1
            let dn_dt = gb[i][3] - gb[i][0]; // dNi/dL4 - dNi/dL1

            j[(0, 0)] += dn_dr * p.x;
            j[(0, 1)] += dn_ds * p.x;
            j[(0, 2)] += dn_dt * p.x;
            j[(1, 0)] += dn_dr * p.y;
            j[(1, 1)] += dn_ds * p.y;
            j[(1, 2)] += dn_dt * p.y;
            j[(2, 0)] += dn_dr * p.z;
            j[(2, 1)] += dn_ds * p.z;
            j[(2, 2)] += dn_dt * p.z;
        }

        let det_j = j.determinant();
        let j_inv = j.try_inverse()
            .expect("Tet10: singular jacobian — degenerate element");

        let mut grad_phys = DMatrix::zeros(10, 3);
        for i in 0..10 {
            let dn_ref = Vector3::new(
                gb[i][1] - gb[i][0],
                gb[i][2] - gb[i][0],
                gb[i][3] - gb[i][0],
            );
            let gp = j_inv.transpose() * dn_ref;
            grad_phys[(i, 0)] = gp.x;
            grad_phys[(i, 1)] = gp.y;
            grad_phys[(i, 2)] = gp.z;
        }

        (grad_phys, det_j)
    }
}

impl FiniteElement for Tet10 {
    type Geometry = Tet10Geometry;
    type Dof = crate::element::traits::Vec3Dof;

    const N_NODES: usize = 10;

    /// Precompute physical-space gradients at each of the 4 Keast points.
    /// Called once per element in Assembler::new.
    fn precompute(nodes: &[Vector3<f64>]) -> Tet10Geometry {
        assert_eq!(nodes.len(), 10, "Tet10 requires exactly 10 nodes");

        let gauss_pts      = Self::integration_points();
        let mut gauss_grad = Vec::with_capacity(4);
        let mut det_j      = Vec::with_capacity(4);

        for (xi, _w) in &gauss_pts {
            let l2 = xi.x;
            let l3 = xi.y;
            let l4 = xi.z;
            let l1 = 1.0 - l2 - l3 - l4;
            let (grad, dj) = Self::compute_gauss_grad(nodes, l1, l2, l3, l4);
            gauss_grad.push(grad);
            det_j.push(dj);
        }

        Tet10Geometry { gauss_grad, det_j }
    }

    /// Quadratic shape functions at xi = (L2, L3, L4).
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64> {
        let l2 = xi.x; let l3 = xi.y; let l4 = xi.z;
        let l1 = 1.0 - l2 - l3 - l4;
        DVector::from_row_slice(&Self::shape_functions_bary(l1, l2, l3, l4))
    }

    /// Physical-space gradients at Gauss point gauss_index.
    /// Returns a (10 x 3) matrix: row i = [dNi/dx, dNi/dy, dNi/dz].
    fn shape_gradients(geometry: &Tet10Geometry, gauss_index: usize) -> DMatrix<f64> {
        geometry.gauss_grad[gauss_index].clone()
    }

    /// 4-point Keast quadrature. Points are permutations of (b,a,a,a)
    /// expressed as xi = (L2,L3,L4). Weights sum to 1 — volume scaling
    /// is handled by gauss_det_j (divided by 6).
    fn integration_points() -> Vec<(Vector3<f64>, f64)> {
        vec![
            (Vector3::new(KEAST4_B, KEAST4_A, KEAST4_A), KEAST4_W), // L2=b
            (Vector3::new(KEAST4_A, KEAST4_B, KEAST4_A), KEAST4_W), // L3=b
            (Vector3::new(KEAST4_A, KEAST4_A, KEAST4_B), KEAST4_W), // L4=b
            (Vector3::new(KEAST4_A, KEAST4_A, KEAST4_A), KEAST4_W), // L1=b
        ]
    }

    /// Build the 6x30 B matrix. Voigt convention, identical structure to Tet4.
    fn b_matrix(grad_n: &DMatrix<f64>, _f: &Matrix3<f64>) -> DMatrix<f64> {
        let n_nodes = grad_n.nrows();
        let n_dofs  = n_nodes * 3;
        let mut b   = DMatrix::zeros(6, n_dofs);

        for i in 0..n_nodes {
            let bi  = grad_n[(i, 0)];
            let ci  = grad_n[(i, 1)];
            let di  = grad_n[(i, 2)];
            let col = i * 3;

            b[(0, col)]     = bi; // eps_xx
            b[(1, col + 1)] = ci; // eps_yy
            b[(2, col + 2)] = di; // eps_zz

            b[(3, col)]     = ci; // eps_xy: dN/dy
            b[(3, col + 1)] = bi; // eps_xy: dN/dx
            b[(4, col + 1)] = di; // eps_yz: dN/dz
            b[(4, col + 2)] = ci; // eps_yz: dN/dy
            b[(5, col)]     = di; // eps_xz: dN/dz
            b[(5, col + 2)] = bi; // eps_xz: dN/dx
        }
        b
    }

    /// Total element volume = sum_g det_J(g) * w(g) / 6.
    /// The /6 factor accounts for the barycentric coordinate normalization
    /// (Keast weights sum to 1, reference tetrahedron volume = 1/6).
    fn element_volume(geometry: &Tet10Geometry) -> f64 {
        let gauss_pts = Self::integration_points();
        geometry.det_j.iter()
            .zip(gauss_pts.iter())
            .map(|(dj, (_xi, w))| dj * w)
            .sum::<f64>() / 6.0
    }

    /// Jacobian determinant at Gauss point gauss_index, scaled by 1/6.
    /// Consistent with Tet4::gauss_det_j = volume * 6 and weight = 1/6:
    /// both give ke = B^T C B * det_J * w with the same effective volume.
    fn gauss_det_j(geometry: &Tet10Geometry, gauss_index: usize) -> f64 {
        geometry.det_j[gauss_index] / 6.0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::subdivision::tet4_to_tet10;
    use crate::mesh::Tet4Mesh;

    fn unit_tet10_nodes() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0), // 0 corner
            Vector3::new(1.0, 0.0, 0.0), // 1 corner
            Vector3::new(0.0, 1.0, 0.0), // 2 corner
            Vector3::new(0.0, 0.0, 1.0), // 3 corner
            Vector3::new(0.0, 0.5, 0.5), // 4 mid(2,3)
            Vector3::new(0.5, 0.0, 0.5), // 5 mid(1,3)
            Vector3::new(0.5, 0.5, 0.0), // 6 mid(1,2)
            Vector3::new(0.0, 0.0, 0.5), // 7 mid(0,3)
            Vector3::new(0.0, 0.5, 0.0), // 8 mid(0,2)
            Vector3::new(0.5, 0.0, 0.0), // 9 mid(0,1)
        ]
    }

    #[test]
    fn test_shape_functions_partition_of_unity() {
        let test_points = vec![
            Vector3::new(0.25, 0.25, 0.25),
            Vector3::new(0.1,  0.2,  0.3),
            Vector3::new(0.5,  0.1,  0.1),
        ];
        for xi in test_points {
            let n   = Tet10::shape_functions(&xi);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "Partition of unity failed: sum = {}", sum);
        }
    }

    #[test]
    fn test_shape_functions_kronecker_delta() {
        let node_xi = vec![
            Vector3::new(0.0, 0.0, 0.0), // node 0: L1=1
            Vector3::new(1.0, 0.0, 0.0), // node 1: L2=1
            Vector3::new(0.0, 1.0, 0.0), // node 2: L3=1
            Vector3::new(0.0, 0.0, 1.0), // node 3: L4=1
        ];
        for (i, xi) in node_xi.iter().enumerate() {
            let n = Tet10::shape_functions(xi);
            assert!((n[i] - 1.0).abs() < 1e-12,
                "N{}({:?}) should be 1, got {}", i, xi, n[i]);
            for j in 0..10 {
                if j != i {
                    assert!(n[j].abs() < 1e-12,
                        "N{}({:?}) should be 0, got {}", j, xi, n[j]);
                }
            }
        }
    }

    #[test]
    fn test_element_volume() {
        let nodes = unit_tet10_nodes();
        let geo   = Tet10::precompute(&nodes);
        let vol   = Tet10::element_volume(&geo);
        assert!((vol - 1.0 / 6.0).abs() < 1e-10,
            "Volume should be 1/6, got {}", vol);
    }

    #[test]
    fn test_gradient_sum_zero() {
        let nodes = unit_tet10_nodes();
        let geo   = Tet10::precompute(&nodes);
        for g in 0..4 {
            let grad = Tet10::shape_gradients(&geo, g);
            for col in 0..3 {
                let sum: f64 = (0..10).map(|i| grad[(i, col)]).sum();
                assert!(sum.abs() < 1e-10,
                    "Gradient sum col {} at Gauss {} should be 0, got {}", col, g, sum);
            }
        }
    }
}