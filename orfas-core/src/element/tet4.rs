// UTF-8
// orfas-core/src/element/tet4.rs

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement};

// ---------------------------------------------------------------------------
// Tet4Geometry
// ---------------------------------------------------------------------------

/// Precomputed geometry for a linear tetrahedron (Tet4).
/// Gradients are constant within the element (one Gauss point).
/// b, c, d are the shape function gradient components derived
/// analytically from the node positions, already divided by 6V.
pub struct Tet4Geometry {
    /// Shape function gradient components along x for each of the 4 nodes.
    pub b: [f64; 4],
    /// Shape function gradient components along y for each of the 4 nodes.
    pub c: [f64; 4],
    /// Shape function gradient components along z for each of the 4 nodes.
    pub d: [f64; 4],
    /// Element volume (always positive, computed from the determinant).
    pub volume: f64,
}

impl ElementGeometry for Tet4Geometry {}

// ---------------------------------------------------------------------------
// Tet4 element
// ---------------------------------------------------------------------------

/// Linear tetrahedron with 4 nodes (CST 3D — Constant Strain Tetrahedron).
/// One Gauss point at the centroid with weight = 1/6.
/// B matrix and F are constant within the element.
pub struct Tet4;

impl FiniteElement for Tet4 {
    type Geometry = Tet4Geometry;
    type Dof = crate::element::traits::Vec3Dof;

    const N_NODES: usize = 4;

    /// Precompute analytical shape function gradients and volume.
    /// Nodes are expected in order: [n0, n1, n2, n3].
    /// Gradients are divided by 6V to get physical-space gradients.
    fn precompute(nodes: &[Vector3<f64>]) -> Tet4Geometry {
        assert_eq!(nodes.len(), 4, "Tet4 requires exactly 4 nodes");

        let [p0, p1, p2, p3] = [nodes[0], nodes[1], nodes[2], nodes[3]];

        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let e3 = p3 - p0;

        // Signed volume = (1/6) * det([e1, e2, e3])
        let signed_vol = e1.dot(&e2.cross(&e3)) / 6.0;
        let volume     = signed_vol.abs();
        let inv6v      = 1.0 / (6.0 * signed_vol);

        let b0 = ((p1.y - p3.y) * (p2.z - p3.z) - (p2.y - p3.y) * (p1.z - p3.z)) * inv6v;
        let b1 = ((p2.y - p3.y) * (p0.z - p3.z) - (p0.y - p3.y) * (p2.z - p3.z)) * inv6v;
        let b2 = ((p0.y - p3.y) * (p1.z - p3.z) - (p1.y - p3.y) * (p0.z - p3.z)) * inv6v;
        let b3 = -(b0 + b1 + b2);

        let c0 = ((p2.x - p3.x) * (p1.z - p3.z) - (p1.x - p3.x) * (p2.z - p3.z)) * inv6v;
        let c1 = ((p0.x - p3.x) * (p2.z - p3.z) - (p2.x - p3.x) * (p0.z - p3.z)) * inv6v;
        let c2 = ((p1.x - p3.x) * (p0.z - p3.z) - (p0.x - p3.x) * (p1.z - p3.z)) * inv6v;
        let c3 = -(c0 + c1 + c2);

        let d0 = ((p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)) * inv6v;
        let d1 = ((p2.x - p3.x) * (p0.y - p3.y) - (p0.x - p3.x) * (p2.y - p3.y)) * inv6v;
        let d2 = ((p0.x - p3.x) * (p1.y - p3.y) - (p1.x - p3.x) * (p0.y - p3.y)) * inv6v;
        let d3 = -(d0 + d1 + d2);

        Tet4Geometry {
            b: [b0, b1, b2, b3],
            c: [c0, c1, c2, c3],
            d: [d0, d1, d2, d3],
            volume,
        }
    }

    /// Linear shape functions in barycentric coordinates.
    /// N0 = 1 - xi - eta - zeta, N1 = xi, N2 = eta, N3 = zeta.
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64> {
        let mut n = DVector::zeros(4);
        n[0] = 1.0 - xi.x - xi.y - xi.z;
        n[1] = xi.x;
        n[2] = xi.y;
        n[3] = xi.z;
        n
    }

    /// Gradients are constant for Tet4 — gauss_index is ignored.
    /// Returns a 4x3 matrix where row i = [b_i, c_i, d_i].
    fn shape_gradients(geometry: &Tet4Geometry, _gauss_index: usize) -> DMatrix<f64> {
        let mut grad = DMatrix::zeros(4, 3);
        for i in 0..4 {
            grad[(i, 0)] = geometry.b[i];
            grad[(i, 1)] = geometry.c[i];
            grad[(i, 2)] = geometry.d[i];
        }
        grad
    }

    /// One Gauss point at the centroid of the tetrahedron.
    /// xi = (1/4, 1/4, 1/4), weight = 1/6 (volume of reference tetrahedron).
    fn integration_points() -> Vec<(Vector3<f64>, f64)> {
        vec![(Vector3::new(0.25, 0.25, 0.25), 1.0 / 6.0)]
    }

    /// Build the 6x12 B matrix from shape function gradients.
    /// For Tet4, F does not affect B (linear kinematics assumption).
    /// Voigt convention: [eps_xx, eps_yy, eps_zz, eps_xy, eps_yz, eps_xz].
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

    /// Total element volume.
    fn element_volume(geometry: &Tet4Geometry) -> f64 {
        geometry.volume
    }

    /// Jacobian determinant at the single Gauss point.
    /// For Tet4, det_J = 6V in barycentric coordinates — consistent with
    /// the Tet10 convention where det_J is the raw jacobian determinant.
    fn gauss_det_j(geometry: &Tet4Geometry, _gauss_index: usize) -> f64 {
        geometry.volume * 6.0
    }
}