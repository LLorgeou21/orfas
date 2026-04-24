// UTF-8
// assembler/geometry.rs — element geometry, BMatrix trait, deformation gradient.

use nalgebra::{Matrix3, SMatrix, Vector3};

pub type Matrix6x12 = SMatrix<f64, 6, 12>;

/// Reference geometry cached per element.
pub struct ElementGeometry {
    pub volume: f64,
    pub b:      [f64; 4],
    pub c:      [f64; 4],
    pub d:      [f64; 4],
}

/// Volume of a tetrahedron: V = |det([p1-p0, p2-p0, p3-p0])| / 6
pub fn tetra_volume(
    p0: &Vector3<f64>, p1: &Vector3<f64>,
    p2: &Vector3<f64>, p3: &Vector3<f64>,
) -> f64 {
    Matrix3::from_columns(&[p1 - p0, p2 - p0, p3 - p0]).determinant().abs() / 6.0
}

/// Shape function gradients (b, c, d) for a linear tetrahedron.
pub fn tetra_bcd(
    p0: &Vector3<f64>, p1: &Vector3<f64>,
    p2: &Vector3<f64>, p3: &Vector3<f64>,
) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let b1 = (p1.y - p3.y) * (p2.z - p3.z) - (p2.y - p3.y) * (p1.z - p3.z);
    let b2 = (p2.y - p3.y) * (p0.z - p3.z) - (p0.y - p3.y) * (p2.z - p3.z);
    let b3 = (p0.y - p3.y) * (p1.z - p3.z) - (p1.y - p3.y) * (p0.z - p3.z);
    let b4 = -(b1 + b2 + b3);

    let c1 = (p2.x - p3.x) * (p1.z - p3.z) - (p1.x - p3.x) * (p2.z - p3.z);
    let c2 = (p0.x - p3.x) * (p2.z - p3.z) - (p2.x - p3.x) * (p0.z - p3.z);
    let c3 = (p1.x - p3.x) * (p0.z - p3.z) - (p0.x - p3.x) * (p1.z - p3.z);
    let c4 = -(c1 + c2 + c3);

    let d1 = (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
    let d2 = (p2.x - p3.x) * (p0.y - p3.y) - (p0.x - p3.x) * (p2.y - p3.y);
    let d3 = (p0.x - p3.x) * (p1.y - p3.y) - (p1.x - p3.x) * (p0.y - p3.y);
    let d4 = -(d1 + d2 + d3);

    ([b1, b2, b3, b4], [c1, c2, c3, c4], [d1, d2, d3, d4])
}

/// B matrix (6x12) in Voigt notation from shape function gradients and volume.
pub fn tetra_b_matrix(b: &[f64; 4], c: &[f64; 4], d: &[f64; 4], volume: f64) -> Matrix6x12 {
    Matrix6x12::from_row_slice(&[
        b[0], 0.0,  0.0,  b[1], 0.0,  0.0,  b[2], 0.0,  0.0,  b[3], 0.0,  0.0,
        0.0,  c[0], 0.0,  0.0,  c[1], 0.0,  0.0,  c[2], 0.0,  0.0,  c[3], 0.0,
        0.0,  0.0,  d[0], 0.0,  0.0,  d[1], 0.0,  0.0,  d[2], 0.0,  0.0,  d[3],
        c[0], b[0], 0.0,  c[1], b[1], 0.0,  c[2], b[2], 0.0,  c[3], b[3], 0.0,
        0.0,  d[0], c[0], 0.0,  d[1], c[1], 0.0,  d[2], c[2], 0.0,  d[3], c[3],
        d[0], 0.0,  b[0], d[1], 0.0,  b[1], d[2], 0.0,  b[2], d[3], 0.0,  b[3],
    ]).scale(1.0 / (6.0 * volume))
}

/// Deformation gradient F = I + sum_i ui x gradNi for a linear tetrahedron.
pub fn compute_deformation_gradient(
    u0: &Vector3<f64>, u1: &Vector3<f64>,
    u2: &Vector3<f64>, u3: &Vector3<f64>,
    b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
) -> Matrix3<f64> {
    let mut f = Matrix3::identity();
    for i in 0..4 {
        let u = [u0, u1, u2, u3][i];
        let grad_n = Vector3::new(b[i], c[i], d[i]);
        f += u * grad_n.transpose();
    }
    f
}

/// Trait defining the strain-displacement B matrix computation.
pub trait BMatrix {
    fn compute(
        b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
        volume: f64,
        f: &Matrix3<f64>,
    ) -> Matrix6x12;
}

/// Linear B matrix — constant, F ignored.
pub struct LinearBMatrix;

impl BMatrix for LinearBMatrix {
    fn compute(
        b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
        volume: f64,
        _f: &Matrix3<f64>,
    ) -> Matrix6x12 {
        tetra_b_matrix(b, c, d, volume)
    }
}