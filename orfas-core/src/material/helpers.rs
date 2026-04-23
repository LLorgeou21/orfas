// UTF-8
// material/helpers.rs — shared mathematical helpers for material implementations.

use nalgebra::{Matrix3, Matrix6};

// ─── Public helpers ───────────────────────────────────────────────────────────

/// Compute the Lame parameters (lambda, mu) from Young's modulus and Poisson's ratio.
///
/// lambda = E*nu / ((1+nu)*(1-2*nu))
/// mu     = E / (2*(1+nu))
pub fn lame(youngs_modulus: f64, poisson_ratio: f64) -> (f64, f64) {
    let nu = poisson_ratio;
    let e  = youngs_modulus;
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu     = e / (2.0 * (1.0 + nu));
    (lambda, mu)
}

/// Build the isotropic tangent stiffness matrix in Voigt notation (6x6)
/// from Lame parameters lambda and mu.
///
/// C = lambda*(I tensor I) + 2*mu*I4, identical to the linear elastic tangent.
/// For SVK, this matrix is constant (independent of F).
///
/// Voigt ordering: [11, 22, 33, 12, 23, 13]
/// Engineering shear convention — no extra factors, consistent with B-matrix.
pub fn hooke_voigt(lambda: f64, mu: f64) -> Matrix6<f64> {
    let d = lambda + 2.0 * mu;
    Matrix6::new(
             d, lambda, lambda, 0.0, 0.0, 0.0,
        lambda,      d, lambda, 0.0, 0.0, 0.0,
        lambda, lambda,      d, 0.0, 0.0, 0.0,
           0.0,    0.0,    0.0,  mu, 0.0, 0.0,
           0.0,    0.0,    0.0, 0.0,  mu, 0.0,
           0.0,    0.0,    0.0, 0.0, 0.0,  mu,
    )
}

// ─── Internal helpers ─────────────────────────────────────────────────────────

/// Compute the Green-Lagrange strain tensor E = 0.5*(F^T F - I) from F.
pub(super) fn green_lagrange(f: &Matrix3<f64>) -> Matrix3<f64> {
    0.5 * (f.transpose() * f - Matrix3::identity())
}

/// Build a symmetric 4th-order tangent in Voigt notation (6x6) from two
/// scalar contributions expressed in terms of a symmetric tensor A (typically C^{-1}):
///
///   out = a * (A tensor A) + b * (A odot A)
///
/// where:
///   (A tensor A)_IJKL = A_IJ * A_KL               (outer product)
///   (A odot  A)_IJKL  = 0.5*(A_IK*A_JL + A_IL*A_JK)  (symmetrized product)
///
/// Voigt index map: 0->(0,0), 1->(1,1), 2->(2,2), 3->(0,1), 4->(1,2), 5->(0,2)
///
/// Verification at A=I, a=lambda, b=mu gives hooke_voigt(lambda, mu).
pub(super) fn cinv_tangent_voigt(a_mat: &Matrix3<f64>, a: f64, b: f64) -> Matrix6<f64> {
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let mut out = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (p, q) = idx[i];
            let (r, s) = idx[j];
            out[(i, j)] =
                a * a_mat[(p, q)] * a_mat[(r, s)]
              + b * (a_mat[(p, r)] * a_mat[(q, s)] + a_mat[(p, s)] * a_mat[(q, r)]);
        }
    }
    out
}

/// Build (I tensor A + A tensor I) in Voigt notation (6x6).
///
/// (I tensor A)_IJKL = delta_IJ * A_KL
/// (A tensor I)_IJKL = A_IJ * delta_KL
pub(super) fn i_tensor_a_plus_a_tensor_i(a_mat: &Matrix3<f64>) -> Matrix6<f64> {
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let mut out = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (p, q) = idx[i];
            let (r, s) = idx[j];
            let d_ij = if p == q { 1.0 } else { 0.0 };
            let d_kl = if r == s { 1.0 } else { 0.0 };
            out[(i, j)] = d_ij * a_mat[(r, s)] + a_mat[(p, q)] * d_kl;
        }
    }
    out
}

/// Build (I odot I) in Voigt notation (6x6).
///
/// (I odot I)_IJKL = 0.5*(delta_IK*delta_JL + delta_IL*delta_JK)
pub(super) fn i_odot_i_voigt() -> Matrix6<f64> {
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let mut out = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (p, q) = idx[i];
            let (r, s) = idx[j];
            out[(i, j)] = 0.5 * (
                (if p==r && q==s { 1.0 } else { 0.0 })
              + (if p==s && q==r { 1.0 } else { 0.0 })
            );
        }
    }
    out
}

/// Build (A tensor B + B tensor A) in Voigt notation (6x6).
///
/// (A tensor B)_IJKL = A_IJ * B_KL
pub(super) fn a_tensor_b_plus_b_tensor_a(
    a_mat: &Matrix3<f64>,
    b_mat: &Matrix3<f64>,
) -> Matrix6<f64> {
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let mut out = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (p, q) = idx[i];
            let (r, s) = idx[j];
            out[(i, j)] = a_mat[(p,q)] * b_mat[(r,s)] + b_mat[(p,q)] * a_mat[(r,s)];
        }
    }
    out
}

/// Build (A odot I + I odot A) symmetrized in Voigt notation (6x6).
///
/// (A odot I)_IJKL = 0.5*(A_IK*delta_JL + A_IL*delta_JK
///                      + delta_IK*A_JL + delta_IL*A_JK)
pub(super) fn a_odot_i_voigt(a_mat: &Matrix3<f64>) -> Matrix6<f64> {
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let mut out = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (p, q) = idx[i];
            let (r, s) = idx[j];
            out[(i, j)] = 0.5 * (
                a_mat[(p,r)] * (if q==s { 1.0 } else { 0.0 })
              + a_mat[(p,s)] * (if q==r { 1.0 } else { 0.0 })
              + (if p==r { 1.0 } else { 0.0 }) * a_mat[(q,s)]
              + (if p==s { 1.0 } else { 0.0 }) * a_mat[(q,r)]
            );
        }
    }
    out
}