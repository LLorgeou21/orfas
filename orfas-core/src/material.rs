use nalgebra::{Matrix3, Matrix6};

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Compute the Lame parameters from Young's modulus and Poisson's ratio.
fn lame(youngs_modulus: f64, poisson_ratio: f64) -> (f64, f64) {
    let nu = poisson_ratio;
    let e  = youngs_modulus;
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu     = e / (2.0 * (1.0 + nu));
    (lambda, mu)
}

/// Compute the Green-Lagrange strain tensor E = 0.5 * (FtF - I) from F.
fn green_lagrange(f: &Matrix3<f64>) -> Matrix3<f64> {
    0.5 * (f.transpose() * f - Matrix3::identity())
}

/// Build the isotropic tangent stiffness matrix in Voigt notation (6x6)
/// from Lame parameters lambda and mu.
/// This is C = lambda*(I tensor I) + 2*mu*I4, identical to the linear elastic C.
/// For SVK, this matrix is constant (independent of F).
fn hooke_voigt(lambda: f64, mu: f64) -> Matrix6<f64> {
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

/// Build the Neo-Hookean tangent stiffness in Voigt notation (6x6).
///
/// C_tangent = lambda * (C^{-1} tensor C^{-1}) + (mu - lambda*ln(J)) * (C^{-1} odot C^{-1})
///
/// Voigt ordering: [11, 22, 33, 12, 23, 13]
/// Index map: 0=>(0,0), 1=>(1,1), 2=>(2,2), 3=>(0,1), 4=>(1,2), 5=>(0,2)
///
/// (A tensor B)_IJ = A_I * B_J  (outer product in Voigt)
/// (A odot B)_IJKL = 0.5*(A_ik*B_jl + A_il*B_jk) — symmetrized product
///
/// Direct mapping of the 4th-order tensor C_ijkl to Voigt notation.
///
/// Voigt ordering: [11,22,33,12,23,13] -> indices [0,1,2,3,4,5]
/// Index map: 0->(0,0), 1->(1,1), 2->(2,2), 3->(0,1), 4->(1,2), 5->(0,2)
///
/// Convention used in hooke_voigt (and throughout the assembler):
///   C_IJ maps directly to C_ijkl with no additional Voigt factors.
///   dS_I = sum_J C_IJ * dE_J  where dE uses engineering shear (dE_3 = 2*de_12).
///
/// Verification at F=I (C^{-1}=I, coeff=mu):
///   C_0000 = lambda*1 + mu*(1+1) = lambda + 2*mu  (= hooke[0,0]) OK
///   C_0011 = lambda*1*1 + mu*(0+0) = lambda        (= hooke[0,1]) OK
///   C_0101 = 0 + mu*(1*1 + 0*0) = mu              (= hooke[3,3]) OK
fn nh_tangent_voigt(c_inv: &Matrix3<f64>, lambda: f64, coeff: f64) -> Matrix6<f64> {
    let idx = [(0usize,0usize),(1,1),(2,2),(0,1),(1,2),(0,2)];

    let mut c_mat = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            let (a, b) = idx[i];
            let (k, l) = idx[j];

            // C_abkl = lambda*C^{-1}_ab*C^{-1}_kl
            //        + coeff*(C^{-1}_ak*C^{-1}_bl + C^{-1}_al*C^{-1}_bk)
            c_mat[(i, j)] = lambda * c_inv[(a, b)] * c_inv[(k, l)]
                + coeff * (
                    c_inv[(a, k)] * c_inv[(b, l)]
                  + c_inv[(a, l)] * c_inv[(b, k)]
                );
        }
    }
    c_mat
}

// ─── Trait ────────────────────────────────────────────────────────────────────

/// Defines a hyperelastic material law in the Lagrangian frame.
///
/// All methods take the deformation gradient F (3x3) as input.
/// Each implementation derives what it needs internally (E, C, J, ...).
///
/// Convention:
/// - pk2_stress returns S (2nd Piola-Kirchhoff, symmetric, reference frame)
/// - The assembler computes P = F*S (1st Piola-Kirchhoff) internally for assembly
/// - tangent_stiffness returns dS/dE in Voigt notation (6x6)
pub trait MaterialLaw {
    /// Mass density (kg/m3) of the material in the reference configuration.
    fn density(&self) -> f64;

    /// Strain energy density W(F).
    fn strain_energy(&self, f: &Matrix3<f64>) -> f64;

    /// 2nd Piola-Kirchhoff stress tensor S = dW/dE (3x3, symmetric).
    fn pk2_stress(&self, f: &Matrix3<f64>) -> Matrix3<f64>;

    /// Material tangent stiffness C = dS/dE in Voigt notation (6x6).
    /// For linear-like materials (LinearElastic, SVK) this is constant.
    /// For Neo-Hookean and beyond this depends on F.
    fn tangent_stiffness(&self, f: &Matrix3<f64>) -> Matrix6<f64>;
}

// ─── Saint Venant-Kirchhoff ───────────────────────────────────────────────────

/// Saint Venant-Kirchhoff hyperelastic material.
///
/// Extends linear elasticity to large deformations by applying the same
/// Hooke law to the Green-Lagrange strain tensor E instead of the small
/// strain tensor epsilon:
///
///   S = lambda*tr(E)*I + 2*mu*E
///
/// where E = 0.5*(FtF - I).
///
/// Properties:
/// - Reduces exactly to linear elasticity for small deformations (F ~ I)
/// - Tangent stiffness C = dS/dE is constant (same as linear elastic C)
/// - Not suitable for very large compressive deformations (S can become non-physical)
///
/// Replaces LinearElastic from v0.3 — all previous tests remain numerically valid
/// since SVK == linear elasticity for small deformations.
pub struct SaintVenantKirchhoff {
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    pub density:        f64,
}

impl MaterialLaw for SaintVenantKirchhoff {

    fn density(&self) -> f64 {
        self.density
    }

    /// W = lambda/2 * tr(E)^2 + mu * E:E
    fn strain_energy(&self, f: &Matrix3<f64>) -> f64 {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let e = green_lagrange(f);
        let trace_e = e.trace();
        lambda / 2.0 * trace_e * trace_e + mu * e.norm_squared()
    }

    /// S = lambda*tr(E)*I + 2*mu*E
    fn pk2_stress(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let e = green_lagrange(f);
        lambda * e.trace() * Matrix3::identity() + 2.0 * mu * e
    }

    /// C = lambda*(I tensor I) + 2*mu*I4 — constant, independent of F.
    fn tangent_stiffness(&self, _f: &Matrix3<f64>) -> Matrix6<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        hooke_voigt(lambda, mu)
    }
}

// ─── Neo-Hookean ──────────────────────────────────────────────────────────────

/// Compressible Neo-Hookean hyperelastic material.
///
/// Strain energy density:
///   W = mu/2 * (I1 - 3) - mu*ln(J) + lambda/2 * (ln J)^2
///
/// where:
///   I1 = tr(C) = tr(F^T F)
///   J  = det(F)
///   C  = F^T F  (right Cauchy-Green tensor)
///
/// 2nd Piola-Kirchhoff stress:
///   S = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1}
///
/// Material tangent stiffness (dS/dE in Voigt, 6x6):
///   C_tangent = lambda * (C^{-1} tensor C^{-1})
///             + (mu - lambda*ln(J)) * (C^{-1} odot C^{-1})
///
/// where odot is the symmetrized tensor product:
///   (A odot B)_ijkl = 0.5*(A_ik B_jl + A_il B_jk)
///
/// Properties:
/// - Physically consistent for large deformations and large compressions
/// - S -> 0 and C_tangent -> Hooke as F -> I (consistent with SVK/linear)
/// - Requires J > 0 — assembler skips degenerate elements (volume < 1e-10)
/// - C_tangent depends on F -> K must be reassembled at each Newton iteration
pub struct NeoHookean {
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    pub density:        f64,
}

impl MaterialLaw for NeoHookean {

    fn density(&self) -> f64 {
        self.density
    }

    /// W = mu/2*(I1-3) - mu*ln(J) + lambda/2*(ln J)^2
    fn strain_energy(&self, f: &Matrix3<f64>) -> f64 {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        let ln_j = j.ln();
        let c = f.transpose() * f;
        let i1 = c.trace();
        mu / 2.0 * (i1 - 3.0) - mu * ln_j + lambda / 2.0 * ln_j * ln_j
    }

    /// S = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1}
    fn pk2_stress(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }
        let ln_j = j.ln();
        let c = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None => return Matrix3::zeros(),
        };
        mu * (Matrix3::identity() - c_inv) + lambda * ln_j * c_inv
    }

    /// C_tangent = lambda*(C^{-1} tensor C^{-1}) + (mu - lambda*ln(J))*(C^{-1} odot C^{-1})
    fn tangent_stiffness(&self, f: &Matrix3<f64>) -> Matrix6<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let j = f.determinant();
        if j <= 0.0 { return hooke_voigt(lambda, mu); }
        let ln_j = j.ln();
        let c = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None => return hooke_voigt(lambda, mu),
        };
        // coeff of the odot term: mu - lambda*ln(J)
        let coeff = mu - lambda * ln_j;
        nh_tangent_voigt(&c_inv, lambda, coeff)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn svk_steel() -> SaintVenantKirchhoff {
        SaintVenantKirchhoff { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1000.0 }
    }

    fn nh_steel() -> NeoHookean {
        NeoHookean { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1000.0 }
    }

    // ─── SVK tests (inchanges) ────────────────────────────────────────────────

    #[test]
    fn test_tangent_stiffness_at_identity() {
        let mat = svk_steel();
        let c = mat.tangent_stiffness(&Matrix3::identity());
        let expected = 1000.0 * 0.7 / (1.3 * 0.4);
        assert!((c[(0, 0)] - expected).abs() < 1e-3, "C[0,0] = {}, expected {}", c[(0,0)], expected);
    }

    #[test]
    fn test_tangent_stiffness_is_constant() {
        let mat = svk_steel();
        let f1 = Matrix3::identity();
        let f2 = Matrix3::new(1.1, 0.05, 0.0, 0.0, 0.95, 0.0, 0.0, 0.0, 1.0);
        let c1 = mat.tangent_stiffness(&f1);
        let c2 = mat.tangent_stiffness(&f2);
        assert!((c1 - c2).norm() < 1e-10, "tangent_stiffness must be constant for SVK");
    }

    #[test]
    fn test_pk2_stress_zero_at_identity() {
        let mat = svk_steel();
        let s = mat.pk2_stress(&Matrix3::identity());
        assert!(s.norm() < 1e-10, "S at F=I must be zero, got norm {}", s.norm());
    }

    #[test]
    fn test_strain_energy_zero_at_identity() {
        let mat = svk_steel();
        let w = mat.strain_energy(&Matrix3::identity());
        assert!(w.abs() < 1e-10, "W at F=I must be zero, got {}", w);
    }

    #[test]
    fn test_strain_energy_non_negative() {
        let mat = svk_steel();
        let f = Matrix3::new(1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.0, 0.0, 1.0);
        let w = mat.strain_energy(&f);
        assert!(w >= 0.0, "Strain energy must be non-negative, got {}", w);
    }

    #[test]
    fn test_svk_converges_to_linear_for_small_deformations() {
        let mat = svk_steel();
        let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
        let eps_scale = 1e-4;
        let grad_u = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
        let f = Matrix3::identity() + grad_u;
        let s_svk = mat.pk2_stress(&f);
        let eps_lin = 0.5 * (grad_u + grad_u.transpose());
        let s_lin = lambda * eps_lin.trace() * Matrix3::identity() + 2.0 * mu * eps_lin;
        let error = (s_svk - s_lin).norm() / s_lin.norm();
        assert!(error < 1e-3, "SVK should converge to linear for small deformations, error = {:.2e}", error);
    }

    // ─── Neo-Hookean tests ────────────────────────────────────────────────────

    /// S at F=I must be zero.
    /// S = mu*(I - C^{-1}) + lambda*ln(J)*C^{-1}
    /// At F=I : C=I, C^{-1}=I, J=1, ln(J)=0 -> S = mu*(I-I) + 0 = 0.
    #[test]
    fn test_nh_pk2_zero_at_identity() {
        let mat = nh_steel();
        let s = mat.pk2_stress(&Matrix3::identity());
        assert!(s.norm() < 1e-10, "NH S at F=I must be zero, norm = {:.2e}", s.norm());
    }

    /// W at F=I must be zero.
    /// I1 = tr(I) = 3, J = 1, ln(J) = 0 -> W = mu/2*(3-3) - 0 + 0 = 0.
    #[test]
    fn test_nh_strain_energy_zero_at_identity() {
        let mat = nh_steel();
        let w = mat.strain_energy(&Matrix3::identity());
        assert!(w.abs() < 1e-10, "NH W at F=I must be zero, got {:.2e}", w);
    }

    /// W must be non-negative for a physically admissible F (det > 0).
    #[test]
    fn test_nh_strain_energy_non_negative() {
        let mat = nh_steel();
        let f = Matrix3::new(1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.0, 0.0, 1.0);
        assert!(f.determinant() > 0.0);
        let w = mat.strain_energy(&f);
        assert!(w >= 0.0, "NH strain energy must be non-negative, got {}", w);
    }

    /// C_tangent at F=I must equal the Hooke matrix (same as SVK at rest).
    /// At F=I : C^{-1}=I, ln(J)=0, coeff = mu - 0 = mu.
    /// C_tangent = lambda*(I tensor I) + mu*(I odot I) = hooke_voigt(lambda, mu).
    #[test]
    fn test_nh_tangent_at_identity_matches_hooke() {
        let mat = nh_steel();
        let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
        let c_nh = mat.tangent_stiffness(&Matrix3::identity());
        let c_hooke = hooke_voigt(lambda, mu);
        let diff = (c_nh - c_hooke).norm();
        assert!(diff < 1e-8, "NH C_tangent at F=I must match Hooke, diff = {:.2e}", diff);
    }

    /// C_tangent must be symmetric for any F.
    #[test]
    fn test_nh_tangent_symmetric() {
        let mat = nh_steel();
        let f = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);
        let c = mat.tangent_stiffness(&f);
        let diff = (&c - c.transpose()).abs().max();
        assert!(diff < 1e-10, "NH C_tangent must be symmetric, diff = {:.2e}", diff);
    }

    /// For small deformations, NH pk2_stress must converge to linear elastic stress.
    /// Both SVK and NH reduce to linear elasticity for F ~ I.
    #[test]
    fn test_nh_converges_to_linear_for_small_deformations() {
        let mat = nh_steel();
        let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
        let eps_scale = 1e-4;
        let grad_u = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
        let f = Matrix3::identity() + grad_u;
        let s_nh = mat.pk2_stress(&f);
        let eps_lin = 0.5 * (grad_u + grad_u.transpose());
        let s_lin = lambda * eps_lin.trace() * Matrix3::identity() + 2.0 * mu * eps_lin;
        let error = (s_nh - s_lin).norm() / s_lin.norm();
        assert!(error < 1e-3, "NH should converge to linear for small deformations, error = {:.2e}", error);
    }

    /// NH and SVK must give close results for small deformations.
    #[test]
    fn test_nh_matches_svk_small_strain() {
        let nh  = nh_steel();
        let svk = svk_steel();
        let eps_scale = 1e-5;
        let grad_u = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
        let f = Matrix3::identity() + grad_u;
        let s_nh  = nh.pk2_stress(&f);
        let s_svk = svk.pk2_stress(&f);
        let error = (s_nh - s_svk).norm() / s_svk.norm().max(1e-14);
        assert!(error < 1e-3, "NH and SVK must agree for small strains, error = {:.2e}", error);
    }

    /// Numerical tangent check via symmetric perturbation of E.
    ///
    /// Strategy: perturb E directly by adding h*e_I (a unit Voigt direction),
    /// then recover a valid F by computing F = chol(I + 2*(E0 + h*e_I))^T.
    /// This ensures dE is exactly h*e_I and the finite difference is clean.
    ///
    /// We verify: C_analytical[:,I] ~ (S(E+h*e_I) - S(E-h*e_I)) / (2*h)
    /// for each Voigt direction I. Tolerance 1e-4 is appropriate for h=1e-7.
    #[test]
    fn test_nh_tangent_numerical_consistency() {
        use nalgebra::linalg::Cholesky;

        let mat = nh_steel();
        let h = 1e-7;

        // Voigt index -> (row, col)
        let idx = [(0usize,0usize),(1,1),(2,2),(0,1),(1,2),(0,2)];

        // Base E — moderate deformation, SPD
        let e0 = Matrix3::new(
            0.05, 0.02, 0.01,
            0.02, 0.04, 0.01,
            0.01, 0.01, 0.03,
        );
        // F0 from E0: C = I + 2*E, F = chol(C)^T
        let c0: Matrix3<f64> = Matrix3::identity() + 2.0 * e0;
        let f0: Matrix3<f64> = Cholesky::new(c0).unwrap().l().transpose();

        let c_analytical = mat.tangent_stiffness(&f0);

        for i in 0..6 {
            let (r, c) = idx[i];
            let is_shear = r != c;

            // dE_mat such that the Voigt component dE_I = h exactly.
            // Voigt E_I = E_rr for normal, E_I = 2*E_rc for shear (engineering).
            // So for shear: dE_rc = dE_cr = h/2 -> dE_I_voigt = 2*(h/2) = h.
            let amp = if is_shear { h / 2.0 } else { h };
            let mut de_mat = Matrix3::zeros();
            de_mat[(r, c)] = amp;
            if is_shear { de_mat[(c, r)] = amp; }

            // F_plus from E0 + dE, F_minus from E0 - dE
            let c_plus:  Matrix3<f64> = Matrix3::identity() + 2.0 * (e0 + de_mat);
            let c_minus: Matrix3<f64> = Matrix3::identity() + 2.0 * (e0 - de_mat);
            let f_plus:  Matrix3<f64> = Cholesky::new(c_plus).unwrap().l().transpose();
            let f_minus: Matrix3<f64> = Cholesky::new(c_minus).unwrap().l().transpose();

            let s_plus  = mat.pk2_stress(&f_plus);
            let s_minus = mat.pk2_stress(&f_minus);

            // dS/dE_I via central difference — Voigt
            let ds_fd: Vec<f64> = idx.iter().map(|&(ri, ci)| {
                (s_plus[(ri, ci)] - s_minus[(ri, ci)]) / (2.0 * h)
            }).collect();

            // Analytical column i of C_tangent
            let ds_analytical: Vec<f64> = (0..6).map(|j| c_analytical[(j, i)]).collect();

            let norm_fd: f64  = ds_fd.iter().map(|x| x*x).sum::<f64>().sqrt();
            let norm_err: f64 = ds_fd.iter().zip(ds_analytical.iter())
                .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();

            if norm_fd > 1e-10 {
                let rel = norm_err / norm_fd;
                assert!(rel < 1e-4,
                    "NH tangent column {} : rel error = {:.4e}", i, rel);
            }
        }
    }
}