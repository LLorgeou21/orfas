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

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn svk_steel() -> SaintVenantKirchhoff {
        SaintVenantKirchhoff { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1000.0 }
    }

    /// tangent_stiffness at F=I must match the former LinearElastic C matrix.
    /// C[0,0] = lambda + 2*mu = E*(1-nu)/((1+nu)*(1-2*nu))
    #[test]
    fn test_tangent_stiffness_at_identity() {
        let mat = svk_steel();
        let c = mat.tangent_stiffness(&Matrix3::identity());
        let expected = 1000.0 * 0.7 / (1.3 * 0.4); // E*(1-nu)/((1+nu)*(1-2nu))
        assert!((c[(0, 0)] - expected).abs() < 1e-3, "C[0,0] = {}, expected {}", c[(0,0)], expected);
    }

    /// tangent_stiffness must be independent of F (constant for SVK).
    #[test]
    fn test_tangent_stiffness_is_constant() {
        let mat = svk_steel();
        let f1 = Matrix3::identity();
        let f2 = Matrix3::new(
            1.1, 0.05, 0.0,
            0.0, 0.95, 0.0,
            0.0, 0.0,  1.0,
        );
        let c1 = mat.tangent_stiffness(&f1);
        let c2 = mat.tangent_stiffness(&f2);
        assert!((c1 - c2).norm() < 1e-10, "tangent_stiffness must be constant for SVK");
    }

    /// pk2_stress at F=I must be zero (no deformation -> no stress).
    #[test]
    fn test_pk2_stress_zero_at_identity() {
        let mat = svk_steel();
        let s = mat.pk2_stress(&Matrix3::identity());
        assert!(s.norm() < 1e-10, "S at F=I must be zero, got norm {}", s.norm());
    }

    /// strain_energy at F=I must be zero.
    #[test]
    fn test_strain_energy_zero_at_identity() {
        let mat = svk_steel();
        let w = mat.strain_energy(&Matrix3::identity());
        assert!(w.abs() < 1e-10, "W at F=I must be zero, got {}", w);
    }

    /// strain_energy must be non-negative for any F with det(F) > 0.
    #[test]
    fn test_strain_energy_non_negative() {
        let mat = svk_steel();
        let f = Matrix3::new(
            1.2, 0.1, 0.0,
            0.0, 0.9, 0.05,
            0.0, 0.0, 1.0,
        );
        let w = mat.strain_energy(&f);
        assert!(w >= 0.0, "Strain energy must be non-negative, got {}", w);
    }

    /// For small deformations (F = I + eps*grad_u), SVK pk2_stress must
    /// converge to the linear elastic stress as eps -> 0.
    /// We check that ||S_svk - C:eps|| / ||C:eps|| < tolerance.
    #[test]
    fn test_svk_converges_to_linear_for_small_deformations() {
        let mat = svk_steel();
        let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);

        // Small displacement gradient
        let eps_scale = 1e-4;
        let grad_u = Matrix3::new(
            0.3, 0.1, 0.05,
            0.1, 0.2, 0.0,
            0.05, 0.0, 0.1,
        ) * eps_scale;

        let f = Matrix3::identity() + grad_u;

        // SVK stress
        let s_svk = mat.pk2_stress(&f);

        // Linear elastic stress: sigma = lambda*tr(eps)*I + 2*mu*eps
        // where eps = 0.5*(grad_u + grad_u^T) (small strain)
        let eps_lin = 0.5 * (grad_u + grad_u.transpose());
        let s_lin = lambda * eps_lin.trace() * Matrix3::identity() + 2.0 * mu * eps_lin;

        let error = (s_svk - s_lin).norm() / s_lin.norm();
        assert!(error < 1e-3, "SVK should converge to linear for small deformations, error = {:.2e}", error);
    }
}