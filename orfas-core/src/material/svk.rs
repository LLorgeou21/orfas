// UTF-8
// material/svk.rs — Saint Venant-Kirchhoff hyperelastic material.

use nalgebra::{Matrix3, Matrix6};
use super::traits::MaterialLaw;
use crate::material::MaterialContext;
use super::helpers::{lame, hooke_voigt, green_lagrange};

// ─── SaintVenantKirchhoff ─────────────────────────────────────────────────────

/// Saint Venant-Kirchhoff hyperelastic material.
///
/// Extends linear elasticity to large deformations by applying Hooke's law
/// to the Green-Lagrange strain tensor E instead of the infinitesimal strain:
///
///   W = lambda/2 * tr(E)^2 + mu * E:E
///   S = lambda*tr(E)*I + 2*mu*E    where E = 0.5*(F^T F - I)
///
/// Properties:
/// - Tangent stiffness C = dS/dE is constant (independent of F).
/// - Reduces exactly to linear elasticity for small deformations (F ~ I).
/// - Not suitable for large compressive deformations (S can become non-physical).
/// - Does NOT use the isochoric/volumetric decomposition — implemented
///   directly as `MaterialLaw` without `CompressibleMaterial<I,V>`.
///
/// Parameters:
///   youngs_modulus — Young's modulus (Pa), must be > 0
///   poisson_ratio  — Poisson's ratio (dimensionless), must be in (0, 0.5)
///   density        — mass density (kg/m^3), must be > 0
#[derive(serde::Serialize, serde::Deserialize)]
pub struct SaintVenantKirchhoff {
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    pub density:        f64,
}

impl SaintVenantKirchhoff {
    /// Safe constructor — validates E > 0, 0 < nu < 0.5, density > 0.
    pub fn new(youngs_modulus: f64, poisson_ratio: f64, density: f64) -> Result<Self, String> {
        if youngs_modulus <= 0.0 {
            return Err(format!("SVK: youngs_modulus must be > 0, got {}", youngs_modulus));
        }
        if poisson_ratio <= 0.0 || poisson_ratio >= 0.5 {
            return Err(format!(
                "SVK: poisson_ratio must be in (0, 0.5), got {}", poisson_ratio
            ));
        }
        if density <= 0.0 {
            return Err(format!("SVK: density must be > 0, got {}", density));
        }
        Ok(Self { youngs_modulus, poisson_ratio, density })
    }
}

impl MaterialLaw for SaintVenantKirchhoff {

    fn density(&self) -> f64 { self.density }

    /// W = lambda/2 * tr(E)^2 + mu * E:E
    fn strain_energy(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> f64 {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let e = green_lagrange(f);
        let trace_e = e.trace();
        lambda / 2.0 * trace_e * trace_e + mu * e.norm_squared()
    }

    /// S = lambda*tr(E)*I + 2*mu*E
    fn pk2_stress(&self, f: &Matrix3<f64>, _ctx: &mut MaterialContext) -> Matrix3<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let e = green_lagrange(f);
        lambda * e.trace() * Matrix3::identity() + 2.0 * mu * e
    }

    /// C = lambda*(I tensor I) + 2*mu*I4 — constant, independent of F.
    fn tangent_stiffness(&self, _f: &Matrix3<f64>, _ctx: &MaterialContext) -> Matrix6<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        hooke_voigt(lambda, mu)
    }
}