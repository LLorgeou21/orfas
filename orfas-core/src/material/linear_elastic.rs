// UTF-8
// material/linear_elastic.rs — Infinitesimal strain linear elastic material.
// Uses small-strain kinematics: eps = sym(grad(u)) = sym(F - I).
// Intended for patch tests and small-deformation validation only.
// For large deformations, use SaintVenantKirchhoff instead.

use nalgebra::{Matrix3, Matrix6};
use super::traits::MaterialLaw;
use crate::material::MaterialContext;
use super::helpers::{lame, hooke_voigt};

/// Infinitesimal strain linear elastic material.
/// sigma = lambda * tr(eps) * I + 2 * mu * eps
/// where eps = sym(F - I) = sym(grad(u)).
///
/// The PK2 stress is approximated as sigma (valid only for small strains).
/// The tangent stiffness is identical to SVK — constant Hooke tensor.
/// Passes the patch test exactly for linear displacement fields.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct LinearElastic {
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    pub density:        f64,
}

impl MaterialLaw for LinearElastic {

    fn density(&self) -> f64 { self.density }

    /// W = lambda/2 * tr(eps)^2 + mu * eps:eps
    fn strain_energy(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> f64 {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        // Small-strain tensor: eps = sym(F - I)
        let eps = Self::small_strain(f);
        let trace_eps = eps.trace();
        lambda / 2.0 * trace_eps * trace_eps + mu * eps.norm_squared()
    }

    /// S = sigma = lambda * tr(eps) * I + 2 * mu * eps
    fn pk2_stress(&self, f: &Matrix3<f64>, _ctx: &mut MaterialContext) -> Matrix3<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        let eps = Self::small_strain(f);
        lambda * eps.trace() * Matrix3::identity() + 2.0 * mu * eps
    }

    /// C = lambda*(I tensor I) + 2*mu*I4 — identical to SVK tangent.
    fn tangent_stiffness(&self, _f: &Matrix3<f64>, _ctx: &MaterialContext) -> Matrix6<f64> {
        let (lambda, mu) = lame(self.youngs_modulus, self.poisson_ratio);
        hooke_voigt(lambda, mu)
    }
}

impl LinearElastic {
    /// Small-strain tensor: eps = sym(F - I) = (F + F^T)/2 - I
    fn small_strain(f: &Matrix3<f64>) -> Matrix3<f64> {
        let grad_u = f - Matrix3::identity();
        (grad_u + grad_u.transpose()) * 0.5
    }
}