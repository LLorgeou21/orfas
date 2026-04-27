// UTF-8
// material/anisotropic/no_anisotropy.rs — null anisotropic part for isotropic viscoelastic materials.

use nalgebra::{Matrix3, Matrix6};
use crate::material::traits::AnisotropicPart;
use crate::material::context::MaterialContext;

// ─── NoAnisotropy ─────────────────────────────────────────────────────────────

/// Null implementation of `AnisotropicPart` — returns zero for all contributions.
///
/// Used with `ViscoelasticMaterial<I, NoAnisotropy, V>` when no fiber reinforcement
/// is needed. The compiler inlines all calls and eliminates the zero additions,
/// so there is no runtime cost compared to a purely isotropic material.
#[derive(serde::Serialize, serde::Deserialize)]
pub struct NoAnisotropy;

impl AnisotropicPart for NoAnisotropy {
    /// W_aniso = 0 — no fiber contribution.
    fn strain_energy_aniso(&self, _f: &Matrix3<f64>, _ctx: &MaterialContext) -> f64 {
        0.0
    }

    /// S_aniso = 0 — no fiber stress contribution.
    fn pk2_stress_aniso(&self, _f: &Matrix3<f64>, _ctx: &mut MaterialContext) -> Matrix3<f64> {
        Matrix3::zeros()
    }

    /// C_aniso = 0 — no fiber tangent contribution.
    fn tangent_aniso(&self, _f: &Matrix3<f64>, _ctx: &MaterialContext) -> Matrix6<f64> {
        Matrix6::zeros()
    }
}