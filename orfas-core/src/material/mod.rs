// UTF-8
// material/mod.rs — public API of the material module.
// Re-exports all public types and helpers. External code only needs to import from here.

mod traits;
pub mod helpers;
mod volumetric;
mod compressible;
mod svk;
pub mod isochoric;
pub mod fiber_fields;
pub mod context;
pub mod anisotropic;
pub mod internal_variables;
pub mod viscoelastic;
pub mod consistency;


#[cfg(test)]
mod tests;

// ─── Public re-exports ────────────────────────────────────────────────────────

pub use traits::{MaterialLaw, IsochoricPart, VolumetricPart};
pub use helpers::{lame, hooke_voigt};
pub use volumetric::{VolumetricLnJ, VolumetricQuad};
pub use compressible::{CompressibleMaterial,CompressibleAnisotropicMaterial};
pub use svk::SaintVenantKirchhoff;
pub use isochoric::{NeoHookeanIso, MooneyRivlinIso, OgdenIso};
pub use anisotropic::holzapfel_ogden::HolzapfelOgden;
pub use context::{MaterialContext, SimulationContext};
pub use internal_variables::{InternalVariables, ElementInternalVars};
pub use anisotropic::NoAnisotropy;
pub use viscoelastic::ViscoelasticMaterial;
pub use consistency::check_thermodynamic_consistency;