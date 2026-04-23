// UTF-8
// material/mod.rs — public API of the material module.
// Re-exports all public types and helpers. External code only needs to import from here.

mod traits;
mod helpers;
mod volumetric;
mod compressible;
mod svk;
pub mod isochoric;

#[cfg(test)]
mod tests;

// ─── Public re-exports ────────────────────────────────────────────────────────

pub use traits::{MaterialLaw, IsochoricPart, VolumetricPart};
pub use helpers::{lame, hooke_voigt};
pub use volumetric::{VolumetricLnJ, VolumetricQuad};
pub use compressible::CompressibleMaterial;
pub use svk::SaintVenantKirchhoff;
pub use isochoric::{NeoHookeanIso, MooneyRivlinIso, OgdenIso};