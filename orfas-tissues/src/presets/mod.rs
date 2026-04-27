// UTF-8
// orfas-tissues/src/presets/mod.rs — re-exports all tissue presets.

mod liver;
mod brain;
mod cardiac;
mod arterial;
mod tendon;
mod ligament;
mod skin;
mod kidney;
mod prostate;

pub use liver::LiverNeoHookean;
pub use brain::{BrainGreyMatter, BrainWhiteMatter};
pub use cardiac::CardiacMyocardium;
pub use arterial::ArterialWallMedia;
pub use tendon::TendonGroundMatrix;
pub use ligament::LigamentMCL;
pub use skin::SkinMooneyRivlin;
pub use kidney::KidneyNeoHookean;
pub use prostate::ProstateNeoHookean;

use crate::metadata::TissuePreset;

/// Returns all built-in tissue presets as boxed trait objects.
///
/// Useful for iterating over all available presets in the viewer
/// dropdown or running batch consistency checks.
///
/// Order matches the declaration order above — do not rely on it
/// for anything other than display.
pub fn all_presets() -> Vec<Box<dyn TissuePreset>> {
    vec![
        Box::new(LiverNeoHookean::default()),
        Box::new(BrainGreyMatter::default()),
        Box::new(BrainWhiteMatter::default()),
        Box::new(CardiacMyocardium::default()),
        Box::new(ArterialWallMedia::default()),
        Box::new(TendonGroundMatrix::default()),
        Box::new(LigamentMCL::default()),
        Box::new(SkinMooneyRivlin::default()),
        Box::new(KidneyNeoHookean::default()),
        Box::new(ProstateNeoHookean::default()),
    ]
}