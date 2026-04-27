// UTF-8
// orfas-tissues/src/lib.rs — public API of the orfas-tissues crate.
//
// Provides:
//   - Hardcoded tissue presets (10 tissues, calibrated from literature)
//   - TissuePreset trait and TissueMetadata struct
//   - JSON loader for custom presets at runtime
//   - all_presets() convenience function

pub mod metadata;
pub mod loader;
pub mod presets;

// ─── Public re-exports ────────────────────────────────────────────────────────

// Metadata and trait
pub use metadata::{TissueMetadata, ConfidenceInterval, TissuePreset};

// JSON loader
pub use loader::{
    TissueMetadataOwned,
    LoadedPreset,
    LoadError,
    load_preset_from_file,
    load_preset_from_str,
    save_preset_to_str,
    save_preset_to_file,
};

// Individual presets
pub use presets::{
    LiverNeoHookean,
    BrainGreyMatter,
    BrainWhiteMatter,
    CardiacMyocardium,
    ArterialWallMedia,
    TendonGroundMatrix,
    LigamentMCL,
    SkinMooneyRivlin,
    KidneyNeoHookean,
    ProstateNeoHookean,
};

// All presets as a collection
pub use presets::all_presets;



#[cfg(test)]
mod tests;