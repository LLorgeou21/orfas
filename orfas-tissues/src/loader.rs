// UTF-8
// orfas-tissues/src/loader.rs — JSON serialization and deserialization for tissue presets.
//
// A preset JSON file contains:
//   - metadata fields (name, model, reference, doi, protocol, notes)
//   - "parameters" — flat map of parameter name -> f64 value
//   - "confidence_intervals" — flat map of parameter name -> { min, max }
//
// The "model" field drives dispatch: the loader constructs the correct
// MaterialLaw implementation from the parameter map.
//
// Owned types (String instead of &'static str) are used for runtime
// deserialization. They are intentionally separate from the static
// TissueMetadata used by hardcoded presets.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use serde::{Deserialize, Serialize};
use orfas_core::{
    MaterialLaw,
    NeoHookeanIso, MooneyRivlinIso, VolumetricLnJ,
    CompressibleMaterial, SaintVenantKirchhoff,
};
use orfas_core::{HolzapfelOgden, NoAnisotropy, CompressibleAnisotropicMaterial};
use crate::metadata::ConfidenceInterval;

// ─── Owned metadata (runtime / JSON) ─────────────────────────────────────────

/// Runtime-owned version of TissueMetadata for JSON deserialization.
///
/// Uses `String` instead of `&'static str` — cannot be converted to the
/// static variant. Used exclusively by the JSON loader; hardcoded presets
/// use `TissueMetadata` directly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TissueMetadataOwned {
    /// Human-readable tissue name.
    pub name: String,

    /// Material model identifier (e.g. "neo_hookean", "mooney_rivlin").
    pub model: String,

    /// Full citation string.
    pub reference: String,

    /// DOI of the source paper, if available.
    pub doi: Option<String>,

    /// Experimental protocol used to identify the parameters.
    pub protocol: String,

    /// Nominal parameter values keyed by parameter name.
    pub parameters: HashMap<String, f64>,

    /// Confidence intervals keyed by parameter name.
    #[serde(default)]
    pub confidence_intervals: HashMap<String, ConfidenceInterval>,

    /// Optional free-text notes.
    #[serde(default)]
    pub notes: String,
}

// ─── Loaded preset ────────────────────────────────────────────────────────────

/// A tissue preset loaded from a JSON file at runtime.
///
/// Holds the owned metadata and the constructed material law.
/// Unlike hardcoded presets, the material type is erased — only
/// `Box<dyn MaterialLaw>` is available.
pub struct LoadedPreset {
    /// Owned metadata deserialized from JSON.
    pub metadata: TissueMetadataOwned,
    /// Constructed material law, type-erased.
    pub material: Box<dyn MaterialLaw>,
}

// ─── Load errors ──────────────────────────────────────────────────────────────

/// Errors that can occur when loading a tissue preset from JSON.
#[derive(Debug)]
pub enum LoadError {
    /// File could not be read.
    Io(std::io::Error),
    /// JSON could not be parsed.
    Json(serde_json::Error),
    /// A required parameter is missing from the "parameters" map.
    MissingParameter(String),
    /// The "model" field names an unsupported material model.
    UnknownModel(String),
}

impl std::fmt::Display for LoadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LoadError::Io(e)               => write!(f, "IO error: {}", e),
            LoadError::Json(e)             => write!(f, "JSON parse error: {}", e),
            LoadError::MissingParameter(p) => write!(f, "Missing parameter: '{}'", p),
            LoadError::UnknownModel(m)     => write!(f, "Unknown model: '{}'", m),
        }
    }
}

impl From<std::io::Error>    for LoadError { fn from(e: std::io::Error)    -> Self { LoadError::Io(e)   } }
impl From<serde_json::Error> for LoadError { fn from(e: serde_json::Error) -> Self { LoadError::Json(e) } }

// ─── Parameter extraction helpers ────────────────────────────────────────────

/// Extract a required f64 parameter from the parameters map.
///
/// Returns `LoadError::MissingParameter` if the key is absent.
fn get(params: &HashMap<String, f64>, key: &str) -> Result<f64, LoadError> {
    params.get(key)
        .copied()
        .ok_or_else(|| LoadError::MissingParameter(key.to_string()))
}

/// Extract an optional f64 parameter, returning a default if absent.
fn get_or(params: &HashMap<String, f64>, key: &str, default: f64) -> f64 {
    params.get(key).copied().unwrap_or(default)
}

// ─── Model dispatch ───────────────────────────────────────────────────────────

/// Construct a `Box<dyn MaterialLaw>` from a model name and parameter map.
///
/// Supported models:
/// - `"neo_hookean"`          — CompressibleMaterial<NeoHookeanIso, VolumetricLnJ>
///                              parameters: mu, kappa, density
/// - `"mooney_rivlin"`        — CompressibleMaterial<MooneyRivlinIso, VolumetricLnJ>
///                              parameters: c1, c2, kappa, density
/// - `"holzapfel_ogden"`      — CompressibleAnisotropicMaterial<NeoHookeanIso, HolzapfelOgden, VolumetricLnJ>
///                              parameters: mu, k1, k2, kappa, density
/// - `"saint_venant_kirchhoff"` — SaintVenantKirchhoff
///                              parameters: youngs_modulus, poisson_ratio, density
fn build_material(
    model:  &str,
    params: &HashMap<String, f64>,
) -> Result<Box<dyn MaterialLaw>, LoadError> {
    match model {
        "neo_hookean" => {
            let mu      = get(params, "mu")?;
            let kappa   = get(params, "kappa")?;
            let density = get_or(params, "density", 1000.0);
            Ok(Box::new(CompressibleMaterial {
                iso:     NeoHookeanIso { mu },
                vol:     VolumetricLnJ { kappa },
                density,
            }))
        }

        "mooney_rivlin" => {
            let c1      = get(params, "c1")?;
            let c2      = get(params, "c2")?;
            let kappa   = get(params, "kappa")?;
            let density = get_or(params, "density", 1000.0);
            Ok(Box::new(CompressibleMaterial {
                iso:     MooneyRivlinIso { c1, c2 },
                vol:     VolumetricLnJ { kappa },
                density,
            }))
        }

        "holzapfel_ogden" => {
            let mu      = get(params, "mu")?;
            let k1      = get(params, "k1")?;
            let k2      = get(params, "k2")?;
            let kappa   = get(params, "kappa")?;
            let density = get_or(params, "density", 1000.0);
            Ok(Box::new(CompressibleAnisotropicMaterial {
                iso:     NeoHookeanIso  { mu },
                aniso:   HolzapfelOgden { k1, k2 },
                vol:     VolumetricLnJ  { kappa },
                density,
            }))
        }

        "saint_venant_kirchhoff" => {
            let youngs_modulus = get(params, "youngs_modulus")?;
            let poisson_ratio  = get(params, "poisson_ratio")?;
            let density        = get_or(params, "density", 1000.0);
            Ok(Box::new(SaintVenantKirchhoff { youngs_modulus, poisson_ratio, density }))
        }

        other => Err(LoadError::UnknownModel(other.to_string())),
    }
}

// ─── Public API ───────────────────────────────────────────────────────────────

/// Load a tissue preset from a JSON file on disk.
///
/// The JSON must follow the orfas-tissues preset schema:
/// ```json
/// {
///     "name": "Liver",
///     "model": "neo_hookean",
///     "reference": "Nava et al. (2008), J. Biomech.",
///     "doi": "10.1016/j.jbiomech.2007.09.017",
///     "protocol": "ex vivo indentation",
///     "parameters": { "mu": 2100.0, "kappa": 50000.0, "density": 1060.0 },
///     "confidence_intervals": { "mu": { "min": 1500.0, "max": 3000.0 } },
///     "notes": "Porcine liver, room temperature"
/// }
/// ```
///
/// Returns a `LoadedPreset` with the constructed material and owned metadata,
/// or a `LoadError` describing what went wrong.
pub fn load_preset_from_file(path: &Path) -> Result<LoadedPreset, LoadError> {
    let contents = fs::read_to_string(path)?;
    load_preset_from_str(&contents)
}

/// Load a tissue preset from a JSON string.
///
/// Same schema as `load_preset_from_file`. Useful for embedding JSON
/// directly in tests or loading from a network source.
pub fn load_preset_from_str(json: &str) -> Result<LoadedPreset, LoadError> {
    let metadata: TissueMetadataOwned = serde_json::from_str(json)?;
    let material = build_material(&metadata.model, &metadata.parameters)?;
    Ok(LoadedPreset { metadata, material })
}

/// Serialize a `TissueMetadataOwned` to a JSON string.
///
/// Useful for saving a modified or custom preset back to disk.
pub fn save_preset_to_str(metadata: &TissueMetadataOwned) -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(metadata)
}

/// Save a `TissueMetadataOwned` to a JSON file on disk.
pub fn save_preset_to_file(
    metadata: &TissueMetadataOwned,
    path:     &Path,
) -> Result<(), LoadError> {
    let contents = save_preset_to_str(metadata)?;
    fs::write(path, contents)?;
    Ok(())
}