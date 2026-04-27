// UTF-8
// orfas-tissues/src/metadata.rs — tissue metadata, confidence intervals, and TissuePreset trait.

use std::collections::HashMap;
use serde::{Deserialize, Serialize};
use orfas_core::MaterialLaw;

// ─── ConfidenceInterval ───────────────────────────────────────────────────────

/// Confidence interval for a single rheological parameter.
///
/// Both bounds are expressed in the same unit as the parameter itself.
/// Source: literature range reported in the reference paper (typically
/// mean +/- standard deviation, or min/max across subjects).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfidenceInterval {
    /// Lower bound of the confidence interval.
    pub min: f64,
    /// Upper bound of the confidence interval.
    pub max: f64,
}

impl ConfidenceInterval {
    /// Create a new confidence interval.
    pub fn new(min: f64, max: f64) -> Self {
        debug_assert!(min <= max, "ConfidenceInterval: min must be <= max");
        ConfidenceInterval { min, max }
    }

    /// Returns true if the given value falls within [min, max].
    pub fn contains(&self, value: f64) -> bool {
        value >= self.min && value <= self.max
    }
}

// ─── TissueMetadata ───────────────────────────────────────────────────────────

/// Scientific metadata attached to a tissue material preset.
///
/// All string fields use `&'static str` — they are compile-time constants
/// embedded in the binary, with zero runtime overhead.
///
/// `confidence_intervals` maps parameter names (e.g. "mu", "k1", "kappa")
/// to their literature confidence intervals. Keys match the field names
/// of the underlying material struct for clarity.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TissueMetadata {
    /// Human-readable tissue name (e.g. "Liver", "Brain white matter").
    pub name: &'static str,

    /// Material model used (e.g. "Neo-Hookean", "Mooney-Rivlin", "Holzapfel-Ogden").
    pub model: &'static str,

    /// Full citation: "Author et al. (year), Journal, DOI".
    pub reference: &'static str,

    /// DOI of the source paper, if available.
    pub doi: Option<&'static str>,

    /// Experimental protocol used to identify the parameters
    /// (e.g. "ex vivo indentation", "in vivo MRE", "uniaxial tensile test").
    pub protocol: &'static str,

    /// Confidence intervals per rheological parameter.
    /// Key = parameter name as used in the material struct (e.g. "mu", "c1", "k1").
    /// Value = literature range in the same unit as the parameter.
    pub confidence_intervals: HashMap<&'static str, ConfidenceInterval>,

    /// Optional free-text notes (species, tissue condition, age range, etc.).
    pub notes: &'static str,
}

impl TissueMetadata {
    /// Returns true if all nominal parameter values fall within their
    /// reported confidence intervals.
    ///
    /// `values` maps parameter names to their nominal values —
    /// the same keys as `confidence_intervals`.
    pub fn all_within_confidence(
        &self,
        values: &HashMap<&'static str, f64>,
    ) -> bool {
        self.confidence_intervals.iter().all(|(key, ci)| {
            values.get(key).map_or(true, |&v| ci.contains(v))
        })
    }
}

// ─── TissuePreset trait ───────────────────────────────────────────────────────

/// A tissue material preset combining a calibrated MaterialLaw with
/// its scientific metadata.
///
/// Implementors provide:
/// - `metadata()` — bibliographic source, confidence intervals, protocol
/// - `material()` — the ready-to-use MaterialLaw instance
///
/// The material is returned as `Box<dyn MaterialLaw>` to allow uniform
/// handling in the viewer dropdown and JSON loader without knowing the
/// concrete type.
///
/// Example usage:
/// ```rust,ignore
/// let preset = LiverNeoHookean::default();
/// let meta   = preset.metadata();
/// println!("Tissue: {} — source: {}", meta.name, meta.reference);
/// let mat = preset.material(); // Box<dyn MaterialLaw>
/// ```
pub trait TissuePreset {
    /// Returns a reference to the scientific metadata for this preset.
    fn metadata(&self) -> &TissueMetadata;

    /// Returns the calibrated material law as a boxed trait object.
    fn material(&self) -> Box<dyn MaterialLaw>;
}