// UTF-8
// orfas-tissues/src/presets/prostate.rs — prostate tissue preset.
//
// Reference: Phipps et al. (2005), "A computational model for needle insertion
// into soft tissue", Journal of Biomechanical Engineering.
// DOI: 10.1115/1.1993268
//
// Model: Neo-Hookean isochoric + VolumetricLnJ.
// Parameters from indentation and needle insertion experiments on
// ex vivo human prostate tissue.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    NeoHookeanIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Prostate tissue preset — Neo-Hookean model.
///
/// Nominal values (Phipps et al. 2005):
///   mu      = 900 Pa     (shear modulus)
///   kappa   = 20000 Pa   (bulk modulus)
///   density = 1040 kg/m^3
pub struct ProstateNeoHookean {
    metadata: TissueMetadata,
}

impl Default for ProstateNeoHookean {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(400.0, 2000.0));
        ci.insert("kappa", ConfidenceInterval::new(10000.0, 50000.0));

        ProstateNeoHookean {
            metadata: TissueMetadata {
                name:                 "Prostate",
                model:                "Neo-Hookean",
                reference:            "Phipps et al. (2005), J. Biomech. Eng., DOI:10.1115/1.1993268",
                doi:                  Some("10.1115/1.1993268"),
                protocol:             "indentation and needle insertion, ex vivo human prostate",
                confidence_intervals: ci,
                notes:                "Whole gland average. Large inter-subject variability reported.",
            },
        }
    }
}

impl TissuePreset for ProstateNeoHookean {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 900.0 },
            vol:     VolumetricLnJ { kappa: 20000.0 },
            density: 1040.0,
        })
    }
}