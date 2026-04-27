// UTF-8
// orfas-tissues/src/presets/liver.rs — liver tissue preset.
//
// Reference: Nava et al. (2008), "In vivo characterization of the mechanical
// response of the human liver", Medical Image Analysis.
// DOI: 10.1016/j.media.2007.09.001
//
// Model: Neo-Hookean isochoric + logarithmic volumetric (VolumetricLnJ).
// Parameters identified from ex vivo indentation tests on porcine liver.
// Density from literature average for soft hepatic tissue.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    NeoHookeanIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Liver tissue preset — Neo-Hookean model.
///
/// Nominal values:
///   mu      = 2100 Pa   (shear modulus)
///   kappa   = 50000 Pa  (bulk modulus)
///   density = 1060 kg/m^3
pub struct LiverNeoHookean {
    metadata: TissueMetadata,
}

impl Default for LiverNeoHookean {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(1500.0, 3000.0));
        ci.insert("kappa", ConfidenceInterval::new(30000.0, 80000.0));

        LiverNeoHookean {
            metadata: TissueMetadata {
                name:                 "Liver",
                model:                "Neo-Hookean",
                reference:            "Nava et al. (2008), Med. Image Anal., DOI:10.1016/j.media.2007.09.001",
                doi:                  Some("10.1016/j.media.2007.09.001"),
                protocol:             "ex vivo indentation, porcine liver",
                confidence_intervals: ci,
                notes:                "Room temperature, fresh tissue, porcine model",
            },
        }
    }
}

impl TissuePreset for LiverNeoHookean {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 2100.0 },
            vol:     VolumetricLnJ { kappa: 50000.0 },
            density: 1060.0,
        })
    }
}