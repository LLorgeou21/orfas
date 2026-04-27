// UTF-8
// orfas-tissues/src/presets/brain.rs — brain tissue presets (grey and white matter).
//
// Reference: Budday et al. (2017), "Mechanical characterization of human brain tissue",
// Acta Biomaterialia.
// DOI: 10.1016/j.actbio.2017.06.024
//
// Model: Mooney-Rivlin isochoric + logarithmic volumetric (VolumetricLnJ).
// Parameters identified from quasi-static compression and tension tests
// on human post-mortem brain tissue.
// Grey and white matter have distinct mechanical properties — two presets provided.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    MooneyRivlinIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Brain grey matter preset — Mooney-Rivlin model.
///
/// Nominal values (Budday et al. 2017, Table 2):
///   c1      = 310 Pa
///   c2      = 80 Pa
///   kappa   = 25000 Pa
///   density = 1040 kg/m^3
pub struct BrainGreyMatter {
    metadata: TissueMetadata,
}

impl Default for BrainGreyMatter {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("c1",    ConfidenceInterval::new(200.0, 450.0));
        ci.insert("c2",    ConfidenceInterval::new(40.0,  130.0));
        ci.insert("kappa", ConfidenceInterval::new(15000.0, 40000.0));

        BrainGreyMatter {
            metadata: TissueMetadata {
                name:                 "Brain (grey matter)",
                model:                "Mooney-Rivlin",
                reference:            "Budday et al. (2017), Acta Biomater., DOI:10.1016/j.actbio.2017.06.024",
                doi:                  Some("10.1016/j.actbio.2017.06.024"),
                protocol:             "quasi-static compression and tension, human post-mortem",
                confidence_intervals: ci,
                notes:                "Human cerebral cortex, age 50-70, tested within 24h post-mortem",
            },
        }
    }
}

impl TissuePreset for BrainGreyMatter {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     MooneyRivlinIso { c1: 310.0, c2: 80.0 },
            vol:     VolumetricLnJ  { kappa: 25000.0 },
            density: 1040.0,
        })
    }
}

/// Brain white matter preset — Mooney-Rivlin model.
///
/// Nominal values (Budday et al. 2017, Table 2):
///   c1      = 160 Pa
///   c2      = 40 Pa
///   kappa   = 25000 Pa
///   density = 1040 kg/m^3
pub struct BrainWhiteMatter {
    metadata: TissueMetadata,
}

impl Default for BrainWhiteMatter {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("c1",    ConfidenceInterval::new(100.0, 250.0));
        ci.insert("c2",    ConfidenceInterval::new(20.0,  80.0));
        ci.insert("kappa", ConfidenceInterval::new(15000.0, 40000.0));

        BrainWhiteMatter {
            metadata: TissueMetadata {
                name:                 "Brain (white matter)",
                model:                "Mooney-Rivlin",
                reference:            "Budday et al. (2017), Acta Biomater., DOI:10.1016/j.actbio.2017.06.024",
                doi:                  Some("10.1016/j.actbio.2017.06.024"),
                protocol:             "quasi-static compression and tension, human post-mortem",
                confidence_intervals: ci,
                notes:                "Human corona radiata, age 50-70, tested within 24h post-mortem",
            },
        }
    }
}

impl TissuePreset for BrainWhiteMatter {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     MooneyRivlinIso { c1: 160.0, c2: 40.0 },
            vol:     VolumetricLnJ  { kappa: 25000.0 },
            density: 1040.0,
        })
    }
}