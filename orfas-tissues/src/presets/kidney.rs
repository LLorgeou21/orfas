// UTF-8
// orfas-tissues/src/presets/kidney.rs — kidney tissue preset.
//
// Reference: Nasseri et al. (2002), "Viscoelastic properties of pig kidney
// in shear, experimental results and modelling",
// Rheologica Acta.
// DOI: 10.1007/s00397-002-0233-2
//
// Model: Neo-Hookean isochoric + VolumetricLnJ.
// Parameters from dynamic shear tests on porcine kidney cortex.
// The elastic (quasi-static) shear modulus is used here;
// viscoelastic effects are not included in this preset.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    NeoHookeanIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Kidney tissue preset — Neo-Hookean model.
///
/// Nominal values (Nasseri et al. 2002):
///   mu      = 1800 Pa    (quasi-static shear modulus, cortex)
///   kappa   = 40000 Pa   (bulk modulus)
///   density = 1050 kg/m^3
pub struct KidneyNeoHookean {
    metadata: TissueMetadata,
}

impl Default for KidneyNeoHookean {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(1000.0, 3000.0));
        ci.insert("kappa", ConfidenceInterval::new(20000.0, 70000.0));

        KidneyNeoHookean {
            metadata: TissueMetadata {
                name:                 "Kidney",
                model:                "Neo-Hookean",
                reference:            "Nasseri et al. (2002), Rheol. Acta, DOI:10.1007/s00397-002-0233-2",
                doi:                  Some("10.1007/s00397-002-0233-2"),
                protocol:             "dynamic shear tests, porcine kidney cortex, quasi-static limit",
                confidence_intervals: ci,
                notes:                "Cortex layer, quasi-static elastic response only. Viscoelastic effects not included.",
            },
        }
    }
}

impl TissuePreset for KidneyNeoHookean {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 1800.0 },
            vol:     VolumetricLnJ { kappa: 40000.0 },
            density: 1050.0,
        })
    }
}