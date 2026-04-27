// UTF-8
// orfas-tissues/src/presets/cardiac.rs — myocardium tissue preset.
//
// Reference: Holzapfel & Ogden (2009), "Constitutive modelling of passive
// myocardium: a structurally based framework for material characterization",
// Philosophical Transactions of the Royal Society A.
// DOI: 10.1098/rsta.2009.0091
//
// Model: Holzapfel-Ogden anisotropic + Neo-Hookean isochoric + VolumetricLnJ.
// Parameters from biaxial mechanical testing of passive ventricular myocardium.
// The myocardium has a pronounced fiber architecture — two fiber families
// (fiber and sheet directions) are accounted for in the HGO model.
//
// Note: fiber directions must be provided via SimulationContext at runtime.
// This preset provides material parameters only — not fiber geometry.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleAnisotropicMaterial,
    NeoHookeanIso, VolumetricLnJ,
    HolzapfelOgden,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Passive myocardium preset — Holzapfel-Ogden model.
///
/// Nominal values (Holzapfel & Ogden 2009, Table 1):
///   mu      = 59 Pa      (Neo-Hookean isochoric shear modulus)
///   k1      = 18472 Pa   (fiber stiffness)
///   k2      = 16.026     (fiber nonlinearity, dimensionless)
///   kappa   = 350000 Pa  (bulk modulus, near-incompressible)
///   density = 1060 kg/m^3
///
/// Fiber directions must be set in SimulationContext — not encoded here.
pub struct CardiacMyocardium {
    metadata: TissueMetadata,
}

impl Default for CardiacMyocardium {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(30.0,   120.0));
        ci.insert("k1",    ConfidenceInterval::new(5000.0, 40000.0));
        ci.insert("k2",    ConfidenceInterval::new(5.0,    35.0));
        ci.insert("kappa", ConfidenceInterval::new(200000.0, 500000.0));

        CardiacMyocardium {
            metadata: TissueMetadata {
                name:                 "Myocardium (passive)",
                model:                "Holzapfel-Ogden",
                reference:            "Holzapfel & Ogden (2009), Phil. Trans. R. Soc. A, DOI:10.1098/rsta.2009.0091",
                doi:                  Some("10.1098/rsta.2009.0091"),
                protocol:             "biaxial mechanical testing, passive ventricular myocardium, porcine",
                confidence_intervals: ci,
                notes:                "Passive myocardium only. Fiber directions must be provided via SimulationContext.",
            },
        }
    }
}

impl TissuePreset for CardiacMyocardium {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleAnisotropicMaterial {
            iso:     NeoHookeanIso  { mu: 59.0 },
            aniso:   HolzapfelOgden { k1: 18472.0, k2: 16.026 },
            vol:     VolumetricLnJ  { kappa: 350000.0 },
            density: 1060.0,
        })
    }
}