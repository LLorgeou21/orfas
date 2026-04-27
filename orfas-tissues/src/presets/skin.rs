// UTF-8
// orfas-tissues/src/presets/skin.rs — skin tissue preset.
//
// Reference: Groves et al. (2013), "An anisotropic, hyperelastic model for
// skin: experimental measurements, finite element modelling and identification
// of parameters for human and murine skin",
// Journal of the Mechanical Behavior of Biomedical Materials.
// DOI: 10.1016/j.jmbbm.2012.05.015
//
// Model: Mooney-Rivlin isochoric + VolumetricLnJ.
// Parameters from biaxial tensile tests on human forearm skin.
// Skin exhibits strong anisotropy due to Langer lines; this preset
// captures the isotropic average response.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    MooneyRivlinIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Skin tissue preset — Mooney-Rivlin model (isotropic average).
///
/// Nominal values (Groves et al. 2013):
///   c1      = 5000 Pa
///   c2      = 1000 Pa
///   kappa   = 150000 Pa
///   density = 1090 kg/m^3
pub struct SkinMooneyRivlin {
    metadata: TissueMetadata,
}

impl Default for SkinMooneyRivlin {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("c1",    ConfidenceInterval::new(2000.0,  12000.0));
        ci.insert("c2",    ConfidenceInterval::new(300.0,   3000.0));
        ci.insert("kappa", ConfidenceInterval::new(80000.0, 300000.0));

        SkinMooneyRivlin {
            metadata: TissueMetadata {
                name:                 "Skin",
                model:                "Mooney-Rivlin",
                reference:            "Groves et al. (2013), J. Mech. Behav. Biomed. Mater., DOI:10.1016/j.jmbbm.2012.05.015",
                doi:                  Some("10.1016/j.jmbbm.2012.05.015"),
                protocol:             "biaxial tensile test, human forearm skin in vivo",
                confidence_intervals: ci,
                notes:                "Isotropic average — Langer line anisotropy not included. Dermis layer.",
            },
        }
    }
}

impl TissuePreset for SkinMooneyRivlin {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     MooneyRivlinIso { c1: 5000.0, c2: 1000.0 },
            vol:     VolumetricLnJ  { kappa: 150000.0 },
            density: 1090.0,
        })
    }
}