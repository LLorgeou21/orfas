// UTF-8
// orfas-tissues/src/presets/arterial.rs — arterial wall tissue preset.
//
// Reference: Holzapfel et al. (2000), "A new constitutive framework for
// arterial wall mechanics and a comparative study with other models",
// Journal of Elasticity.
// DOI: 10.1023/A:1010835316564
//
// Model: Holzapfel-Ogden anisotropic + Neo-Hookean isochoric + VolumetricLnJ.
// Parameters from uniaxial and biaxial tests on human iliac arteries.
// Two helical fiber families at +/- angle to the vessel axis are typical
// for arterial tissue — fiber directions must be set in SimulationContext.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleAnisotropicMaterial,
    NeoHookeanIso, VolumetricLnJ,
    HolzapfelOgden,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Arterial wall preset (media layer) — Holzapfel-Ogden model.
///
/// Nominal values (Holzapfel et al. 2000, Table 1 — human iliac artery, media):
///   mu      = 3000 Pa     (Neo-Hookean isochoric shear modulus)
///   k1      = 2363 Pa     (fiber stiffness)
///   k2      = 0.8393      (fiber nonlinearity, dimensionless)
///   kappa   = 300000 Pa   (bulk modulus, near-incompressible)
///   density = 1050 kg/m^3
///
/// Fiber directions (helix angle ~29 deg) must be set in SimulationContext.
pub struct ArterialWallMedia {
    metadata: TissueMetadata,
}

impl Default for ArterialWallMedia {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(1500.0,  5000.0));
        ci.insert("k1",    ConfidenceInterval::new(500.0,   6000.0));
        ci.insert("k2",    ConfidenceInterval::new(0.1,     2.5));
        ci.insert("kappa", ConfidenceInterval::new(150000.0, 500000.0));

        ArterialWallMedia {
            metadata: TissueMetadata {
                name:                 "Arterial wall (media)",
                model:                "Holzapfel-Ogden",
                reference:            "Holzapfel et al. (2000), J. Elasticity, DOI:10.1023/A:1010835316564",
                doi:                  Some("10.1023/A:1010835316564"),
                protocol:             "uniaxial and biaxial tests, human iliac artery, media layer",
                confidence_intervals: ci,
                notes:                "Media layer only. Two helical fiber families at +/-29 deg. Fiber directions via SimulationContext.",
            },
        }
    }
}

impl TissuePreset for ArterialWallMedia {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleAnisotropicMaterial {
            iso:     NeoHookeanIso  { mu: 3000.0 },
            aniso:   HolzapfelOgden { k1: 2363.0, k2: 0.8393 },
            vol:     VolumetricLnJ  { kappa: 300000.0 },
            density: 1050.0,
        })
    }
}