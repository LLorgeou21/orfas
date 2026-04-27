// UTF-8
// orfas-tissues/src/presets/ligament.rs — ligament tissue preset.
//
// Reference: Weiss et al. (1996), "Finite element implementation of
// incompressible, transversely isotropic hyperelasticity",
// Computer Methods in Applied Mechanics and Engineering.
// DOI: 10.1016/0045-7825(95)00931-0
//
// Model: Holzapfel-Ogden anisotropic + Neo-Hookean isochoric + VolumetricLnJ.
// Ligaments share structural similarities with tendons but exhibit more
// pronounced toe-region nonlinearity. Parameters from tests on human
// medial collateral ligament (MCL).

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleAnisotropicMaterial,
    NeoHookeanIso, VolumetricLnJ,
    HolzapfelOgden,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Medial collateral ligament (MCL) preset — Holzapfel-Ogden model.
///
/// Nominal values:
///   mu      = 1400 Pa     (Neo-Hookean isochoric shear modulus, ground matrix)
///   k1      = 570000 Pa   (fiber stiffness)
///   k2      = 48.0        (fiber nonlinearity, dimensionless)
///   kappa   = 1.0e8 Pa    (bulk modulus, near-incompressible)
///   density = 1200 kg/m^3
///
/// Fiber directions (along ligament axis) must be set in SimulationContext.
pub struct LigamentMCL {
    metadata: TissueMetadata,
}

impl Default for LigamentMCL {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(500.0,    3000.0));
        ci.insert("k1",    ConfidenceInterval::new(200000.0, 1200000.0));
        ci.insert("k2",    ConfidenceInterval::new(20.0,     80.0));
        ci.insert("kappa", ConfidenceInterval::new(5.0e7,    2.0e8));

        LigamentMCL {
            metadata: TissueMetadata {
                name:                 "Ligament (MCL)",
                model:                "Holzapfel-Ogden",
                reference:            "Weiss et al. (1996), Comput. Methods Appl. Mech. Eng., DOI:10.1016/0045-7825(95)00931-0",
                doi:                  Some("10.1016/0045-7825(95)00931-0"),
                protocol:             "uniaxial tensile test, human medial collateral ligament",
                confidence_intervals: ci,
                notes:                "Fiber directions along ligament axis. Fiber dirs via SimulationContext.",
            },
        }
    }
}

impl TissuePreset for LigamentMCL {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleAnisotropicMaterial {
            iso:     NeoHookeanIso  { mu: 1400.0 },
            aniso:   HolzapfelOgden { k1: 570000.0, k2: 48.0 },
            vol:     VolumetricLnJ  { kappa: 1.0e8 },
            density: 1200.0,
        })
    }
}