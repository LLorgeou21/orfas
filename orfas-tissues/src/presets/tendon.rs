// UTF-8
// orfas-tissues/src/presets/tendon.rs — tendon tissue preset.
//
// Reference: Weiss et al. (1996), "Finite element implementation of
// incompressible, transversely isotropic hyperelasticity",
// Computer Methods in Applied Mechanics and Engineering.
// DOI: 10.1016/0045-7825(95)00931-0
//
// Model: Neo-Hookean isochoric + VolumetricLnJ.
// Tendons are highly anisotropic and nearly incompressible. The Neo-Hookean
// model captures the ground matrix response. For full tendon anisotropy,
// use HolzapfelOgden with fiber directions aligned to the tendon axis.
// Parameters from tensile tests on human patellar tendon.

use std::collections::HashMap;
use orfas_core::{
    MaterialLaw, CompressibleMaterial,
    NeoHookeanIso, VolumetricLnJ,
};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

/// Tendon tissue preset — Neo-Hookean model (ground matrix).
///
/// Nominal values:
///   mu      = 1.2e6 Pa    (shear modulus of ground matrix)
///   kappa   = 1.0e8 Pa    (bulk modulus, near-incompressible)
///   density = 1200 kg/m^3
///
/// Note: this captures only the isotropic ground matrix response.
/// For fiber-reinforced tendon, use a HolzapfelOgden preset with
/// fiber directions aligned along the tendon axis.
pub struct TendonGroundMatrix {
    metadata: TissueMetadata,
}

impl Default for TendonGroundMatrix {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(0.5e6, 2.5e6));
        ci.insert("kappa", ConfidenceInterval::new(5.0e7, 2.0e8));

        TendonGroundMatrix {
            metadata: TissueMetadata {
                name:                 "Tendon (ground matrix)",
                model:                "Neo-Hookean",
                reference:            "Weiss et al. (1996), Comput. Methods Appl. Mech. Eng., DOI:10.1016/0045-7825(95)00931-0",
                doi:                  Some("10.1016/0045-7825(95)00931-0"),
                protocol:             "uniaxial tensile test, human patellar tendon",
                confidence_intervals: ci,
                notes:                "Ground matrix only — isotropic component. Collagen fiber anisotropy not included.",
            },
        }
    }
}

impl TissuePreset for TendonGroundMatrix {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 1.2e6 },
            vol:     VolumetricLnJ { kappa: 1.0e8 },
            density: 1200.0,
        })
    }
}