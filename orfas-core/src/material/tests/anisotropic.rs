// UTF-8
// material/tests/anisotropic.rs — tests for anisotropic material laws.

use crate::material::MaterialContext;
use crate::material::helpers::lame;
use crate::material::volumetric::VolumetricLnJ;
use crate::material::compressible::CompressibleAnisotropicMaterial;
use crate::material::isochoric::NeoHookeanIso;
use crate::material::anisotropic::HolzapfelOgden;
use super::helpers::{run_standard_material_tests, run_anisotropic_part_tests};

fn hgo_material() -> CompressibleAnisotropicMaterial<NeoHookeanIso, HolzapfelOgden, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    CompressibleAnisotropicMaterial {
        iso:     NeoHookeanIso  { mu },
        aniso:   HolzapfelOgden { k1: 2.3632e3, k2: 0.8393 },
        vol:     VolumetricLnJ  { kappa },
        density: 1000.0,
    }
}

#[test]
fn test_holzapfel_ogden_aniso_suite() {
    let mat = hgo_material();
    run_anisotropic_part_tests(&mat.aniso, &mut MaterialContext::default());
}

#[test]
fn test_hgo_full_standard_suite() {
    let mat        = hgo_material();
    let mu_eff     = mat.iso.mu;
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
}