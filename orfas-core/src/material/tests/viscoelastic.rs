// UTF-8
// material/tests/viscoelastic.rs — tests for viscoelastic material laws.

use crate::material::MaterialContext;
use crate::material::helpers::lame;
use crate::material::volumetric::VolumetricLnJ;
use crate::material::isochoric::NeoHookeanIso;
use crate::material::anisotropic::{HolzapfelOgden, NoAnisotropy};
use crate::material::viscoelastic::ViscoelasticMaterial;
use super::helpers::{run_standard_material_tests, run_viscoelastic_tests};

fn nh_visco_material() -> ViscoelasticMaterial<NeoHookeanIso, NoAnisotropy, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    ViscoelasticMaterial {
        iso:        NeoHookeanIso { mu },
        aniso:      NoAnisotropy,
        vol:        VolumetricLnJ { kappa },
        density:    1000.0,
        tau_iso:    vec![1.0],
        beta_iso:   vec![0.3],
        tau_aniso:  vec![],
        beta_aniso: vec![],
    }
}

fn hgo_visco_material() -> ViscoelasticMaterial<NeoHookeanIso, HolzapfelOgden, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    ViscoelasticMaterial {
        iso:        NeoHookeanIso  { mu },
        aniso:      HolzapfelOgden { k1: 2363.2, k2: 0.8393 },
        vol:        VolumetricLnJ  { kappa },
        density:    1000.0,
        tau_iso:    vec![1.0],
        beta_iso:   vec![0.3],
        tau_aniso:  vec![0.5],
        beta_aniso: vec![0.2],
    }
}

#[test]
fn test_nh_viscoelastic_suite() {
    run_viscoelastic_tests(&nh_visco_material());
}

#[test]
fn test_hgo_viscoelastic_suite() {
    run_viscoelastic_tests(&hgo_visco_material());
}

#[test]
fn test_nh_viscoelastic_standard_suite() {
    let mat        = nh_visco_material();
    let mu_eff     = mat.iso.mu;
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
}