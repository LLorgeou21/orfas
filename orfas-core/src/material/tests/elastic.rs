// UTF-8
// material/tests/elastic.rs — tests for elastic material laws.

use nalgebra::Matrix3;
use crate::material::helpers::lame;
use crate::material::MaterialContext;
use crate::material::volumetric::{VolumetricLnJ, VolumetricQuad};
use crate::material::compressible::CompressibleMaterial;
use crate::material::svk::SaintVenantKirchhoff;
use crate::material::isochoric::{NeoHookeanIso, MooneyRivlinIso, OgdenIso};
use crate::material::traits::MaterialLaw;
use super::helpers::{run_standard_material_tests, run_numerical_tangent_check};

fn nh_material() -> CompressibleMaterial<NeoHookeanIso, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    CompressibleMaterial { iso: NeoHookeanIso { mu }, vol: VolumetricLnJ { kappa }, density: 1000.0 }
}

fn mr_material() -> CompressibleMaterial<MooneyRivlinIso, VolumetricLnJ> {
    let mu_mr = 2.0 * (300.0 + 100.0);
    let e_eff = 2.0 * mu_mr * 1.3;
    let (lambda, _) = lame(e_eff, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu_mr;
    CompressibleMaterial { iso: MooneyRivlinIso { c1: 300.0, c2: 100.0 }, vol: VolumetricLnJ { kappa }, density: 1000.0 }
}

fn ogden_nh_material() -> CompressibleMaterial<OgdenIso, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    CompressibleMaterial { iso: OgdenIso::new(vec![mu], vec![2.0]).unwrap(), vol: VolumetricLnJ { kappa }, density: 1000.0 }
}

fn svk_material() -> SaintVenantKirchhoff {
    SaintVenantKirchhoff { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1000.0 }
}

#[test]
fn test_svk_standard_suite() {
    let mat = svk_material();
    let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
    run_standard_material_tests(&mat, lambda, mu, &mut MaterialContext::default());
}

#[test]
fn test_neo_hookean_standard_suite() {
    let mat        = nh_material();
    let mu_eff     = mat.iso.mu;
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
}

#[test]
fn test_mooney_rivlin_standard_suite() {
    let mat        = mr_material();
    let mu_eff     = 2.0 * (mat.iso.c1 + mat.iso.c2);
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
}

#[test]
fn test_ogden_standard_suite() {
    let mat        = ogden_nh_material();
    let mu_eff     = mat.iso.mu[0];
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
}

#[test]
fn test_svk_tangent_is_constant() {
    let mat = svk_material();
    let f1  = Matrix3::identity();
    let f2  = Matrix3::new(1.1, 0.05, 0.0, 0.0, 0.95, 0.0, 0.0, 0.0, 1.0);
    let c1  = mat.tangent_stiffness(&f1, &MaterialContext::default());
    let c2  = mat.tangent_stiffness(&f2, &MaterialContext::default());
    assert!((c1 - c2).norm() < 1e-10, "SVK tangent must be constant");
}

#[test]
fn test_svk_converges_to_linear_for_small_deformations() {
    let mat          = svk_material();
    let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
    let eps_scale    = 1e-4;
    let grad_u       = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f            = Matrix3::identity() + grad_u;
    let s_svk        = mat.pk2_stress(&f, &mut MaterialContext::default());
    let eps_lin      = 0.5 * (grad_u + grad_u.transpose());
    let s_lin        = lambda * eps_lin.trace() * Matrix3::identity() + 2.0 * mu * eps_lin;
    let error        = (s_svk - s_lin).norm() / s_lin.norm();
    assert!(error < 1e-3, "SVK must converge to linear, error = {:.2e}", error);
}

#[test]
fn test_nh_matches_svk_small_strain() {
    let nh        = nh_material();
    let svk       = svk_material();
    let eps_scale = 1e-5;
    let grad_u    = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f         = Matrix3::identity() + grad_u;
    let s_nh      = nh.pk2_stress(&f, &mut MaterialContext::default());
    let s_svk     = svk.pk2_stress(&f, &mut MaterialContext::default());
    let error     = (s_nh - s_svk).norm() / s_svk.norm().max(1e-14);
    assert!(error < 1e-3, "NH and SVK must agree for small strains, error = {:.2e}", error);
}

#[test]
fn test_mr_c2_zero_matches_nh() {
    let c1    = 500.0;
    let mu    = 2.0 * c1;
    let (lambda, _) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    let mr = CompressibleMaterial { iso: MooneyRivlinIso { c1, c2: 0.0 }, vol: VolumetricLnJ { kappa }, density: 1000.0 };
    let nh = CompressibleMaterial { iso: NeoHookeanIso { mu }, vol: VolumetricLnJ { kappa }, density: 1000.0 };
    let f  = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);

    let w_mr = mr.strain_energy(&f, &MaterialContext::default());
    let w_nh = nh.strain_energy(&f, &MaterialContext::default());
    assert!((w_mr - w_nh).abs() / w_nh.abs() < 1e-8, "MR(c2=0) W must match NH");

    let s_mr = mr.pk2_stress(&f, &mut MaterialContext::default());
    let s_nh = nh.pk2_stress(&f, &mut MaterialContext::default());
    assert!((s_mr - s_nh).norm() / s_nh.norm().max(1e-14) < 1e-8, "MR(c2=0) S must match NH");
}

#[test]
fn test_ogden_matches_nh() {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    let ogden = CompressibleMaterial { iso: OgdenIso::new(vec![mu], vec![2.0]).unwrap(), vol: VolumetricLnJ { kappa }, density: 1000.0 };
    let nh    = CompressibleMaterial { iso: NeoHookeanIso { mu }, vol: VolumetricLnJ { kappa }, density: 1000.0 };
    let f     = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);

    let w_diff = (ogden.strain_energy(&f, &MaterialContext::default()) - nh.strain_energy(&f, &MaterialContext::default())).abs();
    assert!(w_diff < 1e-8, "Ogden(N=1,a=2) W must match NH");
    let s_diff = (ogden.pk2_stress(&f, &mut MaterialContext::default()) - nh.pk2_stress(&f, &mut MaterialContext::default())).norm();
    assert!(s_diff < 1e-8, "Ogden(N=1,a=2) S must match NH");
    let c_diff = (ogden.tangent_stiffness(&f, &MaterialContext::default()) - nh.tangent_stiffness(&f, &MaterialContext::default())).norm();
    assert!(c_diff < 1e-6, "Ogden(N=1,a=2) C must match NH");
}

#[test]
fn test_neo_hookean_iso_invalid_mu() {
    assert!(NeoHookeanIso::new(-1.0).is_err());
    assert!(NeoHookeanIso::new(0.0).is_err());
    assert!(NeoHookeanIso::new(1.0).is_ok());
}

#[test]
fn test_volumetric_lnj_invalid_kappa() {
    assert!(VolumetricLnJ::new(-1.0).is_err());
    assert!(VolumetricLnJ::new(0.0).is_err());
    assert!(VolumetricLnJ::new(1.0).is_ok());
}

#[test]
fn test_volumetric_quad_invalid_kappa() {
    assert!(VolumetricQuad::new(-1.0).is_err());
    assert!(VolumetricQuad::new(0.0).is_err());
    assert!(VolumetricQuad::new(1.0).is_ok());
}

#[test]
fn test_mooney_rivlin_iso_invalid_params() {
    assert!(MooneyRivlinIso::new(-1.0,  0.0).is_err());
    assert!(MooneyRivlinIso::new( 1.0, -2.0).is_err());
    assert!(MooneyRivlinIso::new( 1.0,  0.0).is_ok());
    assert!(MooneyRivlinIso::new( 1.0, -0.5).is_ok());
}

#[test]
fn test_ogden_iso_invalid_params() {
    assert!(OgdenIso::new(vec![1.0], vec![2.0, 3.0]).is_err());
    assert!(OgdenIso::new(vec![], vec![]).is_err());
    assert!(OgdenIso::new(vec![1.0], vec![-1.0]).is_err());
    assert!(OgdenIso::new(vec![1.0, 0.5], vec![2.0, 4.0]).is_ok());
}

#[test]
fn test_svk_new_invalid_params() {
    assert!(SaintVenantKirchhoff::new(-1.0,   0.3, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.6, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.0, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.3, 1000.0).is_ok());
}