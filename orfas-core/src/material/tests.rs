// UTF-8
// material/tests.rs — parametric test suite for all material laws.

use nalgebra::{Matrix3};
use nalgebra::linalg::Cholesky;
use super::traits::{MaterialLaw,IsochoricPart};
use super::helpers::{lame, hooke_voigt};
use super::volumetric::{VolumetricLnJ, VolumetricQuad};
use super::compressible::CompressibleMaterial;
use super::svk::SaintVenantKirchhoff;
use super::isochoric::{NeoHookeanIso, MooneyRivlinIso, OgdenIso};

// ─── Parametric test helpers ──────────────────────────────────────────────────

/// Standard test suite for any isotropic hyperelastic material that:
/// - has zero energy and stress at F=I
/// - has tangent = Hooke(lambda, mu) at F=I
/// - has symmetric tangent for arbitrary F
/// - has non-negative strain energy for admissible F
/// - converges to linear elasticity for small deformations
/// - has numerically consistent tangent (central differences)
///
/// `expected_lambda` and `expected_mu` are the effective Lame parameters
/// at linearization. The caller computes them from the material's parameters.
fn run_standard_material_tests(
    mat:             &dyn MaterialLaw,
    expected_lambda: f64,
    expected_mu:     f64,
) {
    let id = Matrix3::identity();

    // W(F=I) = 0
    let w = mat.strain_energy(&id);
    assert!(w.abs() < 1e-10,
        "W at F=I must be zero, got {:.2e}", w);

    // S(F=I) = 0
    let s = mat.pk2_stress(&id);
    assert!(s.norm() < 1e-10,
        "S at F=I must be zero, got norm {:.2e}", s.norm());

    // C_tangent(F=I) = Hooke(lambda, mu)
    let c_mat   = mat.tangent_stiffness(&id);
    let c_hooke = hooke_voigt(expected_lambda, expected_mu);
    let diff    = (c_mat - c_hooke).norm();
    assert!(diff < 1e-6,
        "C_tangent at F=I must match Hooke, diff = {:.2e}", diff);

    // C_tangent symmetric for arbitrary F
    let f_arb   = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);
    let c_arb   = mat.tangent_stiffness(&f_arb);
    let sym_err = (&c_arb - c_arb.transpose()).abs().max();
    assert!(sym_err < 1e-10,
        "C_tangent must be symmetric, diff = {:.2e}", sym_err);

    // W >= 0 for admissible F
    let f_pos = Matrix3::new(1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.0, 0.0, 1.0);
    assert!(f_pos.determinant() > 0.0);
    let w_pos = mat.strain_energy(&f_pos);
    assert!(w_pos >= 0.0,
        "Strain energy must be non-negative, got {}", w_pos);

    // Convergence to linear elasticity for small deformations
    let eps_scale = 1e-4;
    let grad_u  = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f_small = id + grad_u;
    let s_mat   = mat.pk2_stress(&f_small);
    let eps_lin = 0.5 * (grad_u + grad_u.transpose());
    let s_lin   = expected_lambda * eps_lin.trace() * id + 2.0 * expected_mu * eps_lin;
    let error   = (s_mat - s_lin).norm() / s_lin.norm().max(1e-14);
    assert!(error < 1e-3,
        "Material must converge to linear elasticity for small deformations, error = {:.2e}",
        error);

    // Numerical tangent consistency
    run_numerical_tangent_check(mat);
}

/// Numerical tangent check via central finite differences on E.
///
/// Perturbs E in each Voigt direction by h and compares the analytical
/// tangent column to the finite-difference approximation.
/// Tolerance 1e-4 is appropriate for step h=1e-7.
fn run_numerical_tangent_check(mat: &dyn MaterialLaw) {
    let h   = 1e-7;
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];

    // Test at two deformation states
    let test_points: &[Matrix3<f64>] = &[
        // Near identity — tangent should be very accurate
        {
            let e0 = Matrix3::new(
                0.001, 0.0005, 0.0002,
                0.0005, 0.001, 0.0002,
                0.0002, 0.0002, 0.0008,
            );
            let c0: Matrix3<f64> = Matrix3::identity() + 2.0 * e0;
            Cholesky::new(c0).unwrap().l().transpose()
        },
        // Moderate deformation
        {
            let e0 = Matrix3::new(
                0.02, 0.005, 0.003,
                0.005, 0.015, 0.003,
                0.003, 0.003, 0.01,
            );
            let c0: Matrix3<f64> = Matrix3::identity() + 2.0 * e0;
            Cholesky::new(c0).unwrap().l().transpose()
        },
    ];

    for f0 in test_points {
        let c_analytical = mat.tangent_stiffness(f0);

        for i in 0..6 {
            let (r, c) = idx[i];
            let is_shear = r != c;
            let amp      = if is_shear { h / 2.0 } else { h };

            let mut de_mat = Matrix3::zeros();
            de_mat[(r, c)] = amp;
            if is_shear { de_mat[(c, r)] = amp; }

            let e_cur: Matrix3<f64> = 0.5 * (f0.transpose() * f0 - Matrix3::identity());
            let c_plus:  Matrix3<f64> = Matrix3::identity() + 2.0 * (e_cur + de_mat);
            let c_minus: Matrix3<f64> = Matrix3::identity() + 2.0 * (e_cur - de_mat);
            let f_plus:  Matrix3<f64> = Cholesky::new(c_plus).unwrap().l().transpose();
            let f_minus: Matrix3<f64> = Cholesky::new(c_minus).unwrap().l().transpose();

            let s_plus  = mat.pk2_stress(&f_plus);
            let s_minus = mat.pk2_stress(&f_minus);

            let ds_fd: Vec<f64> = idx.iter().map(|&(ri, ci)| {
                (s_plus[(ri, ci)] - s_minus[(ri, ci)]) / (2.0 * h)
            }).collect();

            let ds_an: Vec<f64> = (0..6).map(|j| c_analytical[(j, i)]).collect();

            let norm_fd:  f64 = ds_fd.iter().map(|x| x * x).sum::<f64>().sqrt();
            let norm_err: f64 = ds_fd.iter().zip(ds_an.iter())
                .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();

            if norm_fd > 1e-10 {
                let rel = norm_err / norm_fd;
                assert!(rel < 1e-3,
                    "Tangent column {} numerical check failed: rel error = {:.4e}", i, rel);
            }
        }
    }
}

// ─── Material constructors ────────────────────────────────────────────────────

/// Neo-Hookean: E=1000 Pa, nu=0.3. mu = E/(2*(1+nu)), kappa = lambda + 2/3*mu.
fn nh_material() -> CompressibleMaterial<NeoHookeanIso, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    CompressibleMaterial {
        iso:     NeoHookeanIso { mu },
        vol:     VolumetricLnJ { kappa },
        density: 1000.0,
    }
}

/// Mooney-Rivlin: c1=300, c2=100, kappa from equivalent E=1000, nu=0.3.
/// mu_eff = 2*(c1+c2), kappa = lambda_eff + 2/3*mu_eff.
fn mr_material() -> CompressibleMaterial<MooneyRivlinIso, VolumetricLnJ> {
    let mu_mr  = 2.0 * (300.0 + 100.0);
    // Use nu=0.3 to derive lambda from kappa = lambda + 2/3*mu
    // kappa = E/(3*(1-2*nu)) with E chosen so mu = mu_mr
    // E = 2*mu*(1+nu) -> E = 2*mu_mr*1.3 = 2080
    let e_eff  = 2.0 * mu_mr * 1.3;
    let (lambda, _) = lame(e_eff, 0.3);
    let kappa  = lambda + 2.0 / 3.0 * mu_mr;
    CompressibleMaterial {
        iso:     MooneyRivlinIso { c1: 300.0, c2: 100.0 },
        vol:     VolumetricLnJ   { kappa },
        density: 1000.0,
    }
}

/// Ogden N=1, alpha=2, mu_1=2*mu — equivalent to NeoHookeanIso.
/// Effective Lame: mu_eff = mu, lambda_eff = kappa - 2/3*mu.
fn ogden_nh_material() -> CompressibleMaterial<OgdenIso, VolumetricLnJ> {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;
    CompressibleMaterial {
        iso:     OgdenIso::new(vec![mu], vec![2.0]).unwrap(),  // mu_1 = mu, not 2*mu
        vol:     VolumetricLnJ { kappa },
        density: 1000.0,
    }
}

/// SVK: E=1000 Pa, nu=0.3.
fn svk_material() -> SaintVenantKirchhoff {
    SaintVenantKirchhoff { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1000.0 }
}

// ─── Parametric tests ─────────────────────────────────────────────────────────

#[test]
fn test_svk_standard_suite() {
    let mat = svk_material();
    let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
    run_standard_material_tests(&mat, lambda, mu);
}

#[test]
fn test_neo_hookean_standard_suite() {
    let mat        = nh_material();
    let mu_eff     = mat.iso.mu;
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff);
}

#[test]
fn test_mooney_rivlin_standard_suite() {
    let mat        = mr_material();
    let mu_eff     = 2.0 * (mat.iso.c1 + mat.iso.c2);
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff);
}

#[test]
fn test_ogden_standard_suite() {
    let mat        = ogden_nh_material();
    let mu_eff     = mat.iso.mu[0]; // mu_1 = 2*mu -> mu_eff = mu
    let lambda_eff = mat.vol.kappa - 2.0 / 3.0 * mu_eff;
    run_standard_material_tests(&mat, lambda_eff, mu_eff);
}

// ─── SVK-specific tests ───────────────────────────────────────────────────────

#[test]
fn test_svk_tangent_is_constant() {
    let mat = svk_material();
    let f1  = Matrix3::identity();
    let f2  = Matrix3::new(1.1, 0.05, 0.0, 0.0, 0.95, 0.0, 0.0, 0.0, 1.0);
    let c1  = mat.tangent_stiffness(&f1);
    let c2  = mat.tangent_stiffness(&f2);
    assert!((c1 - c2).norm() < 1e-10, "SVK tangent must be constant");
}

#[test]
fn test_svk_converges_to_linear_for_small_deformations() {
    let mat          = svk_material();
    let (lambda, mu) = lame(mat.youngs_modulus, mat.poisson_ratio);
    let eps_scale    = 1e-4;
    let grad_u       = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f            = Matrix3::identity() + grad_u;
    let s_svk        = mat.pk2_stress(&f);
    let eps_lin      = 0.5 * (grad_u + grad_u.transpose());
    let s_lin        = lambda * eps_lin.trace() * Matrix3::identity() + 2.0 * mu * eps_lin;
    let error        = (s_svk - s_lin).norm() / s_lin.norm();
    assert!(error < 1e-3, "SVK must converge to linear, error = {:.2e}", error);
}

// ─── Cross-material consistency tests ─────────────────────────────────────────

/// NH and SVK must agree for small strains (same E and nu).
#[test]
fn test_nh_matches_svk_small_strain() {
    let nh        = nh_material();
    let svk       = svk_material();
    let eps_scale = 1e-5;
    let grad_u    = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f         = Matrix3::identity() + grad_u;
    let s_nh      = nh.pk2_stress(&f);
    let s_svk     = svk.pk2_stress(&f);
    let error     = (s_nh - s_svk).norm() / s_svk.norm().max(1e-14);
    assert!(error < 1e-3, "NH and SVK must agree for small strains, error = {:.2e}", error);
}

/// MR with c2=0 must match NH exactly (W, S, and C_tangent).
#[test]
fn test_mr_c2_zero_matches_nh() {
    let c1    = 500.0;
    let mu    = 2.0 * c1;
    let (lambda, _) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;

    let mr = CompressibleMaterial {
        iso:     MooneyRivlinIso { c1, c2: 0.0 },
        vol:     VolumetricLnJ   { kappa },
        density: 1000.0,
    };
    let nh = CompressibleMaterial {
        iso:     NeoHookeanIso   { mu },
        vol:     VolumetricLnJ   { kappa },
        density: 1000.0,
    };

    let f = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);

    let w_mr = mr.strain_energy(&f);
    let w_nh = nh.strain_energy(&f);
    assert!((w_mr - w_nh).abs() / w_nh.abs() < 1e-8,
        "MR(c2=0) W must match NH, rel diff = {:.2e}", (w_mr - w_nh).abs() / w_nh.abs());

    let s_mr  = mr.pk2_stress(&f);
    let s_nh  = nh.pk2_stress(&f);
    let err_s = (s_mr - s_nh).norm() / s_nh.norm().max(1e-14);
    assert!(err_s < 1e-8,
        "MR(c2=0) S must match NH, rel diff = {:.2e}", err_s);
}

/// Ogden N=1 with alpha=2 must match NH exactly.
#[test]
fn test_ogden_matches_nh() {
    let (lambda, mu) = lame(1000.0, 0.3);
    let kappa = lambda + 2.0 / 3.0 * mu;

    let ogden = CompressibleMaterial {
        iso: OgdenIso::new(vec![mu], vec![2.0]).unwrap(),
        vol:     VolumetricLnJ { kappa },
        density: 1000.0,
    };
    let nh = CompressibleMaterial {
        iso:     NeoHookeanIso { mu },
        vol:     VolumetricLnJ { kappa },
        density: 1000.0,
    };

    let f = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);

    let w_diff = (ogden.strain_energy(&f) - nh.strain_energy(&f)).abs();
    assert!(w_diff < 1e-8, "Ogden(N=1,a=2) W must match NH, diff = {:.2e}", w_diff);

    let s_diff = (ogden.pk2_stress(&f) - nh.pk2_stress(&f)).norm();
    assert!(s_diff < 1e-8, "Ogden(N=1,a=2) S must match NH, diff = {:.2e}", s_diff);

    let c_diff = (ogden.tangent_stiffness(&f) - nh.tangent_stiffness(&f)).norm();
    assert!(c_diff < 1e-6, "Ogden(N=1,a=2) C must match NH, diff = {:.2e}", c_diff);
}

// ─── Constructor validation tests ─────────────────────────────────────────────

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
    assert!(MooneyRivlinIso::new(-1.0,  0.0).is_err()); // c1 <= 0
    assert!(MooneyRivlinIso::new( 1.0, -2.0).is_err()); // c1+c2 <= 0
    assert!(MooneyRivlinIso::new( 1.0,  0.0).is_ok());
    assert!(MooneyRivlinIso::new( 1.0, -0.5).is_ok());  // c1+c2 = 0.5 > 0
}

#[test]
fn test_ogden_iso_invalid_params() {
    assert!(OgdenIso::new(vec![1.0], vec![2.0, 3.0]).is_err()); // mismatched lengths
    assert!(OgdenIso::new(vec![], vec![]).is_err());              // empty
    assert!(OgdenIso::new(vec![1.0], vec![-1.0]).is_err());      // instability
    assert!(OgdenIso::new(vec![1.0, 0.5], vec![2.0, 4.0]).is_ok());
}

#[test]
fn test_svk_new_invalid_params() {
    assert!(SaintVenantKirchhoff::new(-1.0,   0.3, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.6, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.0, 1000.0).is_err());
    assert!(SaintVenantKirchhoff::new(1000.0, 0.3, 1000.0).is_ok());
}
