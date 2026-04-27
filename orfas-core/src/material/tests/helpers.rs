// UTF-8
// material/tests/helpers.rs — shared parametric test helpers.
//
// `run_standard_material_tests` calls the same `check_*` functions from
// `material::consistency` that power `check_thermodynamic_consistency`,
// but panics via assert! instead of accumulating errors.
// The numerical tangent check and viscoelastic suite remain test-only.

use nalgebra::{Matrix3, Vector3, DVector};
use nalgebra::linalg::Cholesky;
use crate::material::traits::{MaterialLaw, AnisotropicPart, IsochoricPart, VolumetricPart};
use crate::material::helpers::{lame, hooke_voigt};
use crate::material::MaterialContext;
use crate::material::fiber_fields::FiberField;
use crate::material::internal_variables::ElementInternalVars;
use crate::material::viscoelastic::ViscoelasticMaterial;
use crate::material::consistency::{
    check_zero_energy_at_rest,
    check_zero_stress_at_rest,
    check_hooke_at_rest,
    check_tangent_symmetry,
    check_positive_strain_energy,
    check_small_strain_linearization,
    check_tangent_positive_definite,
    check_objectivity,
};

// ─── Standard suite ───────────────────────────────────────────────────────────

/// Run the standard thermodynamic consistency suite for a material law.
///
/// Checks (in order):
///   1. W(F=I) = 0
///   2. S(F=I) = 0
///   3. C(F=I) = Hooke  (requires expected_lambda, expected_mu)
///   4. C symmetric at arbitrary F
///   5. W(F) >= 0
///   6. S -> sigma_lin for small deformations
///   7. C(F=I) positive definite
///   8. W(QF) = W(F) for 4 canonical rotations
///   9. Numerical tangent consistency (finite differences)
///
/// Panics with a descriptive message on the first failed check (assert! behavior).
/// `ctx` is passed through to allow viscoelastic materials to carry their state.
pub fn run_standard_material_tests(
    mat:             &dyn MaterialLaw,
    expected_lambda: f64,
    expected_mu:     f64,
    ctx:             &mut MaterialContext,
) {
    let lame_params = Some((expected_lambda, expected_mu));

    // Checks 1-8 — reuse consistency.rs logic, panic on failure
    if let Some(e) = check_zero_energy_at_rest(mat)                    { panic!("{}", e); }
    if let Some(e) = check_zero_stress_at_rest(mat)                    { panic!("{}", e); }
    if let Some(e) = check_hooke_at_rest(mat, lame_params)             { panic!("{}", e); }
    if let Some(e) = check_tangent_symmetry(mat)                       { panic!("{}", e); }
    if let Some(e) = check_positive_strain_energy(mat)                 { panic!("{}", e); }
    if let Some(e) = check_small_strain_linearization(mat, lame_params){ panic!("{}", e); }
    if let Some(e) = check_tangent_positive_definite(mat)              { panic!("{}", e); }
    if let Some(e) = check_objectivity(mat)                            { panic!("{}", e); }

    // Check 9 — numerical tangent (test-only, not in consistency.rs)
    run_numerical_tangent_check(mat, ctx);
}

// ─── Numerical tangent check ──────────────────────────────────────────────────

/// Verify the analytical tangent stiffness against a finite-difference approximation.
///
/// For each Voigt column i, perturbs the Green-Lagrange strain tensor E by h
/// in the corresponding direction, computes dS/dE by central differences, and
/// compares against the analytical column of C_tangent.
///
/// Two representative test points are used (small and moderate deformations).
/// Relative error tolerance: 1e-3.
pub fn run_numerical_tangent_check(mat: &dyn MaterialLaw, ctx: &mut MaterialContext) {
    let h   = 1e-7;
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];

    let test_points: &[Matrix3<f64>] = &[
        {
            let e0 = Matrix3::new(0.001, 0.0005, 0.0002, 0.0005, 0.001, 0.0002, 0.0002, 0.0002, 0.0008);
            let c0: Matrix3<f64> = Matrix3::identity() + 2.0 * e0;
            Cholesky::new(c0).unwrap().l().transpose()
        },
        {
            let e0 = Matrix3::new(0.02, 0.005, 0.003, 0.005, 0.015, 0.003, 0.003, 0.003, 0.01);
            let c0: Matrix3<f64> = Matrix3::identity() + 2.0 * e0;
            Cholesky::new(c0).unwrap().l().transpose()
        },
    ];

    for f0 in test_points {
        let c_analytical = mat.tangent_stiffness(f0, &ctx);

        for i in 0..6 {
            let (r, c) = idx[i];
            let is_shear = r != c;
            let amp      = if is_shear { h / 2.0 } else { h };

            let mut de_mat = Matrix3::zeros();
            de_mat[(r, c)] = amp;
            if is_shear { de_mat[(c, r)] = amp; }

            let e_cur   = 0.5 * (f0.transpose() * f0 - Matrix3::identity());
            let c_plus  = Matrix3::identity() + 2.0 * (e_cur + de_mat);
            let c_minus = Matrix3::identity() + 2.0 * (e_cur - de_mat);
            let f_plus  = Cholesky::new(c_plus).unwrap().l().transpose();
            let f_minus = Cholesky::new(c_minus).unwrap().l().transpose();

            let s_plus  = mat.pk2_stress(&f_plus,  ctx);
            let s_minus = mat.pk2_stress(&f_minus, ctx);

            let ds_fd: Vec<f64> = idx.iter().map(|&(ri, ci)| {
                (s_plus[(ri, ci)] - s_minus[(ri, ci)]) / (2.0 * h)
            }).collect();
            let ds_an: Vec<f64> = (0..6).map(|j| c_analytical[(j, i)]).collect();

            let norm_fd  = ds_fd.iter().map(|x| x * x).sum::<f64>().sqrt();
            let norm_err = ds_fd.iter().zip(ds_an.iter())
                .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();

            if norm_fd > 1e-10 {
                let rel = norm_err / norm_fd;
                assert!(rel < 1e-3,
                    "Tangent column {} numerical check failed: rel error = {:.4e}", i, rel);
            }
        }
    }
}

// ─── Anisotropic part suite ───────────────────────────────────────────────────

/// Run standard checks for an AnisotropicPart implementation.
///
/// Checks:
///   - W_aniso(F=I) = 0
///   - S_aniso(F=I) = 0
///   - S_aniso = 0 under fiber compression (no tension-compression symmetry for fibers)
///   - Numerical tangent consistency for the anisotropic contribution
pub fn run_anisotropic_part_tests(mat: &dyn AnisotropicPart, ctx: &mut MaterialContext) {
    let f_id = Matrix3::identity();

    let w = mat.strain_energy_aniso(&f_id, &ctx);
    assert!(w.abs() < 1e-10, "W_aniso at F=I should be 0, got {:.2e}", w);

    let s = mat.pk2_stress_aniso(&f_id, ctx);
    assert!(s.norm() < 1e-10, "S_aniso at F=I should be 0, got norm {:.2e}", s.norm());

    let mut f_comp = Matrix3::identity();
    f_comp[(0, 0)] = 0.8;
    let ff_comp      = FiberField::uniform(1, &[Vector3::new(1.0, 0.0, 0.0)]);
    let mut ctx_comp = MaterialContext::from_fiber_field(&ff_comp, 0);
    let s_comp       = mat.pk2_stress_aniso(&f_comp, &mut ctx_comp);
    assert!(s_comp.norm() < 1e-10,
        "S_aniso should be 0 under fiber compression, got norm {:.2e}", s_comp.norm());

    let h   = 1e-7;
    let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
    let e0:  Matrix3<f64>  = Matrix3::new(0.08, 0.02, 0.01, 0.02, 0.04, 0.01, 0.01, 0.01, 0.03);
    let c0:  Matrix3<f64>  = Matrix3::identity() + 2.0 * e0;
    let f_def:  Matrix3<f64>        = Cholesky::new(c0).unwrap().l().transpose();
    let c_analytical = mat.tangent_aniso(&f_def, ctx);

    for i in 0..6 {
        let (r, c) = idx[i];
        let is_shear = r != c;
        let amp      = if is_shear { h / 2.0 } else { h };

        let mut de_mat = Matrix3::zeros();
        de_mat[(r, c)] = amp;
        if is_shear { de_mat[(c, r)] = amp; }

        let e_cur:  Matrix3<f64>  = 0.5 * (f_def.transpose() * f_def - Matrix3::identity());
        let c_plus:  Matrix3<f64> = Matrix3::identity() + 2.0 * (e_cur + de_mat);
        let c_minus:  Matrix3<f64>= Matrix3::identity() + 2.0 * (e_cur - de_mat);
        let f_plus:  Matrix3<f64> = Cholesky::new(c_plus).unwrap().l().transpose();
        let f_minus:  Matrix3<f64>= Cholesky::new(c_minus).unwrap().l().transpose();

        let s_plus  = mat.pk2_stress_aniso(&f_plus,  ctx);
        let s_minus = mat.pk2_stress_aniso(&f_minus, ctx);

        let ds_fd: Vec<f64> = idx.iter().map(|&(ri, ci)| {
            (s_plus[(ri, ci)] - s_minus[(ri, ci)]) / (2.0 * h)
        }).collect();
        let ds_an: Vec<f64> = (0..6).map(|j| c_analytical[(j, i)]).collect();

        let norm_fd  = ds_fd.iter().map(|x| x * x).sum::<f64>().sqrt();
        let norm_err = ds_fd.iter().zip(ds_an.iter())
            .map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();

        if norm_fd > 1e-10 {
            let rel = norm_err / norm_fd;
            assert!(rel < 5e-2, "C_aniso FD mismatch col {}: rel error = {:.4e}", i, rel);
        }
    }
}

// ─── Viscoelastic suite ───────────────────────────────────────────────────────

/// Run the viscoelastic-specific test suite.
///
/// Checks (in order):
///   1. Elastic fallback — S(F=I) = 0 without internal variables
///   2. Zero stress at F=I with initialized internal variables
///   3. Stress relaxation — S must converge to S_inf (elastic equilibrium)
///      under constant deformation after many time steps
///   4. Algorithmic tangent scaling — verifies (1 + delta_iso) factor
///      and consistency with finite differences on the elastic tangent
pub fn run_viscoelastic_tests<I, A, V>(mat: &ViscoelasticMaterial<I, A, V>)
where
    I: IsochoricPart,
    A: AnisotropicPart,
    V: VolumetricPart,
    ViscoelasticMaterial<I, A, V>: MaterialLaw,
{
    let f_id    = Matrix3::identity();
    let m_iso   = mat.m_iso();
    let m_aniso = mat.m_aniso();

    // Test 1 — Elastic fallback
    {
        let mut ctx = MaterialContext::default();
        let s = mat.pk2_stress(&f_id, &mut ctx);
        assert!(s.norm() < 1e-10,
            "Viscoelastic fallback: S at F=I must be zero, got {:.2e}", s.norm());
    }

    // Test 2 — Zero stress at F=I with iv
    {
        let iv      = ElementInternalVars::zeros(m_iso, m_aniso);
        let mut ctx = MaterialContext { dt: 0.1, fiber_dirs: &[], iv_ref: Some(&iv), iv: None };
        let s       = mat.pk2_stress(&f_id, &mut ctx);
        assert!(s.norm() < 1e-10,
            "Viscoelastic: S at F=I with iv must be zero, got {:.2e}", s.norm());
    }

    // Test 3 — Relaxation under constant deformation.
    // S must converge toward S_inf = S_eq (elastic equilibrium) after many steps.
    {
        let e0    = Matrix3::new(0.05, 0.01, 0.0, 0.01, 0.03, 0.0, 0.0, 0.0, 0.02);
        let f_def = Matrix3::identity() + e0;
        let s_inf = {
            let mut ctx = MaterialContext::default();
            mat.pk2_stress(&f_def, &mut ctx)
        };

        let dt         = 0.05;
        let mut iv     = ElementInternalVars::zeros(m_iso, m_aniso);
        let mut s_last = Matrix3::zeros();

        for _ in 0..200 {
            s_last = {
                let mut ctx = MaterialContext { dt, fiber_dirs: &[], iv_ref: Some(&iv), iv: None };
                mat.pk2_stress(&f_def, &mut ctx)
            };
            {
                let mut ctx = MaterialContext { dt, fiber_dirs: &[], iv_ref: None, iv: Some(&mut iv) };
                mat.update_state(&f_def, &mut ctx);
            }
        }

        let rel_err = (s_last - s_inf).norm() / s_inf.norm().max(1e-14);
        assert!(rel_err < 0.01,
            "Viscoelastic relaxation: S must converge to S_inf, rel_err = {:.4e}", rel_err);
    }

    // Test 4 — Tangent consistency.
    // The FD measures dS_eq/dF (sum_Q constant w.r.t. F).
    // Verify elastic tangent via FD, and algorithmic scaling (1 + delta_iso).
    {
        let e0    = Matrix3::new(0.05, 0.01, 0.0, 0.01, 0.03, 0.0, 0.0, 0.0, 0.02);
        let f_def = Matrix3::identity() + e0;
        let dt    = 0.1;
        let h     = 1e-7;
        let idx   = [(0usize,0usize),(1,1),(2,2),(0,1),(1,2),(0,2)];

        let c_elastic = {
            let ctx = MaterialContext::default();
            mat.tangent_stiffness(&f_def, &ctx)
        };

        let c_algorithmic = {
            let ctx = MaterialContext::from_dt(dt);
            mat.tangent_stiffness(&f_def, &ctx)
        };

        let delta_iso = mat.beta_iso.iter().zip(mat.tau_iso.iter())
            .map(|(b, t): (&f64, &f64)| b * (-dt / (2.0 * t)).exp())
            .sum::<f64>();

        let j     = f_def.determinant();
        let c     = f_def.transpose() * f_def;
        let c_inv = c.try_inverse().unwrap();
        let dw    = mat.vol.dw_vol(j);
        let d2w   = mat.vol.d2w_vol(j);
        let a_vol =  j * (j * d2w + dw);
        let b_vol = -j * dw;
        use crate::material::helpers::cinv_tangent_voigt;
        let c_vol  = cinv_tangent_voigt(&c_inv, a_vol, b_vol);
        let c_iso  = c_elastic - c_vol;
        let c_diff = c_algorithmic - c_elastic;

        let ratio = c_diff.norm() / c_iso.norm();
        assert!((ratio - delta_iso).abs() < 1e-6,
            "Algorithmic scaling incorrect: expected {:.6}, got {:.6}", delta_iso, ratio);

        // FD — elastic tangent check (sum_Q cancels in finite difference)
        for i in 0..6 {
            let (r, c) = idx[i];
            let is_shear = r != c;
            let amp = if is_shear { h / 2.0 } else { h };
            let mut de = Matrix3::zeros();
            de[(r,c)] = amp;
            if is_shear { de[(c,r)] = amp; }

            let e_cur = 0.5 * (f_def.transpose() * f_def - Matrix3::identity());
            let f_p   = Cholesky::new(Matrix3::identity() + 2.0*(e_cur+de)).unwrap().l().transpose();
            let f_m   = Cholesky::new(Matrix3::identity() + 2.0*(e_cur-de)).unwrap().l().transpose();

            let s_p = { let mut ctx = MaterialContext::default(); mat.pk2_stress(&f_p, &mut ctx) };
            let s_m = { let mut ctx = MaterialContext::default(); mat.pk2_stress(&f_m, &mut ctx) };

            let ds_fd: Vec<f64> = idx.iter().map(|&(ri,ci)|
                (s_p[(ri,ci)] - s_m[(ri,ci)]) / (2.0*h)
            ).collect();
            let ds_an: Vec<f64> = (0..6).map(|j| c_elastic[(j,i)]).collect();

            let norm_fd  = ds_fd.iter().map(|x| x*x).sum::<f64>().sqrt();
            let norm_err = ds_fd.iter().zip(ds_an.iter())
                .map(|(a,b)| (a-b).powi(2)).sum::<f64>().sqrt();

            if norm_fd > 1e-10 {
                let rel = norm_err / norm_fd;
                assert!(rel < 1e-3,
                    "Elastic tangent FD col {}: rel error = {:.4e}", i, rel);
            }
        }

        if !mat.tau_aniso.is_empty() {
            assert!(
                (c_algorithmic - c_elastic).norm() > 1e-10,
                "Algorithmic tangent must differ from elastic tangent when dt > 0"
            );
        }
    }
}