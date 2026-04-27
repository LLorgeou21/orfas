// UTF-8
// material/consistency.rs — runtime thermodynamic consistency checks for MaterialLaw.
//
// Provides `check_thermodynamic_consistency`, callable at runtime (e.g. after loading
// a material from JSON) to verify that a MaterialLaw implementation satisfies the
// necessary thermodynamic conditions.
//
// Each individual check is a `pub(crate)` function returning `Option<String>`:
//   - None  => check passed
//   - Some  => check failed, message describes the violation
//
// This allows `run_standard_material_tests` (tests/helpers.rs) to reuse the same
// logic via assert!, while this module accumulates errors into a Vec<String>.

use nalgebra::{Matrix3, Matrix6};
use super::traits::MaterialLaw;
use super::context::MaterialContext;
use super::helpers::hooke_voigt;

// ─── Individual checks ────────────────────────────────────────────────────────

/// Check 1 — Zero strain energy at rest: W(F=I) = 0.
///
/// A material must store no energy in its reference configuration.
pub(crate) fn check_zero_energy_at_rest(mat: &dyn MaterialLaw) -> Option<String> {
    let ctx = MaterialContext::default();
    let w   = mat.strain_energy(&Matrix3::identity(), &ctx);
    if w.abs() >= 1e-10 {
        Some(format!("W(F=I) must be zero, got {:.2e}", w))
    } else {
        None
    }
}

/// Check 2 — Zero stress at rest: S(F=I) = 0.
///
/// No stress in the reference configuration (natural state condition).
pub(crate) fn check_zero_stress_at_rest(mat: &dyn MaterialLaw) -> Option<String> {
    let mut ctx = MaterialContext::default();
    let s       = mat.pk2_stress(&Matrix3::identity(), &mut ctx);
    let norm    = s.norm();
    if norm >= 1e-10 {
        Some(format!("S(F=I) must be zero, got norm {:.2e}", norm))
    } else {
        None
    }
}

/// Check 3 — Tangent matches Hooke at rest: C(F=I) = lambda*(I tensor I) + 2*mu*I4.
///
/// The linearized response at F=I must recover classical isotropic linear elasticity.
/// Skipped if `lame` is None — caller does not know the equivalent Lame parameters.
pub(crate) fn check_hooke_at_rest(
    mat:  &dyn MaterialLaw,
    lame: Option<(f64, f64)>,
) -> Option<String> {
    let (lambda, mu) = lame?;
    let ctx          = MaterialContext::default();
    let c_mat        = mat.tangent_stiffness(&Matrix3::identity(), &ctx);
    let c_hooke      = hooke_voigt(lambda, mu);
    let diff         = (c_mat - c_hooke).norm();
    if diff >= 1e-6 {
        Some(format!(
            "C_tangent(F=I) must match Hooke (lambda={:.3e}, mu={:.3e}), diff = {:.2e}",
            lambda, mu, diff
        ))
    } else {
        None
    }
}

/// Check 4 — Tangent symmetry: C(F) = C(F)^T for an arbitrary deformation.
///
/// The material tangent stiffness must be symmetric for any admissible F.
/// Tested on one representative deformation gradient.
pub(crate) fn check_tangent_symmetry(mat: &dyn MaterialLaw) -> Option<String> {
    let ctx   = MaterialContext::default();
    let f_arb = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);
    let c     = mat.tangent_stiffness(&f_arb, &ctx);
    let abs_err = (&c - c.transpose()).abs().max();
    let rel_err = abs_err / c.norm().max(1e-14);
    if rel_err >= 1e-12 {
        Some(format!(
            "C_tangent must be symmetric, rel asymmetry = {:.2e} (abs = {:.2e})",
            rel_err, abs_err
        ))
    } else {
        None
    }
}

/// Check 5 — Non-negative strain energy: W(F) >= 0 for admissible F (J > 0).
///
/// Tested on one representative deformation gradient with positive determinant.
pub(crate) fn check_positive_strain_energy(mat: &dyn MaterialLaw) -> Option<String> {
    let ctx   = MaterialContext::default();
    let f_pos = Matrix3::new(1.2, 0.1, 0.0, 0.0, 0.9, 0.05, 0.0, 0.0, 1.0);
    debug_assert!(f_pos.determinant() > 0.0);
    let w = mat.strain_energy(&f_pos, &ctx);
    if w < 0.0 {
        Some(format!("Strain energy must be non-negative, got {:.4e}", w))
    } else {
        None
    }
}

/// Check 6 — Small-strain linearization: S(F) -> sigma_lin for small grad_u.
///
/// For infinitesimal deformations, the PK2 stress must converge to the linear
/// elastic stress: sigma = lambda * tr(eps) * I + 2 * mu * eps.
/// Skipped if `lame` is None.
pub(crate) fn check_small_strain_linearization(
    mat:  &dyn MaterialLaw,
    lame: Option<(f64, f64)>,
) -> Option<String> {
    let (lambda, mu) = lame?;
    let id           = Matrix3::identity();
    let mut ctx      = MaterialContext::default();
    let eps_scale    = 1e-4;
    let grad_u       = Matrix3::new(0.3, 0.1, 0.05, 0.1, 0.2, 0.0, 0.05, 0.0, 0.1) * eps_scale;
    let f_small      = id + grad_u;
    let s_mat        = mat.pk2_stress(&f_small, &mut ctx);
    let eps_lin      = 0.5 * (grad_u + grad_u.transpose());
    let s_lin        = lambda * eps_lin.trace() * id + 2.0 * mu * eps_lin;
    let error        = (s_mat - s_lin).norm() / s_lin.norm().max(1e-14);
    if error >= 1e-3 {
        Some(format!(
            "Material must linearize to Hooke for small deformations, rel error = {:.2e}",
            error
        ))
    } else {
        None
    }
}

/// Check 7 — Positive definiteness of C at rest: all eigenvalues of C(F=I) > 0.
///
/// A positive definite tangent at F=I guarantees local material stability
/// and is necessary for the FEM stiffness matrix to be positive definite.
/// Tested at F=I only — large-deformation stability is not assessed here.
pub(crate) fn check_tangent_positive_definite(mat: &dyn MaterialLaw) -> Option<String> {
    let ctx        = MaterialContext::default();
    let c: Matrix6<f64> = mat.tangent_stiffness(&Matrix3::identity(), &ctx);
    // Symmetrize to avoid spurious imaginary parts from floating point asymmetry
    let c_sym      = 0.5 * (&c + c.transpose());
    let eigenvalues = c_sym.symmetric_eigenvalues();
    let min_eig    = eigenvalues.min();
    if min_eig <= 0.0 {
        Some(format!(
            "C_tangent(F=I) must be positive definite, min eigenvalue = {:.4e}",
            min_eig
        ))
    } else {
        None
    }
}

/// Check 8 — Frame objectivity: W(Q*F) = W(F) for rigid rotations Q.
///
/// A hyperelastic material energy must be invariant under rigid body rotations
/// applied to the current configuration. Tested on 4 rotations:
///   - 90 deg around X
///   - 90 deg around Y
///   - 90 deg around Z
///   - 45 deg around the normalized (1, 1, 1) axis (non-canonical case)
pub(crate) fn check_objectivity(mat: &dyn MaterialLaw) -> Option<String> {
    let f_ref = Matrix3::new(1.2, 0.1, 0.05, 0.0, 0.9, 0.02, 0.0, 0.0, 1.1);
    let ctx   = MaterialContext::default();
    let w_ref = mat.strain_energy(&f_ref, &ctx);

    // Build a rotation matrix from axis (unit vector) and angle (radians)
    // using the Rodrigues formula: R = I*cos + (1-cos)*(a tensor a) + sin*[a]_x
    let rotation = |ax: [f64; 3], angle: f64| -> Matrix3<f64> {
        let (s, c) = angle.sin_cos();
        let [x, y, z] = ax;
        Matrix3::new(
            c + x*x*(1.0-c),   x*y*(1.0-c) - z*s, x*z*(1.0-c) + y*s,
            y*x*(1.0-c) + z*s, c + y*y*(1.0-c),   y*z*(1.0-c) - x*s,
            z*x*(1.0-c) - y*s, z*y*(1.0-c) + x*s, c + z*z*(1.0-c),
        )
    };

    let half_pi = std::f64::consts::FRAC_PI_2;
    let sqrt3   = 3.0_f64.sqrt();

    let rotations: &[(Matrix3<f64>, &str)] = &[
        (rotation([1.0, 0.0, 0.0], half_pi),                   "90 deg around X"),
        (rotation([0.0, 1.0, 0.0], half_pi),                   "90 deg around Y"),
        (rotation([0.0, 0.0, 1.0], half_pi),                   "90 deg around Z"),
        (rotation([1.0/sqrt3, 1.0/sqrt3, 1.0/sqrt3], half_pi), "90 deg around (1,1,1)"),
    ];

    for (q, label) in rotations {
        let f_rotated = q * f_ref;
        let w_rotated = mat.strain_energy(&f_rotated, &ctx);
        let err       = (w_rotated - w_ref).abs() / w_ref.abs().max(1e-14);
        if err >= 1e-8 {
            return Some(format!(
                "Objectivity violated for rotation {}: W(QF) = {:.6e}, W(F) = {:.6e}, rel err = {:.2e}",
                label, w_rotated, w_ref, err
            ));
        }
    }
    None
}

// ─── Public API ───────────────────────────────────────────────────────────────

/// Check that a material law satisfies the necessary thermodynamic conditions.
///
/// Runs 8 checks:
///   1. W(F=I) = 0          — zero energy at rest
///   2. S(F=I) = 0          — zero stress at rest
///   3. C(F=I) = Hooke      — linearized tangent matches linear elasticity (skipped if lame=None)
///   4. C symmetric         — tangent symmetry at an arbitrary deformation
///   5. W(F) >= 0           — non-negative strain energy
///   6. S -> sigma_lin      — small-strain linearization (skipped if lame=None)
///   7. C positive definite — all eigenvalues of C(F=I) > 0
///   8. W(QF) = W(F)        — frame objectivity for 4 canonical rotations
///
/// # Parameters
/// - `mat`  — material law to check
/// - `lame` — equivalent Lame parameters (lambda, mu) at F=I, or None to skip
///            checks 3 and 6
/// - `ctx`  — evaluation context; use `MaterialContext::default()` for isotropic
///            elastic materials
///
/// # Returns
/// - `Ok(())` if all checks pass
/// - `Err(errors)` where `errors` is the list of all failed checks
pub fn check_thermodynamic_consistency(
    mat:  &dyn MaterialLaw,
    lame: Option<(f64, f64)>,
    _ctx: &mut MaterialContext,
) -> Result<(), Vec<String>> {
    let mut errors = Vec::new();

    if let Some(e) = check_zero_energy_at_rest(mat)            { errors.push(e); }
    if let Some(e) = check_zero_stress_at_rest(mat)            { errors.push(e); }
    if let Some(e) = check_hooke_at_rest(mat, lame)            { errors.push(e); }
    if let Some(e) = check_tangent_symmetry(mat)               { errors.push(e); }
    if let Some(e) = check_positive_strain_energy(mat)         { errors.push(e); }
    if let Some(e) = check_small_strain_linearization(mat, lame) { errors.push(e); }
    if let Some(e) = check_tangent_positive_definite(mat)      { errors.push(e); }
    if let Some(e) = check_objectivity(mat)                    { errors.push(e); }

    if errors.is_empty() { Ok(()) } else { Err(errors) }
}