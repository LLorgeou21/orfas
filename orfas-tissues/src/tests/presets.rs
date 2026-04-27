// UTF-8
// orfas-tissues/src/tests/presets.rs — tests for all hardcoded tissue presets.
//
// For each preset:
//   - thermodynamic consistency check via check_thermodynamic_consistency
//   - nominal parameter values fall within their reported confidence intervals
//   - metadata fields are non-empty
//
// Anisotropic presets (HGO-based) skip the Hooke check (lame = None)
// because their effective Lame parameters depend on fiber directions
// which are not set here.

use orfas_core::{check_thermodynamic_consistency, MaterialContext, lame};
use crate::presets::{
    LiverNeoHookean,
    BrainGreyMatter, BrainWhiteMatter,
    CardiacMyocardium,
    ArterialWallMedia,
    TendonGroundMatrix,
    LigamentMCL,
    SkinMooneyRivlin,
    KidneyNeoHookean,
    ProstateNeoHookean,
};
use crate::metadata::TissuePreset;
use crate::presets::all_presets;

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Assert that check_thermodynamic_consistency passes for a material.
///
/// Panics with the full list of violations if any check fails.
fn assert_thermodynamic_consistency(
    preset:     &dyn TissuePreset,
    lame_params: Option<(f64, f64)>,
) {
    let mat = preset.material();
    let mut ctx = MaterialContext::default();
    let result = check_thermodynamic_consistency(mat.as_ref(), lame_params, &mut ctx);
    if let Err(errors) = result {
        panic!(
            "Thermodynamic consistency failed for '{}':\n  - {}",
            preset.metadata().name,
            errors.join("\n  - ")
        );
    }
}

/// Assert that all nominal parameter values fall within their confidence intervals.
///
/// Each key in confidence_intervals must be present in `nominal_values`.
fn assert_within_confidence(
    preset:         &dyn TissuePreset,
    nominal_values: &[(&'static str, f64)],
) {
    let meta = preset.metadata();
    let map: std::collections::HashMap<&str, f64> = nominal_values.iter().copied().collect();

    for (param, ci) in &meta.confidence_intervals {
        let value = map.get(param).copied().unwrap_or_else(|| {
            panic!(
                "Preset '{}': nominal value missing for parameter '{}'",
                meta.name, param
            )
        });
        assert!(
            ci.contains(value),
            "Preset '{}': parameter '{}' = {:.4e} outside CI [{:.4e}, {:.4e}]",
            meta.name, param, value, ci.min, ci.max
        );
    }
}

/// Assert that all metadata string fields are non-empty.
fn assert_metadata_non_empty(preset: &dyn TissuePreset) {
    let meta = preset.metadata();
    assert!(!meta.name.is_empty(),      "Preset name is empty");
    assert!(!meta.model.is_empty(),     "Preset model is empty");
    assert!(!meta.reference.is_empty(), "Preset reference is empty");
    assert!(!meta.protocol.is_empty(),  "Preset protocol is empty");
}

// ─── Liver ────────────────────────────────────────────────────────────────────

#[test]
fn test_liver_thermodynamic_consistency() {
    let preset = LiverNeoHookean::default();
    // mu=2100, kappa=50000 -> lambda = kappa - 2/3*mu = 50000 - 1400 = 48600
    let mu     = 2100.0_f64;
    let kappa  = 50000.0_f64;
    let lambda = kappa - 2.0 / 3.0 * mu;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu)));
}

#[test]
fn test_liver_confidence_intervals() {
    let preset = LiverNeoHookean::default();
    assert_within_confidence(&preset, &[("mu", 2100.0), ("kappa", 50000.0)]);
}

#[test]
fn test_liver_metadata() {
    assert_metadata_non_empty(&LiverNeoHookean::default());
}

// ─── Brain grey matter ────────────────────────────────────────────────────────

#[test]
fn test_brain_grey_thermodynamic_consistency() {
    let preset  = BrainGreyMatter::default();
    let c1      = 310.0_f64;
    let c2      = 80.0_f64;
    let kappa   = 25000.0_f64;
    let mu_eff  = 2.0 * (c1 + c2);
    let lambda  = kappa - 2.0 / 3.0 * mu_eff;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu_eff)));
}

#[test]
fn test_brain_grey_confidence_intervals() {
    let preset = BrainGreyMatter::default();
    assert_within_confidence(&preset, &[("c1", 310.0), ("c2", 80.0), ("kappa", 25000.0)]);
}

#[test]
fn test_brain_grey_metadata() {
    assert_metadata_non_empty(&BrainGreyMatter::default());
}

// ─── Brain white matter ───────────────────────────────────────────────────────

#[test]
fn test_brain_white_thermodynamic_consistency() {
    let preset  = BrainWhiteMatter::default();
    let c1      = 160.0_f64;
    let c2      = 40.0_f64;
    let kappa   = 25000.0_f64;
    let mu_eff  = 2.0 * (c1 + c2);
    let lambda  = kappa - 2.0 / 3.0 * mu_eff;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu_eff)));
}

#[test]
fn test_brain_white_confidence_intervals() {
    let preset = BrainWhiteMatter::default();
    assert_within_confidence(&preset, &[("c1", 160.0), ("c2", 40.0), ("kappa", 25000.0)]);
}

#[test]
fn test_brain_white_metadata() {
    assert_metadata_non_empty(&BrainWhiteMatter::default());
}

// ─── Cardiac ─────────────────────────────────────────────────────────────────

#[test]
fn test_cardiac_thermodynamic_consistency() {
    // HGO — fiber directions not set, skip Hooke check
    let preset = CardiacMyocardium::default();
    assert_thermodynamic_consistency(&preset, None);
}

#[test]
fn test_cardiac_confidence_intervals() {
    let preset = CardiacMyocardium::default();
    assert_within_confidence(&preset, &[
        ("mu",    59.0),
        ("k1",    18472.0),
        ("k2",    16.026),
        ("kappa", 350000.0),
    ]);
}

#[test]
fn test_cardiac_metadata() {
    assert_metadata_non_empty(&CardiacMyocardium::default());
}

// ─── Arterial ────────────────────────────────────────────────────────────────

#[test]
fn test_arterial_thermodynamic_consistency() {
    // HGO — fiber directions not set, skip Hooke check
    let preset = ArterialWallMedia::default();
    assert_thermodynamic_consistency(&preset, None);
}

#[test]
fn test_arterial_confidence_intervals() {
    let preset = ArterialWallMedia::default();
    assert_within_confidence(&preset, &[
        ("mu",    3000.0),
        ("k1",    2363.0),
        ("k2",    0.8393),
        ("kappa", 300000.0),
    ]);
}

#[test]
fn test_arterial_metadata() {
    assert_metadata_non_empty(&ArterialWallMedia::default());
}

// ─── Tendon ───────────────────────────────────────────────────────────────────

#[test]
fn test_tendon_thermodynamic_consistency() {
    let preset = TendonGroundMatrix::default();
    let mu     = 1.2e6_f64;
    let kappa  = 1.0e8_f64;
    let lambda = kappa - 2.0 / 3.0 * mu;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu)));
}

#[test]
fn test_tendon_confidence_intervals() {
    let preset = TendonGroundMatrix::default();
    assert_within_confidence(&preset, &[("mu", 1.2e6), ("kappa", 1.0e8)]);
}

#[test]
fn test_tendon_metadata() {
    assert_metadata_non_empty(&TendonGroundMatrix::default());
}

// ─── Ligament ────────────────────────────────────────────────────────────────

#[test]
fn test_ligament_thermodynamic_consistency() {
    // HGO — fiber directions not set, skip Hooke check
    let preset = LigamentMCL::default();
    assert_thermodynamic_consistency(&preset, None);
}

#[test]
fn test_ligament_confidence_intervals() {
    let preset = LigamentMCL::default();
    assert_within_confidence(&preset, &[
        ("mu",    1400.0),
        ("k1",    570000.0),
        ("k2",    48.0),
        ("kappa", 1.0e8),
    ]);
}

#[test]
fn test_ligament_metadata() {
    assert_metadata_non_empty(&LigamentMCL::default());
}

// ─── Skin ─────────────────────────────────────────────────────────────────────

#[test]
fn test_skin_thermodynamic_consistency() {
    let preset  = SkinMooneyRivlin::default();
    let c1      = 5000.0_f64;
    let c2      = 1000.0_f64;
    let kappa   = 150000.0_f64;
    let mu_eff  = 2.0 * (c1 + c2);
    let lambda  = kappa - 2.0 / 3.0 * mu_eff;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu_eff)));
}

#[test]
fn test_skin_confidence_intervals() {
    let preset = SkinMooneyRivlin::default();
    assert_within_confidence(&preset, &[("c1", 5000.0), ("c2", 1000.0), ("kappa", 150000.0)]);
}

#[test]
fn test_skin_metadata() {
    assert_metadata_non_empty(&SkinMooneyRivlin::default());
}

// ─── Kidney ───────────────────────────────────────────────────────────────────

#[test]
fn test_kidney_thermodynamic_consistency() {
    let preset = KidneyNeoHookean::default();
    let mu     = 1800.0_f64;
    let kappa  = 40000.0_f64;
    let lambda = kappa - 2.0 / 3.0 * mu;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu)));
}

#[test]
fn test_kidney_confidence_intervals() {
    let preset = KidneyNeoHookean::default();
    assert_within_confidence(&preset, &[("mu", 1800.0), ("kappa", 40000.0)]);
}

#[test]
fn test_kidney_metadata() {
    assert_metadata_non_empty(&KidneyNeoHookean::default());
}

// ─── Prostate ─────────────────────────────────────────────────────────────────

#[test]
fn test_prostate_thermodynamic_consistency() {
    let preset = ProstateNeoHookean::default();
    let mu     = 900.0_f64;
    let kappa  = 20000.0_f64;
    let lambda = kappa - 2.0 / 3.0 * mu;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu)));
}

#[test]
fn test_prostate_confidence_intervals() {
    let preset = ProstateNeoHookean::default();
    assert_within_confidence(&preset, &[("mu", 900.0), ("kappa", 20000.0)]);
}

#[test]
fn test_prostate_metadata() {
    assert_metadata_non_empty(&ProstateNeoHookean::default());
}

// ─── all_presets ──────────────────────────────────────────────────────────────

#[test]
fn test_all_presets_count() {
    assert_eq!(all_presets().len(), 10, "Expected 10 built-in presets");
}

#[test]
fn test_all_presets_metadata_non_empty() {
    for preset in all_presets() {
        assert_metadata_non_empty(preset.as_ref());
    }
}

#[test]
fn test_all_presets_thermodynamic_consistency_no_lame() {
    // Quick smoke test — skip Hooke check for all (lame = None).
    // Per-preset tests above cover the full check with correct lame params.
    for preset in all_presets() {
        let mat     = preset.material();
        let mut ctx = MaterialContext::default();
        let result  = check_thermodynamic_consistency(mat.as_ref(), None, &mut ctx);
        if let Err(errors) = result {
            panic!(
                "Thermodynamic consistency failed for '{}':\n  - {}",
                preset.metadata().name,
                errors.join("\n  - ")
            );
        }
    }
}