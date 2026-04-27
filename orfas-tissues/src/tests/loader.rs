// UTF-8
// orfas-tissues/src/tests/loader.rs — tests for the JSON preset loader.

use crate::loader::{load_preset_from_str, save_preset_to_str, LoadError};

// ─── Valid JSON fixtures ───────────────────────────────────────────────────────

const LIVER_JSON: &str = r#"
{
    "name": "Liver",
    "model": "neo_hookean",
    "reference": "Nava et al. (2008), Med. Image Anal.",
    "doi": "10.1016/j.media.2007.09.001",
    "protocol": "ex vivo indentation",
    "parameters": {
        "mu": 2100.0,
        "kappa": 50000.0,
        "density": 1060.0
    },
    "confidence_intervals": {
        "mu": { "min": 1500.0, "max": 3000.0 }
    },
    "notes": "Porcine liver"
}
"#;

const BRAIN_JSON: &str = r#"
{
    "name": "Brain grey matter",
    "model": "mooney_rivlin",
    "reference": "Budday et al. (2017), Acta Biomater.",
    "doi": null,
    "protocol": "quasi-static compression",
    "parameters": {
        "c1": 310.0,
        "c2": 80.0,
        "kappa": 25000.0,
        "density": 1040.0
    },
    "confidence_intervals": {}
}
"#;

const CARDIAC_JSON: &str = r#"
{
    "name": "Myocardium",
    "model": "holzapfel_ogden",
    "reference": "Holzapfel & Ogden (2009)",
    "doi": null,
    "protocol": "biaxial testing",
    "parameters": {
        "mu": 59.0,
        "k1": 18472.0,
        "k2": 16.026,
        "kappa": 350000.0,
        "density": 1060.0
    },
    "confidence_intervals": {}
}
"#;

const SVK_JSON: &str = r#"
{
    "name": "Bone (cortical, linearized)",
    "model": "saint_venant_kirchhoff",
    "reference": "Generic cortical bone",
    "doi": null,
    "protocol": "compression test",
    "parameters": {
        "youngs_modulus": 17000000000.0,
        "poisson_ratio": 0.3,
        "density": 1900.0
    },
    "confidence_intervals": {}
}
"#;

// ─── Load success tests ───────────────────────────────────────────────────────

#[test]
fn test_load_neo_hookean_from_str() {
    let preset = load_preset_from_str(LIVER_JSON).expect("Failed to load liver JSON");
    assert_eq!(preset.metadata.name,  "Liver");
    assert_eq!(preset.metadata.model, "neo_hookean");
    assert_eq!(preset.metadata.doi,   Some("10.1016/j.media.2007.09.001".to_string()));
    assert_eq!(*preset.metadata.parameters.get("mu").unwrap(),      2100.0);
    assert_eq!(*preset.metadata.parameters.get("kappa").unwrap(),   50000.0);
    assert_eq!(*preset.metadata.parameters.get("density").unwrap(), 1060.0);
    // Material must be constructible — density accessible via trait
    assert!((preset.material.density() - 1060.0).abs() < 1e-10);
}

#[test]
fn test_load_mooney_rivlin_from_str() {
    let preset = load_preset_from_str(BRAIN_JSON).expect("Failed to load brain JSON");
    assert_eq!(preset.metadata.name,  "Brain grey matter");
    assert_eq!(preset.metadata.model, "mooney_rivlin");
    assert!(preset.metadata.doi.is_none());
    assert!((preset.material.density() - 1040.0).abs() < 1e-10);
}

#[test]
fn test_load_holzapfel_ogden_from_str() {
    let preset = load_preset_from_str(CARDIAC_JSON).expect("Failed to load cardiac JSON");
    assert_eq!(preset.metadata.name,  "Myocardium");
    assert_eq!(preset.metadata.model, "holzapfel_ogden");
    assert!((preset.material.density() - 1060.0).abs() < 1e-10);
}

#[test]
fn test_load_saint_venant_kirchhoff_from_str() {
    let preset = load_preset_from_str(SVK_JSON).expect("Failed to load SVK JSON");
    assert_eq!(preset.metadata.name,  "Bone (cortical, linearized)");
    assert_eq!(preset.metadata.model, "saint_venant_kirchhoff");
    assert!((preset.material.density() - 1900.0).abs() < 1e-10);
}

#[test]
fn test_load_default_density_when_absent() {
    // No "density" key in parameters — should default to 1000.0
    let json = r#"
    {
        "name": "Test",
        "model": "neo_hookean",
        "reference": "test",
        "doi": null,
        "protocol": "test",
        "parameters": { "mu": 1000.0, "kappa": 10000.0 },
        "confidence_intervals": {}
    }"#;
    let preset = load_preset_from_str(json).expect("Should load with default density");
    assert!((preset.material.density() - 1000.0).abs() < 1e-10);
}

#[test]
fn test_load_notes_default_empty_when_absent() {
    // "notes" field is optional — should default to empty string
    let preset = load_preset_from_str(BRAIN_JSON).expect("Failed to load brain JSON");
    assert_eq!(preset.metadata.notes, "");
}

// ─── Load error tests ─────────────────────────────────────────────────────────

#[test]
fn test_load_unknown_model_returns_error() {
    let json = r#"
    {
        "name": "Test",
        "model": "unknown_model",
        "reference": "test",
        "doi": null,
        "protocol": "test",
        "parameters": { "mu": 1000.0 },
        "confidence_intervals": {}
    }"#;
    let result = load_preset_from_str(json);
    assert!(matches!(result, Err(LoadError::UnknownModel(_))),
        "Expected UnknownModel error");
}

#[test]
fn test_load_missing_parameter_returns_error() {
    // neo_hookean requires "mu" and "kappa" — omit "kappa"
    let json = r#"
    {
        "name": "Test",
        "model": "neo_hookean",
        "reference": "test",
        "doi": null,
        "protocol": "test",
        "parameters": { "mu": 1000.0 },
        "confidence_intervals": {}
    }"#;
    let result = load_preset_from_str(json);
    assert!(matches!(result, Err(LoadError::MissingParameter(_))),
        "Expected MissingParameter error");
}

#[test]
fn test_load_invalid_json_returns_error() {
    let result = load_preset_from_str("{ not valid json }");
    assert!(matches!(result, Err(LoadError::Json(_))),
        "Expected Json parse error");
}

// ─── Round-trip test ──────────────────────────────────────────────────────────

#[test]
fn test_round_trip_serialize_deserialize() {
    // Load from JSON string, serialize back, load again — metadata must match
    let preset_a = load_preset_from_str(LIVER_JSON).expect("First load failed");
    let json_out = save_preset_to_str(&preset_a.metadata).expect("Serialize failed");
    let preset_b = load_preset_from_str(&json_out).expect("Second load failed");

    assert_eq!(preset_a.metadata.name,    preset_b.metadata.name);
    assert_eq!(preset_a.metadata.model,   preset_b.metadata.model);
    assert_eq!(preset_a.metadata.reference, preset_b.metadata.reference);
    assert!(
        (preset_a.metadata.parameters["mu"] - preset_b.metadata.parameters["mu"]).abs() < 1e-10,
        "mu parameter must survive round-trip"
    );
    assert!(
        (preset_a.metadata.parameters["kappa"] - preset_b.metadata.parameters["kappa"]).abs() < 1e-10,
        "kappa parameter must survive round-trip"
    );
}