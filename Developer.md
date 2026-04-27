# DEVELOPER.md

Technical reference for contributors. Covers codebase architecture, conventions,
and step-by-step guides for the most common contribution tasks.

---

## Workspace structure

```
orfas/
├── orfas-core/        # Core library — mesh, materials, assembly, solvers, integrators
├── orfas-io/          # Mesh I/O — .vtk reader
├── orfas-tissues/     # Tissue preset library — calibrated materials from literature
└── orfas-viewer/      # Interactive egui viewer
```

### orfas-core

```
orfas-core/src/
├── lib.rs                      # Public re-exports
├── mesh.rs                     # Node, Tetrahedron, Mesh
├── mechanical_state.rs         # MechanicalState (dynamic DOF vector)
├── boundary.rs                 # BoundaryConditions, Constraint, Load, PenaltyMethod, EliminationMethod
├── damping.rs                  # RayleighDamping
├── integrator.rs               # ImplicitEulerIntegrator
├── sparse.rs                   # CsrMatrix, CgSolver, ILU(0)
├── solver.rs                   # DenseSolver, DirectSolver, NewtonRaphson, NewtonRaphsonCachedK
├── assembler/
│   ├── mod.rs                  # Assembler, public API
│   ├── geometry.rs             # BMatrix trait, LinearBMatrix, tetra geometry helpers
│   ├── pattern.rs              # CSR sparsity pattern pre-computation
│   ├── assembly.rs             # assemble_tangent, assemble_internal_forces, assemble_mass
│   └── tests.rs                # Assembler tests
└── material/
    ├── mod.rs                  # Public re-exports
    ├── traits.rs               # MaterialLaw, IsochoricPart, VolumetricPart, AnisotropicPart
    ├── context.rs              # MaterialContext, SimulationContext
    ├── helpers.rs              # lame(), hooke_voigt(), cinv_tangent_voigt(), ...
    ├── consistency.rs          # check_thermodynamic_consistency, check_* functions
    ├── svk.rs                  # SaintVenantKirchhoff
    ├── compressible.rs         # CompressibleMaterial<I,V>, CompressibleAnisotropicMaterial<I,A,V>
    ├── isochoric/              # NeoHookeanIso, MooneyRivlinIso, OgdenIso
    ├── volumetric.rs           # VolumetricLnJ, VolumetricQuad
    ├── anisotropic/            # HolzapfelOgden, NoAnisotropy
    ├── viscoelastic.rs         # ViscoelasticMaterial<I,A,V>
    ├── fiber_fields.rs         # FiberField
    ├── internal_variables.rs   # InternalVariables, ElementInternalVars
    └── tests/                  # run_standard_material_tests, run_numerical_tangent_check, ...
```

### orfas-tissues

```
orfas-tissues/src/
├── lib.rs              # Public re-exports
├── metadata.rs         # TissueMetadata, ConfidenceInterval, TissuePreset trait
├── loader.rs           # TissueMetadataOwned, load_preset_from_file, load_preset_from_str
├── presets/
│   ├── mod.rs          # all_presets()
│   ├── liver.rs
│   ├── brain.rs        # BrainGreyMatter, BrainWhiteMatter
│   ├── cardiac.rs
│   ├── arterial.rs
│   ├── tendon.rs
│   ├── ligament.rs
│   ├── skin.rs
│   ├── kidney.rs
│   └── prostate.rs
└── tests/
    ├── mod.rs
    ├── presets.rs      # Thermodynamic consistency + CI validation for all presets
    └── loader.rs       # JSON load/save round-trip tests
```

### orfas-viewer

```
orfas-viewer/src/
├── main.rs             # Entry point
├── app.rs              # MyEguiApp, UI panels (mesh, material, solver, BC, node inspector)
├── state.rs            # AppState, MaterialConfig, MaterialParams, make_material
├── simulation.rs       # build_simulation, run_simulation_static, init_simulation_dynamic
└── render.rs           # depth_sorted_nodes, screen_to_node
```

---

## Key design principles

### Lagrangian frame, PK2 stress

All material computations are in the reference configuration. `pk2_stress(F)` returns S
(2nd Piola-Kirchhoff). The assembler computes P = F·S (1st Piola-Kirchhoff) internally.
`tangent_stiffness(F)` returns dS/dE in Voigt notation (6×6).

Voigt ordering: [11, 22, 33, 12, 23, 13] — engineering shear convention (no extra factors
in C; they live in B). Do not change this convention — it is consistent throughout the codebase.

### Isochoric/volumetric split

Hyperelastic materials follow the Flory decoupled form W = W_iso(C̄) + W_vol(J).
`CompressibleMaterial<I, V>` composes any `IsochoricPart` and `VolumetricPart`.
`CompressibleAnisotropicMaterial<I, A, V>` adds an `AnisotropicPart`.
`ViscoelasticMaterial<I, A, V>` wraps any composition with Prony series dissipation.

### Geometry cache

`Assembler::new(&mesh)` pre-computes all element volumes and shape function gradients
(b, c, d coefficients). The mesh is assumed immutable during simulation until v0.16.0.
Do not add implicit assumptions about mesh mutability.

### Reduced vs full DOF space

`MechanicalState` lives in the reduced DOF space (free DOFs only). The assembler works in
the full space. `bc_result.reconstruct` maps from reduced to full. This separation is
fundamental — respect it in all new code.

### Generic hot paths

Performance-critical functions (`assemble_tangent::<B>`) are generic over `BMatrix` to
allow monomorphization. Less critical paths (`make_material`) return `Box<dyn MaterialLaw>`.
Do not change this principle.

---

## Code conventions

- All source files must start with `// UTF-8`
- All code and comments in English — no French
- Every public item must have a `///` doc comment
- Run `cargo fmt` and `cargo clippy` before submitting
- Tests live in `tests/` subdirectories, not inline in source files (except assembler tests)
- Use `debug_assert!` for invariants that should never fire in release builds

---

## How to add a new material law

A material law can be added at three levels depending on its structure.

### Level 1 — New isochoric model

Implement `IsochoricPart` in a new file under `orfas-core/src/material/isochoric/`:

```rust
use nalgebra::{Matrix3, Matrix6};
use crate::material::traits::IsochoricPart;

pub struct MyIsochoricModel {
    pub mu: f64,
}

impl IsochoricPart for MyIsochoricModel {
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64 { ... }
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64> { ... }
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64> { ... }
}
```

Then compose it:
```rust
let mat = CompressibleMaterial {
    iso:     MyIsochoricModel { mu: 1000.0 },
    vol:     VolumetricLnJ { kappa: 50000.0 },
    density: 1000.0,
};
```

### Level 2 — New full material law

Implement `MaterialLaw` directly (e.g. for coupled formulations like SVK):

```rust
impl MaterialLaw for MyMaterial {
    fn density(&self) -> f64 { ... }
    fn strain_energy(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> f64 { ... }
    fn pk2_stress(&self, f: &Matrix3<f64>, ctx: &mut MaterialContext) -> Matrix3<f64> { ... }
    fn tangent_stiffness(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> Matrix6<f64> { ... }
    // update_state has a default no-op — override only for history-dependent materials
}
```

### Level 3 — Viscoelastic material

Wrap any `CompressibleMaterial` or `CompressibleAnisotropicMaterial` with `ViscoelasticMaterial`:

```rust
let mat = ViscoelasticMaterial {
    iso:        NeoHookeanIso { mu },
    aniso:      NoAnisotropy,
    vol:        VolumetricLnJ { kappa },
    density:    1000.0,
    tau_iso:    vec![1.0],        // relaxation times (s)
    beta_iso:   vec![0.3],        // relative moduli
    tau_aniso:  vec![],
    beta_aniso: vec![],
};
```

### Validation

After implementing, run the standard test suite in your test file:

```rust
use crate::material::tests::helpers::run_standard_material_tests;

#[test]
fn test_my_material() {
    let mat = MyMaterial { ... };
    run_standard_material_tests(&mat, lambda, mu, &mut MaterialContext::default());
}
```

`run_standard_material_tests` checks: W(I)=0, S(I)=0, C(I)=Hooke, tangent symmetry,
non-negative energy, small-strain linearization, positive definiteness, and objectivity.
It also runs a finite-difference tangent check.

You can also call `check_thermodynamic_consistency` at runtime for materials loaded from JSON.

---

## How to add a tissue preset

### Step 1 — Create the preset file

Add a new file in `orfas-tissues/src/presets/`, e.g. `muscle.rs`:

```rust
// UTF-8
// orfas-tissues/src/presets/muscle.rs — skeletal muscle tissue preset.
//
// Reference: Author et al. (year), Journal. DOI: ...
// Model: ...
// Parameters from: ...

use std::collections::HashMap;
use orfas_core::{MaterialLaw, CompressibleMaterial, NeoHookeanIso, VolumetricLnJ};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

pub struct SkeletalMuscle {
    metadata: TissueMetadata,
}

impl Default for SkeletalMuscle {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(1000.0, 5000.0));
        ci.insert("kappa", ConfidenceInterval::new(20000.0, 100000.0));

        SkeletalMuscle {
            metadata: TissueMetadata {
                name:                 "Skeletal muscle",
                model:                "Neo-Hookean",
                reference:            "Author et al. (year), Journal, DOI:...",
                doi:                  Some("10.xxxx/..."),
                protocol:             "description of experimental protocol",
                confidence_intervals: ci,
                notes:                "any relevant notes",
            },
        }
    }
}

impl TissuePreset for SkeletalMuscle {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }

    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 2500.0 },      // nominal value from paper
            vol:     VolumetricLnJ { kappa: 50000.0 },  // nominal value from paper
            density: 1060.0,
        })
    }
}
```

### Step 2 — Register the preset

In `orfas-tissues/src/presets/mod.rs`:

```rust
mod muscle;
pub use muscle::SkeletalMuscle;
```

And add to `all_presets()`:
```rust
Box::new(SkeletalMuscle::default()),
```

### Step 3 — Add tests

In `orfas-tissues/src/tests/presets.rs`, add:

```rust
#[test]
fn test_muscle_thermodynamic_consistency() {
    let preset = SkeletalMuscle::default();
    let mu     = 2500.0_f64;
    let kappa  = 50000.0_f64;
    let lambda = kappa - 2.0 / 3.0 * mu;
    assert_thermodynamic_consistency(&preset, Some((lambda, mu)));
}

#[test]
fn test_muscle_confidence_intervals() {
    let preset = SkeletalMuscle::default();
    assert_within_confidence(&preset, &[("mu", 2500.0), ("kappa", 50000.0)]);
}

#[test]
fn test_muscle_metadata() {
    assert_metadata_non_empty(&SkeletalMuscle::default());
}
```

Update `test_all_presets_count` to reflect the new total.

### Step 4 — Add JSON support (optional)

If the preset uses a model already supported by `loader.rs` (neo_hookean, mooney_rivlin,
holzapfel_ogden, saint_venant_kirchhoff), no changes are needed — it can be loaded from JSON
immediately. For a new model, add a new match arm in `build_material` in `loader.rs`.

### Notes on parameter values

- Use the nominal values from the paper directly — do not use the midpoint of CI.
- The `material()` method must return parameters consistent with `confidence_intervals`.
  `assert_within_confidence` in the test suite will catch any mismatch.
- For anisotropic materials (HGO), fiber directions are not encoded in the preset —
  they are provided at runtime via `SimulationContext`. Document this in the `notes` field.
- For near-incompressible tissues, `kappa/mu` should typically be >= 50.
  Lower ratios may cause Newton divergence for large deformations.

---

## How to add a solver

Implement `NonlinearSolver` or `NonlinearSparseSolver` in `orfas-core/src/solver.rs`:

```rust
pub struct MySolver { ... }

impl NonlinearSolver for MySolver {
    fn solve<B: BMatrix>(
        &self,
        assembler:  &Assembler,
        mesh:       &Mesh,
        material:   &dyn MaterialLaw,
        bc_result:  &BoundaryConditionResult,
        lin_solver: &dyn DenseSolver,
        sim_ctx:    &SimulationContext,
    ) -> Result<DVector<f64>, SolverError> { ... }
}
```

Add it to `SolverChoice` in `orfas-viewer/src/state.rs` and expose it in `run_simulation_static`.

---

## Running tests

```bash
# All tests
cargo test

# Specific crate
cargo test -p orfas-core
cargo test -p orfas-tissues

# Specific test
cargo test -p orfas-tissues test_liver_thermodynamic_consistency

# With output
cargo test -- --nocapture
```

---

## Validation

Numerical validation results (test suites, benchmarks, convergence curves) are documented
in [VALIDATION.md](VALIDATION.md).