# DEVELOPER.md

Technical reference for contributors. Covers codebase architecture, conventions,
and step-by-step guides for the most common contribution tasks.

---

## Workspace structure

```
orfas/
├── orfas-core/        # Core library — mesh, elements, materials, assembly, solvers, integrators
├── orfas-io/          # Mesh I/O — .vtk reader (Tet4 only; Tet10 planned for v0.9.0)
├── orfas-tissues/     # Tissue preset library — calibrated materials from literature
└── orfas-viewer/      # Interactive egui viewer
```

### orfas-core

```
orfas-core/src/
├── lib.rs                      # Public re-exports
├── mesh.rs                     # Node, Mesh<N>, Tet4Mesh, Tet10Mesh, Hex8Mesh, FemMesh trait
├── mechanical_state.rs         # MechanicalState (dynamic DOF vector)
├── boundary.rs                 # BoundaryConditions, Constraint, Load, FixedNode,
│                               # PenaltyMethod, EliminationMethod (supports non-zero Dirichlet)
├── damping.rs                  # RayleighDamping
├── integrator.rs               # ImplicitEulerIntegrator
├── sparse.rs                   # CsrMatrix, CgSolver, ILU(0)
├── solver.rs                   # DenseSolver, DirectSolver, NewtonRaphson, NewtonRaphsonCachedK
├── element/
│   ├── mod.rs                  # Public re-exports
│   ├── traits.rs               # FiniteElement, ElementGeometry, FemMesh, DofType (Vec3Dof, Vec6Dof)
│   ├── tet4.rs                 # Tet4, Tet4Geometry — linear tetrahedron
│   ├── tet10.rs                # Tet10, Tet10Geometry — quadratic tetrahedron (basix convention)
│   ├── hex8.rs                 # Hex8, Hex8Geometry — trilinear hexahedron, 2×2×2 Gauss
│   ├── beam2.rs                # Beam2, Beam2Geometry — Euler-Bernoulli 3D beam, analytical stiffness
│   ├── shell4.rs               # Shell4, Shell4Geometry — MITC4 shell (Bathe & Dvorkin 1986)
│   ├── subdivision.rs          # tet4_to_tet10() — temporary, replaced by VTK loader in v0.9.0
│   └── tests.rs                # Generic element tests (unit checks + bending benchmark)
├── assembler/
│   ├── mod.rs                  # Assembler<E: FiniteElement>, public API
│   ├── geometry.rs             # (kept for reference, no longer used — logic moved to element/)
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
    ├── linear_elastic.rs       # LinearElastic — small-strain, for patch tests and validation
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

### Generic element abstraction

`Assembler<E: FiniteElement>` is generic over element type. All element-specific logic
(shape functions, gradients, Gauss points, B matrix, jacobian) lives in the `FiniteElement`
implementation. The assembler knows nothing about the internals of the element.

`FemMesh` is a trait over `Mesh<N>` that provides a uniform interface (`nodes()`,
`connectivity()`, `n_nodes()`, `n_elements()`) regardless of the number of nodes per element.
Use `Tet4Mesh = Mesh<4>`, `Tet10Mesh = Mesh<10>`, `Hex8Mesh = Mesh<8>` at call sites.

### DofType trait

`DofType` encodes the number of DOFs per node at compile time — following the same pattern
as SOFA's `DataTypes` template. Two implementations:

- `Vec3Dof` (N_DOF = 3) — all volumetric elements: `Tet4`, `Tet10`, `Hex8`
- `Vec6Dof` (N_DOF = 6) — structural elements: `Beam2`, `Shell4`

Every downstream type (`MechanicalState`, `Assembler`, `BoundaryConditions`) uses
`E::Dof::N_DOF` for index arithmetic. Never hardcode `3 *` for DOF offsets — always use
`ndof` derived from `<E::Dof as DofType>::N_DOF`.

### Analytical stiffness path

Structural elements (`Beam2`, `Shell4`) bypass the Gauss-point loop entirely by overriding
`element_stiffness()` in `FiniteElement`. The assembler detects `Some(ke)` and scatters
the matrix directly into `K`, skipping `b_matrix`, `integration_points`, and `gauss_det_j`.

Volumetric elements return `None` (default) and follow the standard `B^T C B` path.
Do not mix the two paths for the same element type.

### Geometry cache

`Assembler::<E>::new(&mesh)` pre-computes all element geometry (`E::Geometry`) once.
For Tet4: analytical gradients and volume. For Tet10/Hex8: physical-space gradients and
jacobian determinants at each Gauss point. For Beam2: element length. For Shell4: director
vectors and rotation bases at each node. The mesh is assumed immutable during simulation.

### Reduced vs full DOF space

`MechanicalState` lives in the reduced DOF space (free DOFs only). The assembler works in
the full space. `bc_result.reconstruct` maps from reduced to full, including prescribed
non-zero Dirichlet values. This separation is fundamental — respect it in all new code.

### Generic hot paths

`Assembler<E>` is generic over `E: FiniteElement` — monomorphized at compile time, zero cost.
Less critical paths (`make_material`) return `Box<dyn MaterialLaw>`.

---

## Code conventions

- All source files must start with `// UTF-8`
- All code and comments in English — no French
- Every public item must have a `///` doc comment
- Run `cargo fmt` and `cargo clippy` before submitting
- Tests live in `tests/` subdirectories or dedicated test modules, not inline in source
- Use `debug_assert!` for invariants that should never fire in release builds

---

## How to add a new element type

Implement `FiniteElement` in a new file under `orfas-core/src/element/`:

### Step 1 — Define the geometry struct

```rust
// orfas-core/src/element/hex8.rs
use crate::element::traits::ElementGeometry;

pub struct Hex8Geometry {
    // precomputed data — whatever the element needs
    pub grad: Vec<nalgebra::DMatrix<f64>>, // gradients at Gauss points
    pub det_j: Vec<f64>,
}

impl ElementGeometry for Hex8Geometry {}
```

### Step 2 — Implement FiniteElement

```rust
use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use crate::element::traits::{ElementGeometry, FiniteElement, Vec3Dof};

pub struct MyElement;

impl FiniteElement for MyElement {
    type Geometry = MyElementGeometry;
    type Dof      = Vec3Dof; // Vec3Dof for volumetric, Vec6Dof for structural

    const N_NODES: usize = 8;

    fn precompute(nodes: &[Vector3<f64>]) -> MyElementGeometry { ... }
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64> { ... }
    fn shape_gradients(geo: &MyElementGeometry, gauss_index: usize) -> DMatrix<f64> { ... }
    fn integration_points() -> Vec<(Vector3<f64>, f64)> { ... }
    fn b_matrix(grad_n: &DMatrix<f64>, f: &Matrix3<f64>) -> DMatrix<f64> { ... }
    fn element_volume(geo: &MyElementGeometry) -> f64 { ... }
    fn gauss_det_j(geo: &MyElementGeometry, gauss_index: usize) -> f64 { ... }

    // Override for analytical stiffness (Beam2, Shell4 style) — optional:
    // fn element_stiffness(nodes: &[Vector3<f64>], material: &dyn MaterialLaw)
    //     -> Option<DMatrix<f64>> { Some(ke) }
}
```

### Step 3 — Implement FemMesh for Mesh<8>

In `mesh.rs`, add:

```rust
impl FemMesh for Mesh<8> {
    fn nodes(&self) -> &[Node] { &self.nodes }
    fn connectivity(&self) -> Vec<Vec<usize>> {
        self.elements.iter().map(|e| e.to_vec()).collect()
    }
}
pub type Hex8Mesh = Mesh<8>;
```

### Step 4 — Register in element/mod.rs

```rust
pub mod hex8;
pub use hex8::{Hex8, Hex8Geometry};
```

### Step 5 — Add unit tests

In `element/tests.rs`, add `tet4_partition_of_unity`-style tests for Hex8:
- Partition of unity
- Gradient sum zero
- Volume correctness
- Bending benchmark

### Quadrature convention

`gauss_det_j` must be consistent with `integration_points` weights so that:
```
sum_g gauss_det_j(g) * w(g) == element_volume
```
For elements using barycentric coordinates (Tet4, Tet10), divide `det_J` by the
volume of the reference simplex (6 for tetrahedra). See `tet4.rs` and `tet10.rs`
for the exact convention used.

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
    tau_iso:    vec![1.0],
    beta_iso:   vec![0.3],
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

---

## How to add a tissue preset

### Step 1 — Create the preset file

Add a new file in `orfas-tissues/src/presets/`, e.g. `muscle.rs`:

```rust
// UTF-8
// orfas-tissues/src/presets/muscle.rs — skeletal muscle tissue preset.

use std::collections::HashMap;
use orfas_core::{MaterialLaw, CompressibleMaterial, NeoHookeanIso, VolumetricLnJ};
use crate::metadata::{ConfidenceInterval, TissueMetadata, TissuePreset};

pub struct SkeletalMuscle { metadata: TissueMetadata }

impl Default for SkeletalMuscle {
    fn default() -> Self {
        let mut ci = HashMap::new();
        ci.insert("mu",    ConfidenceInterval::new(1000.0, 5000.0));
        ci.insert("kappa", ConfidenceInterval::new(20000.0, 100000.0));
        SkeletalMuscle {
            metadata: TissueMetadata {
                name: "Skeletal muscle", model: "Neo-Hookean",
                reference: "Author et al. (year), Journal, DOI:...",
                doi: Some("10.xxxx/..."),
                protocol: "description of experimental protocol",
                confidence_intervals: ci,
                notes: "any relevant notes",
            },
        }
    }
}

impl TissuePreset for SkeletalMuscle {
    fn metadata(&self) -> &TissueMetadata { &self.metadata }
    fn material(&self) -> Box<dyn MaterialLaw> {
        Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: 2500.0 },
            vol:     VolumetricLnJ { kappa: 50000.0 },
            density: 1060.0,
        })
    }
}
```

### Step 2 — Register, test, and add JSON support

See the full guide in the v0.7.2 section of the git history — steps are unchanged from v0.7.2.

---

## How to add a solver

Implement `NonlinearSolver` in `orfas-core/src/solver.rs`.
Solvers no longer take a `BMatrix` generic — the element type is carried by `Assembler<E>`:

```rust
pub struct MySolver { ... }

impl NonlinearSolver for MySolver {
    fn solve(
        &self,
        assembler:  &Assembler<impl FiniteElement>,
        mesh:       &impl FemMesh,
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
cargo test -p orfas-core element::tests::bending -- --nocapture

# With output
cargo test -- --nocapture
```

---

## Known limitations (v0.8.1)

### Tet10 on structured meshes

Tet10 shows ~20% error on the cantilever bending benchmark when using the structured mesh
generated by `Mesh::generate`. This is a known behavior of structured tetrahedral meshes
with alternating element orientation — not a bug in the Tet10 implementation.

Tet10 gradients have been validated against basix (FEniCS) to machine precision.
Element energy is consistent with Tet4 on the same geometry.

Full validation on unstructured meshes is planned for v0.9.0 when the VTK Tet10 loader
is available. For production use of Tet10, import meshes from Gmsh or other mesh generators.

### Shell4 — physics validation pending

`Shell4` (MITC4) passes structural tests (symmetry, drilling stabilisation, partition of unity,
director construction). Full physics benchmarks (Morley skew plate, pinched cylinder,
Scordelis-Lo shell) are planned for v0.9.0. The B-matrix formulation for the bending
block uses a linearised approximation suitable for small-strain analysis; geometric
nonlinearity for shells is deferred to a future version.

### Patch test and large-deformation formulation

The classical FEM patch test (exact reproduction of linear displacement fields) is not
directly applicable to our large-deformation assembler formulation. The assembler uses
the full Lagrangian kinematics (F, PK2, PK1), which introduces nonlinear terms even for
linear displacement fields. The patch test passes in the infinitesimal strain limit
(u_scale → 0) as expected.

---

## Validation

Numerical validation results (test suites, benchmarks, convergence curves) are documented
in [VALIDATION.md](VALIDATION.md).