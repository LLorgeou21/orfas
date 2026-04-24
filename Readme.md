# ORFAS
### Open Rust Framework Architecture for Simulation

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Status: Pre-Alpha](https://img.shields.io/badge/Status-Pre--Alpha-red.svg)]()
[![Rust](https://img.shields.io/badge/Rust-2021%20Edition-orange.svg)](https://www.rust-lang.org/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)]()

> *A generic, extensible Finite Element Method framework written in Rust, with a primary focus on medical and biomechanical simulation.*

---
![ORFAS dynamic simulation](simulation.gif)
---

## Abstract

ORFAS is an open-source framework for physics-based simulation using the Finite Element Method (FEM),
implemented in Rust. It provides a modular, composable architecture designed to support a wide range
of material laws, element types, and numerical solvers. While the framework is domain-agnostic, its
primary design motivation is soft tissue simulation for medical and surgical applications, in the
tradition of projects such as SOFA.

ORFAS is built on three foundational principles: **genericity** (material laws, solvers, and element
types are interchangeable via traits), **extensibility** (each component can be replaced or augmented
without modifying the core), and **performance** (Rust's ownership model and zero-cost abstractions
enable safe, high-performance numerical computation without a garbage collector).

This project is in active early development. The current milestone targets a dynamic 3D FEM solver
with nonlinear hyperelastic materials on tetrahedral meshes, fiber-reinforced anisotropic models,
viscoelastic dissipation via Prony series, and implicit Euler time integration with Newton-Raphson.
Contributions, feedback, and collaborations are warmly welcomed.

---

## Motivation

### Why a new simulation framework?

Existing FEM frameworks for medical simulation are predominantly written in C++ (SOFA, FEBio, deal.II).
While mature and capable, these frameworks carry significant adoption friction: complex CMake build
systems, heavy dependency chains, and architectures designed before modern concurrency and memory
safety were first-class concerns.

Rust offers a compelling alternative: memory safety without a garbage collector, a world-class
package manager (Cargo), and a trait system that enables generic programming with compile-time
guarantees. These properties are particularly valuable in simulation code, where subtle memory
bugs and performance regressions are costly.

### Why not extend SOFA or Fenris?

[SOFA](https://www.sofa-framework.org/) is the reference framework for medical simulation and a
direct source of inspiration for ORFAS. However, its C++ codebase and plugin architecture present
a high integration barrier for projects seeking a modern, memory-safe foundation.

[Fenris](https://github.com/InteractiveComputerGraphics/fenris) is the most mature FEM library in
the Rust ecosystem and an excellent piece of work. Its focus is solid mechanics for computer graphics,
and it is no longer actively maintained as of 2022. ORFAS targets a different domain (medical
simulation) and a different design philosophy (full replaceability of every numerical component
via traits).

### Long-term vision

ORFAS aims to become a reference implementation for FEM-based medical simulation in Rust, providing:
- Native Rust performance and safety guarantees
- Python bindings for the scientific computing community (via PyO3)
- WebAssembly support for browser-based simulation and visualization
- A path toward real-time surgical simulation and haptic feedback

---

## Architecture

ORFAS is organized as a Cargo workspace with three crates:

```
orfas/
├── orfas-core/     # Core library: mesh, materials, assembly, solvers, integrators
├── orfas-io/       # Mesh I/O: .vtk file format
└── orfas-viewer/   # Interactive viewer (egui)
```

### orfas-viewer structure

The viewer is organized as a set of focused modules:

```
orfas-viewer/src/
├── main.rs        # Entry point only
├── app.rs         # MyEguiApp, impl eframe::App, all UI logic
├── state.rs       # AppState, Camera, enums, make_material
├── simulation.rs  # build_simulation, run_simulation_static, init_simulation_dynamic
└── render.rs      # project, screen_to_node, depth_sorted_nodes
```

### Core abstractions (`orfas-core`)

The framework is built around several central traits:

**`MaterialLaw`** — Defines the constitutive relationship between strain and stress in the Lagrangian
frame. Any hyperelastic material model implements this trait via four methods: `pk2_stress(F, ctx)`,
`tangent_stiffness(F, ctx)`, `strain_energy(F, ctx)`, and `update_state(F, ctx)`. All methods take
a `MaterialContext` grouping fiber directions, time step, and internal variables. Swapping material
laws requires no changes to the assembly or solver.

**`IsochoricPart`** — Defines the isochoric (volume-preserving) contribution of a hyperelastic
material. Implemented by `NeoHookeanIso`, `MooneyRivlinIso`, and `OgdenIso`. Each provides
`strain_energy_iso(F)`, `pk2_stress_iso(F)`, and `tangent_iso(F)`.

**`VolumetricPart`** — Defines the volumetric (pressure) contribution of a hyperelastic material.
Implemented by `VolumetricLnJ` (W = κ/2·(ln J)²) and `VolumetricQuad` (W = κ/2·(J−1)²).
Each provides `dw_vol(J)` and `d2w_vol(J)`.

**`AnisotropicPart`** — Defines the anisotropic fiber contribution of a hyperelastic material.
Implemented by `HolzapfelOgden` and `NoAnisotropy` (zero contribution, zero runtime cost).
Each provides `strain_energy_aniso(F, ctx)`, `pk2_stress_aniso(F, ctx)`, and `tangent_aniso(F, ctx)`.
Fiber directions are passed via `MaterialContext` — stored globally in `FiberField` and borrowed
per element at assembly time.

**`CompressibleMaterial<I, V>`** — Generic struct composing any `IsochoricPart` and `VolumetricPart`
into a full `MaterialLaw`. The isochoric/volumetric split follows the Flory decomposition and enables
the same volumetric penalty to be reused across all isochoric models.

**`CompressibleAnisotropicMaterial<I, A, V>`** — Generic struct composing any `IsochoricPart`,
`AnisotropicPart`, and `VolumetricPart` into a full `MaterialLaw`. Natural extension of
`CompressibleMaterial` for fiber-reinforced biological tissues.

**`ViscoelasticMaterial<I, A, V>`** — Generic viscoelastic wrapper using a Prony series on the
isochoric and anisotropic contributions. The volumetric part remains purely elastic. Implements
the algorithmic update from Holzapfel & Gasser (2001), Box 1. Internal variables (Prony tensors
Q_α, previous stresses, precomputed ΣQ) are stored per element in `InternalVariables`.
`pk2_stress` reads ΣQ in O(1); `update_state` performs the full Prony update once per time step
after Newton convergence — following the SOFA/FEBio pattern.

**`MaterialContext<'a>`** — Per-element evaluation context passed to all material law methods.
Lightweight stack-allocated struct borrowing from `SimulationContext`. Contains: time step `dt`,
fiber directions `fiber_dirs` (borrowed from `FiberField`), read-only internal variables `iv_ref`
(used by `pk2_stress`), and mutable internal variables `iv` (used by `update_state`).

**`SimulationContext`** — Global owned context for a simulation step. Holds `FiberField`, `dt`,
and `Option<InternalVariables>`. The assembler constructs a `MaterialContext` per element via
`material_context_for` (read-only) or `material_context_for_mut` (mutable). Extensible — future
fields (temperature, damage) are added here without changing assembler signatures.

**`FiberField`** — Per-element fiber direction storage. Constructed once (`uniform`, `helix`,
`helix_two_families`) and stored in `SimulationContext`. The assembler borrows slices per element
at assembly time — zero allocation per call.

**`InternalVariables`** / **`ElementInternalVars`** — Per-element internal variable storage for
viscoelastic materials. Layout per element: `[Q_iso × m_iso][Q_aniso × m_aniso][S_iso_prev][S_aniso_prev][sum_Q]`,
total `6*(m_iso + m_aniso + 3)` scalars. `sum_Q` is precomputed by `update_state` and read in O(1)
by `pk2_stress` at every Newton iteration.

**`BoundaryConditionMethod`** — Defines how Dirichlet boundary conditions are applied to the system.
Current implementations: penalty method and elimination method.

**`DenseSolver`** — Solves the assembled dense linear system `K·u = f`.
Current implementation: `DirectSolver` — direct LU decomposition via nalgebra.

**`SparseSolver`** — Solves the assembled sparse linear system `K·u = f` (via `CsrMatrix`).
Current implementation: `CgSolver` — preconditioned conjugate gradient with `Identity` or `ILU(0)` preconditioner.

**`NonlinearSolver`** — Solves the nonlinear static problem `R(u) = f_int(u) - f_ext = 0` using dense assembly.
Current implementations: `NewtonRaphson` (reassembles K at each iteration) and
`NewtonRaphsonCachedK` (factorizes K once, for materials with constant tangent stiffness such as SVK).

**`NonlinearSparseSolver`** — Solves the same nonlinear static problem using sparse assembly (`CsrMatrix`).
Current implementation: `NewtonRaphsonSparse` — preferred for large meshes (> 1000 nodes).

**`DampingModel`** — Defines how the damping matrix `C` is computed from the mass and stiffness matrices.
Current implementation: Rayleigh damping `C = α·M + β·K`.

**`IntegratorMethod`** — Defines a time integration scheme advancing the `MechanicalState` by one step `dt`.
Current implementation: implicit Euler with internal Newton-Raphson loop, supporting nonlinear materials.

**`MechanicalState`** — Holds the dynamic state of the simulated object: position, velocity, and acceleration vectors.
Exposes vector operations (`v_op`, `add_mv`) inspired by SOFA's `MechanicalState` abstraction.

### Material library

ORFAS implements a growing library of hyperelastic material laws, all based on the isochoric/volumetric
split architecture. Every isochoric model is independently composable with any volumetric penalty,
and any combination can be wrapped in `ViscoelasticMaterial` for time-dependent dissipation.

| Model | Type | Energy function |
|-------|------|-----------------|
| `SaintVenantKirchhoff` | Coupled | W = λ/2·(tr E)² + μ·tr(E²) |
| `NeoHookeanIso` | Isochoric | W = μ/2·(Ī₁ − 3) |
| `MooneyRivlinIso` | Isochoric | W = c₁·(Ī₁ − 3) + c₂·(Ī₂ − 3) |
| `OgdenIso` | Isochoric | W = Σᵢ μᵢ/αᵢ·(λ̄₁^αᵢ + λ̄₂^αᵢ + λ̄₃^αᵢ − 3) |
| `HolzapfelOgden` | Anisotropic | W = k₁/(2k₂)·Σᵢ[exp(k₂·(Ī₄ᵢ−1)²)−1] |
| `NoAnisotropy` | Anisotropic | W = 0 (zero-cost null implementation) |
| `VolumetricLnJ` | Volumetric | U = κ/2·(ln J)² |
| `VolumetricQuad` | Volumetric | U = κ/2·(J − 1)² |
| `ViscoelasticMaterial<I,A,V>` | Viscoelastic wrapper | S = S_eq + ΣQ_α (Prony series) |

The `OgdenIso` tangent stiffness is derived from the spectral decomposition of C and implemented
following Connolly et al. (2019), with L'Hôpital's rule applied when principal stretches are
numerically equal or similar (tolerance 10⁻⁶). The `NeoHookeanIso`, `MooneyRivlinIso`, and
`HolzapfelOgden` tangent stiffness tensors are derived analytically following Cheng & Zhang (2018).
The `ViscoelasticMaterial` algorithmic tangent follows Holzapfel & Gasser (2001), Box 1.

---

## How to add a material

ORFAS supports several ways to implement a new hyperelastic material law.

### Option A — Isochoric model (recommended)

If your model admits an isochoric/volumetric split, implement the `IsochoricPart` trait and
compose it with an existing `VolumetricPart` via `CompressibleMaterial`.

```rust
use nalgebra::{Matrix3, Matrix6};
use orfas_core::material::{IsochoricPart, CompressibleMaterial, VolumetricLnJ};

pub struct MyIso {
    pub mu: f64,
}

impl IsochoricPart for MyIso {
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64 { todo!() }
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64> { todo!() }
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64> { todo!() }
}
```

Then compose it with a volumetric penalty:

```rust
let mat = CompressibleMaterial {
    iso:     MyIso { mu: 1000.0 },
    vol:     VolumetricLnJ { kappa: 10000.0 },
    density: 1000.0,
};
```

### Option B — Anisotropic model

If your model has fiber reinforcement, implement `AnisotropicPart` and compose it via
`CompressibleAnisotropicMaterial`. Fiber directions are passed via `MaterialContext` — no
changes to the assembler are needed.

```rust
use orfas_core::material::{
    AnisotropicPart, CompressibleAnisotropicMaterial,
    NeoHookeanIso, VolumetricLnJ, MaterialContext,
};

pub struct MyFibers { pub k1: f64, pub k2: f64 }

impl AnisotropicPart for MyFibers {
    fn strain_energy_aniso(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> f64 { todo!() }
    fn pk2_stress_aniso(&self, f: &Matrix3<f64>, ctx: &mut MaterialContext) -> Matrix3<f64> { todo!() }
    fn tangent_aniso(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> Matrix6<f64> { todo!() }
}

let mat = CompressibleAnisotropicMaterial {
    iso:     NeoHookeanIso { mu: 500.0 },
    aniso:   MyFibers { k1: 2363.0, k2: 0.84 },
    vol:     VolumetricLnJ { kappa: 10000.0 },
    density: 1000.0,
};
```

### Option C — Viscoelastic model

Any `CompressibleMaterial` or `CompressibleAnisotropicMaterial` can be wrapped in
`ViscoelasticMaterial` to add Prony series dissipation. Use `NoAnisotropy` for isotropic
viscoelastic materials.

```rust
use orfas_core::material::{
    ViscoelasticMaterial, NeoHookeanIso, NoAnisotropy, VolumetricLnJ,
    InternalVariables, SimulationContext,
};

let mat = ViscoelasticMaterial {
    iso:        NeoHookeanIso { mu: 500.0 },
    aniso:      NoAnisotropy,
    vol:        VolumetricLnJ { kappa: 10000.0 },
    density:    1000.0,
    tau_iso:    vec![1.0, 10.0],   // two relaxation times
    beta_iso:   vec![0.3, 0.1],    // Prony factors
    tau_aniso:  vec![],
    beta_aniso: vec![],
};

// Initialize internal variables and simulation context
let n_elements = mesh.elements.len();
let iv = InternalVariables::new(n_elements, mat.m_iso(), mat.m_aniso());
let mut sim_ctx = SimulationContext::isotropic_viscoelastic(n_elements, dt, iv);

// After Newton convergence at each time step:
assembler.update_internal_variables(&mesh, &mat, &u, &mut sim_ctx);
```

### Option D — Fully coupled model

If your model does not admit an isochoric/volumetric split, implement `MaterialLaw` directly:

```rust
use nalgebra::{Matrix3, Matrix6};
use orfas_core::material::{MaterialLaw, MaterialContext};

pub struct MyMaterial { pub param: f64, pub density: f64 }

impl MaterialLaw for MyMaterial {
    fn density(&self) -> f64 { self.density }
    fn strain_energy(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> f64 { todo!() }
    fn pk2_stress(&self, f: &Matrix3<f64>, _ctx: &mut MaterialContext) -> Matrix3<f64> { todo!() }
    fn tangent_stiffness(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> Matrix6<f64> { todo!() }
}
```

### Validation

Validate your implementation using the standard test suite:

```rust
run_standard_material_tests(&mat, lambda_eff, mu_eff, &mut MaterialContext::default());
```

This checks that at `F = I`: `W = 0`, `S = 0`, `C = C_Hooke`, and that the analytical tangent
matches central finite differences at moderate deformation.

---

## Changelog

### v0.7.1

#### Materials
- **`HolzapfelOgden`** — anisotropic isochoric fiber model (Holzapfel-Gasser-Ogden); strain energy
  W = k₁/(2k₂)·Σᵢ[exp(k₂·(Ī₄ᵢ−1)²)−1]; fibers only contribute under tension (Ī₄ᵢ > 1);
  analytical PK2 stress and tangent derived following Cheng & Zhang (2018) eqs. (49), (56)
- **`NoAnisotropy`** — null implementation of `AnisotropicPart`; returns zero for all contributions;
  zero runtime cost (compiler inlines and eliminates all additions)
- **`CompressibleAnisotropicMaterial<I, A, V>`** — generic struct composing `IsochoricPart`,
  `AnisotropicPart`, and `VolumetricPart` into a full `MaterialLaw`; natural extension of
  `CompressibleMaterial` for fiber-reinforced tissues
- **`ViscoelasticMaterial<I, A, V>`** — generic viscoelastic wrapper using Prony series on the
  isochoric and anisotropic contributions; volumetric part remains purely elastic; implements the
  algorithmic update from Holzapfel & Gasser (2001) Box 1; composable with any `IsochoricPart`,
  `AnisotropicPart`, and `VolumetricPart`
  - `pk2_stress` — reads precomputed `sum_Q` from `iv_ref` in O(1); read-only on internal variables
  - `update_state` — computes Q_α_n+1 via Prony recurrence, writes Q_α, S_prev, sum_Q; called once
    per time step after Newton convergence (SOFA/FEBio pattern)
  - `tangent_stiffness` — algorithmic tangent C = C_vol + (1+δ_iso)·C_iso + (1+δ_aniso)·C_aniso
    where δ_a = Σ_α β_αa·exp(−Δt/2τ_αa)

#### Architecture
- **`MaterialContext<'a>`** — new `iv_ref: Option<&'a ElementInternalVars>` field for read-only
  access to internal variables in `pk2_stress`; `iv: Option<&'a mut ElementInternalVars>` retained
  for write access in `update_state`; separates read and write paths cleanly
- **`SimulationContext`** — added `iv: Option<InternalVariables>` field; new constructors
  `viscoelastic` and `isotropic_viscoelastic`; `material_context_for` now populates `iv_ref`;
  `material_context_for_mut` populates `iv` for `update_internal_variables`
- **`FiberField`** — fiber direction storage decoupled from `MaterialContext`; stored once in
  `SimulationContext`, borrowed per element at assembly time via `directions_for(elem_idx)`
- **`InternalVariables`** / **`ElementInternalVars`** — new module `material/internal_variables.rs`;
  flat `DVector<f64>` storage per element with layout `[Q_iso×m_iso][Q_aniso×m_aniso][S_iso_prev][S_aniso_prev][sum_Q]`,
  total `6*(m_iso + m_aniso + 3)` scalars; `sum_Q` precomputed for O(1) reads

#### Assembler
- **`update_internal_variables`** — new method on `Assembler`; calls `material.update_state` for
  each element after Newton convergence; no-op for elastic materials (`iv.is_none()`); takes
  `&mut SimulationContext` — the only assembler method that mutates the context
- **`assemble_internal_forces`** — now passes `iv_ref` via `material_context_for`; reads `sum_Q`
  correctly for viscoelastic materials during Newton iterations without corrupting internal state
- **Module split** — `assembler.rs` split into `assembler/mod.rs`, `assembler/geometry.rs`,
  `assembler/pattern.rs`, `assembler/assembly.rs`, `assembler/tests.rs`

#### Tests
- **`material/tests/`** — `tests.rs` split into `mod.rs`, `helpers.rs`, `elastic.rs`,
  `anisotropic.rs`, `viscoelastic.rs`
- **`run_viscoelastic_tests`** — parametric suite: elastic fallback, zero stress at F=I,
  relaxation convergence, algorithmic tangent scaling verification
- **`test_viscoelastic_relaxation`** — integration test in `lib.rs`; full pipeline with assembler,
  `update_internal_variables`, and reaction force convergence to elastic equilibrium

#### References
- Holzapfel, G.A., & Gasser, T.C. (2001). A viscoelastic model for fiber-reinforced composites at
  finite strains. *Computer Methods in Applied Mechanics and Engineering*, 190, 4379–4403.

### v0.7.0
- **Isochoric/volumetric split architecture** — new `IsochoricPart` and `VolumetricPart` traits replacing the monolithic `NeoHookean` struct; `CompressibleMaterial<I, V>` composes any pair into a full `MaterialLaw`
- **`NeoHookeanIso`** — isochoric Neo-Hookean: W = μ/2·(Ī₁ − 3); analytical tangent derived following Cheng & Zhang (2018) eq. (39) using modified right Cauchy-Green tensor C̄⁻¹ = J^{2/3}·C⁻¹
- **`MooneyRivlinIso`** — isochoric Mooney-Rivlin: W = c₁·(Ī₁ − 3) + c₂·(Ī₂ − 3); analytical tangent derived following Cheng & Zhang (2018) eq. (25)
- **`OgdenIso`** — isochoric Ogden (N-term): W = Σᵢ μᵢ/αᵢ·(λ̄₁^αᵢ + λ̄₂^αᵢ + λ̄₃^αᵢ − 3); tangent from spectral decomposition of C via SymmetricEigen; L'Hôpital's rule for equal/similar eigenvalues; following Connolly et al. (2019) eqs. (11), (12), (22), (35)
- **`VolumetricLnJ`** — volumetric penalty U = κ/2·(ln J)²; analytical tangent C_vol = κ·(C⁻¹ ⊗ C⁻¹) − 2κ·ln(J)·(C⁻¹ ⊙ C⁻¹)
- **`VolumetricQuad`** — volumetric penalty U = κ/2·(J − 1)²
- **Breaking change** — `NeoHookean { youngs_modulus, poisson_ratio }` removed; replaced by `CompressibleMaterial { iso: NeoHookeanIso { mu }, vol: VolumetricLnJ { kappa }, density }`
- **Viewer** — `MooneyRivlin` and `Ogden` material choices added; material parameters derived from Young's modulus and Poisson's ratio for consistent small-strain behavior

### v0.6.1
- **`SparseAssemblyStrategy` trait** — generic assembly strategy for `NewtonRaphsonSparse`; two implementations: `Sequential` (existing) and `Parallel` (new)
- **`NewtonRaphsonSparse<S: SparseAssemblyStrategy>`** — refactored to be generic over the assembly strategy; zero code duplication between sequential and parallel variants
- **`assemble_tangent_sparse_parallel`** — parallel sparse assembly via rayon + atomic f64 additions (`fetch_update`); ~10x speedup on meshes > 8000 nodes (measured on 20-core machine)
- **Pre-built CSR pattern** in `Assembler` — sparsity pattern computed once at `Assembler::new` via `build_csr_pattern`; reused at every assembly call; eliminates COO→CSR sorting overhead
- **`entry_map`** — `HashMap<(i,j), idx>` cached at `Assembler::new`; enables O(1) direct write into CSR values array per element
- **`build_element_colors`** — greedy graph coloring of elements (available, unused); stored as `Option<Vec<Vec<usize>>>` in `Assembler` for future use
- **Viewer** — `SolverChoice::NewtonSparseParallel` added; selectable as "Newton (sparse CG parallel)" in the UI

### v0.6.0
- **Sparse solvers** (`orfas-core/src/sparse.rs`) — new module providing sparse linear and nonlinear solvers
- **`SparseSolver` trait** — mirror of `DenseSolver` operating on `CsrMatrix<f64>` (nalgebra-sparse)
- **`CgSolver`** — preconditioned conjugate gradient solver; supports `Preconditioner::Identity` and `Preconditioner::Ilu(k)`
- **`ILU(0)` preconditioner** — incomplete LU factorization with zero fill-in
- **`assemble_tangent_sparse`** — new sparse assembly method on `Assembler`
- **`NonlinearSparseSolver` trait** and **`NewtonRaphsonSparse`** — sparse equivalent of `NewtonRaphson`
- **`DenseSolver` rename** — `Solver` trait renamed to `DenseSolver` for clarity

### v0.5
- **Neo-Hookean hyperelastic material** (`NeoHookean`) — compressible coupled formulation
- **`NewtonRaphsonCachedK`** — Newton-Raphson variant that factorizes `K_tangent` once at `u=0`
- **Viewer refactored** into 5 focused modules
- **Improved camera** — orbit, pan, zoom, depth-sorted node rendering

### v0.4
- Nonlinear hyperelastic material: `SaintVenantKirchhoff`
- Fully refactored `MaterialLaw` trait
- `NonlinearSolver` trait and `NewtonRaphson` implementation
- Nonlinear implicit Euler integrator

### v0.3
- Dynamic simulation: `MechanicalState`
- Lumped mass assembly
- Rayleigh damping
- Implicit Euler integrator

### v0.2
- `read_vtk` in `orfas-io`
- `EliminationMethod`
- `FixedNode` with per-axis flags
- Node inspector in viewer

### v0.1
- Linear tetrahedral elements (CST 3D)
- Structured mesh generation
- Penalty method for boundary conditions
- Direct LU solver
- Interactive egui viewer

---

## Roadmap

| Version | Status | Description |
|---------|--------|-------------|
| **v0.1** | ✅ Done | Static 3D FEM, linear elasticity, tetrahedral mesh, direct solver, basic egui viewer |
| **v0.2** | ✅ Done | VTK mesh loading, elimination boundary conditions, per-node forces, interactive node inspector |
| **v0.3** | ✅ Done | Dynamic simulation, implicit Euler time integration, Rayleigh damping, MechanicalState |
| **v0.4** | ✅ Done | Nonlinear materials (SVK), Newton-Raphson solver, nonlinear implicit Euler, refactored MaterialLaw |
| **v0.5** | ✅ Done | Neo-Hookean material, NewtonRaphsonCachedK, viewer refactor, improved camera |
| **v0.6.0** | ✅ Done | **Performance** — sparse solvers (CsrMatrix), conjugate gradient with ILU(0) preconditioner, NewtonRaphsonSparse |
| **v0.6.1** | ✅ Done | **Performance** — parallel sparse assembly (rayon, ~10x speedup), pre-built CSR pattern, SparseAssemblyStrategy trait |
| **v0.7.0** | ✅ Done | **Materials** — isochoric/volumetric split architecture, NeoHookeanIso, MooneyRivlinIso, OgdenIso, VolumetricLnJ, VolumetricQuad |
| **v0.7.1** | ✅ Done | **Materials** — HolzapfelOgden anisotropic fibers, ViscoelasticMaterial (Prony series), MaterialContext/SimulationContext architecture, FiberField, InternalVariables |
| **v0.7.2** | ⬜ Planned | **Materials** — `orfas-tissues` preset library with nominal values, confidence intervals and literature references; automatic thermodynamic consistency checks |
| **v0.8.0** | ⬜ Planned | **Elements** — `FiniteElement` trait abstraction, Tet10 (quadratic, reduces shear locking) |
| **v0.8.1** | ⬜ Planned | **Elements** — Hex8, shells, beams |
| **v0.9.0** | ⬜ Planned | **I/O** — VTU/VTK export (ParaView), OBJ/STL/MSH import |
| **v0.9.1** | ⬜ Planned | **I/O** — `.orfas` simulation file format (mesh, material, BCs, solver params, ORFAS version hash, results checksum) |
| **v0.10.0** | ⬜ Planned | **Scientific** — SI unit system via `uom` crate, compile-time unit safety |
| **v0.10.1** | ⬜ Planned | **Scientific** — standardized validation metrics (Hausdorff distance, nodal RMSE, force-displacement curve error) |
| **v0.10.2** | ⬜ Planned | **Scientific** — automated analytical benchmarks with pass/fail reports and convergence curves |
| **v0.11.0** | ⬜ Planned | **Architecture** — stable public trait audit, `SimulationPlugin` trait, static and dynamic plugin loading |
| **v0.11.1** | ⬜ Planned | **Architecture** — multi-object scene graph, `BarycentricMapping`, `RigidMapping` |
| **v0.11.2** | ⬜ Planned | **Architecture** — `.orfas` format extended for full scene graph serialization |
| **v0.12.0** | ⬜ Planned | **Rendering** — OpenGL/wgpu renderer, separate visual mesh from simulation mesh |
| **v0.12.1** | ⬜ Planned | **Rendering** — shaders, transparency, colormap pipeline, clinician-friendly viewer |
| **v0.13.0** | ⬜ Planned | **Multi-body** — rigid bodies (6 DOF), FEM-rigid coupling |
| **v0.13.1** | ⬜ Planned | **Multi-body** — articulated constraints (revolute, ball-and-socket) |
| **v0.14.0** | ⬜ Planned | **Contact** — broad-phase collision detection (BVH/AABB) |
| **v0.14.1** | ⬜ Planned | **Contact** — LCP formulation, Projected Gauss-Seidel solver, Coulomb friction |
| **v0.14.2** | ⬜ Planned | **Contact** — self-contact (hollow organs, folding soft tissue, brain under compression) |
| **v0.14.3** | ⬜ Planned | **Contact** — tool/mesh contact, surgical instrument interaction, haptic-ready force feedback pipeline |
| **v0.15.0** | ⬜ Planned | **Reduction** — modal decomposition, model order reduction (ROM) |
| **v0.15.1** | ⬜ Planned | **Reduction** — reduced basis methods; prerequisite for real-time and haptic on realistic meshes |
| **v0.16.0** | ⬜ Planned | **Topology** — infrastructure for runtime mesh modification (add/remove elements) |
| **v0.16.1** | ⬜ Planned | **Topology** — real-time cutting with local remeshing, suture simulation |
| **v0.17.0** | ⬜ Planned | **Inverse** — gradient descent and Nelder-Mead parameter identification |
| **v0.17.1** | ⬜ Planned | **Inverse** — Bayesian uncertainty propagation, Monte Carlo, confidence interval maps |
| **v0.18.0** | ⬜ Planned | **Interfaces** — Python bindings via PyO3, numpy-compatible API |
| **v0.18.1** | ⬜ Planned | **Interfaces** — WebAssembly via wasm-bindgen, browser-based simulation |
| **v0.18.2** | ⬜ Planned | **Interfaces** — DICOM/NIfTI import, automated patient-specific meshing from CT/MRI |
| **v0.19.0** | ⬜ Planned | **Real-time** — hard timing constraints (<10ms/frame), adaptive time stepping |
| **v0.19.1** | ⬜ Planned | **Real-time** — 1000Hz haptic loop, force rendering, integration with Geomagic/OpenHaptics |
| **v1.0** | ⬜ Planned | **Release** — stable C ABI, C/C++ FFI bindings, integration as a SOFA plugin, clinical-grade API documentation |
| **v1.x** | ⬜ Planned | **Clinical validation** — benchmarks against experimental datasets, patient-specific case studies |

---

## Getting Started

> **Note:** ORFAS is in pre-alpha. The API is unstable and subject to change between versions.

### Prerequisites

- Rust 1.75 or later (`rustup update stable`)
- Cargo (included with Rust)

### Build

```bash
git clone https://github.com/LLorgeou21/orfas.git
cd orfas
cargo build --release
```

### Run the viewer

```bash
cargo run -p orfas-viewer
```

### Run the tests

```bash
cargo test
```

Numerical validation results are documented in [VALIDATION.md](VALIDATION.md).

---

## References

The analytical derivations of stress and elasticity tensors implemented in ORFAS are based on the
following publications:

- Cheng, J., & Zhang, L. T. (2018). A general approach to derive stress and elasticity tensors for hyperelastic isotropic and anisotropic biomaterials. *International Journal of Computational Methods*, 15(04), 1850028. https://doi.org/10.1142/S0219876218500287

- Connolly, S. J., Mackenzie, D., & Gorash, Y. (2019). Isotropic hyperelasticity in principal stretches: explicit elasticity tensors and numerical implementation. *Computational Mechanics*, 64(5), 1273–1288. https://doi.org/10.1007/s00466-019-01707-1

- Holzapfel, G. A., & Gasser, T. C. (2001). A viscoelastic model for fiber-reinforced composites at finite strains: Continuum basis, computational aspects and applications. *Computer Methods in Applied Mechanics and Engineering*, 190, 4379–4403. https://doi.org/10.1016/S0045-7825(00)00323-6

---

## Contributing

ORFAS is a young project and contributions of all kinds are welcome.

Whether you are a numerical methods researcher, a Rust developer, a biomedical engineer, or simply
curious about simulation — there is a place for you here.

### Ways to contribute

- **Report bugs** — open an issue with a minimal reproducible case
- **Suggest features** — open a discussion before submitting a PR for large changes
- **Implement a new material law** — implement the `IsochoricPart` or `MaterialLaw` trait and open a PR
- **Improve documentation** — doc comments, examples, and guides are always needed
- **Validate results** — comparisons against SOFA or analytical solutions are invaluable

### Good first issues

- Additional boundary condition types
- New element types (Tet10, Hex8)
- Export to `.vtu` for ParaView visualization
- Benchmarks against known analytical solutions

### Code style

ORFAS follows standard Rust conventions. All public items must be documented with `///` comments.
Run `cargo fmt` and `cargo clippy` before submitting a pull request.

---

## Related Projects

| Project | Language | Focus | Notes |
|---------|----------|-------|-------|
| [SOFA](https://www.sofa-framework.org/) | C++ | Medical simulation | Primary inspiration for ORFAS |
| [Fenris](https://github.com/InteractiveComputerGraphics/fenris) | Rust | Solid mechanics / graphics | Most mature Rust FEM library, inactive since 2022 |
| [FEBio](https://febio.org/) | C++ | Biomechanics | Strong focus on biological tissues |
| [deal.II](https://www.dealii.org/) | C++ | General FEM | Reference academic FEM library |

---

## License

ORFAS is distributed under the terms of the [Apache License, Version 2.0](LICENSE).

---

## Author

Developed by [LLorgeou21](https://github.com/LLorgeou21).

*Contributions and collaborations welcome — see the Contributing section above.*