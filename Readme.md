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
with nonlinear hyperelastic materials on tetrahedral meshes and implicit Euler time integration with
Newton-Raphson. Contributions, feedback, and collaborations are warmly welcomed.

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
frame. Any hyperelastic material model (Saint Venant-Kirchhoff, Neo-Hookean, ...) implements this
trait via three methods: `pk2_stress(F)`, `tangent_stiffness(F)`, and `strain_energy(F)`. Swapping
material laws requires no changes to the assembly or solver.

**`BoundaryConditionMethod`** — Defines how Dirichlet boundary conditions are applied to the system.
Current implementations: penalty method and elimination method.

**`Solver`** — Solves the assembled linear system `K·u = f`.
Current implementation: direct LU decomposition via nalgebra.

**`NonlinearSolver`** — Solves the nonlinear static problem `R(u) = f_int(u) - f_ext = 0`.
Current implementations: `NewtonRaphson` (reassembles K at each iteration) and
`NewtonRaphsonCachedK` (factorizes K once, for materials with constant tangent stiffness such as SVK).

**`DampingModel`** — Defines how the damping matrix `C` is computed from the mass and stiffness matrices.
Current implementation: Rayleigh damping `C = α·M + β·K`.

**`IntegratorMethod`** — Defines a time integration scheme advancing the `MechanicalState` by one step `dt`.
Current implementation: implicit Euler with internal Newton-Raphson loop, supporting nonlinear materials.

**`MechanicalState`** — Holds the dynamic state of the simulated object: position, velocity, and acceleration vectors.
Exposes vector operations (`v_op`, `add_mv`) inspired by SOFA's `MechanicalState` abstraction.

---

## Changelog

### v0.5
- **Neo-Hookean hyperelastic material** (`NeoHookean`) — compressible formulation:
  `W = μ/2·(I₁−3) − μ·ln(J) + λ/2·(ln J)²`, with `S = μ(I − C⁻¹) + λ·ln(J)·C⁻¹`
  and full analytical tangent stiffness `C_tangent = dS/dE` depending on `F` — requires K reassembly at each Newton iteration
- **`NewtonRaphsonCachedK`** — Newton-Raphson variant that factorizes `K_tangent` once at `u=0`
  and reuses the LU factorization across all Newton iterations; valid for SVK where `C_tangent` is constant;
  reduces cost from `N×O(n³)` to `O(n³) + N×O(n²)` for `N` Newton iterations
- **Viewer refactored** into 5 focused modules (`app.rs`, `state.rs`, `simulation.rs`, `render.rs`, `main.rs`)
- **Improved camera** — orbit with explicit target point, Shift+drag pan, multiplicative scroll zoom,
  pitch clamp to avoid gimbal flip, auto-focus on mesh bounding box at load/simulate,
  proper view matrix via `right/up/forward` vectors, depth-sorted node rendering (painter's algorithm)
- Neo-Hookean and Newton (cached K) selectable in the viewer UI

### v0.4
- Nonlinear hyperelastic material: `SaintVenantKirchhoff` — replaces `LinearElastic`, reduces to linear elasticity for small deformations
- Fully refactored `MaterialLaw` trait: `pk2_stress(F)`, `tangent_stiffness(F)`, `strain_energy(F)` — all taking the deformation gradient `F` as input
- `NonlinearSolver` trait and `NewtonRaphson` implementation with SOFA-style normalized convergence criteria
- `assemble_tangent` and `assemble_internal_forces` on `Assembler` — geometry cached at init, forces computed from `F = I + Σ uᵢ ⊗ ∇Nᵢ`
- `ImplicitEulerIntegrator` extended with internal Newton-Raphson loop — nonlinear dynamic simulation
- `BoundaryConditionResult::reconstruct_ref` — non-consuming reconstruction for use inside iterative loops
- `restrict_matrix` / `restrict_vector` — shared helpers for free-DOF extraction, used by both solver and integrator
- Viewer updated: `SaintVenantKirchhoff` material, `Newton-Raphson` solver selectable in UI

### v0.3
- Dynamic simulation via implicit Euler time integration
- Lumped mass matrix assembly (`assemble_mass` on `Assembler`)
- Rayleigh damping (`DampingModel` trait, `RayleighDamping` implementation)
- `MechanicalState` struct with position, velocity, acceleration and vector operations
- `IntegratorMethod` trait and `ImplicitEulerIntegrator` implementation
- Refactored boundary conditions: `Constraint` and `Load` as independent scene components (removed `fixed` from `Node`)
- Real-time dynamic simulation in the viewer: Initialize / Play / Pause / Reset
- Static and dynamic simulation modes selectable in the UI
- Scrollable left panel in the viewer
- `density` added to `MaterialLaw` trait and `LinearElastic`
- Rayleigh coefficients α, β and time step `dt` exposed in the UI

### v0.2
- Load external tetrahedral meshes from `.vtk` files (`orfas-io`)
- Per-node force assignment via interactive viewer
- Per-direction boundary conditions (`FixedNode` with x/y/z bools)
- Elimination method for boundary conditions (reduces system size, more numerically stable than penalty)
- Interactive node inspector: click a node to view its position, toggle fixed, set applied force
- Fixed nodes and force-loaded nodes highlighted in the viewer

### v0.1
- Static 3D FEM solver for linear elastic materials
- Linear tetrahedral elements (CST 3D) with constant strain-displacement matrix B
- Structured mesh generation (`Mesh::generate`)
- Penalty method for zero-displacement boundary conditions
- Direct LU solver
- Interactive egui viewer with 3D perspective projection, rotation, zoom
- Deformed shape visualization with displacement colormap

---

## Roadmap

| Version | Status | Description |
|---------|--------|-------------|
| **v0.1** | ✅ Done | Static 3D FEM, linear elasticity, tetrahedral mesh, direct solver, basic egui viewer |
| **v0.2** | ✅ Done | VTK mesh loading, elimination boundary conditions, per-node forces, interactive node inspector |
| **v0.3** | ✅ Done | Dynamic simulation, implicit Euler time integration, Rayleigh damping, MechanicalState |
| **v0.4** | ✅ Done | Nonlinear materials (SVK), Newton-Raphson solver, nonlinear implicit Euler, refactored MaterialLaw |
| **v0.5** | ✅ Done | Neo-Hookean material, NewtonRaphsonCachedK, viewer refactor, improved camera |
| **v0.6.0** | ⬜ Planned | **Performance** — sparse solvers (nalgebra-sparse), conjugate gradient with ILU preconditioner |
| **v0.6.1** | ⬜ Planned | **Performance** — multi-threaded assembly (rayon), parallel force assembly |
| **v0.7.0** | ⬜ Planned | **Materials** — Mooney-Rivlin, Ogden |
| **v0.7.1** | ⬜ Planned | **Materials** — Holzapfel-Ogden (anisotropic fibers), viscoelastic (Maxwell, Kelvin-Voigt) |
| **v0.7.2** | ⬜ Planned | **Materials** — `orfas-tissues` preset library with nominal values, confidence intervals and literature references; automatic thermodynamic consistency checks (non-negative energy, objectivity, positive-definiteness of tangent) |
| **v0.8.0** | ⬜ Planned | **Elements** — `FiniteElement` trait abstraction, Tet10 (quadratic, reduces shear locking) |
| **v0.8.1** | ⬜ Planned | **Elements** — Hex8, shells, beams |
| **v0.9.0** | ⬜ Planned | **I/O** — VTU/VTK export (ParaView), OBJ/STL/MSH import |
| **v0.9.1** | ⬜ Planned | **I/O** — `.orfas` simulation file format (mesh, material, BCs, solver params, ORFAS version hash, results checksum — designed for scientific reproducibility and publication); per-element field export (stress, strain, energy), force-displacement CSV, structured per-timestep logging |
| **v0.10.0** | ⬜ Planned | **Scientific** — SI unit system via `uom` crate, compile-time unit safety (Pa vs kPa, mm vs m), documented numerical precision per solver |
| **v0.10.1** | ⬜ Planned | **Scientific** — standardized validation metrics (Hausdorff distance, nodal RMSE, force-displacement curve error), parametric sensitivity analysis |
| **v0.10.2** | ⬜ Planned | **Scientific** — automated analytical benchmarks (traction, bending, torsion, internal pressure) with pass/fail reports and convergence curves; numerical regression snapshot testing |
| **v0.11.0** | ⬜ Planned | **Architecture** — stable public trait audit, `SimulationPlugin` trait, static and dynamic plugin loading (`libloading`), `orfas-plugin-template` crate, API versioning with semver |
| **v0.11.1** | ⬜ Planned | **Architecture** — multi-object scene graph, `BarycentricMapping`, `RigidMapping`; prerequisite for multi-body and visual rendering |
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
| **v0.17.0** | ⬜ Planned | **Inverse** — gradient descent and Nelder-Mead parameter identification, identify E/ν/μ from experimental data, patient-specific calibration |
| **v0.17.1** | ⬜ Planned | **Inverse** — Bayesian uncertainty propagation, Monte Carlo and polynomial chaos expansion, confidence interval maps on simulation results |
| **v0.18.0** | ⬜ Planned | **Interfaces** — Python bindings via PyO3, numpy-compatible API, Jupyter-friendly scripting |
| **v0.18.1** | ⬜ Planned | **Interfaces** — WebAssembly via wasm-bindgen, browser-based simulation, JavaScript/TypeScript API; enables surgical training platforms and shareable preoperative planning |
| **v0.18.2** | ⬜ Planned | **Interfaces** — DICOM/NIfTI import, ITK/SimpleITK integration via Python bindings, automated patient-specific meshing from CT/MRI |
| **v0.19.0** | ⬜ Planned | **Real-time** — hard timing constraints (<10ms/frame), adaptive time stepping, simulation thread decoupled from render thread |
| **v0.19.1** | ⬜ Planned | **Real-time** — 1000Hz haptic loop, force rendering, integration with Geomagic/OpenHaptics |
| **v1.0** | ⬜ Planned | **Release** — stable C ABI, C/C++ FFI bindings, integration as a SOFA plugin, clinical-grade API documentation |
| **v1.x** | ⬜ Planned | **Clinical validation** — benchmarks against experimental datasets, patient-specific case studies, publication of validation results |

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

## Contributing

ORFAS is a young project and contributions of all kinds are welcome.

Whether you are a numerical methods researcher, a Rust developer, a biomedical engineer, or simply
curious about simulation — there is a place for you here.

### Ways to contribute

- **Report bugs** — open an issue with a minimal reproducible case
- **Suggest features** — open a discussion before submitting a PR for large changes
- **Implement a new material law** — implement the `MaterialLaw` trait and open a PR
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