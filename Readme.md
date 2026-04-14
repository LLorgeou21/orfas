# ORFAS
### Open Rust Framework Architecture for Simulation

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Status: Pre-Alpha](https://img.shields.io/badge/Status-Pre--Alpha-red.svg)]()
[![Rust](https://img.shields.io/badge/Rust-2021%20Edition-orange.svg)](https://www.rust-lang.org/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)]()

> *A generic, extensible Finite Element Method framework written in Rust, with a primary focus on medical and biomechanical simulation.*

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

This project is in active early development. The current milestone targets a static 3D FEM solver
with linear elasticity on tetrahedral meshes. Contributions, feedback, and collaborations are
warmly welcomed.

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
├── orfas-core/     # Core library: mesh, materials, assembly, solvers
├── orfas-io/       # Mesh I/O: .vtk, .obj file formats
└── orfas-viewer/   # Interactive viewer for development and testing (egui)
```

### Core abstractions (`orfas-core`)

The framework is built around three central traits:

**`MaterialLaw`** — Defines the constitutive relationship between strain and stress.
Any material model (linear elasticity, Neo-Hookean, Saint Venant-Kirchhoff, ...) implements
this trait. Swapping material laws requires no changes to the assembly or solver.

```rust
pub trait MaterialLaw {
    fn stress(&self, deformation_gradient: &Matrix3) -> Matrix3;
    fn tangent_moduli(&self, deformation_gradient: &Matrix3) -> Tensor4;
}
```

**`Element`** — Defines how a finite element contributes to the global stiffness matrix.
The first implementation targets linear tetrahedral elements (Tet4).

```rust
pub trait Element {
    fn stiffness_matrix(&self, nodes: &[Node], material: &dyn MaterialLaw) -> DMatrix;
    fn shape_functions(&self, xi: &Point3) -> DVector;
}
```

**`Solver`** — Solves the assembled linear or nonlinear system.
Static and dynamic solvers implement this trait independently.

```rust
pub trait Solver {
    fn solve(&self, system: &AssembledSystem) -> Result<DVec, SolverError>;
}
```

---

## Roadmap

| Version | Status | Description |
|---------|--------|-------------|
| **v0.1** | 🔴 In Progress | Static 3D FEM, linear elasticity, tetrahedral mesh, direct solver, egui viewer |
| **v0.2** | ⬜ Planned | External mesh loading (`.vtk`, `.obj`), configurable boundary conditions |
| **v0.3** | ⬜ Planned | Dynamic simulation, implicit Euler time integration |
| **v0.4** | ⬜ Planned | Nonlinear material laws: Neo-Hookean, Saint Venant-Kirchhoff |
| **v0.5** | ⬜ Planned | Python bindings via PyO3, numpy-compatible API |
| **v0.6** | ⬜ Planned | Mesh/mesh contact and collision detection |
| **v0.7** | ⬜ Planned | Tool/mesh contact for surgical interaction |
| **v0.8** | ⬜ Planned | Sparse solvers, parallelism, performance optimization |
| **v0.9** | ⬜ Planned | Real-time simulation with hard timing constraints |
| **v0.10** | ⬜ Planned | Haptic feedback interface |
| **v1.0** | ⬜ Planned | C/C++ FFI bindings — integration as a SOFA plugin |

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

### Run the viewer (v0.1)

```bash
cargo run -p orfas-viewer -- --width 4 --height 4 --depth 4 --youngs-modulus 1e6 --poisson 0.45
```

*Full usage documentation will be added as v0.1 stabilizes.*

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

The following areas are good entry points for new contributors:

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