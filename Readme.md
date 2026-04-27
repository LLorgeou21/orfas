# ORFAS
### Open Rust Framework for Articulated Simulation

[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Status: Pre-Alpha](https://img.shields.io/badge/Status-Pre--Alpha-red.svg)]()
[![Rust](https://img.shields.io/badge/Rust-2024%20Edition-orange.svg)](https://www.rust-lang.org/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)]()

> *An open-source Finite Element Method framework in Rust, designed for medical and biomechanical simulation.*

---

![ORFAS dynamic simulation](simulation.gif)

---

## What is ORFAS?

ORFAS is a modular FEM framework targeting soft tissue simulation for medical and surgical
applications. It is written entirely in Rust and designed around trait-based composability:
material laws, element types, solvers, and integrators are all interchangeable without
modifying the core assembly pipeline.

The primary domain is biomechanics — passive soft tissue, fiber-reinforced biological
structures, and viscoelastic organs — in the tradition of frameworks such as
[SOFA](https://www.sofa-framework.org/) and [FEBio](https://febio.org/), with a modern
memory-safe foundation.

---

## Motivation

### Why Rust for FEM?

Existing FEM frameworks for medical simulation are predominantly written in C++ (SOFA, FEBio,
deal.II). While mature and capable, these frameworks carry significant adoption friction:
complex build systems, heavy dependency chains, and architectures designed before modern
concurrency and memory safety were first-class language concerns.

Rust offers: memory safety without a garbage collector, zero-cost abstractions, a
world-class package manager (Cargo), and a trait system that enables generic programming
with compile-time guarantees. In simulation code, where subtle memory bugs and performance
regressions are costly, these properties translate directly into correctness and
maintainability.

### Why not extend SOFA or Fenris?

[SOFA](https://www.sofa-framework.org/) is the reference framework for medical simulation
and the direct inspiration for ORFAS. Its C++ plugin architecture presents a high integration
barrier for projects seeking a modern, memory-safe foundation and a simpler contribution model.

[Fenris](https://github.com/InteractiveComputerGraphics/fenris) is the most mature FEM library
in the Rust ecosystem. Its focus is solid mechanics for computer graphics, and it has been
inactive since 2022. ORFAS targets medical simulation and a different design philosophy: full
replaceability of every numerical component via traits, and a calibrated tissue library
grounded in the biomechanics literature.

---

## Current Capabilities (v0.7.2)

### Simulation

- Static and dynamic (implicit Euler) nonlinear FEM on tetrahedral meshes
- Newton-Raphson and Newton-Raphson with cached tangent
- Sparse solvers with ILU(0)-preconditioned conjugate gradient
- Parallel sparse assembly via Rayon (~10x speedup on large meshes)
- Rayleigh damping

### Material library (`orfas-core`)

All hyperelastic materials are implemented in the Lagrangian frame with analytical
PK2 stress and consistent tangent stiffness tensors.

| Model | Type | Notes |
|---|---|---|
| `SaintVenantKirchhoff` | Hyperelastic | Linear elastic tangent, nonlinear at large strains |
| `NeoHookeanIso` | Isochoric | Flory decoupled, W = μ/2·(Ī₁ − 3) |
| `MooneyRivlinIso` | Isochoric | W = c₁·(Ī₁ − 3) + c₂·(Ī₂ − 3) |
| `OgdenIso` | Isochoric | Principal stretch formulation, arbitrary N |
| `VolumetricLnJ` | Volumetric | W = κ/2·(ln J)² |
| `VolumetricQuad` | Volumetric | W = κ/2·(J − 1)² |
| `HolzapfelOgden` | Anisotropic | Fiber-reinforced, exponential strain energy |
| `ViscoelasticMaterial<I,A,V>` | Viscoelastic | Prony series, Holzapfel & Gasser (2001) |

Any isochoric model composes with any volumetric model via `CompressibleMaterial<I, V>`.
Anisotropic fiber contributions are added via `CompressibleAnisotropicMaterial<I, A, V>`.
Viscoelastic dissipation wraps any composition via `ViscoelasticMaterial<I, A, V>`.

### Tissue preset library (`orfas-tissues`)

Calibrated material presets for 10 biological tissues, each with literature reference,
experimental protocol, and per-parameter confidence intervals.

| Tissue | Model | Reference |
|---|---|---|
| Liver | Neo-Hookean | Nava et al. (2008) |
| Brain grey matter | Mooney-Rivlin | Budday et al. (2017) |
| Brain white matter | Mooney-Rivlin | Budday et al. (2017) |
| Myocardium (passive) | Holzapfel-Ogden | Holzapfel & Ogden (2009) |
| Arterial wall (media) | Holzapfel-Ogden | Holzapfel et al. (2000) |
| Tendon (ground matrix) | Neo-Hookean | Weiss et al. (1996) |
| Ligament (MCL) | Holzapfel-Ogden | Weiss et al. (1996) |
| Skin | Mooney-Rivlin | Groves et al. (2013) |
| Kidney | Neo-Hookean | Nasseri et al. (2002) |
| Prostate | Neo-Hookean | Phipps et al. (2005) |

Presets can be loaded programmatically or from JSON files at runtime.
Each preset is validated against thermodynamic consistency checks at compile time.

### Thermodynamic consistency checks (`orfas-core`)

`check_thermodynamic_consistency` verifies any `MaterialLaw` implementation against
8 necessary conditions: zero energy at rest, zero stress at rest, Hooke linearization,
tangent symmetry, non-negative strain energy, small-strain convergence, positive
definiteness, and frame objectivity. Available at runtime for custom materials loaded
from JSON.

### Viewer (`orfas-viewer`)

Interactive egui viewer with orbit camera, node picking, displacement colormap,
material selection (manual parameters or tissue presets), boundary condition editor,
and static/dynamic simulation controls.

---

## Getting Started

> ORFAS is in pre-alpha. The public API is unstable and subject to change between versions.

### Prerequisites

- Rust 2024 edition (`rustup update stable`)
- Cargo (included with Rust)

### Build

```bash
git clone https://github.com/LLorgeou21/orfas.git
cd orfas
cargo build --release
```

### Run the viewer

```bash
cargo run --release -p orfas-viewer
```

### Run the tests

```bash
cargo test
```

Numerical validation results are documented in [VALIDATION.md](VALIDATION.md).

---

## Roadmap

| Version | Status | Description |
|---|---|---|
| **v0.1** | ✅ Done | Static 3D FEM, linear elasticity, tetrahedral mesh, direct solver, basic egui viewer |
| **v0.2** | ✅ Done | VTK mesh loading, elimination boundary conditions, per-node forces, interactive node inspector |
| **v0.3** | ✅ Done | Dynamic simulation, implicit Euler time integration, Rayleigh damping, MechanicalState |
| **v0.4** | ✅ Done | Nonlinear materials (SVK), Newton-Raphson solver, nonlinear implicit Euler, refactored MaterialLaw |
| **v0.5** | ✅ Done | Neo-Hookean material, NewtonRaphsonCachedK, viewer refactor, improved camera |
| **v0.6.0** | ✅ Done | **Performance** — sparse solvers (CsrMatrix), conjugate gradient with ILU(0) preconditioner, NewtonRaphsonSparse |
| **v0.6.1** | ✅ Done | **Performance** — parallel sparse assembly (rayon, ~10x speedup), pre-built CSR pattern, SparseAssemblyStrategy trait |
| **v0.7.0** | ✅ Done | **Materials** — isochoric/volumetric split architecture, NeoHookeanIso, MooneyRivlinIso, OgdenIso, VolumetricLnJ, VolumetricQuad |
| **v0.7.1** | ✅ Done | **Materials** — HolzapfelOgden anisotropic fibers, ViscoelasticMaterial (Prony series), MaterialContext/SimulationContext architecture, FiberField, InternalVariables |
| **v0.7.2** | ✅ Done | **Materials** — `orfas-tissues` preset library with confidence intervals and literature references; thermodynamic consistency checks; viewer UI reorganization |
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

## Related Projects

| Project | Language | Focus | Notes |
|---|---|---|---|
| [SOFA](https://www.sofa-framework.org/) | C++ | Medical simulation | Primary inspiration for ORFAS |
| [FEBio](https://febio.org/) | C++ | Biomechanics | Strong focus on biological tissues |
| [Fenris](https://github.com/InteractiveComputerGraphics/fenris) | Rust | Solid mechanics / graphics | Most mature Rust FEM library, inactive since 2022 |
| [deal.II](https://www.dealii.org/) | C++ | General FEM | Reference academic FEM library |

---

## References

### Material law derivations

- Cheng, J., & Zhang, L. T. (2018). A general approach to derive stress and elasticity tensors for hyperelastic isotropic and anisotropic biomaterials. *International Journal of Computational Methods*, 15(04). https://doi.org/10.1142/S0219876218500287

- Connolly, S. J., Mackenzie, D., & Gorash, Y. (2019). Isotropic hyperelasticity in principal stretches: explicit elasticity tensors and numerical implementation. *Computational Mechanics*, 64(5), 1273–1288. https://doi.org/10.1007/s00466-019-01707-1

- Holzapfel, G. A., & Gasser, T. C. (2001). A viscoelastic model for fiber-reinforced composites at finite strains. *Computer Methods in Applied Mechanics and Engineering*, 190, 4379–4403. https://doi.org/10.1016/S0045-7825(00)00323-6

### Tissue presets (orfas-tissues)

- **Liver** — Nava, A., Mazza, E., Furrer, M., Villiger, P., & Reinhart, W. H. (2008). In vivo mechanical characterization of human liver. *Medical Image Analysis*, 12(2), 203–216. https://doi.org/10.1016/j.media.2007.09.001

- **Brain (grey and white matter)** — Budday, S., Nay, R., de Rooij, R., Steinmann, P., Wyrobek, T., Ovaert, T. C., & Kuhl, E. (2017). Mechanical properties of gray and white matter brain tissue by indentation. *Acta Biomaterialia*, 48, 319–330. https://doi.org/10.1016/j.actbio.2016.10.036

- **Myocardium (passive)** — Holzapfel, G. A., & Ogden, R. W. (2009). Constitutive modelling of passive myocardium: a structurally based framework for material characterization. *Philosophical Transactions of the Royal Society A*, 367, 3445–3475. https://doi.org/10.1098/rsta.2009.0091

- **Arterial wall (media)** — Holzapfel, G. A., Gasser, T. C., & Ogden, R. W. (2000). A new constitutive framework for arterial wall mechanics and a comparative study with other models. *Journal of Elasticity*, 61, 1–48. https://doi.org/10.1023/A:1010835316564

- **Tendon and ligament (MCL)** — Weiss, J. A., Maker, B. N., & Govindjee, S. (1996). Finite element implementation of incompressible, transversely isotropic hyperelasticity. *Computer Methods in Applied Mechanics and Engineering*, 135(1–2), 107–128. https://doi.org/10.1016/0045-7825(95)00931-0

- **Skin** — Groves, R. B., Coulman, S. A., Birchall, J. C., & Evans, S. L. (2013). An anisotropic, hyperelastic model for skin: experimental measurements, finite element modelling and identification of parameters for human and murine skin. *Journal of the Mechanical Behavior of Biomedical Materials*, 18, 167–180. https://doi.org/10.1016/j.jmbbm.2012.10.021

- **Kidney** — Nasseri, S., Bilston, L. E., & Phan-Thien, N. (2002). Viscoelastic properties of pig kidney in shear, experimental results and modelling. *Rheologica Acta*, 41(1–2), 180–192. https://doi.org/10.1007/s00397-002-0233-2

- **Prostate** — Phipps, S., Yang, T. H. J., Habib, F. K., Reuben, R. L., & McNeill, S. A. (2005). Measurement of the mechanical properties of the prostate. *Journal of Biomechanics*, 38(8), 1733–1741. https://doi.org/10.1016/j.jbiomech.2004.07.028

---

## License

ORFAS is distributed under the terms of the [Apache License, Version 2.0](LICENSE).

## Contributing

See [DEVELOPER.md](DEVELOPER.md) for architecture details, code conventions, and contribution guides.

## Author

Developed by [LLorgeou21](https://github.com/LLorgeou21).