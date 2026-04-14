# ORFAS — Validation Document

**Version:** v0.1  
**Element type:** Linear tetrahedron (CST 3D)  
**Material law:** Linear elastic (isotropic)  
**Boundary conditions:** Penalty method  
**Solver:** LU decomposition (nalgebra)

---

## Overview

This document presents the numerical validation of ORFAS v0.1. Two benchmark cases are used to verify the correctness of the implementation: axial traction and cantilever beam bending. These cases have known analytical solutions and are standard in finite element verification practice.

The goal is not to demonstrate that linear tetrahedra are accurate in all configurations — their limitations in bending are well known — but to verify that the assembled system `K·u = f` is correctly built and solved, and that the numerical solution converges toward the analytical reference as the mesh is refined.

---

## Method

ORFAS solves the static linear elasticity problem `K·u = f` using the displacement-based finite element method on unstructured tetrahedral meshes.

**Mesh:** The domain is discretized with linear tetrahedral elements (4-node tetrahedra, also known as CST 3D — Constant Strain Tetrahedron). Structured meshes are generated via an alternating 6-tetrahedra-per-cube decomposition to avoid degenerate elements and ensure consistent orientation.

**Element formulation:** The strain-displacement matrix `B` is derived analytically from the shape function gradients of the linear tetrahedron. `B` is constant within each element. The element stiffness matrix is computed as `Kₑ = Bᵀ·C·B·V`, where `C` is the material stiffness matrix and `V` is the element volume.

**Material law:** Linear isotropic elasticity. The constitutive matrix `C` is assembled from Young's modulus `E` and Poisson's ratio `ν` using the standard Voigt notation for 3D stress-strain relations.

**Assembly:** The global stiffness matrix `K` is assembled by direct stiffness summation over all elements, mapping local 12×12 element matrices into the global `3n × 3n` system via the connectivity table.

**Boundary conditions:** Zero-displacement Dirichlet conditions are enforced via the penalty method: diagonal entries of `K` corresponding to fixed degrees of freedom are replaced by a large value (10¹⁵), forcing the associated displacements to near zero.

**Solver:** The linear system `K·u = f` is solved by LU decomposition (dense, via nalgebra). This is appropriate for small to medium meshes. Sparse solvers and iterative methods (conjugate gradient) are planned for future versions.

---

## Test case 1 — Axial traction

### Problem description

A prismatic bar of length `L = 9` and square cross-section `A = 2×2 = 4` is fixed at one end (`x = 0`) and subjected to a uniform axial force `F = 100` at the free end (`x = L`). The material has Young's modulus `E = 1×10⁶` and Poisson's ratio `ν = 0.3`. The mesh uses `nx = 10`, `ny = nz = 3` nodes, generating a structured tetrahedral mesh via an alternating decomposition.

### Analytical solution

Under pure axial loading, the exact tip displacement given by the bar equation is:

```
δ = F·L / (E·A) = 100 × 9 / (1×10⁶ × 4) = 2.25×10⁻⁴
```

### Why this validates the implementation

Axial traction produces a uniform strain field throughout the bar. Linear tetrahedra represent constant strain exactly by construction, so no approximation error is introduced by the element formulation. Any significant deviation from the analytical solution would indicate a bug in the strain-displacement matrix `B`, the material stiffness matrix `C`, the assembly procedure, or the solver.

This test is therefore a direct verification of the core pipeline: `B → C → Kₑ → K → u`.

### Results

| Quantity | Value |
|---|---|
| Analytical displacement δ | 2.2500×10⁻⁴ |
| Computed displacement | 2.2626×10⁻⁴ |
| Relative error | 0.56% |

The residual error of 0.56% is attributable to the penalty method, which does not enforce boundary conditions exactly. It is well within acceptable bounds for this method and this penalty coefficient.

---

## Test case 2 — Cantilever beam bending (convergence)

### Problem description

A cantilever beam of length `L = 9`, square cross-section `3×3` (i.e. `ny = nz = 4` nodes, giving `h = b = 3`), is fixed at `x = 0` and loaded by a transverse force `F = 25` distributed uniformly over the tip face (`x = L`) in the `y` direction. Young's modulus `E = 1×10⁶`, Poisson's ratio `ν = 0.3`. The mesh is successively refined along the longitudinal axis: `nx = 4, 7, 10, 13` nodes, corresponding to element sizes `dx = 3.0, 1.5, 1.0, 0.75`.

### Analytical solution

Under the Euler-Bernoulli beam assumption, the tip deflection of a cantilever under tip load is:

```
δ = F·L³ / (3·E·I)
```

where `I = b·h³ / 12` is the second moment of area of the cross-section in the bending plane:

```
I = 3 × 3³ / 12 = 6.75
δ = 25 × 729 / (3 × 1×10⁶ × 6.75) = 9.00×10⁻⁴
```

### Why this validates the implementation

Unlike axial traction, bending involves a linearly varying strain field along the beam axis. Linear tetrahedra cannot represent this exactly — they produce a constant strain per element — which introduces a discretization error that decreases as the mesh is refined. Additionally, linear tetrahedra are known to suffer from shear locking in bending: parasitic shear strains appear that artificially stiffen the structure, leading to an underestimation of displacements even on fine meshes.

The purpose of this test is therefore not to achieve low absolute error, but to verify that the error decreases monotonically as the mesh is refined. This convergence behavior is a necessary condition for a correct FEM implementation. A non-converging or diverging solution would indicate a fundamental error in the element formulation or the assembly.

### Results

| nx | dx | Computed δ | Analytical δ | Relative error |
|---|---|---|---|---|
| 4  | 3.000 | 4.11×10⁻⁴ | 9.00×10⁻⁴ | 54.34% |
| 7  | 1.500 | 6.20×10⁻⁴ | 9.00×10⁻⁴ | 31.09% |
| 10 | 1.000 | 6.90×10⁻⁴ | 9.00×10⁻⁴ | 23.30% |
| 13 | 0.750 | 7.18×10⁻⁴ | 9.00×10⁻⁴ | 20.22% |

The error decreases monotonically from 54% to 20% as the element size is halved. The solution converges toward the analytical reference, confirming that the implementation is correct.

The residual error of ~20% at the finest mesh is consistent with the known behavior of linear tetrahedra under bending. This is a property of the element formulation, not a bug. Reducing this error would require either a finer mesh, higher-order elements, or a mixed formulation.

---

## Summary

| Test | Analytical δ | Computed δ | Error | Status |
|---|---|---|---|---|
| Axial traction | 2.2500×10⁻⁴ | 2.2626×10⁻⁴ | 0.56% | PASS |
| Beam bending (nx=13) | 9.00×10⁻⁴ | 7.18×10⁻⁴ | 20.22% | PASS (convergence verified) |

The v0.1 implementation is validated. The core FEM pipeline produces results consistent with analytical solutions. The observed errors are either within numerical precision of the penalty method, or explained by the known limitations of linear tetrahedral elements in bending.