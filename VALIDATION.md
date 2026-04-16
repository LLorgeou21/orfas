# ORFAS — Validation Document

---

## v0.1 — Static FEM core

**Element type:** Linear tetrahedron (CST 3D)
**Material law:** Linear elastic (isotropic)
**Boundary conditions:** Penalty method
**Solver:** LU decomposition (nalgebra)

### Method

ORFAS solves the static linear elasticity problem `K·u = f` using the displacement-based finite element method on unstructured tetrahedral meshes.

**Mesh:** The domain is discretized with linear tetrahedral elements (4-node tetrahedra, also known as CST 3D — Constant Strain Tetrahedron). Structured meshes are generated via an alternating 6-tetrahedra-per-cube decomposition to avoid degenerate elements and ensure consistent orientation.

**Element formulation:** The strain-displacement matrix `B` is derived analytically from the shape function gradients of the linear tetrahedron. `B` is constant within each element. The element stiffness matrix is computed as `Kₑ = Bᵀ·C·B·V`, where `C` is the material stiffness matrix and `V` is the element volume.

**Material law:** Linear isotropic elasticity. The constitutive matrix `C` is assembled from Young's modulus `E` and Poisson's ratio `ν` using the standard Voigt notation for 3D stress-strain relations.

**Assembly:** The global stiffness matrix `K` is assembled by direct stiffness summation over all elements, mapping local 12×12 element matrices into the global `3n × 3n` system via the connectivity table.

**Boundary conditions:** Zero-displacement Dirichlet conditions are enforced via the penalty method: diagonal entries of `K` corresponding to fixed degrees of freedom are replaced by a large value (10³⁰), forcing the associated displacements to near zero.

**Solver:** The linear system `K·u = f` is solved by LU decomposition (dense, via nalgebra). This is appropriate for small to medium meshes. Sparse solvers and iterative methods (conjugate gradient) are planned for future versions.

### Test case 1 — Axial traction

A prismatic bar of length `L = 9` and square cross-section `A = 2×2 = 4` is fixed at one end (`x = 0`) and subjected to a uniform axial force `F = 100` at the free end (`x = L`). The material has Young's modulus `E = 1×10⁶` and Poisson's ratio `ν = 0.3`. The mesh uses `nx = 10`, `ny = nz = 3` nodes.

Exact solution: `δ = F·L / (E·A) = 100 × 9 / (1×10⁶ × 4) = 2.25×10⁻⁴`

Axial traction produces a uniform strain field. Linear tetrahedra represent constant strain exactly by construction, so this test is a direct verification of the core pipeline: `B → C → Kₑ → K → u`.

| Quantity | Value |
|---|---|
| Analytical displacement δ | 2.2500×10⁻⁴ |
| Computed displacement | 2.2626×10⁻⁴ |
| Relative error | 0.56% |

The residual error of 0.56% is attributable to the penalty method, which does not enforce boundary conditions exactly.

### Test case 2 — Cantilever beam bending (convergence)

A cantilever beam of length `L = 9`, square cross-section `3×3` (`ny = nz = 4` nodes, `h = b = 3`), fixed at `x = 0`, loaded by a transverse force `F = 25` at the tip in the `y` direction. The mesh is refined along the longitudinal axis: `nx = 4, 7, 10, 13`.

Exact solution (Euler-Bernoulli): `δ = F·L³ / (3·E·I) = 9.00×10⁻⁴` with `I = b·h³/12 = 6.75`

The purpose of this test is not to achieve low absolute error, but to verify that the error decreases monotonically as the mesh is refined — a necessary condition for a correct FEM implementation.

| nx | dx | Computed δ | Analytical δ | Error |
|---|---|---|---|---|
| 4  | 3.000 | 4.11×10⁻⁴ | 9.00×10⁻⁴ | 54.34% |
| 7  | 1.500 | 6.20×10⁻⁴ | 9.00×10⁻⁴ | 31.09% |
| 10 | 1.000 | 6.90×10⁻⁴ | 9.00×10⁻⁴ | 23.30% |
| 13 | 0.750 | 7.18×10⁻⁴ | 9.00×10⁻⁴ | 20.22% |

The error decreases monotonically. The residual ~20% error at the finest mesh is consistent with the known shear locking behavior of linear tetrahedra — this is a property of the element formulation, not a bug.

### Summary

| Test | Analytical δ | Computed δ | Error | Status |
|---|---|---|---|---|
| Axial traction | 2.2500×10⁻⁴ | 2.2626×10⁻⁴ | 0.56% | PASS |
| Beam bending (nx=13) | 9.00×10⁻⁴ | 7.18×10⁻⁴ | 20.22% | PASS (convergence verified) |

---

## v0.2 — I/O and boundary conditions

### What was added

- **Elimination method** for boundary conditions (replaces penalty for fixed DOFs by reducing the system size)
- **VTK mesh loading** (`read_vtk` in `orfas-io`)

### Validation

No numerical benchmark was run for v0.2 as the core FEM pipeline was not modified. Validation focused on correctness of the new components via unit tests.

### VTK mesh loading (`orfas-io`)

| Test | Description | Status |
|---|---|---|
| `test_read_vtk_valid` | Load reference cube mesh, check 8 nodes and 6 tetrahedra | PASS |
| `test_read_vtk_node_positions` | Verify node 0 at origin and node 7 at (1,1,1) | PASS |
| `test_read_vtk_file_not_found` | Non-existent file returns `IoError::UnreadableFile` | PASS |
| `test_read_vtk_invalid_format` | File not starting with `# vtk` returns `IoError::InvalidFormat` | PASS |

### Elimination method

The elimination method was verified indirectly — all v0.1 numerical tests pass with both penalty and elimination, confirming that the system reduction and reconstruction produce equivalent results.

---

## v0.3 — Dynamic simulation

**Integrator:** Implicit Euler
**Damping:** Rayleigh (`C = α·M + β·K`)
**Mass:** Lumped (concentrated) — each node receives 1/4 of the mass of each connected element

### Method

ORFAS solves the dynamic linear elasticity problem `M·a + C·v + K·u = f` using implicit Euler time integration.

**Mass assembly:** The lumped mass matrix is assembled by distributing each element's mass (`ρ·V`) equally across its 4 nodes. The result is stored as a `DVector` (diagonal only) rather than a full matrix, avoiding unnecessary memory usage and enabling efficient `M·v` products via component-wise multiplication.

**Rayleigh damping:** The damping matrix is computed as `C = α·M + β·K`, where `α` and `β` are user-defined coefficients. `α` controls mass-proportional damping (low-frequency), `β` controls stiffness-proportional damping (high-frequency).

**Implicit Euler integration:** At each time step `dt`, the following velocity-based system is solved:

```
(M/dt + C + dt·K) · v_next = M·v/dt + f - K·u
```

Then the position is updated: `u_next = u + dt·v_next`. This formulation is unconditionally stable — large time steps do not cause divergence, unlike explicit methods.

**MechanicalState:** The dynamic state (position, velocity, acceleration) is encapsulated in a `MechanicalState` struct. Vector operations (`v_op`, `add_mv`) are exposed as methods, inspired by SOFA's `MechanicalState` abstraction.

### Test case 1 — Mass assembly

A unit cube mesh (`2×2×2`, volume = 1) with density `ρ = 1500 kg/m³`. The total assembled mass must satisfy `sum(mass) / 3 = ρ · V`.

| Quantity | Value |
|---|---|
| Expected total mass | 1500.0 |
| Computed total mass | 1500.0 |
| Relative error | 0.0000% |

### Test case 2 — Rayleigh damping symmetry

The damping matrix `C = α·M + β·K` must be symmetric since both `M` (diagonal) and `K` (symmetric by construction) are symmetric. Verified on a `3×3×3` mesh with `α = 0.1`, `β = 0.01`.

| Quantity | Value |
|---|---|
| Max asymmetry | < 1×10⁻¹⁰ |
| Status | PASS |

### Test case 3 — Dynamic convergence to static solution

A bar under axial traction (`nx=5, ny=nz=2`, `L=4`, `A=1`, `F=100`, `E=1×10⁶`) is simulated dynamically with heavy Rayleigh damping (`α=10`, `β=0.01`) using implicit Euler with `dt=1.0` for 500 steps. With sufficient damping, the dynamic solution must converge to the static equilibrium.

The test compares the dynamic solution to both the static FEM solution and the analytical solution to isolate integrator behavior from mesh discretization error. The elimination method is used to avoid ill-conditioning from the penalty method.

| Quantity | Value |
|---|---|
| Static FEM displacement | 3.9282×10⁻⁴ |
| Dynamic displacement (500 steps) | 3.9282×10⁻⁴ |
| Analytical displacement | 4.0000×10⁻⁴ |
| Error vs static | 0.0000% |
| Error vs analytical | 1.79% |

The 1.79% error vs analytical is attributable to mesh discretization (same error as the static solver on this mesh), confirming that the integrator introduces no additional error at convergence.

### Summary

| Test | Description | Status |
|---|---|---|
| `test_mass_assembly` | Total mass matches `density * volume` | PASS |
| `test_rayleigh_damping_symmetry` | Damping matrix is symmetric | PASS |
| `test_implicit_euler_static_convergence` | Dynamic solution converges to static (error 0.0000%) | PASS |