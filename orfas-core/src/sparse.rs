// ─── Sparse solvers ───────────────────────────────────────────────────────────
//
// This module provides sparse linear and nonlinear solvers for use with large
// meshes where dense LU factorization (DirectSolver) is too slow or too memory
// intensive. All solvers here operate on CsrMatrix<f64> instead of DMatrix<f64>.
//
// Intended usage:
//   - CgSolver         : drop-in replacement for DirectSolver on large systems
//   - NewtonRaphsonSparse : drop-in replacement for NewtonRaphson using sparse assembly

use crate::assembler::{Assembler, BMatrix};
use crate::boundary::BoundaryConditionResult;
use crate::material::MaterialLaw;
use crate::material::SimulationContext;
use crate::mesh::Mesh;
use crate::solver::{restrict_vector, SolverError};
use nalgebra::{DVector};
use nalgebra_sparse::{CooMatrix, CsrMatrix};

// ─── Preconditioner ───────────────────────────────────────────────────────────

/// Preconditioner for the conjugate gradient solver.
#[derive(Debug, Clone)]
pub enum Preconditioner {
    /// Identity — no preconditioning.
    /// Simple but slow to converge on ill-conditioned systems.
    Identity,
    /// Incomplete LU factorization at fill level k.
    /// ILU(0) (k=0) is the most common choice — zero fill-in, good balance
    /// between setup cost and convergence improvement.
    Ilu(usize),
}

// ─── SparseSolver trait ───────────────────────────────────────────────────────

/// Sparse linear solver trait.
/// Mirror of DenseSolver but operates on CsrMatrix<f64>.
/// Implement this trait to provide alternative sparse solvers (e.g. GMRES, MINRES).
pub trait SparseSolver {
    fn solve(&self, k: &CsrMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError>;
}

// ─── CgSolver ─────────────────────────────────────────────────────────────────

/// Preconditioned conjugate gradient solver.
///
/// Solves the sparse linear system K*x = f iteratively.
/// Requires K to be symmetric positive definite — guaranteed after proper
/// boundary condition application (EliminationMethod or PenaltyMethod).
///
/// The `precond` field is reserved for future ILU(k) preconditioning.
/// Currently only `Preconditioner::Identity` is active.
pub struct CgSolver {
    pub max_iter:  usize,
    pub tolerance: f64,
    pub precond:   Preconditioner,
}

impl Default for CgSolver {
    fn default() -> Self {
        CgSolver {
            max_iter:  1000,
            tolerance: 1e-8,
            precond:   Preconditioner::Identity,
        }
    }
}

impl SparseSolver for CgSolver {
    fn solve(&self, k: &CsrMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError> {
        let n = f.len();
        let mut x = DVector::zeros(n);

        // Build preconditioner once before the iteration loop
        let precond: Option<Ilu0> = match &self.precond {
            Preconditioner::Identity => None,
            Preconditioner::Ilu(_)   => Some(Ilu0::new(k)),
        };

        let apply = |r: &DVector<f64>| -> DVector<f64> {
            match &precond {
                None    => r.clone(),
                Some(m) => m.apply(r),
            }
        };

        // r = f - K*x, with x=0 => r = f
        let mut r = f.clone();
        let mut z = apply(&r);
        let mut p = z.clone();
        let mut rz = r.dot(&z);

        for _ in 0..self.max_iter {
            if r.norm() < self.tolerance {
                return Ok(x);
            }

            // Sparse matrix-vector product K*p
            let kp = k * &p;

            let alpha = rz / p.dot(&kp);
            x += alpha * &p;
            r -= alpha * &kp;

            z = apply(&r);
            let rz_new = r.dot(&z);
            let beta = rz_new / rz;
            p = &z + beta * &p;
            rz = rz_new;
        }

        if r.norm() < self.tolerance {
            Ok(x)
        } else {
            Err(SolverError::MaxIterationsReached {
                residual:   r.norm(),
                correction: 0.0,
            })
        }
    }
}


/// Strategy trait for sparse tangent assembly.
/// Allows NewtonRaphsonSparse to be generic over sequential vs parallel assembly.
pub trait SparseAssemblyStrategy {
    fn assemble_tangent<B: BMatrix>(
        assembler:   &Assembler,
        mesh:        &Mesh,
        material:    &dyn MaterialLaw,
        u:           &DVector<f64>,
        sim_ctx: &SimulationContext,
    ) -> CsrMatrix<f64>;
}

/// Sequential sparse assembly strategy.
pub struct Sequential;

/// Parallel sparse assembly strategy (uses atomic f64 additions via rayon).
/// Recommended for meshes with more than ~1000 nodes.
pub struct Parallel;

impl SparseAssemblyStrategy for Sequential {
    fn assemble_tangent<B: BMatrix>(
        assembler: &Assembler,
        mesh:      &Mesh,
        material:  &dyn MaterialLaw,
        u:         &DVector<f64>,
        sim_ctx: &SimulationContext,
    ) -> CsrMatrix<f64> {
        assembler.assemble_tangent_sparse::<B>(mesh, material, u, sim_ctx)
    }
}

impl SparseAssemblyStrategy for Parallel {
    fn assemble_tangent<B: BMatrix>(
        assembler: &Assembler,
        mesh:      &Mesh,
        material:  &dyn MaterialLaw,
        u:         &DVector<f64>,
        sim_ctx: &SimulationContext,
    ) -> CsrMatrix<f64> {
        assembler.assemble_tangent_sparse_parallel::<B>(mesh, material, u,sim_ctx)
    }
}


// ─── Sparse matrix restriction ────────────────────────────────────────────────

/// Restricts a sparse square matrix to the free DOFs.
/// Filters triplets (i, j, v) — keeps only entries where both i and j are in free_dofs.
/// For penalty method (free_dofs = None), returns the full matrix.
pub fn restrict_matrix_sparse(
    k: &CsrMatrix<f64>,
    free_dofs: &Option<Vec<usize>>,
) -> CsrMatrix<f64> {
    match free_dofs {
        None => k.clone(),
        Some(free) => {
            let n = free.len();
            let global_to_local: std::collections::HashMap<usize, usize> = free
                .iter()
                .enumerate()
                .map(|(local, &global)| (global, local))
                .collect();

            let mut coo = CooMatrix::zeros(n, n);
            for (i, j, v) in k.triplet_iter() {
                if let (Some(&li), Some(&lj)) = (
                    global_to_local.get(&i),
                    global_to_local.get(&j),
                ) {
                    coo.push(li, lj, *v);
                }
            }
            CsrMatrix::from(&coo)
        }
    }
}

// ─── NonlinearSparseSolver trait ──────────────────────────────────────────────

/// Solves the nonlinear static problem R(u) = f_int(u) - f_ext = 0
/// using sparse assembly and a SparseSolver for the inner linear system.
///
/// Mirror of NonlinearSolver but uses assemble_tangent_sparse internally.
pub trait NonlinearSparseSolver {
    fn solve<B: BMatrix>(
        &self,
        assembler:     &Assembler,
        mesh:          &Mesh,
        material:      &dyn MaterialLaw,
        bc_result:     &BoundaryConditionResult,
        linear_solver: &dyn SparseSolver,
        sim_ctx: &SimulationContext,
    ) -> Result<DVector<f64>, SolverError>;
}

// ─── NewtonRaphsonSparse ──────────────────────────────────────────────────────

/// Sparse Newton-Raphson nonlinear solver.
///
/// At each iteration:
///   1. Assemble K_tangent_sparse(u) and f_int(u) on the full mesh
///   2. Restrict to free DOFs via bc_result.free_dofs
///   3. Compute residual R = f_int_reduced - f_ext_reduced
///   4. Solve K_tangent_sparse_reduced * du = -R  (via SparseSolver)
///   5. Update u_reduced += du
///
/// Convergence is checked on two normalized criteria (SOFA convention):
///   - ||R|| / ||f_ext|| < tol_residual
///   - ||du|| / (||u|| + 1e-14) < tol_correction
///
/// Stops as soon as either criterion is met, or after max_iter iterations.
pub struct NewtonRaphsonSparse<S: SparseAssemblyStrategy = Sequential> {
    pub max_iter:       usize,
    pub tol_residual:   f64,
    pub tol_correction: f64,
    _strategy:          std::marker::PhantomData<S>,
}

impl<S: SparseAssemblyStrategy> Default for NewtonRaphsonSparse<S> {
    fn default() -> Self {
        NewtonRaphsonSparse {
            max_iter:       20,
            tol_residual:   1e-6,
            tol_correction: 1e-6,
            _strategy:      std::marker::PhantomData,
        }
    }
}

impl<S: SparseAssemblyStrategy> NonlinearSparseSolver for NewtonRaphsonSparse<S> {
    fn solve<B: BMatrix>(
        &self,
        assembler:     &Assembler,
        mesh:          &Mesh,
        material:      &dyn MaterialLaw,
        bc_result:     &BoundaryConditionResult,
        linear_solver: &dyn SparseSolver,
        sim_ctx: &SimulationContext,
    ) -> Result<DVector<f64>, SolverError> {
        let n_full     = bc_result.n_dofs;
        let f_ext_red  = &bc_result.f;
        let f_ext_norm = f_ext_red.norm().max(1e-14);

        let mut u_reduced = DVector::zeros(f_ext_red.len());

        for _iter in 0..self.max_iter {
            let u_full = bc_result.reconstruct_ref(&u_reduced, n_full);

            let k_full     = S::assemble_tangent::<B>(assembler, mesh, material, &u_full, sim_ctx);
            let f_int_full = assembler.assemble_internal_forces(mesh, material, &u_full, sim_ctx);

            let k_red     = restrict_matrix_sparse(&k_full,     &bc_result.free_dofs);
            let f_int_red = restrict_vector(&f_int_full, &bc_result.free_dofs);

            let residual      = &f_int_red - f_ext_red;
            let residual_norm = residual.norm() / f_ext_norm;

            if residual_norm < self.tol_residual {
                return Ok(u_reduced);
            }

            let du = linear_solver.solve(&k_red, &(-&residual))?;

            let correction_norm = du.norm() / (u_reduced.norm() + 1e-14);
            u_reduced += &du;

            if correction_norm < self.tol_correction {
                return Ok(u_reduced);
            }
        }

        // Did not converge
        let u_full     = bc_result.reconstruct_ref(&u_reduced, n_full);
        let f_int_full = assembler.assemble_internal_forces(mesh, material, &u_full, sim_ctx);
        let f_int_red  = restrict_vector(&f_int_full, &bc_result.free_dofs);
        let residual   = (&f_int_red - f_ext_red).norm() / f_ext_norm;

        Err(SolverError::MaxIterationsReached { residual, correction: 0.0 })
    }
}

pub struct Ilu0 {
    l: CsrMatrix<f64>,
    u: CsrMatrix<f64>,
}


impl Ilu0 {
    pub fn new(k: &CsrMatrix<f64>) -> Self {
        let (l, u) = ilu0(k);
        Ilu0 { l, u }
    }

    /// Applies the ILU(0) preconditioner: solves M*x = r where M = L*U.
    /// Equivalent to two triangular solves: L*y = r, then U*x = y.
    pub fn apply(&self, r: &DVector<f64>) -> DVector<f64> {
        let y = forward_substitution(&self.l, r);
        backward_substitution(&self.u, &y)
    }
}

/// Performs an ILU(0) factorization of a sparse matrix K.
/// Returns (L, U) where L is lower triangular (unit diagonal) and U is upper triangular.
/// Only entries that are non-zero in K are updated — zero fill-in.
///
/// Uses a dense row buffer of size n during factorization to allow
/// in-place updates without mutating CsrMatrix directly.
pub fn ilu0(k: &CsrMatrix<f64>) -> (CsrMatrix<f64>, CsrMatrix<f64>) {
    let n = k.nrows();

    // Accumulators for L and U entries
    let mut l_coo = CooMatrix::zeros(n, n);
    let mut u_coo = CooMatrix::zeros(n, n);

    // Stores the fully factored rows, needed for updates of subsequent rows.
    // Row i is stored as a sparse vec: (col, value).
    let mut factored_rows: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n);

    for i in 0..n {
        // Load row i of K into a dense buffer
        let mut buf = vec![0.0f64; n];
        if let Some(row) = k.get_row(i) {
            for (&j, &v) in row.col_indices().iter().zip(row.values().iter()) {
                buf[j] = v;
            }
        }

        // Apply updates from previously factored rows k < i
        for k in 0..i {
            let l_ik = buf[k];
            if l_ik.abs() < 1e-30 { continue; }

            // l[i,k] = buf[k] / u[k,k]
            let u_kk = factored_rows[k]
                .iter()
                .find(|&&(j, _)| j == k)
                .map(|&(_, v)| v)
                .unwrap_or(1.0);
            let l_ik = l_ik / u_kk;
            buf[k] = l_ik;

            // buf[j] -= l[i,k] * u[k,j] for j > k
            for &(j, u_kj) in factored_rows[k].iter().filter(|&&(j, _)| j > k) {
                if buf[j].abs() > 1e-30 || u_kj.abs() > 1e-30 {
                    buf[j] -= l_ik * u_kj;
                }
            }
        }

        // Extract L and U entries from buffer, store factored row
        let mut factored = Vec::new();
        for j in 0..n {
            let v = buf[j];
            if v.abs() < 1e-30 { continue; }
            if j < i {
                l_coo.push(i, j, v);
            } else {
                u_coo.push(i, j, v);
                factored.push((j, v));
            }
        }
        // Unit diagonal for L
        l_coo.push(i, i, 1.0);
        factored_rows.push(factored);
    }

    (CsrMatrix::from(&l_coo), CsrMatrix::from(&u_coo))
}

/// Solves L*y = r where L is a lower triangular sparse matrix with unit diagonal.
/// Forward substitution — processes rows top to bottom.
fn forward_substitution(l: &CsrMatrix<f64>, r: &DVector<f64>) -> DVector<f64> {
    let n = r.len();
    let mut y = r.clone();
    for i in 0..n {
        if let Some(row) = l.get_row(i) {
            for (&j, &l_ij) in row.col_indices().iter().zip(row.values().iter()) {
                if j < i {
                    y[i] -= l_ij * y[j];
                }
            }
        }
        // Unit diagonal — no division needed
    }
    y
}

/// Solves U*x = y where U is an upper triangular sparse matrix.
/// Backward substitution — processes rows bottom to top.
fn backward_substitution(u: &CsrMatrix<f64>, y: &DVector<f64>) -> DVector<f64> {
    let n = y.len();
    let mut x = y.clone();
    for i in (0..n).rev() {
        if let Some(row) = u.get_row(i) {
            let mut diag = 1.0;
            for (&j, &u_ij) in row.col_indices().iter().zip(row.values().iter()) {
                if j > i {
                    x[i] -= u_ij * x[j];
                } else if j == i {
                    diag = u_ij;
                }
            }
            x[i] /= diag;
        }
    }
    x
}




// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembler::{Assembler, LinearBMatrix};
    use crate::boundary::{BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load};
    use crate::material::SaintVenantKirchhoff;
    use crate::mesh::Mesh;
    use crate::solver::{DirectSolver, NewtonRaphson, NonlinearSolver};
    use nalgebra::{DVector, Vector3};

    /// NewtonRaphsonSparse + CgSolver must match NewtonRaphson + DirectSolver
    /// within numerical tolerance. Validates both sparse assembly and CG solve.
    #[test]
    fn test_newton_sparse_matches_dense() {
        let nx = 3; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let mat  = SaintVenantKirchhoff { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 };
        let assembler = Assembler::new(&mesh);

        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero,&SimulationContext::isotropic_static(mesh.elements.len()));

        let tips: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();

        let constraint = Constraint {
            list: (0..ny).flat_map(|j| {
                (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
            }).collect(),
        };
        let load = Load { list: tips.clone(), force: Vector3::new(1.0, 0.0, 0.0) };
        let bc = BoundaryConditions::new(constraint, vec![load], Box::new(EliminationMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());

        let u_dense = NewtonRaphson::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &DirectSolver,&SimulationContext::isotropic_static(mesh.elements.len()))
            .unwrap();

        let u_sparse = NewtonRaphsonSparse::<Sequential>::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &CgSolver::default(),&SimulationContext::isotropic_static(mesh.elements.len()))
            .unwrap();

        let diff = (&u_dense - &u_sparse).norm() / u_dense.norm();
        assert!(diff < 1e-6, "sparse != dense, diff = {:.2e}", diff);

        let u_sparse_par = NewtonRaphsonSparse::<Parallel>::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &CgSolver::default(),&SimulationContext::isotropic_static(mesh.elements.len()))
            .unwrap();

        let diff = (&u_dense - &u_sparse_par).norm() / u_dense.norm();
        assert!(diff < 1e-6, "parallel sparse != dense, diff = {:.2e}", diff);
    }


    /// ILU(0) preconditioner must give the same solution as Identity preconditioner.
    /// Validates that the ILU(0) factorization and triangular solves are correct.
    #[test]
    fn test_cg_ilu_matches_identity() {
        let nx = 3; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let mat  = SaintVenantKirchhoff { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 };
        let assembler = Assembler::new(&mesh);

        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero,&SimulationContext::isotropic_static(mesh.elements.len()));

        let tips: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();

        let constraint = Constraint {
            list: (0..ny).flat_map(|j| {
                (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
            }).collect(),
        };
        let load = Load { list: tips.clone(), force: Vector3::new(1.0, 0.0, 0.0) };
        let bc = BoundaryConditions::new(constraint, vec![load], Box::new(EliminationMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());


        let cgsolver : CgSolver = CgSolver { max_iter: 1000, tolerance: 1e-8, precond: Preconditioner::Ilu(0) };
        let u_sparse_ilu0 = NewtonRaphsonSparse::<Sequential>::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &cgsolver,&SimulationContext::isotropic_static(mesh.elements.len()))
            .unwrap();

        let u_sparse = NewtonRaphsonSparse::<Sequential>::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &CgSolver::default(),&SimulationContext::isotropic_static(mesh.elements.len()))
            .unwrap();

        let diff = (&u_sparse_ilu0 - &u_sparse).norm() / u_sparse_ilu0.norm();
        assert!(diff < 1e-6, "sparse != dense, diff = {:.2e}", diff);
    }



}