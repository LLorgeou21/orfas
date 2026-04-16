use nalgebra::{DMatrix, DVector};
use nalgebra::linalg::LU;
use crate::assembler::{Assembler, BMatrix, LinearBMatrix};
use crate::boundary::BoundaryConditionResult;
use crate::material::MaterialLaw;
use crate::mesh::Mesh;

// ─── Errors ───────────────────────────────────────────────────────────────────

/// All errors that can occur during solving.
#[derive(Debug)]
pub enum SolverError {
    /// K is singular and cannot be inverted.
    /// Typically caused by insufficient boundary conditions (rigid body motion)
    /// or a degenerate element with zero/negative volume.
    SingularMatrix,

    /// Newton-Raphson did not converge within max_iter iterations.
    /// Contains the final residual norm and correction norm for diagnostics.
    MaxIterationsReached { residual: f64, correction: f64 },
}

// ─── Linear solver ────────────────────────────────────────────────────────────

/// Solves a linear system K*u = f.
/// U is organized as [ux0, uy0, uz0, ux1, uy1, uz1, ...] —
/// 3 components per node in the same order as the mesh node list.
pub trait Solver {
    fn solve(&self, k: &DMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError>;
}

/// Direct LU decomposition solver (dense, via nalgebra).
/// Suitable for small to medium meshes. O(n^3) in time, O(n^2) in memory.
pub struct DirectSolver;

impl Solver for DirectSolver {
    fn solve(&self, k: &DMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError> {
        let lu = LU::new(k.clone());
        lu.solve(f).ok_or(SolverError::SingularMatrix)
    }
}

// ─── Nonlinear solver ─────────────────────────────────────────────────────────

/// Solves the nonlinear static problem R(u) = f_int(u) - f_ext = 0
/// by iterating on the reduced system produced by BoundaryConditionResult.
///
/// Returns the reduced displacement vector u_reduced on convergence.
/// Use BoundaryConditionResult::reconstruct to get the full displacement vector.
pub trait NonlinearSolver {
    fn solve<B: BMatrix>(
        &self,
        assembler:     &Assembler,
        mesh:          &Mesh,
        material:      &dyn MaterialLaw,
        bc_result:     &BoundaryConditionResult,
        linear_solver: &dyn Solver,
    ) -> Result<DVector<f64>, SolverError>;
}

/// Newton-Raphson nonlinear solver.
///
/// At each iteration:
///   1. Assemble K_tangent(u) and f_int(u) on the full mesh
///   2. Restrict to free DOFs via bc_result.free_dofs
///   3. Compute residual R = f_int_reduced - f_ext_reduced
///   4. Solve K_tangent_reduced * du = -R
///   5. Update u_reduced += du
///
/// Convergence is checked on two normalized criteria (SOFA convention):
///   - ||R|| / ||f_ext|| < tol_residual
///   - ||du|| / (||u|| + 1e-14) < tol_correction
///
/// Stops as soon as either criterion is met, or after max_iter iterations.
pub struct NewtonRaphson {
    pub max_iter:       usize,
    pub tol_residual:   f64,
    pub tol_correction: f64,
}

impl Default for NewtonRaphson {
    fn default() -> Self {
        NewtonRaphson {
            max_iter:       20,
            tol_residual:   1e-6,
            tol_correction: 1e-6,
        }
    }
}

impl NonlinearSolver for NewtonRaphson {
    fn solve<B: BMatrix>(
        &self,
        assembler:     &Assembler,
        mesh:          &Mesh,
        material:      &dyn MaterialLaw,
        bc_result:     &BoundaryConditionResult,
        linear_solver: &dyn Solver,
    ) -> Result<DVector<f64>, SolverError> {

        let n_full     = bc_result.n_dofs;
        let f_ext_red  = &bc_result.f;
        let f_ext_norm = f_ext_red.norm().max(1e-14);

        // Start from zero displacement
        let mut u_reduced = DVector::zeros(f_ext_red.len());

        for _iter in 0..self.max_iter {

            // Reconstruct full displacement to pass to assembler
            let u_full = bc_result.reconstruct_ref(&u_reduced, n_full);

            // Assemble K_tangent and f_int on the full mesh
            let k_full    = assembler.assemble_tangent::<B>(mesh, material, &u_full);
            let f_int_full = assembler.assemble_internal_forces(mesh, material, &u_full);

            // Restrict to free DOFs
            let k_red     = restrict_matrix(&k_full,     &bc_result.free_dofs);
            let f_int_red = restrict_vector(&f_int_full, &bc_result.free_dofs);

            // Residual R = f_int - f_ext (we want R = 0)
            let residual      = &f_int_red - f_ext_red;
            let residual_norm = residual.norm() / f_ext_norm;

            if residual_norm < self.tol_residual {
                return Ok(u_reduced);
            }

            // Solve K * du = -R
            let du = linear_solver.solve(&k_red, &(-&residual))?;

            let correction_norm = du.norm() / (u_reduced.norm() + 1e-14);
            u_reduced += &du;

            if correction_norm < self.tol_correction {
                return Ok(u_reduced);
            }
        }

        // Did not converge — compute final residual for diagnostics
        let u_full     = bc_result.reconstruct_ref(&u_reduced, n_full);
        let f_int_full = assembler.assemble_internal_forces(mesh, material, &u_full);
        let f_int_red  = restrict_vector(&f_int_full, &bc_result.free_dofs);
        let residual   = (&f_int_red - f_ext_red).norm() / f_ext_norm;

        Err(SolverError::MaxIterationsReached { residual, correction: 0.0 })
    }
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Restricts a full square matrix to the free DOFs.
/// Uses row-major indexing (r, s) — consistent with EliminationMethod.
/// For penalty method (free_dofs = None), returns the full matrix.
pub fn restrict_matrix(
    k_full:    &DMatrix<f64>,
    free_dofs: &Option<Vec<usize>>,
) -> DMatrix<f64> {
    match free_dofs {
        None => k_full.clone(),
        Some(free) => {
            let n = free.len();
            DMatrix::from_fn(n, n, |r, s| k_full[(free[r], free[s])])
        }
    }
}

/// Restricts a full vector to the free DOFs.
/// For penalty method (free_dofs = None), returns the full vector.
pub fn restrict_vector(
    f_full:    &DVector<f64>,
    free_dofs: &Option<Vec<usize>>,
) -> DVector<f64> {
    match free_dofs {
        None => f_full.clone(),
        Some(free) => DVector::from_iterator(free.len(), free.iter().map(|&i| f_full[i])),
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembler::{Assembler, LinearBMatrix};
    use crate::boundary::{BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load};
    use crate::material::SaintVenantKirchhoff;
    use crate::mesh::Mesh;
    use nalgebra::Vector3;

    fn make_bar(nx: usize, ny: usize, nz: usize)
        -> (Mesh, SaintVenantKirchhoff, Assembler, Vec<usize>)
    {
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let mat  = SaintVenantKirchhoff { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 };
        let assembler = Assembler::new(&mesh);
        let tips: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();
        (mesh, mat, assembler, tips)
    }

    fn make_bc_elim(
        mesh: &Mesh, nx: usize, ny: usize, nz: usize,
        tips: &[usize], force: Vector3<f64>,
    ) -> BoundaryConditions {
        let constraint = Constraint {
            list: (0..ny).flat_map(|j| {
                (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
            }).collect(),
        };
        let load = Load { list: tips.to_vec(), force };
        BoundaryConditions::new(constraint, vec![load], Box::new(EliminationMethod))
    }

    /// Newton doit converger pour SVK (tangente constante).
    /// SVK n'est pas lineaire — f_int depend de E = 0.5*(F^T*F - I).
    /// Newton converge vers la vraie solution SVK, pas la solution lineaire exacte.
    #[test]
    fn test_newton_converges_for_svk() {
        let nx = 3; let ny = 2; let nz = 2;
        let (mesh, mat, assembler, tips) = make_bar(nx, ny, nz);
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let bc = make_bc_elim(&mesh, nx, ny, nz, &tips, Vector3::new(1.0, 0.0, 0.0));
        let bc_result = bc.apply(&k, mesh.nodes.len());

        let newton = NewtonRaphson { max_iter: 20, tol_residual: 1e-8, tol_correction: 1e-8 };
        let result = newton.solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &DirectSolver);
        assert!(result.is_ok(), "Newton doit converger pour SVK : {:?}", result);
    }

    /// Newton doit satisfaire le residu f_int(u) - f_ext = 0 a la convergence.
    /// On verifie que la solution Newton satisfait bien l'equilibre.
    #[test]
    fn test_newton_satisfies_residual() {
        let nx = 3; let ny = 2; let nz = 2;
        let (mesh, mat, assembler, tips) = make_bar(nx, ny, nz);
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let bc = make_bc_elim(&mesh, nx, ny, nz, &tips, Vector3::new(1.0, 0.0, 0.0));
        let bc_result = bc.apply(&k, mesh.nodes.len());
        let n_full = bc_result.n_dofs;
        let f_ext = bc_result.f.clone();

        let u_newton_red = NewtonRaphson::default()
            .solve::<LinearBMatrix>(&assembler, &mesh, &mat, &bc_result, &DirectSolver)
            .unwrap();

        // Verifier que f_int(u) ~ f_ext
        let u_full = bc_result.reconstruct_ref(&u_newton_red, n_full);
        let f_int_full = assembler.assemble_internal_forces(&mesh, &mat, &u_full);
        let f_int_red = restrict_vector(&f_int_full, &bc_result.free_dofs);
        let residual_norm = (&f_int_red - &f_ext).norm() / f_ext.norm();
        println!("Newton residual : {:.2e}", residual_norm);
        assert!(residual_norm < 1e-6, "Residu Newton trop grand : {:.2e}", residual_norm);
    }
}