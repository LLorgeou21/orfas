use nalgebra::{DMatrix, DVector, linalg::Cholesky};


/// Defines the interface for linear system solvers.
/// Takes the global stiffness matrix K and force vector F as input.
/// Returns a Result containing the displacement vector U on success,
/// or a SolverError on failure.
/// U is organized as [ux0, uy0, uz0, ux1, uy1, uz1, ...] —
/// 3 components per node in the same order as the mesh node list.
pub trait Solver {
    fn solve(&self, k: &DMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError>;
}




/// Groups all errors that can occur during solving.
pub enum SolverError {
    /// K is singular and cannot be inverted.
    /// This typically happens when boundary conditions are insufficient
    /// (not enough fixed nodes), allowing rigid body motion.
    /// Can also occur if a tetrahedron has zero or negative volume
    /// due to a degenerate mesh.
    SingularMatrix,
}


pub struct DirectSolver;

impl Solver for DirectSolver {
    fn solve(&self, k: &DMatrix<f64>, f: &DVector<f64>) -> Result<DVector<f64>, SolverError> {
        let chol = Cholesky::new(k.clone())
            .ok_or(SolverError::SingularMatrix)?;
        Ok(chol.solve(f))
    }
}