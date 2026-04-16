use nalgebra::{DVector, DMatrix};
use crate::mechanical_state::MechanicalState;
use crate::solver::{Solver, SolverError};

/// Trait defining a time integration method.
/// Each integrator advances the MechanicalState by one time step dt.
/// The solver is passed as a parameter to decouple the integration scheme
/// from the linear system resolution — consistent with SOFA's architecture.
pub trait IntegratorMethod {
    fn step(
        &self,
        state: &mut MechanicalState,
        mass: &DVector<f64>,
        k: &DMatrix<f64>,
        c: &DMatrix<f64>,
        f: &DVector<f64>,
        dt: f64,
        solver: &dyn Solver,
    ) -> Result<(), SolverError>;
}

pub struct ImplicitEulerIntegrator;

impl IntegratorMethod for ImplicitEulerIntegrator {

    /// Advances the state by one time step using implicit Euler integration.
    /// Solves the velocity-based system:
    /// (M/dt + C + dt*K) * v_next = M*v/dt + f - K*u
    /// Then updates position: u_next = u + dt*v_next
    fn step(
        &self,
        state: &mut MechanicalState,
        mass: &DVector<f64>,
        k: &DMatrix<f64>,
        c: &DMatrix<f64>,
        f: &DVector<f64>,
        dt: f64,
        solver: &dyn Solver,
    ) -> Result<(), SolverError> {

        let n = mass.len();

        // Build the full mass matrix from the lumped mass vector (diagonal)
        let mut m_mat = DMatrix::zeros(n, n);
        for i in 0..n {
            m_mat[(i, i)] = mass[i];
        }

        // System matrix: A = M/dt + C + dt*K
        let a = m_mat.scale(1.0 / dt) + c + k.scale(dt);

        // Right-hand side: rhs = M*v/dt + f - K*u
        let mut rhs = DVector::zeros(n);
        MechanicalState::add_mv(&mut rhs, 1.0 / dt, mass, &state.velocity);
        rhs += f - k * &state.position;

        // Solve A * v_next = rhs
        let v_next = match solver.solve(&a, &rhs) {
            Ok(v)  => v,
            Err(e) => return Err(e),
        };

        // Update state: u_next = u + dt*v_next
        state.position = &state.position + dt * &v_next;
        state.velocity = v_next;

        Ok(())
    }
}