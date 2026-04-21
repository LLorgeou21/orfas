use nalgebra::{DMatrix, DVector};
use crate::assembler::{Assembler, BMatrix, LinearBMatrix};
use crate::boundary::BoundaryConditionResult;
use crate::material::MaterialLaw;
use crate::mechanical_state::MechanicalState;
use crate::mesh::Mesh;
use crate::solver::{restrict_matrix, restrict_vector, DenseSolver, SolverError};

// ─── Trait ────────────────────────────────────────────────────────────────────

/// Defines a time integration method advancing MechanicalState by one step dt.
///
/// In v0.4 the integrator owns the Newton-Raphson loop internally —
/// consistent with SOFA's EulerImplicitSolver which embeds its own
/// nonlinear solve rather than delegating to an external Newton DenseSolver.
///
/// The MechanicalState (position, velocity) lives in the reduced space
/// of free DOFs produced by BoundaryConditionResult.
///
/// Arguments:
/// - state        : current position and velocity in reduced space
/// - mass         : lumped mass vector in reduced space (diagonal)
/// - c            : damping matrix in reduced space
/// - f_ext        : external force vector in reduced space
/// - dt           : time step
/// - assembler    : used to reassemble K_tangent and f_int each Newton iter
/// - mesh         : reference mesh (undeformed)
/// - material     : material law — provides pk2_stress and tangent_stiffness
/// - bc_result    : boundary condition result — used to map reduced <-> full
/// - linear_solver: solves the linear system at each Newton iteration
pub trait IntegratorMethod {
    fn step<B: BMatrix>(
        &self,
        state:          &mut MechanicalState,
        mass:           &DVector<f64>,
        c:              &DMatrix<f64>,
        f_ext:          &DVector<f64>,
        dt:             f64,
        assembler:      &Assembler,
        mesh:           &Mesh,
        material:       &dyn MaterialLaw,
        bc_result:      &BoundaryConditionResult,
        linear_solver:  &dyn DenseSolver,
    ) -> Result<(), SolverError>;
}

// ─── Implicit Euler ───────────────────────────────────────────────────────────

/// Implicit Euler time integrator with internal Newton-Raphson loop.
///
/// Solves at each time step:
///   R(v_next) = M*(v_next - v)/dt + C*v_next + f_int(u + dt*v_next) - f_ext = 0
///
/// Newton-Raphson linearizes R around the current v_next estimate:
///   A * dv = -R(v_next)
///   A = M/dt + C + dt*K_tangent(u_next)
///
/// Then updates:
///   v_next <- v_next + dv
///   u_next  = u + dt*v_next
///
/// The dt factor in dt*K_tangent comes from the chain rule:
///   d(f_int(u_next))/d(v_next) = K_tangent * d(u_next)/d(v_next) = K_tangent * dt
///
/// Convergence criteria (normalized, SOFA convention):
///   ||R|| / ||f_eff|| < tol_residual
///   ||dv|| / (||v_next|| + 1e-14) < tol_correction
///
/// where f_eff = M*v/dt + f_ext is the effective RHS norm for normalization.
pub struct ImplicitEulerIntegrator {
    pub max_iter:       usize,
    pub tol_residual:   f64,
    pub tol_correction: f64,
}

impl Default for ImplicitEulerIntegrator {
    fn default() -> Self {
        ImplicitEulerIntegrator {
            max_iter:       20,
            tol_residual:   1e-6,
            tol_correction: 1e-6,
        }
    }
}

impl IntegratorMethod for ImplicitEulerIntegrator {

    fn step<B: BMatrix>(
        &self,
        state:         &mut MechanicalState,
        mass:          &DVector<f64>,
        c:             &DMatrix<f64>,
        f_ext:         &DVector<f64>,
        dt:            f64,
        assembler:     &Assembler,
        mesh:          &Mesh,
        material:      &dyn MaterialLaw,
        bc_result:     &BoundaryConditionResult,
        linear_solver: &dyn DenseSolver,
    ) -> Result<(), SolverError> {

        let n = mass.len();
        let n_full = bc_result.n_dofs;

        // Build diagonal mass matrix from lumped mass vector
        let mut m_mat = DMatrix::zeros(n, n);
        for i in 0..n {
            m_mat[(i, i)] = mass[i];
        }

        // Effective RHS norm for residual normalization:
        // f_eff = M*v/dt + f_ext
        let mut f_eff = DVector::zeros(n);
        MechanicalState::add_mv(&mut f_eff, 1.0 / dt, mass, &state.velocity);
        f_eff += f_ext;
        let f_eff_norm = f_eff.norm().max(1e-14);

        // Initial guess: v_next = v (explicit predictor)
        let mut v_next = state.velocity.clone();

        for _iter in 0..self.max_iter {

            // u_next = u + dt*v_next  (in reduced space)
            let u_next_red = &state.position + dt * &v_next;

            // Reconstruct full displacement to pass to assembler
            let u_next_full = bc_result.reconstruct_ref(&u_next_red, n_full);

            // Assemble K_tangent and f_int on the full mesh
            let k_full     = assembler.assemble_tangent::<B>(mesh, material, &u_next_full);
            let f_int_full = assembler.assemble_internal_forces(mesh, material, &u_next_full);

            // Restrict K_tangent and f_int to free DOFs
            let k_red     = restrict_matrix(&k_full,     &bc_result.free_dofs);
            let f_int_red = restrict_vector(&f_int_full, &bc_result.free_dofs);

            // Residual: R = M*(v_next - v)/dt + C*v_next + f_int - f_ext
            let residual =
                m_mat.scale(1.0 / dt) * (&v_next - &state.velocity)
                + c * &v_next
                + &f_int_red
                - f_ext;

            let residual_norm = residual.norm() / f_eff_norm;
            if residual_norm < self.tol_residual {
                // Converged — commit state
                state.position = u_next_red;
                state.velocity = v_next;
                return Ok(());
            }

            // System matrix: A = M/dt + C + dt*K_tangent
            let a = m_mat.scale(1.0 / dt) + c + k_red.scale(dt);

            // Solve A * dv = -R
            let dv = linear_solver.solve(&a, &(-&residual))?;

            let correction_norm = dv.norm() / (v_next.norm() + 1e-14);
            v_next += &dv;

            if correction_norm < self.tol_correction {
                // Converged — commit state
                let u_next_red = &state.position + dt * &v_next;
                state.position = u_next_red;
                state.velocity = v_next;
                return Ok(());
            }
        }

        Err(SolverError::MaxIterationsReached {
            residual:   f64::NAN,
            correction: f64::NAN,
        })
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembler::Assembler;
    use crate::boundary::{BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load};
    use crate::damping::RayleighDamping;
    use crate::damping::DampingModel;
    use crate::material::SaintVenantKirchhoff;
    use crate::mesh::Mesh;
    use nalgebra::Vector3;

    /// With strong damping, the dynamic simulation must converge
    /// to the static solution (same benchmark as v0.3 but with SVK).
    #[test]
    fn test_implicit_euler_converges_to_static() {
        let nx = 5; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let mat  = SaintVenantKirchhoff {
            youngs_modulus: 1e6,
            poisson_ratio:  0.3,
            density:        1000.0,
        };
        let assembler = Assembler::new(&mesh);

        let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();
        let nb_tip = tip_nodes.len() as f64;
        let f_total = 100.0;

        let make_bc = || {
            let constraint = Constraint {
                list: (0..ny).flat_map(|j| {
                    (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
                }).collect(),
            };
            let load = Load {
                list:  tip_nodes.clone(),
                force: Vector3::new(f_total / nb_tip, 0.0, 0.0),
            };
            BoundaryConditions::new(constraint, vec![load], Box::new(EliminationMethod))
        };

        // Reference statique via Newton
        let bc_static = make_bc();
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k_ref = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let bc_result_static = bc_static.apply(&k_ref, mesh.nodes.len());
        let u_static_red = crate::solver::DirectSolver
            .solve(&bc_result_static.k, &bc_result_static.f).unwrap();
        let u_static = bc_result_static.reconstruct(u_static_red);
        let mean_static = tip_nodes.iter()
            .map(|&i| u_static[3 * i]).sum::<f64>() / nb_tip;

        // Simulation dynamique
        let k_full = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let mass   = assembler.assemble_mass(&mesh, &mat);
        let damping = RayleighDamping { alpha: 10.0, beta: 0.01 };
        let c_full  = damping.compute(&mass, &k_full);

        let bc_dyn = make_bc();
        let bc_result_dyn = bc_dyn.apply(&k_full, mesh.nodes.len());

        let f_ext = bc_result_dyn.f.clone();
        let (mass_red, c_red) = match &bc_result_dyn.free_dofs {
            Some(free) => {
                let m = DVector::from_iterator(free.len(), free.iter().map(|&i| mass[i]));
                let c = DMatrix::from_fn(free.len(), free.len(), |r, s| {
                    c_full[(free[r], free[s])]
                });
                (m, c)
            }
            None => (mass.clone(), c_full.clone()),
        };

        let n_red = mass_red.len();
        let mut state = MechanicalState::new(n_red / 3);
        let integrator = ImplicitEulerIntegrator::default();
        let solver = crate::solver::DirectSolver;

        for _ in 0..500 {
            integrator.step::<LinearBMatrix>(
                &mut state, &mass_red, &c_red, &f_ext, 1.0,
                &assembler, &mesh, &mat, &bc_result_dyn, &solver,
            ).unwrap();
        }

        let mean_dyn = match &bc_result_dyn.free_dofs {
            Some(free) => tip_nodes.iter().map(|&idx| {
                free.iter().position(|&d| d == 3 * idx)
                    .map(|j| state.position[j])
                    .unwrap_or(0.0)
            }).sum::<f64>() / nb_tip,
            None => tip_nodes.iter().map(|&i| state.position[3 * i]).sum::<f64>() / nb_tip,
        };

        let error = (mean_dyn - mean_static).abs() / mean_static;
        println!(
            "dynamic={:.8} static={:.8} error={:.4}%",
            mean_dyn, mean_static, error * 100.0
        );
        assert!(error < 0.01, "Dynamic simulation must converge to static solution: {:.4}%", error * 100.0);
    }
}