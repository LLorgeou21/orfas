// UTF-8
// orfas-viewer/src/simulation.rs — static and dynamic simulation entry points.

use nalgebra::{DMatrix, DVector};
use orfas_core::{
    assembler::Assembler,
    boundary::{
        BoundaryConditionResult, BoundaryConditions, Constraint,
        EliminationMethod, PenaltyMethod,
    },
    damping::{DampingModel, RayleighDamping},
    element::Tet4,
    material::SimulationContext,
    mechanical_state::MechanicalState,
    mesh::Tet4Mesh,
    solver::{DirectSolver, NewtonRaphson, NewtonRaphsonCachedK, NonlinearSolver, DenseSolver},
    sparse::{CgSolver, NewtonRaphsonSparse, NonlinearSparseSolver, Parallel, Sequential},
};

use orfas_core::element::FiniteElement;
use orfas_core::element::traits::DofType;
use crate::state::{AppState, BoundaryChoice, SolverChoice, make_material};

// ─── Internal build helper ────────────────────────────────────────────────────

/// Build the mesh, assemble initial K and mass, apply boundary conditions.
/// Returns (mesh, mass_red, c_red, bc_result) or None if mesh is unavailable.
pub fn build_simulation(
    state: &mut AppState,
) -> Option<(Tet4Mesh, DVector<f64>, DMatrix<f64>, BoundaryConditionResult)> {
    let mesh = match state.mesh.take() {
        Some(m) => m,
        None => {
            let m = Tet4Mesh::generate(
                state.nx, state.ny, state.nz,
                state.dx, state.dy, state.dz,
            );
            state.constraint    = Constraint { list: Vec::new() };
            state.loads         = Vec::new();
            state.selected_load = None;
            m
        }
    };

    state.camera.focus_on_mesh(&mesh.nodes);

    let material  = make_material(state);
    let assembler = Assembler::<Tet4>::new(&mesh);
    let u_zero = DVector::zeros(<Tet4 as FiniteElement>::Dof::N_DOF * mesh.nodes.len());
    let sim_ctx   = SimulationContext::isotropic_static(mesh.elements.len());

    let k_full = assembler.assemble_tangent(&mesh, material.as_ref(), &u_zero, &sim_ctx);
    let mass   = assembler.assemble_mass(&mesh, material.as_ref());
    let c_full = RayleighDamping { alpha: state.alpha, beta: state.beta }
        .compute(&mass, &k_full);

    let bc = match state.boundary {
        BoundaryChoice::Penalty => BoundaryConditions::new(
            state.constraint.clone(),
            state.loads.clone(),
            Box::new(PenaltyMethod),
        ),
        BoundaryChoice::Elimination => BoundaryConditions::new(
            state.constraint.clone(),
            state.loads.clone(),
            Box::new(EliminationMethod),
        ),
    };
    let bc_result = bc.apply(&k_full, mesh.nodes.len(), <Tet4 as FiniteElement>::Dof::N_DOF);

    let (mass_red, c_red) = match &bc_result.free_dofs {
        Some(free) => (
            DVector::from_iterator(free.len(), free.iter().map(|&i| mass[i])),
            DMatrix::from_fn(free.len(), free.len(), |r, s| c_full[(free[r], free[s])]),
        ),
        None => (mass, c_full),
    };

    Some((mesh, mass_red, c_red, bc_result))
}

// ─── Static simulation ────────────────────────────────────────────────────────

/// Run a static FEM simulation and store the displacement field in state.
pub fn run_simulation_static(state: &mut AppState) {
    if let Some((mesh, _mass, _c, bc_result)) = build_simulation(state) {
        let material  = make_material(state);
        let assembler = Assembler::<Tet4>::new(&mesh);
        let sim_ctx   = SimulationContext::isotropic_static(mesh.elements.len());

        let result = match state.solver {
            SolverChoice::Direct => DirectSolver
                .solve(&bc_result.k, &bc_result.f)
                .map(|u_red| bc_result.reconstruct(u_red)),

            SolverChoice::Newton => NewtonRaphson::default()
                .solve(
                    &assembler, &mesh, material.as_ref(), &bc_result, &DirectSolver,
                    &sim_ctx,
                )
                .map(|u_red| bc_result.reconstruct(u_red)),

            SolverChoice::NewtonCachedK => NewtonRaphsonCachedK::default()
                .solve(
                    &assembler, &mesh, material.as_ref(), &bc_result, &DirectSolver,
                    &sim_ctx,
                )
                .map(|u_red| bc_result.reconstruct(u_red)),

            SolverChoice::NewtonSparse => NewtonRaphsonSparse::<Sequential>::default()
                .solve(
                    &assembler, &mesh, material.as_ref(), &bc_result, &CgSolver::default(),
                    &sim_ctx,
                )
                .map(|u_red| bc_result.reconstruct(u_red)),

            SolverChoice::NewtonSparseParallel => NewtonRaphsonSparse::<Parallel>::default()
                .solve(
                    &assembler, &mesh, material.as_ref(), &bc_result, &CgSolver::default(),
                    &sim_ctx,
                )
                .map(|u_red| bc_result.reconstruct(u_red)),
        };

        match result {
            Ok(u)  => state.displacements = Some(u),
            Err(e) => println!("Simulation error: {:?}", e),
        }
        state.mesh = Some(mesh);
    }
}

// ─── Dynamic simulation ───────────────────────────────────────────────────────

/// Initialize dynamic simulation state (mass, damping, BCs, mechanical state).
/// Called once before the time-stepping loop.
pub fn init_simulation_dynamic(state: &mut AppState) {
    if let Some((mesh, mass_red, c_red, bc_result)) = build_simulation(state) {
        let n_red        = mass_red.len();
        let n_dofs_total = <Tet4 as FiniteElement>::Dof::N_DOF * mesh.nodes.len();
        let mech         = MechanicalState::new(<Tet4 as FiniteElement>::Dof::N_DOF, n_red / <Tet4 as FiniteElement>::Dof::N_DOF);
        let init_pos     = mech.position.clone();
        let f_ext        = bc_result.f.clone();

        state.initial_positions = Some(init_pos);
        state.n_dofs_total      = n_dofs_total;
        state.mechanical_state  = Some(mech);
        state.mass              = Some(mass_red);
        state.c_cached          = Some(c_red);
        state.f_cached          = Some(f_ext);
        state.bc_result_cached  = Some(bc_result);
        state.displacements     = None;
        state.mesh              = Some(mesh);
    }
}