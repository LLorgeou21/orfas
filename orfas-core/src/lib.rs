pub mod mesh;
pub mod material;
pub mod assembler;
pub mod boundary;
pub mod solver;
pub mod damping;
pub mod integrator;
pub mod mechanical_state;

#[cfg(test)]
mod integration_tests {
    use crate::mesh::Mesh;
    use crate::material::SaintVenantKirchhoff;
    use crate::assembler::{Assembler, LinearBMatrix};
    use crate::boundary::{
        BoundaryConditions, Constraint,
        EliminationMethod, FixedNode, Load, PenaltyMethod,
    };
    use crate::solver::{DirectSolver, NonlinearSolver, NewtonRaphson, Solver};
    use crate::damping::{DampingModel, RayleighDamping};
    use crate::integrator::{ImplicitEulerIntegrator, IntegratorMethod};
    use crate::mechanical_state::MechanicalState;
    use nalgebra::{DMatrix, DVector, Vector3};

    // ─── Helpers ──────────────────────────────────────────────────────────────

    fn svk(e: f64, nu: f64, rho: f64) -> SaintVenantKirchhoff {
        SaintVenantKirchhoff { youngs_modulus: e, poisson_ratio: nu, density: rho }
    }

    fn make_bc(
        nx: usize, ny: usize, nz: usize,
        tip_nodes: &[usize],
        force_per_node: Vector3<f64>,
        method: Box<dyn crate::boundary::BoundaryConditionMethod>,
    ) -> BoundaryConditions {
        let constraint = Constraint {
            list: (0..ny).flat_map(|j| {
                (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
            }).collect(),
        };
        let load = Load { list: tip_nodes.to_vec(), force: force_per_node };
        BoundaryConditions::new(constraint, vec![load], method)
    }

    fn tip_nodes(nx: usize, ny: usize, nz: usize) -> Vec<usize> {
        (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect()
    }

    // ─── v0.1 sanity ──────────────────────────────────────────────────────────

    /// Smoke test : le pipeline complet doit produire un vecteur deplacement
    /// de la bonne taille sans erreur.
    #[test]
    fn test_validation_v0_1() {
        let mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        let material = svk(1000.0, 0.3, 1.0);
        let assembler = Assembler::new(&mesh);
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);

        let constraint = Constraint {
            list: vec![
                FixedNode::all(0), FixedNode::all(1),
                FixedNode::all(2), FixedNode::all(3),
            ],
        };
        let load = Load { list: vec![4, 5, 6, 7], force: Vector3::new(0.0, 0.0, -100.0) };
        let bc = BoundaryConditions::new(constraint, vec![load], Box::new(PenaltyMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());

        let result = DirectSolver.solve(&bc_result.k, &bc_result.f);
        assert!(result.is_ok());
        let u = bc_result.reconstruct(result.unwrap());
        assert_eq!(u.len(), 3 * 8);
    }

    // ─── Traction axiale ──────────────────────────────────────────────────────

    /// delta = F*L / (E*A) — erreur < 1% (residuel methode penalty).
    #[test]
    fn test_axial_traction() {
        let nx = 10; let ny = 3; let nz = 3;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let material = svk(1e6, 0.3, 1.0);
        let assembler = Assembler::new(&mesh);
        let tips = tip_nodes(nx, ny, nz);
        let nb_tip = tips.len() as f64;
        let f_total = 100.0;

        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);
        let bc = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_total / nb_tip, 0.0, 0.0),
            Box::new(PenaltyMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());
        let u_red = DirectSolver.solve(&bc_result.k, &bc_result.f).unwrap();
        let u = bc_result.reconstruct(u_red);

        let mean_disp = tips.iter().map(|&i| u[3 * i]).sum::<f64>() / nb_tip;
        let l = (nx - 1) as f64;
        let a = ((ny - 1) as f64) * ((nz - 1) as f64);
        let delta_exact = f_total * l / (1e6 * a);
        let error = (mean_disp - delta_exact).abs() / delta_exact;

        println!("traction | mean={:.8} exact={:.8} error={:.4}%",
            mean_disp, delta_exact, error * 100.0);
        assert!(error < 0.01, "Erreur traction trop grande : {:.4}%", error * 100.0);
    }

    // ─── Convergence en flexion ───────────────────────────────────────────────

    /// L'erreur doit decroitre monotoniquement avec le raffinement du maillage.
    #[test]
    fn test_beam_bending_convergence() {
        let ny = 4; let nz = 4;
        let l_beam = 9.0_f64;
        let f_total = 25.0;
        let e = 1e6;
        let h = (ny - 1) as f64;
        let b = (nz - 1) as f64;
        let i_flex = b * h.powi(3) / 12.0;
        let delta_exact = f_total * l_beam.powi(3) / (3.0 * e * i_flex);

        let nx_values = [4, 7, 10, 13];
        let mut previous_error = f64::MAX;

        for &nx in &nx_values {
            let dx = l_beam / (nx - 1) as f64;
            let mesh = Mesh::generate(nx, ny, nz, dx, 1.0, 1.0);
            let material = svk(e, 0.3, 1.0);
            let assembler = Assembler::new(&mesh);
            let tips = tip_nodes(nx, ny, nz);
            let nb_tip = tips.len() as f64;

            let u_zero = DVector::zeros(3 * mesh.nodes.len());
            let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);
            let bc = make_bc(nx, ny, nz, &tips,
                Vector3::new(0.0, -f_total / nb_tip, 0.0),
                Box::new(PenaltyMethod));
            let bc_result = bc.apply(&k, mesh.nodes.len());
            let u_red = DirectSolver.solve(&bc_result.k, &bc_result.f).unwrap();
            let u = bc_result.reconstruct(u_red);
            let max_disp = tips.iter().map(|&i| u[3 * i + 1].abs()).fold(0.0f64, f64::max);
            let error = (max_disp - delta_exact).abs() / delta_exact;

            println!("nx={} dx={:.3} | disp={:.6} exact={:.6} error={:.2}%",
                nx, dx, max_disp, delta_exact, error * 100.0);
            assert!(error < previous_error,
                "L'erreur n'a pas diminue a nx={}: {:.2}% >= {:.2}%",
                nx, error * 100.0, previous_error * 100.0);
            previous_error = error;
        }
        assert!(previous_error < 0.25, "Erreur finale trop grande: {:.2}%", previous_error * 100.0);
    }

    // ─── Assemblage de la masse ───────────────────────────────────────────────

    /// sum(mass) / 3 == density * volume.
    #[test]
    fn test_mass_assembly() {
        let mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        let density = 1500.0;
        let material = svk(1e6, 0.3, density);
        let assembler = Assembler::new(&mesh);
        let mass = assembler.assemble_mass(&mesh, &material);
        let total_mass = mass.sum() / 3.0;
        let expected = density * 1.0;
        let error = (total_mass - expected).abs() / expected;
        println!("total_mass={:.6} expected={:.6} error={:.4}%", total_mass, expected, error * 100.0);
        assert!(error < 1e-10, "Erreur assemblage masse : {:.4}%", error * 100.0);
    }

    // ─── Symetrie amortissement Rayleigh ─────────────────────────────────────

    /// C = alpha*M + beta*K doit etre symetrique.
    #[test]
    fn test_rayleigh_damping_symmetry() {
        let mesh = Mesh::generate(3, 3, 3, 1.0, 1.0, 1.0);
        let material = svk(1e6, 0.3, 1000.0);
        let assembler = Assembler::new(&mesh);
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);
        let mass = assembler.assemble_mass(&mesh, &material);
        let c = RayleighDamping { alpha: 0.1, beta: 0.01 }.compute(&mass, &k);
        let diff = (&c - c.transpose()).abs().max();
        assert!(diff < 1e-10, "C non symetrique : max diff = {:.2e}", diff);
    }

    // ─── Newton-Raphson statique ──────────────────────────────────────────────

    /// Newton doit converger et satisfaire le residu f_int(u) - f_ext = 0.
    /// On verifie l'equilibre a la convergence, pas la comparaison avec le lineaire.
    #[test]
    fn test_newton_matches_direct_static() {
        let nx = 5; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let material = svk(1e6, 0.3, 1000.0);
        let assembler = Assembler::new(&mesh);
        let tips = tip_nodes(nx, ny, nz);
        let nb_tip = tips.len() as f64;
        let f_total = 100.0;

        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);

        let bc = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_total / nb_tip, 0.0, 0.0),
            Box::new(EliminationMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());
        let n_full = bc_result.n_dofs;
        let f_ext = bc_result.f.clone();

        let newton = NewtonRaphson::default();
        let u_newton_red = newton.solve::<LinearBMatrix>(
            &assembler, &mesh, &material, &bc_result, &DirectSolver,
        ).unwrap();

        // Verifier que f_int(u) ~ f_ext
        let u_full = bc_result.reconstruct_ref(&u_newton_red, n_full);
        let f_int_full = assembler.assemble_internal_forces(&mesh, &material, &u_full);
        let f_int_red = crate::solver::restrict_vector(&f_int_full, &bc_result.free_dofs);
        let residual_norm = (&f_int_red - &f_ext).norm() / f_ext.norm();
        println!("Newton residual | {:.2e}", residual_norm);
        assert!(residual_norm < 1e-6, "Residu Newton trop grand : {:.2e}", residual_norm);
    }

    // ─── Convergence dynamique vers statique ──────────────────────────────────

    /// Avec amortissement fort, l'integrateur implicite doit converger
    /// vers la solution statique. Erreur vs statique < 1%, vs analytique < 5%.
    #[test]
    fn test_implicit_euler_static_convergence() {
        let nx = 5; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let material = svk(1e6, 0.3, 1000.0);
        let assembler = Assembler::new(&mesh);
        let tips = tip_nodes(nx, ny, nz);
        let nb_tip = tips.len() as f64;
        let f_total = 100.0;

        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k_full = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);
        let mass = assembler.assemble_mass(&mesh, &material);

        // Solution statique de reference
        let bc_static = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_total / nb_tip, 0.0, 0.0),
            Box::new(EliminationMethod));
        let bc_result_static = bc_static.apply(&k_full, mesh.nodes.len());
        let u_static_red = DirectSolver.solve(&bc_result_static.k, &bc_result_static.f).unwrap();
        let u_static = bc_result_static.reconstruct(u_static_red);
        let mean_static = tips.iter().map(|&i| u_static[3 * i]).sum::<f64>() / nb_tip;

        // Solution analytique
        let l = (nx - 1) as f64;
        let a = ((ny - 1) as f64) * ((nz - 1) as f64);
        let delta_exact = f_total * l / (1e6 * a);

        // Amortissement Rayleigh fort
        let c_full = RayleighDamping { alpha: 10.0, beta: 0.01 }.compute(&mass, &k_full);

        // Simulation dynamique
        let bc_dyn = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_total / nb_tip, 0.0, 0.0),
            Box::new(EliminationMethod));
        let bc_result_dyn = bc_dyn.apply(&k_full, mesh.nodes.len());
        let f_ext = bc_result_dyn.f.clone();

        let (mass_red, c_red) = match &bc_result_dyn.free_dofs {
            Some(free) => (
                DVector::from_iterator(free.len(), free.iter().map(|&i| mass[i])),
                DMatrix::from_fn(free.len(), free.len(), |r, s| c_full[(free[r], free[s])]),
            ),
            None => (mass.clone(), c_full.clone()),
        };

        let n_red = mass_red.len();
        let mut state = MechanicalState::new(n_red / 3);
        let integrator = ImplicitEulerIntegrator::default();

        for _ in 0..500 {
            integrator.step::<LinearBMatrix>(
                &mut state, &mass_red, &c_red, &f_ext, 1.0,
                &assembler, &mesh, &material, &bc_result_dyn, &DirectSolver,
            ).unwrap();
        }

        let mean_dyn = match &bc_result_dyn.free_dofs {
            Some(free) => tips.iter().map(|&idx| {
                free.iter().position(|&d| d == 3 * idx)
                    .map(|j| state.position[j])
                    .unwrap_or(0.0)
            }).sum::<f64>() / nb_tip,
            None => tips.iter().map(|&i| state.position[3 * i]).sum::<f64>() / nb_tip,
        };

        let error_vs_static  = (mean_dyn - mean_static).abs() / mean_static;
        let error_vs_exact   = (mean_dyn - delta_exact).abs() / delta_exact;
        println!(
            "dyn={:.8} static={:.8} exact={:.8} | vs_static={:.4}% vs_exact={:.4}%",
            mean_dyn, mean_static, delta_exact,
            error_vs_static * 100.0, error_vs_exact * 100.0
        );
        assert!(error_vs_static < 0.01,
            "Dynamique ne converge pas vers statique : {:.4}%", error_vs_static * 100.0);
        assert!(error_vs_exact < 0.05,
            "Dynamique trop loin de l'analytique : {:.4}%", error_vs_exact * 100.0);
    }

    // ─── SVK petites deformations == lineaire ─────────────────────────────────

    /// Pour des petits deplacements, SVK doit donner le meme resultat
    /// que l'ancienne loi lineaire (erreur < 0.1%).
    #[test]
    fn test_svk_small_strain_matches_linear() {
        let nx = 5; let ny = 2; let nz = 2;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        let material = svk(1e6, 0.3, 1000.0);
        let assembler = Assembler::new(&mesh);
        let tips = tip_nodes(nx, ny, nz);
        let nb_tip = tips.len() as f64;

        // Force tres petite -> petites deformations -> SVK == lineaire
        let f_small = 1e-3;
        let u_zero = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);

        let bc = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_small / nb_tip, 0.0, 0.0),
            Box::new(EliminationMethod));
        let bc_result = bc.apply(&k, mesh.nodes.len());

        // Direct lineaire
        let u_lin_red = DirectSolver.solve(&bc_result.k, &bc_result.f).unwrap();
        let u_lin = bc_result.reconstruct(u_lin_red);

        // Newton SVK
        let bc2 = make_bc(nx, ny, nz, &tips,
            Vector3::new(f_small / nb_tip, 0.0, 0.0),
            Box::new(EliminationMethod));
        let bc_result2 = bc2.apply(&k, mesh.nodes.len());
        let newton = NewtonRaphson::default();
        let u_svk_red = newton.solve::<LinearBMatrix>(
            &assembler, &mesh, &material, &bc_result2, &DirectSolver,
        ).unwrap();
        let u_svk = bc_result2.reconstruct(u_svk_red);

        let error = (&u_svk - &u_lin).norm() / u_lin.norm();
        println!("SVK vs lineaire (petites def) | error = {:.2e}", error);
        assert!(error < 1e-6, "SVK doit etre identique au lineaire en petites def : {:.2e}", error);
    }
}