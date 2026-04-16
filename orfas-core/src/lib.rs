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
    use crate::material::LinearElastic;
    use crate::assembler::{Assembler, LinearBMatrix};
    use crate::boundary::{BoundaryConditionResult, BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load, PenaltyMethod};
    use crate::solver::{DirectSolver, Solver};
    use crate::damping::{RayleighDamping, DampingModel};
    use crate::integrator::{ImplicitEulerIntegrator, IntegratorMethod};
    use crate::mechanical_state::MechanicalState;
    use nalgebra::Vector3;

    #[test]
    fn test_validation_v0_1(){
        let mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        let material: LinearElastic = LinearElastic { youngs_modulus: 1000.0, poisson_ratio: 0.3, density: 1.0 };

        let constraint = Constraint {
            list: vec![FixedNode::all(0), FixedNode::all(1), FixedNode::all(2), FixedNode::all(3)]
        };
        let load = Load {
            list: vec![4, 5, 6, 7],
            force: Vector3::new(0.0, 0.0, -100.0),
        };
        let bc = BoundaryConditions::new(constraint, vec![load], Box::new(PenaltyMethod));

        let assembler = Assembler::new(&mesh);
        let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        let bc_result: BoundaryConditionResult = bc.apply(&k, mesh.nodes.len());
        let result = DirectSolver{}.solve(&bc_result.k, &bc_result.f);
        assert!(result.is_ok());
        let u = bc_result.reconstruct(result.unwrap());
        assert_eq!(u.len(), 3 * 8);
    }

    #[test]
    fn test_beam_bending_convergence() {
        let ny = 4;
        let nz = 4;
        let l_beam = 9.0_f64;
        let dy = 1.0;
        let dz = 1.0;
        let f_total = 25.0;
        let e = 1e6;

        let h = (ny - 1) as f64 * dy;
        let b = (nz - 1) as f64 * dz;
        let i_flex = b * h.powi(3) / 12.0;
        let delta_exact = f_total * l_beam.powi(3) / (3.0 * e * i_flex);

        let nx_values = [4, 7, 10, 13];
        let mut previous_error = f64::MAX;

        for &nx in &nx_values {
            let dx = l_beam / (nx - 1) as f64;
            let mesh = Mesh::generate(nx, ny, nz, dx, dy, dz);
            let material = LinearElastic { youngs_modulus: e, poisson_ratio: 0.3, density: 1.0 };

            let constraint = Constraint {
                list: (0..ny).flat_map(|j| {
                    (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
                }).collect()
            };

            let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
                (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
            }).collect();
            let nb_tip = tip_nodes.len() as f64;
            let load = Load {
                list: tip_nodes.clone(),
                force: Vector3::new(0.0, -f_total / nb_tip, 0.0),
            };

            let bc = BoundaryConditions::new(constraint, vec![load], Box::new(PenaltyMethod));
            let assembler = Assembler::new(&mesh);
            let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
            let bc_result: BoundaryConditionResult = bc.apply(&k, mesh.nodes.len());
            let result = DirectSolver{}.solve(&bc_result.k, &bc_result.f).unwrap();
            let u = bc_result.reconstruct(result);

            let max_disp = tip_nodes.iter()
                .map(|&idx| u[3 * idx + 1].abs())
                .fold(0.0f64, f64::max);

            let error = (max_disp - delta_exact).abs() / delta_exact;
            println!(
                "nx={} dx={:.3} | max_disp={:.6} | delta_exact={:.6} | error={:.2}%",
                nx, dx, max_disp, delta_exact, error * 100.0
            );

            assert!(
                error < previous_error,
                "Error did not decrease at nx={}: {:.2}% >= {:.2}%",
                nx, error * 100.0, previous_error * 100.0
            );
            previous_error = error;
        }

        assert!(previous_error < 0.25, "Final error too large: {:.2}%", previous_error * 100.0);
    }



    #[test]
    fn test_axial_traction() {
        let nx = 10;
        let ny = 3;
        let nz = 3;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let mesh = Mesh::generate(nx, ny, nz, dx, dy, dz);
        let material = LinearElastic { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1.0 };

        let constraint = Constraint {
            list: (0..ny).flat_map(|j| {
                (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
            }).collect()
        };

        let f_total = 100.0;
        let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();
        let nb_tip = tip_nodes.len() as f64;
        let load = Load {
            list: tip_nodes.clone(),
            force: Vector3::new(f_total / nb_tip, 0.0, 0.0),
        };

        let bc = BoundaryConditions::new(constraint, vec![load], Box::new(PenaltyMethod));
        let assembler = Assembler::new(&mesh);
        let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        let bc_result: BoundaryConditionResult = bc.apply(&k, mesh.nodes.len());
        let result = DirectSolver{}.solve(&bc_result.k, &bc_result.f).unwrap();
        let u = bc_result.reconstruct(result);

        let mean_disp = tip_nodes.iter()
            .map(|&idx| u[3 * idx])
            .sum::<f64>() / nb_tip;

        let l = (nx - 1) as f64 * dx;
        let a = ((ny - 1) as f64 * dy) * ((nz - 1) as f64 * dz);
        let delta_exact = f_total * l / (1e6 * a);

        let error = (mean_disp - delta_exact).abs() / delta_exact;
        println!(
            "traction | L={} A={} | mean_disp={:.8} | delta_exact={:.8} | error={:.4}%",
            l, a, mean_disp, delta_exact, error * 100.0
        );

        assert!(error < 0.01, "Error too large: {:.4}%", error * 100.0);
    }

    #[test]
    fn test_mass_assembly() {
        // Total assembled mass must equal density * total_volume
        // For a unit cube (1x1x1): volume = 1, total_mass = density * 1
        // sum(mass_vector) / 3 == density * volume (factor 3 because 3 DOFs per node)
        let mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        let density = 1500.0;
        let material = LinearElastic { youngs_modulus: 1e6, poisson_ratio: 0.3, density };

        let assembler = Assembler::new(&mesh);
        let mass = assembler.assemble_mass(&mesh, &material);

        let total_mass = mass.sum() / 3.0;
        let expected_mass = density * 1.0; // density * volume of unit cube

        let error = (total_mass - expected_mass).abs() / expected_mass;
        println!("total_mass={:.6} expected={:.6} error={:.4}%", total_mass, expected_mass, error * 100.0);
        assert!(error < 1e-10, "Mass assembly error too large: {:.4}%", error * 100.0);
    }

    #[test]
    fn test_rayleigh_damping_symmetry() {
        // C = alpha*M + beta*K must be symmetric since both M and K are symmetric
        let mesh = Mesh::generate(3, 3, 3, 1.0, 1.0, 1.0);
        let material = LinearElastic { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 };

        let assembler = Assembler::new(&mesh);
        let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        let mass = assembler.assemble_mass(&mesh, &material);

        let damping = RayleighDamping { alpha: 0.1, beta: 0.01 };
        let c = damping.compute(&mass, &k);

        let diff = (&c - c.transpose()).abs().max();
        assert!(diff < 1e-10, "Damping matrix is not symmetric: max diff = {:.2e}", diff);
    }

    #[test]
    fn test_implicit_euler_static_convergence() {
        // A bar under axial traction, simulated dynamically with heavy damping,
        // must converge to the static solution and to the analytical solution.
        // EliminationMethod is used to avoid ill-conditioning from the penalty method.
        let nx = 5;
        let ny = 2;
        let nz = 2;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let mesh = Mesh::generate(nx, ny, nz, dx, dy, dz);
        let material = LinearElastic { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 };

        let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();
        let nb_tip = tip_nodes.len() as f64;
        let f_total = 100.0;

        let make_bc = || {
            let constraint = Constraint {
                list: (0..ny).flat_map(|j| {
                    (0..nz).map(move |k| FixedNode::all(j * nx + k * nx * ny))
                }).collect()
            };
            let load = Load {
                list: tip_nodes.clone(),
                force: Vector3::new(f_total / nb_tip, 0.0, 0.0),
            };
            BoundaryConditions::new(constraint, vec![load], Box::new(crate::boundary::EliminationMethod))
        };

        let assembler = Assembler::new(&mesh);
        let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        let mass = assembler.assemble_mass(&mesh, &material);

        // Compute static reference
        let bc_static = make_bc();
        let bc_result_static = bc_static.apply(&k, mesh.nodes.len());
        let result_static = DirectSolver{}.solve(&bc_result_static.k, &bc_result_static.f).unwrap();
        let u_static = bc_result_static.reconstruct(result_static);
        let mean_disp_static = tip_nodes.iter()
            .map(|&idx| u_static[3 * idx])
            .sum::<f64>() / nb_tip;

        // Analytical solution: delta = F*L / (E*A)
        let l = (nx - 1) as f64 * dx;
        let a = ((ny - 1) as f64 * dy) * ((nz - 1) as f64 * dz);
        let delta_exact = f_total * l / (1e6 * a);

        // Dynamic simulation with heavy damping
        let damping = RayleighDamping { alpha: 10.0, beta: 0.01 };
        let c = damping.compute(&mass, &k);

        let bc_dynamic = make_bc();
        let bc_result_dynamic = bc_dynamic.apply(&k, mesh.nodes.len());
        let f = bc_result_dynamic.f.clone();

        // Reduce mass and c to match the reduced system from EliminationMethod
        let (mass_reduced, c_reduced) = match &bc_result_dynamic.free_dofs {
            Some(free_dofs) => {
                let m = nalgebra::DVector::from_iterator(
                    free_dofs.len(),
                    free_dofs.iter().map(|&i| mass[i])
                );
                let c = nalgebra::DMatrix::from_fn(free_dofs.len(), free_dofs.len(), |r, s| {
                    c[(free_dofs[r], free_dofs[s])]
                });
                (m, c)
            }
            None => (mass.clone(), c),
        };

        let n_reduced = mass_reduced.len();
        let mut state = MechanicalState::new(n_reduced / 3);
        let integrator = ImplicitEulerIntegrator;
        let solver = DirectSolver {};
        let dt = 1.0;

        for _ in 0..500 {
            integrator.step(&mut state, &mass_reduced, &bc_result_dynamic.k, &c_reduced, &f, dt, &solver).unwrap();
        }

        // Reconstruct mean x displacement on tip face
        let mean_disp_dynamic = match &bc_result_dynamic.free_dofs {
            Some(free_dofs) => {
                tip_nodes.iter().map(|&idx| {
                    free_dofs.iter().position(|&d| d == 3 * idx)
                        .map(|j| state.position[j])
                        .unwrap_or(0.0)
                }).sum::<f64>() / nb_tip
            }
            None => tip_nodes.iter().map(|&idx| state.position[3 * idx]).sum::<f64>() / nb_tip,
        };

        let error_vs_static = (mean_disp_dynamic - mean_disp_static).abs() / mean_disp_static;
        let error_vs_exact  = (mean_disp_dynamic - delta_exact).abs() / delta_exact;

        println!(
            "dynamic convergence | mean_disp_dynamic={:.8} | mean_disp_static={:.8} | delta_exact={:.8} | error_vs_static={:.4}% | error_vs_exact={:.4}%",
            mean_disp_dynamic, mean_disp_static, delta_exact, error_vs_static * 100.0, error_vs_exact * 100.0
        );

        assert!(error_vs_static < 0.01, "Dynamic did not converge to static: {:.4}%", error_vs_static * 100.0);
        assert!(error_vs_exact  < 0.05, "Dynamic solution too far from analytical: {:.4}%", error_vs_exact * 100.0);
    }



}