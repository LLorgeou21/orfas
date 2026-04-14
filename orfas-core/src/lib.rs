pub mod mesh;
pub mod material;
pub mod assembler;
pub mod boundary;
pub mod solver;



#[cfg(test)]
mod integration_tests {
    use crate::mesh::Mesh;
    use crate::material::LinearElastic;
    use crate::assembler::{Assembler, LinearBMatrix};
    use crate::boundary::{BoundaryConditions, PenaltyMethod};
    use crate::solver::{DirectSolver, Solver};
    use nalgebra::DVector;

    #[test]
    fn test_validation_v0_1(){
        let mut mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        let material: LinearElastic = LinearElastic { youngs_modulus: 1000.0, poisson_ratio: 0.3 };
        for i in 0..2 {
            for j in 0..2 {
                let idx = i + j * 2;
                mesh.nodes[idx].fixed = true;
            }
        }
        let fixed_nodes = vec![0, 1, 2, 3];
        let bc = BoundaryConditions::new(fixed_nodes, Box::new(PenaltyMethod));
        let n = 3 * mesh.nodes.len();
        let mut f = DVector::zeros(n);

        let force = -100.0;
        for node_idx in [4, 5, 6, 7] {
            f[3 * node_idx + 2] = force;
        }
        let assembler = Assembler::new(&mesh);
        let mut k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        bc.apply(&mut k, &mut f);
        let solver = DirectSolver{};
        let result = solver.solve(&k, &f);
        assert!(result.is_ok());
    }

    #[test]
    fn test_beam_bending_convergence() {
        // Convergence test for a cantilever beam under tip load (bending in xy plane)
        // Geometry is fixed: L=9, section 3x3 (h=b=3)
        // nx is refined: 4, 7, 10, 13 nodes along x (dx shrinks accordingly)
        // This refines the same problem, so the error must strictly decrease
        // Linear tetrahedra converge in O(h) in bending (with shear locking)
        let ny = 4;
        let nz = 4;
        let l_beam = 9.0_f64;
        let dy = 1.0;
        let dz = 1.0;
        let f_total = 25.0;
        let e = 1e6;

        // Fixed geometry: section height h=(ny-1)*dy, width b=(nz-1)*dz
        let h = (ny - 1) as f64 * dy;
        let b = (nz - 1) as f64 * dz;
        let i_flex = b * h.powi(3) / 12.0;
        let delta_exact = f_total * l_beam.powi(3) / (3.0 * e * i_flex);

        let nx_values = [4, 7, 10, 13];
        let mut previous_error = f64::MAX;

        for &nx in &nx_values {
            // dx is adjusted so that (nx-1)*dx = L always
            let dx = l_beam / (nx - 1) as f64;

            let mesh = Mesh::generate(nx, ny, nz, dx, dy, dz);
            let material = LinearElastic { youngs_modulus: e, poisson_ratio: 0.3 };

            let fixed_nodes: Vec<usize> = (0..ny).flat_map(|j| {
                (0..nz).map(move |k| j * nx + k * nx * ny)
            }).collect();
            let bc = BoundaryConditions::new(fixed_nodes, Box::new(PenaltyMethod));

            let n = 3 * mesh.nodes.len();
            let mut f = DVector::zeros(n);

            let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
                (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
            }).collect();
            let nb_tip = tip_nodes.len() as f64;
            for &idx in &tip_nodes {
                f[3 * idx + 1] = -f_total / nb_tip;
            }

            let assembler = Assembler::new(&mesh);
            let mut k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
            bc.apply(&mut k, &mut f);

            let u = DirectSolver {}.solve(&k, &f).unwrap();

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

        // Final refinement must be within 25% -- realistic for linear tetrahedra with locking
        assert!(previous_error < 0.25, "Final error too large: {:.2}%", previous_error * 100.0);
    }

    #[test]
    fn test_axial_traction() {
        // Pure axial traction in x: no bending, no locking
        // Exact solution: delta = F * L / (E * A)
        // Expected error < 1%
        let nx = 10;
        let ny = 3;
        let nz = 3;
        let dx = 1.0;
        let dy = 1.0;
        let dz = 1.0;

        let mesh = Mesh::generate(nx, ny, nz, dx, dy, dz);
        let material = LinearElastic { youngs_modulus: 1e6, poisson_ratio: 0.3 };

        // Fix all nodes on the x=0 face
        let fixed_nodes: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| j * nx + k * nx * ny)
        }).collect();
        let bc = BoundaryConditions::new(fixed_nodes, Box::new(PenaltyMethod));

        let n = 3 * mesh.nodes.len();
        let mut f = DVector::zeros(n);

        // Total force F=100 in x direction on the tip face (x=L)
        let f_total = 100.0;
        let tip_nodes: Vec<usize> = (0..ny).flat_map(|j| {
            (0..nz).map(move |k| (nx - 1) + j * nx + k * nx * ny)
        }).collect();
        let nb_tip = tip_nodes.len() as f64;
        for &idx in &tip_nodes {
            f[3 * idx] = f_total / nb_tip; // x component = index 0
        }

        let assembler = Assembler::new(&mesh);
        let mut k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
        bc.apply(&mut k, &mut f);

        let solver = DirectSolver {};
        let u = solver.solve(&k, &f).unwrap();

        // Mean x displacement on the tip face
        let mean_disp = tip_nodes.iter()
            .map(|&idx| u[3 * idx])
            .sum::<f64>() / nb_tip;

        // delta = F * L / (E * A)
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

}