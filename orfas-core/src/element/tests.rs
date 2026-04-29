// UTF-8
// orfas-core/src/element/tests.rs
// Generic element validation tests and numerical benchmarks.
// Designed to be reused for any E: FiniteElement added in future versions.
//
// Two categories:
// 1. unit_checks  — mathematical properties every element must satisfy
//                   (partition of unity, volume, gradient consistency)
// 2. bending      — cantilever beam benchmark against Euler-Bernoulli analytical solution
//                   delta = F * L^3 / (3 * E * I)

#[cfg(test)]
mod unit_checks {
    use crate::element::{FiniteElement, Tet4, Tet10};
    use crate::element::subdivision::tet4_to_tet10;
    use crate::mesh::Tet4Mesh;
    use nalgebra::Vector3;

    /// Helper: build a single unit tetrahedron as a Tet4Mesh (1 element).
    fn unit_tet4_mesh() -> Tet4Mesh {
        Tet4Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0)
    }

    /// Check that shape functions sum to 1 at several interior points.
    /// Valid for any element type — tests partition of unity.
    fn check_partition_of_unity<E: FiniteElement>() {
        let test_points = vec![
            Vector3::new(0.25, 0.25, 0.25),
            Vector3::new(0.1,  0.2,  0.1),
            Vector3::new(0.5,  0.1,  0.1),
        ];
        for xi in &test_points {
            let n   = E::shape_functions(xi);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12,
                "{}: partition of unity failed at {:?}: sum = {}",
                std::any::type_name::<E>(), xi, sum);
        }
    }

    /// Check that gradient sum over all nodes is zero at every Gauss point.
    /// Consequence of partition of unity: sum_i grad(Ni) = grad(1) = 0.
    fn check_gradient_sum_zero<E: FiniteElement>(nodes: &[nalgebra::Vector3<f64>]) {
        let geo = E::precompute(nodes);
        let n_gauss = E::integration_points().len();
        for g in 0..n_gauss {
            let grad = E::shape_gradients(&geo, g);
            for col in 0..3 {
                let sum: f64 = (0..E::N_NODES).map(|i| grad[(i, col)]).sum();
                assert!(sum.abs() < 1e-10,
                    "{}: gradient sum col {} at Gauss {} should be 0, got {}",
                    std::any::type_name::<E>(), col, g, sum);
            }
        }
    }

    /// Check that element volume matches the expected value.
    fn check_volume<E: FiniteElement>(nodes: &[Vector3<f64>], expected: f64) {
        let geo = E::precompute(nodes);
        let vol = E::element_volume(&geo);
        assert!((vol - expected).abs() / expected < 1e-10,
            "{}: volume should be {}, got {}",
            std::any::type_name::<E>(), expected, vol);
    }

    // -------------------------------------------------------------------------
    // Tet4 unit checks
    // -------------------------------------------------------------------------

    #[test]
    fn tet4_partition_of_unity() {
        check_partition_of_unity::<Tet4>();
    }

    #[test]
    fn tet4_gradient_sum_zero() {
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        check_gradient_sum_zero::<Tet4>(&nodes);
    }

    #[test]
    fn tet4_volume() {
        let nodes = vec![
            Vector3::new(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        check_volume::<Tet4>(&nodes, 1.0 / 6.0);
    }

    // -------------------------------------------------------------------------
    // Tet10 unit checks
    // -------------------------------------------------------------------------

    fn unit_tet10_nodes() -> Vec<Vector3<f64>> {
        vec![
            Vector3::new(0.0, 0.0, 0.0), // 0 corner
            Vector3::new(1.0, 0.0, 0.0), // 1 corner
            Vector3::new(0.0, 1.0, 0.0), // 2 corner
            Vector3::new(0.0, 0.0, 1.0), // 3 corner
            Vector3::new(0.0, 0.5, 0.5), // 4 mid(2,3)
            Vector3::new(0.5, 0.0, 0.5), // 5 mid(1,3)
            Vector3::new(0.5, 0.5, 0.0), // 6 mid(1,2)
            Vector3::new(0.0, 0.0, 0.5), // 7 mid(0,3)
            Vector3::new(0.0, 0.5, 0.0), // 8 mid(0,2)
            Vector3::new(0.5, 0.0, 0.0), // 9 mid(0,1)
        ]
    }

    #[test]
    fn tet10_partition_of_unity() {
        check_partition_of_unity::<Tet10>();
    }

    #[test]
    fn tet10_gradient_sum_zero() {
        check_gradient_sum_zero::<Tet10>(&unit_tet10_nodes());
    }

    #[test]
    fn tet10_volume() {
        check_volume::<Tet10>(&unit_tet10_nodes(), 1.0 / 6.0);
    }
}

// ---------------------------------------------------------------------------
// Cantilever beam bending benchmark
// ---------------------------------------------------------------------------
// Reference: Euler-Bernoulli beam theory
//   delta = F * L^3 / (3 * E * I)
//   I = b * h^3 / 12  (second moment of area)
//
// Expected behavior:
//   Tet4  — large error (~20-50%) due to shear locking
//   Tet10 — small error (< 5%) due to quadratic displacement field
//
// This test documents the quantitative improvement from Tet4 to Tet10
// and serves as the primary validation for v0.8.0.

#[cfg(test)]
mod bending {
    use crate::assembler::Assembler;
    use crate::boundary::{BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load};
    use crate::element::{FiniteElement, Tet4, Tet10};
    use crate::element::subdivision::tet4_to_tet10;
    use crate::element::traits::FemMesh;
    use crate::material::{SaintVenantKirchhoff, SimulationContext};
    use crate::mesh::Tet4Mesh;
    use crate::solver::{DirectSolver, DenseSolver};
    use nalgebra::{DVector, Vector3};
    use crate::element::traits::DofType;
    const NDOF: usize = <Tet4 as FiniteElement>::Dof::N_DOF;

    /// Run the cantilever beam benchmark for any element type.
    /// Returns the relative error vs the Euler-Bernoulli analytical solution.
    ///
    /// Geometry: beam of length L = (nx-1)*dx, cross-section h x b.
    /// Boundary: clamped at x=0, tip load F in -y direction at x=L.
    /// Reference: delta = F * L^3 / (3 * E * I), I = b * h^3 / 12.
    fn cantilever_error<E: FiniteElement>(
        fem_mesh:  &impl FemMesh,
        nx: usize, ny: usize, nz: usize,
        dx: f64,
        f_total: f64,
        e: f64,
    ) -> (f64, f64, f64)  {
        let assembler = Assembler::<E>::new(fem_mesh);
        let material  = SaintVenantKirchhoff {
            youngs_modulus: e,
            poisson_ratio:  0.3,
            density:        1.0,
        };

        // Tip nodes: x = (nx-1)*dx, all y and z
        let x_tip = (nx - 1) as f64 * dx;
        let tips: Vec<usize> = (0..fem_mesh.n_nodes())
            .filter(|&i| (fem_mesh.nodes()[i].position.x - x_tip).abs() < 1e-10)
            .collect();
        let nb_tip = tips.len() as f64;

        // Map Tet4 tip node indices to the FEM mesh node indices.
        // For Tet4: indices are identical.
        // For Tet10: corner nodes keep their Tet4 indices (subdivision preserves them).
        let u_zero   = DVector::zeros(3 * fem_mesh.n_nodes());
        let sim_ctx  = SimulationContext::isotropic_static(fem_mesh.n_elements());
        let k        = assembler.assemble_tangent(fem_mesh, &material, &u_zero, &sim_ctx);

        // Boundary conditions: clamp x=0 face, apply -y force at tip
        let constraint = Constraint {
            list: (0..fem_mesh.n_nodes())
                .filter(|&i| fem_mesh.nodes()[i].position.x.abs() < 1e-10)
                .map(|i| FixedNode::all(i))
                .collect(),
        };
        let load = Load {
            list:  tips.clone(),
            force: Vector3::new(0.0, -f_total / nb_tip, 0.0),
        };
        let bc        = BoundaryConditions::new(constraint, vec![load], Box::new(EliminationMethod));
        let bc_result = bc.apply(&k, fem_mesh.n_nodes(),NDOF);
        let u_red     = DirectSolver.solve(&bc_result.k, &bc_result.f).unwrap();
        let u_full    = bc_result.reconstruct(u_red);

        // Max tip displacement in y
        let max_disp = tips.iter()
            .map(|&i| u_full[3 * i + 1].abs())
            .fold(0.0_f64, f64::max);

        // Euler-Bernoulli analytical solution
        let l       = (nx - 1) as f64 * dx;
        let h       = (ny - 1) as f64;
        let b       = (nz - 1) as f64;
        let i_flex  = b * h.powi(3) / 12.0;
        let delta   = f_total * l.powi(3) / (3.0 * e * i_flex);

        ((max_disp - delta).abs() / delta, max_disp, delta)
    }

    /// Tet4 must show significant locking error (> 10%) on coarse mesh.
    /// This documents the known limitation that Tet10 is designed to fix.
    #[test]
    fn tet4_bending_has_locking() {
        let nx = 7; let ny = 4; let nz = 4;
        let dx = 1.0; let f_total = 25.0; let e = 1e6;

        let mesh = Tet4Mesh::generate(nx, ny, nz, dx, 1.0, 1.0);
        let (error,a,b) = cantilever_error::<Tet4>(&mesh, nx, ny, nz, dx, f_total, e);

        println!("Tet4 bending error: {:.2}%", error * 100.0);
        assert!(error > 0.10,
            "Tet4 should show locking (>10% error), got {:.2}%", error * 100.0);
    }


    /// Side-by-side comparison printed for documentation.
    #[test]
    fn tet4_vs_tet10_bending_comparison() {
        let nx = 7; let ny = 4; let nz = 4;
        let dx = 1.0; let f_total = 25.0; let e = 1e6;

        let tet4_mesh  = Tet4Mesh::generate(nx, ny, nz, dx, 1.0, 1.0);
        let tet10_mesh = tet4_to_tet10(&tet4_mesh);

        let (err4,  disp4,  delta) = cantilever_error::<Tet4>(&tet4_mesh,  nx, ny, nz, dx, f_total, e);
        let (err10, disp10, _)     = cantilever_error::<Tet10>(&tet10_mesh, nx, ny, nz, dx, f_total, e);

        println!("--- Cantilever bending benchmark (nx={} ny={} nz={}) ---", nx, ny, nz);
        println!("delta_exact:    {:.8}", delta);
        println!("Tet4  max_disp: {:.8}  error: {:.2}%", disp4,  err4  * 100.0);
        println!("Tet10 max_disp: {:.8}  error: {:.2}%", disp10, err10 * 100.0);
    }





}
