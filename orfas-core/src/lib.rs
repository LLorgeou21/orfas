


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
    
}