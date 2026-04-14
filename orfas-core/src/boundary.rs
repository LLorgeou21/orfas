use nalgebra::{DMatrix, DVector};

pub trait BoundaryConditionMethod {
    fn apply(&self, k: &mut DMatrix<f64>, f: &mut DVector<f64>, fixed_nodes: &[usize]);
}


/// Struct for the penalty method.
/// Doesn't store any state — the method only needs K, F and the fixed node indices.
/// Fixes a node by replacing its diagonal entries in K with a very large value (1e30),
/// which forces the displacement at that node to near zero.
/// This is a first-order approximation — it will be replaced by proper elimination in v0.2.
pub struct PenaltyMethod;

pub struct BoundaryConditions {
    pub fixed_nodes: Vec<usize>,
    pub method: Box<dyn BoundaryConditionMethod>,
}

impl BoundaryConditionMethod for PenaltyMethod {
    fn apply(&self, k: &mut DMatrix<f64>, f: &mut DVector<f64>, fixed_nodes: &[usize]) {
        for &node in fixed_nodes {
            for dir in 0..3 {
                let i = 3 * node + dir;
                k[(i,i)] = 1e30;
                f[i] = 0.0;
            }
        }
    }
}

impl BoundaryConditions {
    pub fn apply(&self, k: &mut DMatrix<f64>, f: &mut DVector<f64>) {
        self.method.apply(k, f, &self.fixed_nodes);
    }
    pub fn new(fixed_nodes: Vec<usize>, method: Box<dyn BoundaryConditionMethod>) -> Self {
    BoundaryConditions { fixed_nodes, method }
    }
}

