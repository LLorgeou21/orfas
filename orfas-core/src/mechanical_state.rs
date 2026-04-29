use nalgebra::DVector;

pub struct MechanicalState {
    pub position: DVector<f64>,
    pub velocity: DVector<f64>,
    pub acceleration: DVector<f64>,
}

impl MechanicalState {

    /// Creates a new MechanicalState with all vectors initialized to zero.
    // new() — remplacer la signature et les 3 occurrences de 3 *
    pub fn new(n_dof_per_node: usize, n_nodes: usize) -> MechanicalState {
        let total = n_dof_per_node * n_nodes;
        MechanicalState {
            position:     DVector::zeros(total),
            velocity:     DVector::zeros(total),
            acceleration: DVector::zeros(total),
        }
    }

    /// Vector operation: result = a + scalar * b
    /// Equivalent to SOFA's vOp — modifies result in place to avoid allocation.
    pub fn v_op(result: &mut DVector<f64>, a: &DVector<f64>, scalar: f64, b: &DVector<f64>) {
        *result = a + scalar * b;
    }

    /// Adds scalar * M * v to result, where M is the lumped mass vector (diagonal).
    /// Since M is diagonal, M * v reduces to a component-wise multiplication.
    /// Equivalent to SOFA's addMv — avoids constructing the full mass matrix.
    pub fn add_mv(result: &mut DVector<f64>, scalar: f64, mass: &DVector<f64>, v: &DVector<f64>) {
        *result += scalar * mass.component_mul(v);
    }
}