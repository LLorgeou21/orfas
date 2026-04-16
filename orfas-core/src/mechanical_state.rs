use nalgebra::DVector;

pub struct MechanicalState {
    pub position: DVector<f64>,
    pub velocity: DVector<f64>,
    pub acceleration: DVector<f64>,
}

impl MechanicalState {

    /// Creates a new MechanicalState with all vectors initialized to zero.
    /// n_nodes is the number of nodes in the mesh — vectors are of size 3*n_nodes (x, y, z per node).
    pub fn new(n_nodes: usize) -> MechanicalState {
        MechanicalState {
            position:     DVector::zeros(3 * n_nodes),
            velocity:     DVector::zeros(3 * n_nodes),
            acceleration: DVector::zeros(3 * n_nodes),
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