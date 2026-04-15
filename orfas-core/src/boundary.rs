use nalgebra::{DMatrix, DVector};
use std::collections::HashMap;

// ─── Data structures ──────────────────────────────────────────────────────────

/// Stores a fixed node index and which displacement directions are constrained.
/// Use the constructors (all, only_x, only_y, only_z) for convenience.
pub struct FixedNode {
    pub indice: usize,
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl FixedNode {

    /// Fix all three displacement directions (x, y, z).
    pub fn all(indice: usize) -> FixedNode {
        FixedNode { indice, x: true, y: true, z: true }
    }

    /// Fix only the x displacement direction.
    pub fn only_x(indice: usize) -> FixedNode {
        FixedNode { indice, x: true, y: false, z: false }
    }

    /// Fix only the y displacement direction.
    pub fn only_y(indice: usize) -> FixedNode {
        FixedNode { indice, x: false, y: true, z: false }
    }

    /// Fix only the z displacement direction.
    pub fn only_z(indice: usize) -> FixedNode {
        FixedNode { indice, x: false, y: false, z: true }
    }
}

/// Result of applying boundary conditions to the system K·u = f.
/// Contains the (possibly reduced) stiffness matrix and force vector,
/// along with the information needed to reconstruct the full displacement vector.
pub struct BoundaryConditionResult {
    /// Modified or reduced stiffness matrix
    pub k: DMatrix<f64>,
    /// Modified or reduced force vector
    pub f: DVector<f64>,
    /// Indices of free DOFs in the original system.
    /// None for penalty method (system size unchanged),
    /// Some for elimination method (system is reduced).
    pub free_dofs: Option<Vec<usize>>,
    /// Size of the original (unreduced) system — needed to reconstruct u.
    pub n_dofs: usize,
}

impl BoundaryConditionResult {

    /// Reconstructs the full displacement vector from the solved (possibly reduced) system.
    /// For the penalty method (free_dofs = None), returns u unchanged.
    /// For the elimination method, inserts zeros at fixed DOF positions.
    pub fn reconstruct(self, u_reduced: DVector<f64>) -> DVector<f64> {
        match self.free_dofs {
            Some(free_dofs) => {
                let mut u_complete = DVector::zeros(self.n_dofs);
                for i in 0..free_dofs.len() {
                    u_complete[free_dofs[i]] = u_reduced[i];
                }
                u_complete
            }
            None => u_reduced
        }
    }
}

/// Holds the fixed nodes and the chosen boundary condition method.
pub struct BoundaryConditions {
    pub fixed_nodes: Vec<FixedNode>,
    pub method: Box<dyn BoundaryConditionMethod>,
}

impl BoundaryConditions {

    pub fn new(fixed_nodes: Vec<FixedNode>, method: Box<dyn BoundaryConditionMethod>) -> Self {
        BoundaryConditions { fixed_nodes, method }
    }

    pub fn apply(&self, k: &DMatrix<f64>, f: &DVector<f64>) -> BoundaryConditionResult {
        self.method.apply(k, f, &self.fixed_nodes)
    }
}

// ─── Trait ────────────────────────────────────────────────────────────────────

/// Defines how boundary conditions are applied to the system K·u = f.
/// Returns a BoundaryConditionResult containing the modified system
/// and the information needed to reconstruct the full displacement vector.
pub trait BoundaryConditionMethod {
    fn apply(&self, k: &DMatrix<f64>, f: &DVector<f64>, fixed_nodes: &[FixedNode]) -> BoundaryConditionResult;
}

// ─── Implementations ──────────────────────────────────────────────────────────

/// Penalty method: enforces zero displacement by replacing diagonal entries
/// of K with a very large value (1e30), forcing near-zero displacement.
/// Does not change the size of the system.
/// Simple but introduces numerical ill-conditioning proportional to the penalty value.
pub struct PenaltyMethod;

impl BoundaryConditionMethod for PenaltyMethod {
    fn apply(&self, k: &DMatrix<f64>, f: &DVector<f64>, fixed_nodes: &[FixedNode]) -> BoundaryConditionResult {
        let mut new_k = k.clone();
        let mut new_f = f.clone();
        for node in fixed_nodes {
            let i = 3 * node.indice;
            if node.x { new_k[(i,   i  )] = 1e30; new_f[i  ] = 0.0; }
            if node.y { new_k[(i+1, i+1)] = 1e30; new_f[i+1] = 0.0; }
            if node.z { new_k[(i+2, i+2)] = 1e30; new_f[i+2] = 0.0; }
        }
        BoundaryConditionResult { k: new_k, f: new_f, free_dofs: None, n_dofs: f.len() }
    }
}

/// Elimination method: removes fixed DOF rows and columns from K and f,
/// solving a smaller system. More numerically stable than the penalty method.
/// The full displacement vector is reconstructed after solving via BoundaryConditionResult::reconstruct.
pub struct EliminationMethod;

impl BoundaryConditionMethod for EliminationMethod {
    fn apply(&self, k: &DMatrix<f64>, f: &DVector<f64>, fixed_nodes: &[FixedNode]) -> BoundaryConditionResult {
        // Build a map from node index to blocked directions for O(1) lookup
        let mut blocked: HashMap<usize, (bool, bool, bool)> = HashMap::new();
        for node in fixed_nodes {
            blocked.insert(node.indice, (node.x, node.y, node.z));
        }

        // Collect the indices of all free DOFs
        let mut free_dofs = Vec::new();
        for i in 0..f.len() / 3 {
            match blocked.get(&i) {
                Some((x, y, z)) => {
                    if !*x { free_dofs.push(3 * i);     }
                    if !*y { free_dofs.push(3 * i + 1); }
                    if !*z { free_dofs.push(3 * i + 2); }
                }
                None => {
                    free_dofs.push(3 * i);
                    free_dofs.push(3 * i + 1);
                    free_dofs.push(3 * i + 2);
                }
            }
        }

        // Build the reduced system by extracting free DOF rows and columns
        let f_reduce = DVector::from_iterator(
            free_dofs.len(),
            free_dofs.iter().map(|&i| f[i])
        );
        let k_reduce = DMatrix::from_iterator(
            free_dofs.len(),
            free_dofs.len(),
            free_dofs.iter().flat_map(|&col| free_dofs.iter().map(move |&row| k[(row, col)]))
        );

        BoundaryConditionResult { k: k_reduce, f: f_reduce, free_dofs: Some(free_dofs), n_dofs: f.len() }
    }
}