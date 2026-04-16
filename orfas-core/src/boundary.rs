use nalgebra::{DMatrix, DVector, Vector3};
use std::collections::HashMap;

// ─── Data structures ──────────────────────────────────────────────────────────

/// Stores a fixed node index and which displacement directions are constrained.
/// Use the constructors (all, only_x, only_y, only_z) for convenience.
#[derive(Clone)]
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

/// Result of applying boundary conditions to the system K*u = f.
/// Contains the (possibly reduced) stiffness matrix and force vector,
/// along with the information needed to reconstruct the full displacement vector.
pub struct BoundaryConditionResult {
    /// Modified or reduced stiffness matrix.
    pub k: DMatrix<f64>,
    /// Modified or reduced force vector.
    pub f: DVector<f64>,
    /// Indices of free DOFs in the original system.
    /// None for penalty method (system size unchanged).
    /// Some for elimination method (system is reduced).
    pub free_dofs: Option<Vec<usize>>,
    /// Size of the original (unreduced) system — needed to reconstruct u.
    pub n_dofs: usize,
}

impl BoundaryConditionResult {

    /// Reconstructs the full displacement vector, consuming self.
    /// For the penalty method (free_dofs = None), returns u_reduced unchanged.
    /// For the elimination method, inserts zeros at fixed DOF positions.
    pub fn reconstruct(self, u_reduced: DVector<f64>) -> DVector<f64> {
        match self.free_dofs {
            Some(free_dofs) => {
                let mut u_full = DVector::zeros(self.n_dofs);
                for (i, &dof) in free_dofs.iter().enumerate() {
                    u_full[dof] = u_reduced[i];
                }
                u_full
            }
            None => u_reduced,
        }
    }

    /// Reconstructs the full displacement vector without consuming self.
    /// Used inside Newton-Raphson and implicit Euler loops where bc_result
    /// must remain available across multiple iterations.
    /// n_full is passed explicitly since self.n_dofs is not always set by
    /// the penalty method (which keeps the full system size).
    pub fn reconstruct_ref(&self, u_reduced: &DVector<f64>, n_full: usize) -> DVector<f64> {
        match &self.free_dofs {
            Some(free_dofs) => {
                let mut u_full = DVector::zeros(n_full);
                for (i, &dof) in free_dofs.iter().enumerate() {
                    u_full[dof] = u_reduced[i];
                }
                u_full
            }
            None => u_reduced.clone(),
        }
    }
}

/// Holds the constraint, loads, and the chosen boundary condition method.
pub struct BoundaryConditions {
    pub constraint: Constraint,
    pub loads:      Vec<Load>,
    pub method:     Box<dyn BoundaryConditionMethod>,
}

impl BoundaryConditions {

    pub fn new(
        constraint: Constraint,
        loads:      Vec<Load>,
        method:     Box<dyn BoundaryConditionMethod>,
    ) -> Self {
        BoundaryConditions { constraint, loads, method }
    }

    /// Builds f from loads and applies the boundary condition method.
    /// k is only used by the penalty method — elimination ignores it.
    pub fn apply(&self, k: &DMatrix<f64>, n_nodes: usize) -> BoundaryConditionResult {
        let mut f = DVector::zeros(3 * n_nodes);
        for load in &self.loads {
            for &idx in &load.list {
                f[3 * idx]     += load.force.x;
                f[3 * idx + 1] += load.force.y;
                f[3 * idx + 2] += load.force.z;
            }
        }
        self.method.apply(k, &f, &self.constraint)
    }
}

#[derive(Clone)]
pub struct Constraint {
    pub list: Vec<FixedNode>,
}

#[derive(Clone)]
pub struct Load {
    pub list:  Vec<usize>,
    pub force: Vector3<f64>,
}

// ─── Trait ────────────────────────────────────────────────────────────────────

/// Defines how boundary conditions are applied to the system K*u = f.
/// Returns a BoundaryConditionResult containing the modified system
/// and the information needed to reconstruct the full displacement vector.
pub trait BoundaryConditionMethod {
    fn apply(
        &self,
        k:          &DMatrix<f64>,
        f:          &DVector<f64>,
        constraint: &Constraint,
    ) -> BoundaryConditionResult;
}

// ─── Implementations ──────────────────────────────────────────────────────────

/// Penalty method: enforces zero displacement by replacing diagonal entries
/// of K with a very large value (1e30), forcing near-zero displacement.
/// Does not change the size of the system.
/// Simple but introduces numerical ill-conditioning proportional to the penalty value.
pub struct PenaltyMethod;

impl BoundaryConditionMethod for PenaltyMethod {
    fn apply(
        &self,
        k:          &DMatrix<f64>,
        f:          &DVector<f64>,
        constraint: &Constraint,
    ) -> BoundaryConditionResult {
        let mut new_k = k.clone();
        let mut new_f = f.clone();
        for node in &constraint.list {
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
/// The full displacement vector is reconstructed after solving via reconstruct / reconstruct_ref.
pub struct EliminationMethod;

impl BoundaryConditionMethod for EliminationMethod {
    fn apply(
        &self,
        k:          &DMatrix<f64>,
        f:          &DVector<f64>,
        constraint: &Constraint,
    ) -> BoundaryConditionResult {
        // Build a map from node index to blocked directions for O(1) lookup
        let mut blocked: HashMap<usize, (bool, bool, bool)> = HashMap::new();
        for node in &constraint.list {
            blocked.insert(node.indice, (node.x, node.y, node.z));
        }

        // Collect the indices of all free DOFs
        let mut free_dofs = Vec::new();
        for i in 0..f.len() / 3 {
            match blocked.get(&i) {
                Some((bx, by, bz)) => {
                    if !bx { free_dofs.push(3 * i);     }
                    if !by { free_dofs.push(3 * i + 1); }
                    if !bz { free_dofs.push(3 * i + 2); }
                }
                None => {
                    free_dofs.push(3 * i);
                    free_dofs.push(3 * i + 1);
                    free_dofs.push(3 * i + 2);
                }
            }
        }

        // Build the reduced system by extracting free DOF rows and columns
        let f_reduced = DVector::from_iterator(
            free_dofs.len(),
            free_dofs.iter().map(|&i| f[i]),
        );
        let k_reduced = DMatrix::from_iterator(
            free_dofs.len(),
            free_dofs.len(),
            free_dofs.iter().flat_map(|&col| {
                free_dofs.iter().map(move |&row| k[(row, col)])
            }),
        );

        BoundaryConditionResult {
            k:         k_reduced,
            f:         f_reduced,
            free_dofs: Some(free_dofs),
            n_dofs:    f.len(),
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// reconstruct et reconstruct_ref doivent produire le meme resultat.
    #[test]
    fn test_reconstruct_ref_matches_reconstruct() {
        let free_dofs = vec![0, 1, 3, 4];
        let n_dofs = 6;
        let u_reduced = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);

        // reconstruct_ref via une instance temporaire
        let bc = BoundaryConditionResult {
            k:         DMatrix::zeros(4, 4),
            f:         DVector::zeros(4),
            free_dofs: Some(free_dofs.clone()),
            n_dofs,
        };
        let u_ref = bc.reconstruct_ref(&u_reduced, n_dofs);

        // reconstruct (consume)
        let bc2 = BoundaryConditionResult {
            k:         DMatrix::zeros(4, 4),
            f:         DVector::zeros(4),
            free_dofs: Some(free_dofs),
            n_dofs,
        };
        let u_consume = bc2.reconstruct(u_reduced);

        assert_eq!(u_ref, u_consume);
    }

    /// Les DDL fixes doivent valoir zero apres reconstruction.
    #[test]
    fn test_fixed_dofs_are_zero_after_reconstruction() {
        let free_dofs = vec![0, 2, 4];
        let n_dofs = 6;
        let u_reduced = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let bc = BoundaryConditionResult {
            k:         DMatrix::zeros(3, 3),
            f:         DVector::zeros(3),
            free_dofs: Some(free_dofs),
            n_dofs,
        };
        let u_full = bc.reconstruct_ref(&u_reduced, n_dofs);
        // DOFs 1, 3, 5 doivent etre a zero
        assert_eq!(u_full[1], 0.0);
        assert_eq!(u_full[3], 0.0);
        assert_eq!(u_full[5], 0.0);
        // DOFs libres doivent avoir les bonnes valeurs
        assert_eq!(u_full[0], 1.0);
        assert_eq!(u_full[2], 2.0);
        assert_eq!(u_full[4], 3.0);
    }

    /// La methode penalty ne doit pas changer la taille du systeme.
    #[test]
    fn test_penalty_preserves_system_size() {
        let n = 12;
        let k = DMatrix::zeros(n, n);
        let constraint = Constraint {
            list: vec![FixedNode::all(0)],
        };
        let loads = vec![Load {
            list:  vec![3],
            force: Vector3::new(1.0, 0.0, 0.0),
        }];
        let bc = BoundaryConditions::new(constraint, loads, Box::new(PenaltyMethod));
        let result = bc.apply(&k, n / 3);
        assert_eq!(result.k.nrows(), n);
        assert_eq!(result.f.len(), n);
        assert!(result.free_dofs.is_none());
    }

    /// La methode elimination doit reduire la taille du systeme.
    #[test]
    fn test_elimination_reduces_system_size() {
        let n = 12;
        let k = DMatrix::zeros(n, n);
        let constraint = Constraint {
            list: vec![FixedNode::all(0)], // bloque 3 DDL
        };
        let loads = vec![Load {
            list:  vec![3],
            force: Vector3::new(1.0, 0.0, 0.0),
        }];
        let bc = BoundaryConditions::new(constraint, loads, Box::new(EliminationMethod));
        let result = bc.apply(&k, n / 3);
        assert_eq!(result.k.nrows(), n - 3);
        assert_eq!(result.f.len(), n - 3);
        assert!(result.free_dofs.is_some());
    }
}