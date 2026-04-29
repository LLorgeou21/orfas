// UTF-8
// orfas-core/src/boundary.rs

use nalgebra::{DMatrix, DVector, Vector3};
use std::collections::HashMap;

// ─── Data structures ──────────────────────────────────────────────────────────

/// Stores a fixed node index and which displacement directions are constrained.
/// Use the constructors (all, only_x, only_y, only_z) for convenience.
/// `value` holds the prescribed displacement — zero for standard fixed BCs.
#[derive(Clone)]
pub struct FixedNode {
    pub indice: usize,
    pub x:      bool,
    pub y:      bool,
    pub z:      bool,
    pub rx:     bool,  // rotation about x (theta_x)
    pub ry:     bool,  // rotation about y (theta_y)
    pub rz:     bool,  // rotation about z (theta_z)
    pub value:  Vector3<f64>,
}

impl FixedNode {
    /// Fix all three translational DOFs to zero.
    pub fn all(indice: usize) -> FixedNode {
        FixedNode { indice, x: true, y: true, z: true,
                    rx: false, ry: false, rz: false, value: Vector3::zeros() }
    }

    /// Fix only the x translational DOF to zero.
    pub fn only_x(indice: usize) -> FixedNode {
        FixedNode { indice, x: true, y: false, z: false,
            rx: false, ry: false, rz: false, value: Vector3::zeros() }
    }

    /// Fix only the y translational DOF to zero.
    pub fn only_y(indice: usize) -> FixedNode {
        FixedNode { indice, x: false, y: true, z: false,
            rx: false, ry: false, rz: false, value: Vector3::zeros() }
    }

    /// Fix only the z translational DOF to zero.
    pub fn only_z(indice: usize) -> FixedNode {
        FixedNode { indice, x: false, y: false, z: true,
            rx: false, ry: false, rz: false, value: Vector3::zeros() }
    }

    /// Fix all three translational DOFs with prescribed non-zero displacement.
    pub fn all_with_value(indice: usize, value: Vector3<f64>) -> FixedNode {
        FixedNode { indice, x: true, y: true, z: true,
            rx: false, ry: false, rz: false, value }
    }

    /// Fix all 6 DOFs (3 translations + 3 rotations) to zero.
    /// Use for clamped boundary conditions on Beam2 and Shell4 elements.
    pub fn clamped(indice: usize) -> FixedNode {
        FixedNode { indice, x: true, y: true, z: true,
                    rx: true, ry: true, rz: true, value: Vector3::zeros() }
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
    /// Prescribed DOF indices and their values for non-zero Dirichlet BCs.
    pub prescribed_dofs: Vec<usize>,
    pub prescribed_vals: Vec<f64>,
}

impl BoundaryConditionResult {

    /// Reconstructs the full displacement vector, consuming self.
    /// For the penalty method (free_dofs = None), returns u_reduced unchanged.
    /// For the elimination method, inserts zeros at fixed DOF positions
    /// and prescribed values at constrained DOF positions.
    pub fn reconstruct(self, u_reduced: DVector<f64>) -> DVector<f64> {
        match self.free_dofs {
            Some(free_dofs) => {
                let mut u_full = DVector::zeros(self.n_dofs);
                for (i, &dof) in free_dofs.iter().enumerate() {
                    u_full[dof] = u_reduced[i];
                }
                for (i, &dof) in self.prescribed_dofs.iter().enumerate() {
                    u_full[dof] = self.prescribed_vals[i];
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
                for (i, &dof) in self.prescribed_dofs.iter().enumerate() {
                    u_full[dof] = self.prescribed_vals[i];
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

    /// Creates a new BoundaryConditions instance.
    pub fn new(
        constraint: Constraint,
        loads:      Vec<Load>,
        method:     Box<dyn BoundaryConditionMethod>,
    ) -> Self {
        BoundaryConditions { constraint, loads, method }
    }

    /// Builds f from loads and applies the boundary condition method.
    /// k is only used by the penalty method — elimination ignores it.
    ///
    /// `n_dof_per_node` must match the element DofType::N_DOF used in the simulation:
    /// 3 for volumetric elements (Tet4, Tet10, Hex8), 6 for structural elements
    /// (Beam2, Shell). Pass `E::Dof::N_DOF` at call sites.
    ///
    /// Note: currently only translational DOFs (x, y, z) are constrained via
    /// FixedNode. For Vec6Dof elements the rotational DOFs are left free unless
    /// explicitly added to the constraint list in a future extension.
    pub fn apply(
        &self,
        k:              &DMatrix<f64>,
        n_nodes:        usize,
        n_dof_per_node: usize,
    ) -> BoundaryConditionResult {
        let mut f = DVector::zeros(n_dof_per_node * n_nodes);
        for load in &self.loads {
            for &idx in &load.list {
                f[n_dof_per_node * idx]     += load.force.x;
                f[n_dof_per_node * idx + 1] += load.force.y;
                f[n_dof_per_node * idx + 2] += load.force.z;
            }
        }
        self.method.apply(k, &f, &self.constraint, n_dof_per_node)
    }
}

/// A set of nodal displacement constraints.
#[derive(Clone)]
pub struct Constraint {
    pub list: Vec<FixedNode>,
}

/// A set of nodal forces applied to a list of nodes.
#[derive(Clone)]
pub struct Load {
    pub list:  Vec<usize>,
    pub force: Vector3<f64>,
}

// ─── Trait ────────────────────────────────────────────────────────────────────

/// Defines how boundary conditions are applied to the system K*u = f.
/// Returns a BoundaryConditionResult containing the modified system
/// and the information needed to reconstruct the full displacement vector.
///
/// `n_dof_per_node` is forwarded from `BoundaryConditions::apply` and must
/// match the element DofType — 3 for volumetric, 6 for structural.
pub trait BoundaryConditionMethod {
    fn apply(
        &self,
        k:              &DMatrix<f64>,
        f:              &DVector<f64>,
        constraint:     &Constraint,
        n_dof_per_node: usize,
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
        k:              &DMatrix<f64>,
        f:              &DVector<f64>,
        constraint:     &Constraint,
        n_dof_per_node: usize,
    ) -> BoundaryConditionResult {
        let mut new_k = k.clone();
        let mut new_f = f.clone();
        for node in &constraint.list {
            let i = n_dof_per_node * node.indice;
            if node.x { new_k[(i,   i  )] = 1e30; new_f[i  ] = 0.0; }
            if node.y { new_k[(i+1, i+1)] = 1e30; new_f[i+1] = 0.0; }
            if node.z { new_k[(i+2, i+2)] = 1e30; new_f[i+2] = 0.0; }
            if node.rx { new_k[(i+3, i+3)] = 1e30; new_f[i+3] = 0.0; }
            if node.ry { new_k[(i+4, i+4)] = 1e30; new_f[i+4] = 0.0; }
            if node.rz { new_k[(i+5, i+5)] = 1e30; new_f[i+5] = 0.0; }
        }
        BoundaryConditionResult {
            k:               new_k,
            f:               new_f,
            free_dofs:       None,
            n_dofs:          f.len(),
            prescribed_dofs: vec![],
            prescribed_vals: vec![],
        }
    }
}

/// Elimination method: removes fixed DOF rows and columns from K and f,
/// solving a smaller system. More numerically stable than the penalty method.
/// The full displacement vector is reconstructed after solving via
/// reconstruct / reconstruct_ref.
pub struct EliminationMethod;

impl BoundaryConditionMethod for EliminationMethod {
    fn apply(
        &self,
        k:              &DMatrix<f64>,
        f:              &DVector<f64>,
        constraint:     &Constraint,
        n_dof_per_node: usize,
    ) -> BoundaryConditionResult {
        // Build map: node index -> (bx, by, bz, brx, bry, brz).
        let mut blocked: HashMap<usize, (bool, bool, bool, bool, bool, bool)> = HashMap::new();
        for node in &constraint.list {
            blocked.insert(node.indice, (node.x, node.y, node.z, node.rx, node.ry, node.rz));
        }
 
        let n_nodes             = f.len() / n_dof_per_node;
        let mut free_dofs       = Vec::new();
        let mut prescribed_dofs = Vec::new();
        let mut prescribed_vals = Vec::new();
 
        for i in 0..n_nodes {
            match blocked.get(&i) {
                Some((bx, by, bz, brx, bry, brz)) => {
                    let node = constraint.list.iter().find(|n| n.indice == i).unwrap();
 
                    // Translational DOFs (always present)
                    if *bx {
                        prescribed_dofs.push(n_dof_per_node * i);
                        prescribed_vals.push(node.value.x);
                    } else {
                        free_dofs.push(n_dof_per_node * i);
                    }
                    if *by {
                        prescribed_dofs.push(n_dof_per_node * i + 1);
                        prescribed_vals.push(node.value.y);
                    } else {
                        free_dofs.push(n_dof_per_node * i + 1);
                    }
                    if *bz {
                        prescribed_dofs.push(n_dof_per_node * i + 2);
                        prescribed_vals.push(node.value.z);
                    } else {
                        free_dofs.push(n_dof_per_node * i + 2);
                    }
 
                    // Rotational DOFs — only for Vec6Dof elements (n_dof_per_node > 3)
                    if n_dof_per_node > 3 {
                        if *brx {
                            prescribed_dofs.push(n_dof_per_node * i + 3);
                            prescribed_vals.push(0.0);
                        } else {
                            free_dofs.push(n_dof_per_node * i + 3);
                        }
                        if *bry {
                            prescribed_dofs.push(n_dof_per_node * i + 4);
                            prescribed_vals.push(0.0);
                        } else {
                            free_dofs.push(n_dof_per_node * i + 4);
                        }
                        if *brz {
                            prescribed_dofs.push(n_dof_per_node * i + 5);
                            prescribed_vals.push(0.0);
                        } else {
                            free_dofs.push(n_dof_per_node * i + 5);
                        }
                    }
                }
                None => {
                    for dof_offset in 0..n_dof_per_node {
                        free_dofs.push(n_dof_per_node * i + dof_offset);
                    }
                }
            }
        }
 
        // Build reduced f: subtract coupling terms from prescribed DOFs.
        let mut f_reduced = DVector::from_iterator(
            free_dofs.len(),
            free_dofs.iter().map(|&i| f[i]),
        );
        for (j, &pdof) in prescribed_dofs.iter().enumerate() {
            let u_val = prescribed_vals[j];
            for (i, &rdof) in free_dofs.iter().enumerate() {
                f_reduced[i] -= k[(rdof, pdof)] * u_val;
            }
        }
 
        // Build reduced K from free DOF submatrix.
        let k_reduced = DMatrix::from_iterator(
            free_dofs.len(),
            free_dofs.len(),
            free_dofs.iter().flat_map(|&col| {
                free_dofs.iter().map(move |&row| k[(row, col)])
            }),
        );
 
        BoundaryConditionResult {
            k:               k_reduced,
            f:               f_reduced,
            free_dofs:       Some(free_dofs),
            n_dofs:          f.len(),
            prescribed_dofs,
            prescribed_vals,
        }
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    /// reconstruct and reconstruct_ref must produce identical results.
    #[test]
    fn test_reconstruct_ref_matches_reconstruct() {
        let free_dofs = vec![0, 1, 3, 4];
        let n_dofs    = 6;
        let u_reduced = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);

        let bc = BoundaryConditionResult {
            k:               DMatrix::zeros(4, 4),
            f:               DVector::zeros(4),
            free_dofs:       Some(free_dofs.clone()),
            n_dofs,
            prescribed_dofs: vec![],
            prescribed_vals: vec![],
        };
        let u_ref = bc.reconstruct_ref(&u_reduced, n_dofs);

        let bc2 = BoundaryConditionResult {
            k:               DMatrix::zeros(4, 4),
            f:               DVector::zeros(4),
            free_dofs:       Some(free_dofs),
            n_dofs,
            prescribed_dofs: vec![],
            prescribed_vals: vec![],
        };
        let u_consume = bc2.reconstruct(u_reduced);

        assert_eq!(u_ref, u_consume);
    }

    /// Fixed DOFs must be zero after reconstruction.
    #[test]
    fn test_fixed_dofs_are_zero_after_reconstruction() {
        let free_dofs = vec![0, 2, 4];
        let n_dofs    = 6;
        let u_reduced = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let bc = BoundaryConditionResult {
            k:               DMatrix::zeros(3, 3),
            f:               DVector::zeros(3),
            free_dofs:       Some(free_dofs),
            n_dofs,
            prescribed_dofs: vec![],
            prescribed_vals: vec![],
        };
        let u_full = bc.reconstruct_ref(&u_reduced, n_dofs);
        assert_eq!(u_full[1], 0.0);
        assert_eq!(u_full[3], 0.0);
        assert_eq!(u_full[5], 0.0);
        assert_eq!(u_full[0], 1.0);
        assert_eq!(u_full[2], 2.0);
        assert_eq!(u_full[4], 3.0);
    }

    /// Penalty method must not change the system size.
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
        let result = bc.apply(&k, n / 3, 3);
        assert_eq!(result.k.nrows(), n);
        assert_eq!(result.f.len(), n);
        assert!(result.free_dofs.is_none());
    }

    /// Elimination method must reduce the system size by the number of fixed DOFs.
    #[test]
    fn test_elimination_reduces_system_size() {
        let n = 12;
        let k = DMatrix::zeros(n, n);
        let constraint = Constraint {
            list: vec![FixedNode::all(0)], // blocks 3 DOFs
        };
        let loads = vec![Load {
            list:  vec![3],
            force: Vector3::new(1.0, 0.0, 0.0),
        }];
        let bc = BoundaryConditions::new(constraint, loads, Box::new(EliminationMethod));
        let result = bc.apply(&k, n / 3, 3);
        assert_eq!(result.k.nrows(), n - 3);
        assert_eq!(result.f.len(), n - 3);
        assert!(result.free_dofs.is_some());
    }
}