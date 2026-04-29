// UTF-8
// assembler/mod.rs — Assembler struct, constructors, mass assembly, re-exports.

mod pattern;
mod assembly;

use nalgebra::{DVector, Vector3};
use nalgebra_sparse::CsrMatrix;
use std::collections::HashMap;
use std::marker::PhantomData;
use crate::element::{FiniteElement, FemMesh};
use crate::element::traits::DofType;
use crate::material::MaterialLaw;
use pattern::build_csr_pattern;

// ---------------------------------------------------------------------------
// Assembler<E>
// ---------------------------------------------------------------------------

/// Generic finite element assembler.
/// Parameterized over E: FiniteElement — works with Tet4, Tet10, Hex8,
/// Beam2, and any future element type without code duplication.
/// Precomputes all geometric data at construction time (geometry cache).
pub struct Assembler<E: FiniteElement> {
    /// Node indices per element, stored as Vec<usize> for runtime flexibility.
    pub(crate) connectivity: Vec<Vec<usize>>,
    /// Precomputed geometric data per element (type depends on E).
    pub(crate) geometry:     Vec<E::Geometry>,
    /// Pre-built CSR sparsity pattern (values reset to 0 before each assembly).
    pub(crate) csr_pattern:  CsrMatrix<f64>,
    /// Maps (row, col) DOF pairs to CSR value array indices for O(1) scatter.
    pub(crate) entry_map:    HashMap<(usize, usize), usize>,
    /// Element coloring for conflict-free parallel assembly (optional).
    pub(crate) colors:       Option<Vec<Vec<usize>>>,
    pub(crate) _phantom:     PhantomData<E>,
}

impl<E: FiniteElement> Assembler<E> {
    /// Build a new assembler from any mesh implementing FemMesh.
    /// Precomputes element geometry and CSR sparsity pattern.
    /// Called once before simulation — O(n_elements) construction cost.
    ///
    /// # Arguments
    /// * `mesh` - any mesh implementing FemMesh (Tet4Mesh, Tet10Mesh, Hex8Mesh, ...)
    pub fn new(mesh: &impl FemMesh) -> Assembler<E> {
        let raw_connectivity = mesh.connectivity();

        let geometry: Vec<E::Geometry> = raw_connectivity.iter()
            .map(|elem_nodes| {
                let positions: Vec<Vector3<f64>> = elem_nodes.iter()
                    .map(|&i| mesh.nodes()[i].position)
                    .collect();
                E::precompute(&positions)
            })
            .collect();

        // CSR pattern size uses N_DOF per node — works for Vec3Dof (3) and Vec6Dof (6).
        let ndof = <E::Dof as DofType>::N_DOF;
        let n    = ndof * mesh.n_nodes();
        let (csr_pattern, entry_map) = build_csr_pattern(&raw_connectivity, n);

        Assembler {
            connectivity: raw_connectivity,
            geometry,
            csr_pattern,
            entry_map,
            colors: None,
            _phantom: PhantomData,
        }
    }

    /// Extract displacement vector for node i from the global DOF vector.
    /// Uses N_DOF from the element DofType — works for both Vec3Dof and Vec6Dof.
    /// Only the first 3 components (translations) are returned; rotational DOFs
    /// are handled separately in structural element formulations.
    pub(crate) fn node_displacement(u: &DVector<f64>, i: usize) -> Vector3<f64> {
        let ndof = <E::Dof as DofType>::N_DOF;
        Vector3::new(u[ndof * i], u[ndof * i + 1], u[ndof * i + 2])
    }

    /// Assemble the lumped mass vector.
    /// Each element contributes density * volume / N_NODES to each of its nodes,
    /// distributed equally (lumped mass matrix diagonal).
    /// Vector size is N_DOF * n_nodes to match the global DOF layout.
    ///
    /// # Arguments
    /// * `mesh`     - mesh providing node count and connectivity
    /// * `material` - material providing density
    pub fn assemble_mass(&self, mesh: &impl FemMesh, material: &dyn MaterialLaw) -> DVector<f64> {
        let ndof     = <E::Dof as DofType>::N_DOF;
        let mut mass = DVector::zeros(ndof * mesh.n_nodes());

        for (elem_idx, nodes) in self.connectivity.iter().enumerate() {
            let geo       = &self.geometry[elem_idx];
            let vol       = E::element_volume(geo);
            let elem_mass = material.density() * vol;
            for &node_idx in nodes {
                for j in 0..ndof {
                    mass[ndof * node_idx + j] += elem_mass / (E::N_NODES as f64);
                }
            }
        }
        mass
    }
}

#[cfg(test)]
mod tests;