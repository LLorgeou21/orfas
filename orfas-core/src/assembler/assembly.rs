// UTF-8
// assembler/assembly.rs — tangent and internal forces assembly.

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::CsrMatrix;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use crate::element::{FiniteElement, FemMesh};
use crate::element::traits::DofType;
use crate::material::{MaterialLaw, SimulationContext};
use super::Assembler;
use nalgebra::Vector3;

impl<E: FiniteElement> Assembler<E> {

    /// Assembles the dense tangent stiffness matrix.
    /// K = sum_e sum_g B(g)^T * C(F(g)) * B(g) * det_J(g) * w(g)
    /// Loops over elements then over Gauss points within each element.
    /// Matrix size is (N_DOF * n_nodes) x (N_DOF * n_nodes) where
    /// N_DOF = E::Dof::N_DOF (3 for volumetric, 6 for structural elements).
    pub fn assemble_tangent(
        &self,
        mesh:     &impl FemMesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> DMatrix<f64> {
        let ndof = <E::Dof as DofType>::N_DOF;
        let n            = ndof * mesh.n_nodes();
        let mut k        = DMatrix::zeros(n, n);
        let gauss_points = E::integration_points();

        for (elem_idx, nodes) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];

            // Analytical stiffness path (Beam2, Shell): bypass Gauss loop.
            let node_positions: Vec<Vector3<f64>> = nodes.iter()
                .map(|&i| mesh.nodes()[i].position)
                .collect();
            if let Some(ke) = E::element_stiffness(&node_positions, material) {
                for r in 0..E::N_NODES {
                    for s in 0..E::N_NODES {
                        let gr = nodes[r];
                        let gs = nodes[s];
                        for dr in 0..ndof {
                            for dc in 0..ndof {
                                k[(ndof * gr + dr, ndof * gs + dc)]
                                    += ke[(ndof * r + dr, ndof * s + dc)];
                            }
                        }
                    }
                }
                continue;
            }

            let displacements: Vec<_> = nodes.iter()
                .map(|&i| Self::node_displacement(u, i))
                .collect();

            for (g_idx, (_xi, weight)) in gauss_points.iter().enumerate() {
                let grad_n = E::shape_gradients(geo, g_idx);
                let f_grad = E::deformation_gradient(&grad_n, &displacements);
                let b_mat  = E::b_matrix(&grad_n, &f_grad);
                let det_j  = E::gauss_det_j(geo, g_idx);
                let ctx    = sim_ctx.material_context_for(elem_idx);
                let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
                let ke     = b_mat.transpose() * c_mat * &b_mat * det_j * *weight;

                for r in 0..E::N_NODES {
                    for s in 0..E::N_NODES {
                        let gr = nodes[r];
                        let gs = nodes[s];
                        for dr in 0..ndof {
                            for dc in 0..ndof {
                                k[(ndof * gr + dr, ndof * gs + dc)]
                                    += ke[(ndof * r + dr, ndof * s + dc)];
                            }
                        }
                    }
                }
            }
        }
        k
    }

    /// Assembles the sparse tangent stiffness matrix (sequential).
    /// Scatters element contributions into the pre-built CSR pattern
    /// via entry_map for O(1) index lookup per entry.
    pub fn assemble_tangent_sparse(
        &self,
        mesh:     &impl FemMesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> CsrMatrix<f64> {
        let ndof = <E::Dof as DofType>::N_DOF;
        let mut k      = self.csr_pattern.clone();
        let values     = k.values_mut();
        for v in values.iter_mut() { *v = 0.0; }

        let gauss_points = E::integration_points();

        for (elem_idx, nodes) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];

            // Analytical stiffness path (Beam2, Shell): bypass Gauss loop.
            let node_positions: Vec<Vector3<f64>> = nodes.iter()
                .map(|&i| mesh.nodes()[i].position)
                .collect();
            if let Some(ke) = E::element_stiffness(&node_positions, material) {
                for r in 0..E::N_NODES {
                    for s in 0..E::N_NODES {
                        let gr = nodes[r];
                        let gs = nodes[s];
                        for dr in 0..ndof {
                            for dc in 0..ndof {
                                let i   = ndof * gr + dr;
                                let j   = ndof * gs + dc;
                                let idx = self.entry_map[&(i, j)];
                                values[idx] += ke[(ndof * r + dr, ndof * s + dc)];
                            }
                        }
                    }
                }
                continue;
            }

            let displacements: Vec<_> = nodes.iter()
                .map(|&i| Self::node_displacement(u, i))
                .collect();

            for (g_idx, (_xi, weight)) in gauss_points.iter().enumerate() {
                let grad_n = E::shape_gradients(geo, g_idx);
                let f_grad = E::deformation_gradient(&grad_n, &displacements);
                let b_mat  = E::b_matrix(&grad_n, &f_grad);
                let det_j  = E::gauss_det_j(geo, g_idx);
                let ctx    = sim_ctx.material_context_for(elem_idx);
                let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
                let ke     = b_mat.transpose() * c_mat * &b_mat * det_j * *weight;

                for r in 0..E::N_NODES {
                    for s in 0..E::N_NODES {
                        let gr = nodes[r];
                        let gs = nodes[s];
                        for dr in 0..ndof {
                            for dc in 0..ndof {
                                let i   = ndof * gr + dr;
                                let j   = ndof * gs + dc;
                                let idx = self.entry_map[&(i, j)];
                                values[idx] += ke[(ndof * r + dr, ndof * s + dc)];
                            }
                        }
                    }
                }
            }
        }
        k
    }

    /// Assembles the sparse tangent stiffness matrix (parallel, lock-free atomic).
    /// Each thread accumulates into AtomicU64 via bit-cast f64 addition.
    /// Requires E::Geometry: Sync (guaranteed by ElementGeometry: Sync).
    pub fn assemble_tangent_sparse_parallel(
        &self,
        mesh:     &(impl FemMesh + Sync),
        material: &(dyn MaterialLaw + Sync),
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> CsrMatrix<f64> {
        let ndof = <E::Dof as DofType>::N_DOF;
        let mut k      = self.csr_pattern.clone();
        let n_vals     = k.values().len();

        let atomic_values: Vec<AtomicU64> = (0..n_vals)
            .map(|_| AtomicU64::new(0u64))
            .collect();

        let gauss_points = E::integration_points();

        self.connectivity
            .par_iter()
            .enumerate()
            .for_each(|(elem_idx, nodes)| {
                let geo = &self.geometry[elem_idx];

                // Analytical stiffness path (Beam2, Shell): bypass Gauss loop.
                let node_positions: Vec<Vector3<f64>> = nodes.iter()
                    .map(|&i| mesh.nodes()[i].position)
                    .collect();
                if let Some(ke) = E::element_stiffness(&node_positions, material) {
                    for r in 0..E::N_NODES {
                        for s in 0..E::N_NODES {
                            let gr = nodes[r];
                            let gs = nodes[s];
                            for dr in 0..ndof {
                                for dc in 0..ndof {
                                    let i   = ndof * gr + dr;
                                    let j   = ndof * gs + dc;
                                    let idx = self.entry_map[&(i, j)];
                                    let v   = ke[(ndof * r + dr, ndof * s + dc)];
                                    atomic_values[idx].fetch_update(
                                        Ordering::Relaxed, Ordering::Relaxed,
                                        |old| Some((f64::from_bits(old) + v).to_bits()),
                                    ).unwrap();
                                }
                            }
                        }
                    }
                    return;
                }

                let displacements: Vec<_> = nodes.iter()
                    .map(|&i| Self::node_displacement(u, i))
                    .collect();

                for (g_idx, (_xi, weight)) in gauss_points.iter().enumerate() {
                    let grad_n = E::shape_gradients(geo, g_idx);
                    let f_grad = E::deformation_gradient(&grad_n, &displacements);
                    let b_mat  = E::b_matrix(&grad_n, &f_grad);
                    let det_j  = E::gauss_det_j(geo, g_idx);
                    let ctx    = sim_ctx.material_context_for(elem_idx);
                    let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
                    let ke     = b_mat.transpose() * c_mat * &b_mat * det_j * *weight;

                    for r in 0..E::N_NODES {
                        for s in 0..E::N_NODES {
                            let gr = nodes[r];
                            let gs = nodes[s];
                            for dr in 0..ndof {
                                for dc in 0..ndof {
                                    let i   = ndof * gr + dr;
                                    let j   = ndof * gs + dc;
                                    let idx = self.entry_map[&(i, j)];
                                    let v   = ke[(ndof * r + dr, ndof * s + dc)];
                                    atomic_values[idx].fetch_update(
                                        Ordering::Relaxed, Ordering::Relaxed,
                                        |old| Some((f64::from_bits(old) + v).to_bits()),
                                    ).unwrap();
                                }
                            }
                        }
                    }
                }
            });

        let values = k.values_mut();
        for (i, av) in atomic_values.iter().enumerate() {
            values[i] = f64::from_bits(av.load(Ordering::Relaxed));
        }
        k
    }

    /// Assembles the internal forces vector f_int.
    /// f_int_i = sum_e sum_g P^T * grad(N_i) * det_J(g) * w(g)
    /// where P = F * S is the first Piola-Kirchhoff stress.
    /// Vector size is N_DOF * n_nodes.
    pub fn assemble_internal_forces(
        &self,
        mesh:     &impl FemMesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> DVector<f64> {
        let ndof = <E::Dof as DofType>::N_DOF;
        let mut f_int    = DVector::zeros(ndof * mesh.n_nodes());
        let gauss_points = E::integration_points();

        for (elem_idx, nodes) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];

            let displacements: Vec<_> = nodes.iter()
                .map(|&i| Self::node_displacement(u, i))
                .collect();

            for (g_idx, (_xi, weight)) in gauss_points.iter().enumerate() {
                let grad_n  = E::shape_gradients(geo, g_idx);
                let f_grad  = E::deformation_gradient(&grad_n, &displacements);
                let det_j   = E::gauss_det_j(geo, g_idx);
                let mut ctx = sim_ctx.material_context_for(elem_idx);
                let s       = material.pk2_stress(&f_grad, &mut ctx);
                let p       = f_grad * s;

                for i in 0..E::N_NODES {
                    let grad_ni  = grad_n.row(i).transpose();
                    let f_node   = det_j * weight * p.transpose() * grad_ni;
                    let global_i = nodes[i];
                    // Scatter the 3 translational DOFs — for Vec6Dof elements,
                    // rotational contributions are handled in the B-matrix and
                    // do not appear in f_int via the PK1 path.
                    f_int[ndof * global_i]     += f_node[0];
                    f_int[ndof * global_i + 1] += f_node[1];
                    f_int[ndof * global_i + 2] += f_node[2];
                }
            }
        }
        f_int
    }

    /// Updates internal variables for all elements after Newton convergence.
    /// No-op for elastic materials (sim_ctx.iv is None).
    /// Must be called once per time step after the Newton loop converges.
    /// Uses the first Gauss point only — consistent with single-point Tet4 convention.
    /// Will be extended to per-Gauss-point update when InternalVariables supports it.
    pub fn update_internal_variables(
        &self,
        mesh:     &impl FemMesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &mut SimulationContext,
    ) {
        if sim_ctx.iv.is_none() { return; }

        for (elem_idx, nodes) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];

            let displacements: Vec<_> = nodes.iter()
                .map(|&i| Self::node_displacement(u, i))
                .collect();

            // Use first Gauss point for state update — single-point convention.
            let grad_n  = E::shape_gradients(geo, 0);
            let f_grad  = E::deformation_gradient(&grad_n, &displacements);
            let mut ctx = sim_ctx.material_context_for_mut(elem_idx);
            material.update_state(&f_grad, &mut ctx);
        }

        // Suppress unused variable warning — mesh used for node count in other methods.
        let _ = mesh.n_nodes();
    }
}