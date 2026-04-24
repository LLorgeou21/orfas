// UTF-8
// assembler/assembly.rs — tangent and internal forces assembly.

use nalgebra::{DMatrix, DVector};
use nalgebra_sparse::CsrMatrix;
use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, Ordering};
use crate::material::{MaterialLaw, SimulationContext};
use crate::mesh::Mesh;
use super::Assembler;
use super::geometry::{BMatrix, compute_deformation_gradient};

impl Assembler {
    /// Assembles the dense tangent stiffness matrix K = sum_e Bt * C * B * V.
    pub fn assemble_tangent<B: BMatrix>(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> DMatrix<f64> {
        let n = 3 * mesh.nodes.len();
        let mut k = DMatrix::zeros(n, n);

        for (elem_idx, &[a, b, c, d]) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];
            if geo.volume < 1e-10 { continue; }

            let u0 = Self::node_displacement(u, a);
            let u1 = Self::node_displacement(u, b);
            let u2 = Self::node_displacement(u, c);
            let u3 = Self::node_displacement(u, d);

            let f_grad = compute_deformation_gradient(&u0, &u1, &u2, &u3, &geo.b, &geo.c, &geo.d);
            let b_mat  = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
            let ctx    = sim_ctx.material_context_for(elem_idx);
            let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
            let ke     = b_mat.transpose() * c_mat * b_mat * geo.volume;

            let global_indices = [a, b, c, d];
            for r in 0..4 {
                for s in 0..4 {
                    let mut block = k.fixed_view_mut::<3, 3>(3 * global_indices[r], 3 * global_indices[s]);
                    block += ke.fixed_view::<3, 3>(3 * r, 3 * s);
                }
            }
        }
        k
    }

    /// Assembles the sparse tangent stiffness matrix (sequential).
    pub fn assemble_tangent_sparse<B: BMatrix>(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> CsrMatrix<f64> {
        let mut k = self.csr_pattern.clone();
        let values = k.values_mut();
        for v in values.iter_mut() { *v = 0.0; }

        for (elem_idx, &[a, b, c, d]) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];
            if geo.volume < 1e-10 { continue; }

            let u0 = Self::node_displacement(u, a);
            let u1 = Self::node_displacement(u, b);
            let u2 = Self::node_displacement(u, c);
            let u3 = Self::node_displacement(u, d);

            let f_grad = compute_deformation_gradient(&u0, &u1, &u2, &u3, &geo.b, &geo.c, &geo.d);
            let b_mat  = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
            let ctx    = sim_ctx.material_context_for(elem_idx);
            let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
            let ke     = b_mat.transpose() * c_mat * b_mat * geo.volume;

            let global_indices = [a, b, c, d];
            for r in 0..4 {
                for s in 0..4 {
                    let block = ke.fixed_view::<3, 3>(3 * r, 3 * s);
                    for dr in 0..3 {
                        for dc in 0..3 {
                            let i   = 3 * global_indices[r] + dr;
                            let j   = 3 * global_indices[s] + dc;
                            let idx = self.entry_map[&(i, j)];
                            values[idx] += block[(dr, dc)];
                        }
                    }
                }
            }
        }
        k
    }

    /// Assembles the sparse tangent stiffness matrix (parallel, atomic).
    pub fn assemble_tangent_sparse_parallel<B: BMatrix>(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> CsrMatrix<f64> {
        let mut k = self.csr_pattern.clone();
        let n_vals = k.values().len();

        let atomic_values: Vec<AtomicU64> = (0..n_vals)
            .map(|_| AtomicU64::new(0u64))
            .collect();

        self.connectivity
            .par_iter()
            .enumerate()
            .for_each(|(elem_idx, &[a, b, c, d])| {
                let geo = &self.geometry[elem_idx];
                if geo.volume < 1e-10 { return; }

                let u0 = Self::node_displacement(u, a);
                let u1 = Self::node_displacement(u, b);
                let u2 = Self::node_displacement(u, c);
                let u3 = Self::node_displacement(u, d);

                let f_grad = compute_deformation_gradient(&u0, &u1, &u2, &u3, &geo.b, &geo.c, &geo.d);
                let b_mat  = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
                let ctx    = sim_ctx.material_context_for(elem_idx);
                let c_mat  = material.tangent_stiffness(&f_grad, &ctx);
                let ke     = b_mat.transpose() * c_mat * b_mat * geo.volume;

                let global_indices = [a, b, c, d];
                for r in 0..4 {
                    for s in 0..4 {
                        let block = ke.fixed_view::<3, 3>(3 * r, 3 * s);
                        for dr in 0..3 {
                            for dc in 0..3 {
                                let i   = 3 * global_indices[r] + dr;
                                let j   = 3 * global_indices[s] + dc;
                                let idx = self.entry_map[&(i, j)];
                                let v   = block[(dr, dc)];
                                atomic_values[idx].fetch_update(
                                    Ordering::Relaxed, Ordering::Relaxed,
                                    |old| Some((f64::from_bits(old) + v).to_bits()),
                                ).unwrap();
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
    pub fn assemble_internal_forces(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &SimulationContext,
    ) -> DVector<f64> {
        let n = 3 * mesh.nodes.len();
        let mut f_int = DVector::zeros(n);

        for (elem_idx, &[a, b, c, d]) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];
            if geo.volume < 1e-10 { continue; }

            let u0 = Self::node_displacement(u, a);
            let u1 = Self::node_displacement(u, b);
            let u2 = Self::node_displacement(u, c);
            let u3 = Self::node_displacement(u, d);

            let f_grad = compute_deformation_gradient(&u0, &u1, &u2, &u3, &geo.b, &geo.c, &geo.d);
            let mut ctx = sim_ctx.material_context_for(elem_idx);
            let s = material.pk2_stress(&f_grad, &mut ctx);
            let p = f_grad * s;

            let global_indices = [a, b, c, d];
            for i in 0..4 {
                let grad_n   = nalgebra::Vector3::new(geo.b[i], geo.c[i], geo.d[i]);
                let f_node   = geo.volume * p.transpose() * grad_n;
                let global_i = global_indices[i];
                f_int[3 * global_i]     += f_node[0];
                f_int[3 * global_i + 1] += f_node[1];
                f_int[3 * global_i + 2] += f_node[2];
            }
        }
        f_int
    }


    /// Updates internal variables for all elements after Newton convergence.
    /// Calls `material.update_state` once per element — no-op for elastic materials.
    /// Must be called once per time step after the Newton loop converges.
    pub fn update_internal_variables(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
        sim_ctx:  &mut SimulationContext,
    ) {
        if sim_ctx.iv.is_none() { return; }

        for (elem_idx, &[a, b, c, d]) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];
            if geo.volume < 1e-10 { continue; }

            let u0 = Self::node_displacement(u, a);
            let u1 = Self::node_displacement(u, b);
            let u2 = Self::node_displacement(u, c);
            let u3 = Self::node_displacement(u, d);

            let f_grad = compute_deformation_gradient(
                &u0, &u1, &u2, &u3,
                &geo.b, &geo.c, &geo.d,
            );

            let mut ctx = sim_ctx.material_context_for_mut(elem_idx);
            material.update_state(&f_grad, &mut ctx);
        }
    }
}