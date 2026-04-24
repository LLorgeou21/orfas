// UTF-8
// assembler/mod.rs — Assembler struct, constructors, mass assembly, re-exports.

mod geometry;
mod pattern;
mod assembly;

pub use geometry::{BMatrix, LinearBMatrix};

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use nalgebra_sparse::CsrMatrix;
use std::collections::HashMap;
use crate::material::{MaterialLaw, SimulationContext};
use crate::mesh::Mesh;
use geometry::ElementGeometry;
use pattern::{build_csr_pattern, build_element_colors};

pub struct Assembler {
    pub(crate) connectivity: Vec<[usize; 4]>,
    pub(crate) geometry:     Vec<ElementGeometry>,
    pub(crate) csr_pattern:  CsrMatrix<f64>,
    pub(crate) entry_map:    HashMap<(usize, usize), usize>,
    pub(crate) colors:       Option<Vec<Vec<usize>>>,
}

impl Assembler {
    pub fn new(mesh: &Mesh) -> Assembler {
        let connectivity: Vec<[usize; 4]> = mesh.elements.iter()
            .map(|tetra| tetra.indices)
            .collect();

        let geometry: Vec<ElementGeometry> = mesh.elements.iter()
            .map(|tetra| {
                let [a, b, c, d] = tetra.indices;
                let p0 = &mesh.nodes[a].position;
                let p1 = &mesh.nodes[b].position;
                let p2 = &mesh.nodes[c].position;
                let p3 = &mesh.nodes[d].position;
                let volume = geometry::tetra_volume(p0, p1, p2, p3);
                let (bv, cv, dv) = geometry::tetra_bcd(p0, p1, p2, p3);
                ElementGeometry { volume, b: bv, c: cv, d: dv }
            })
            .collect();

        let n = 3 * mesh.nodes.len();
        let (csr_pattern, entry_map) = build_csr_pattern(&connectivity, n);
        Assembler { connectivity, geometry, csr_pattern, entry_map, colors: None }
    }

    pub(crate) fn node_displacement(u: &DVector<f64>, i: usize) -> Vector3<f64> {
        Vector3::new(u[3 * i], u[3 * i + 1], u[3 * i + 2])
    }

    pub fn assemble_mass(&self, mesh: &Mesh, material: &dyn MaterialLaw) -> DVector<f64> {
        let n_nodes = mesh.nodes.len();
        let mut mass = DVector::zeros(3 * n_nodes);
        for (elem_idx, &[a, b, c, d]) in self.connectivity.iter().enumerate() {
            let geo = &self.geometry[elem_idx];
            let elem_mass = material.density() * geo.volume;
            for &node_idx in &[a, b, c, d] {
                for j in 0..3 {
                    mass[3 * node_idx + j] += elem_mass / 4.0;
                }
            }
        }
        mass
    }
}

#[cfg(test)]
mod tests;