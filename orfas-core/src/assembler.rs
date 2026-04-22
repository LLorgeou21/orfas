use nalgebra::{DMatrix, DVector, Matrix3, SMatrix, Vector3};
use nalgebra_sparse::{CooMatrix,CsrMatrix};
use crate::material::MaterialLaw;
use crate::mesh::Mesh;
use rayon::prelude::*;
use std::collections::{HashMap,HashSet};

type Matrix6x12 = SMatrix<f64, 6, 12>;

// ─── Helpers geometriques ─────────────────────────────────────────────────────

/// Volume d'un tetraedre : V = |det([p1-p0, p2-p0, p3-p0])| / 6
fn tetra_volume(
    p0: &Vector3<f64>, p1: &Vector3<f64>,
    p2: &Vector3<f64>, p3: &Vector3<f64>,
) -> f64 {
    Matrix3::from_columns(&[p1 - p0, p2 - p0, p3 - p0]).determinant().abs() / 6.0
}

/// Gradients des fonctions de forme (b, c, d) pour un tetraedre lineaire.
/// bi = dNi/dx, ci = dNi/dy, di = dNi/dz — constants par element.
/// Le 4eme coefficient est derive de la partition de l'unite : b4 = -(b1+b2+b3).
fn tetra_bcd(
    p0: &Vector3<f64>, p1: &Vector3<f64>,
    p2: &Vector3<f64>, p3: &Vector3<f64>,
) -> ([f64; 4], [f64; 4], [f64; 4]) {
    let b1 = (p1.y - p3.y) * (p2.z - p3.z) - (p2.y - p3.y) * (p1.z - p3.z);
    let b2 = (p2.y - p3.y) * (p0.z - p3.z) - (p0.y - p3.y) * (p2.z - p3.z);
    let b3 = (p0.y - p3.y) * (p1.z - p3.z) - (p1.y - p3.y) * (p0.z - p3.z);
    let b4 = -(b1 + b2 + b3);

    let c1 = (p2.x - p3.x) * (p1.z - p3.z) - (p1.x - p3.x) * (p2.z - p3.z);
    let c2 = (p0.x - p3.x) * (p2.z - p3.z) - (p2.x - p3.x) * (p0.z - p3.z);
    let c3 = (p1.x - p3.x) * (p0.z - p3.z) - (p0.x - p3.x) * (p1.z - p3.z);
    let c4 = -(c1 + c2 + c3);

    let d1 = (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
    let d2 = (p2.x - p3.x) * (p0.y - p3.y) - (p0.x - p3.x) * (p2.y - p3.y);
    let d3 = (p0.x - p3.x) * (p1.y - p3.y) - (p1.x - p3.x) * (p0.y - p3.y);
    let d4 = -(d1 + d2 + d3);

    ([b1, b2, b3, b4], [c1, c2, c3, c4], [d1, d2, d3, d4])
}

/// Matrice B (6x12) en notation de Voigt depuis les gradients et le volume.
fn tetra_b_matrix(b: &[f64; 4], c: &[f64; 4], d: &[f64; 4], volume: f64) -> Matrix6x12 {
    Matrix6x12::from_row_slice(&[
        b[0], 0.0,  0.0,  b[1], 0.0,  0.0,  b[2], 0.0,  0.0,  b[3], 0.0,  0.0,
        0.0,  c[0], 0.0,  0.0,  c[1], 0.0,  0.0,  c[2], 0.0,  0.0,  c[3], 0.0,
        0.0,  0.0,  d[0], 0.0,  0.0,  d[1], 0.0,  0.0,  d[2], 0.0,  0.0,  d[3],
        c[0], b[0], 0.0,  c[1], b[1], 0.0,  c[2], b[2], 0.0,  c[3], b[3], 0.0,
        0.0,  d[0], c[0], 0.0,  d[1], c[1], 0.0,  d[2], c[2], 0.0,  d[3], c[3],
        d[0], 0.0,  b[0], d[1], 0.0,  b[1], d[2], 0.0,  b[2], d[3], 0.0,  b[3],
    ]).scale(1.0 / (6.0 * volume))
}

/// Gradient de deformation F = I + sum_i ui x gradNi pour un tetraedre lineaire.
///
/// F = dx/dX = d(X + u)/dX = I + du/dX
///
/// ui sont les deplacements nodaux (pas les positions deformees).
/// gradNi = [b[i], c[i], d[i]]^T — gradients des fonctions de forme.
///
/// Propriete cle : avec u = 0, F = I.
/// Garantit pk2_stress(I) = 0 et tangent_stiffness(I) = C_hooke,
/// coherence avec le cas lineaire au repos.
fn compute_deformation_gradient(
    u0: &Vector3<f64>, u1: &Vector3<f64>,
    u2: &Vector3<f64>, u3: &Vector3<f64>,
    b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
) -> Matrix3<f64> {
    let mut f = Matrix3::identity();
    for i in 0..4 {
        let u = [u0, u1, u2, u3][i];
        let grad_n = Vector3::new(b[i], c[i], d[i]);
        f += u * grad_n.transpose();
    }
    f
}

// ─── Trait BMatrix ────────────────────────────────────────────────────────────

/// Trait definissant le calcul de la matrice strain-displacement B.
/// Prend les gradients caches et F — LinearBMatrix ignore F.
/// Une future NonlinearBMatrix (Tet10, formulation mixte) pourrait utiliser F.
pub trait BMatrix {
    fn compute(
        b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
        volume: f64,
        f: &Matrix3<f64>,
    ) -> Matrix6x12;
}

/// Implementation lineaire — B constant, F ignore.
pub struct LinearBMatrix;

impl BMatrix for LinearBMatrix {
    fn compute(
        b: &[f64; 4], c: &[f64; 4], d: &[f64; 4],
        volume: f64,
        _f: &Matrix3<f64>,
    ) -> Matrix6x12 {
        tetra_b_matrix(b, c, d, volume)
    }
}

// ─── Assembler ────────────────────────────────────────────────────────────────

/// Donnees geometriques de reference d'un element, calculees une fois a l'init.
/// Inspire de TetrahedronSetGeometryAlgorithms dans SOFA.
struct ElementGeometry {
    volume: f64,
    b:      [f64; 4],
    c:      [f64; 4],
    d:      [f64; 4],
}

/// Computes a greedy element coloring for parallel assembly.
/// Elements of the same color share no nodes — safe to assemble in parallel
/// without atomic operations or locks.
///
/// Algorithm:
///   1. Build node -> elements map (adjacency inverse)
///   2. For each element, collect colors used by neighboring elements
///   3. Assign the smallest available color
///
/// Returns colors[c] = list of element indices with color c.
/// Number of colors is typically 5-10 for tetrahedral meshes.
fn build_element_colors(connectivity: &[[usize; 4]]) -> Vec<Vec<usize>> {
    let n_elems = connectivity.len();

    // Step 1 — node -> elements map
    let mut node_to_elems: HashMap<usize, Vec<usize>> = HashMap::new();
    for (elem_idx, &[a, b, c, d]) in connectivity.iter().enumerate() {
        for &node in &[a, b, c, d] {
            node_to_elems.entry(node).or_insert_with(Vec::new).push(elem_idx);
        }
    }

    // Step 2 & 3 — greedy coloring
    let mut elem_color = vec![usize::MAX; n_elems];

    for elem_idx in 0..n_elems {
        // Collect colors used by neighboring elements
        let mut used_colors: std::collections::HashSet<usize> = HashSet::new();
        for &node in &connectivity[elem_idx] {
            if let Some(neighbors) = node_to_elems.get(&node) {
                for &neighbor in neighbors {
                    if neighbor != elem_idx && elem_color[neighbor] != usize::MAX {
                        used_colors.insert(elem_color[neighbor]);
                    }
                }
            }
        }

        // Assign smallest available color
        let mut color = 0;
        while used_colors.contains(&color) {
            color += 1;
        }
        elem_color[elem_idx] = color;
    }

    // Group elements by color
    let n_colors = *elem_color.iter().max().unwrap_or(&0) + 1;
    let mut colors: Vec<Vec<usize>> = vec![Vec::new(); n_colors];
    for (elem_idx, &color) in elem_color.iter().enumerate() {
        colors[color].push(elem_idx);
    }

    colors
}

/// Pre-computes the CSR sparsity pattern from the mesh connectivity.
/// For each tetrahedron, the 4 nodes are fully coupled — all 4x4x9 = 144
/// (i,j) DOF pairs are non-zero. Pairs are deduplicated via BTreeSet
/// (sorted row-major) and converted directly to CSR data arrays.
/// Also builds entry_map: (i,j) -> flat index in the CSR values Vec.
/// This is called once at Assembler::new and reused at every assembly.
fn build_csr_pattern(
    connectivity: &[[usize; 4]],
    n: usize,
) -> (CsrMatrix<f64>, HashMap<(usize, usize), usize>) {
    // BTreeSet gives sorted row-major order for free — required by CSR format
    let mut pairs: std::collections::BTreeSet<(usize, usize)> = Default::default();
    for &[a, b, c, d] in connectivity {
        let nodes = [a, b, c, d];
        for &r in &nodes {
            for &s in &nodes {
                for dr in 0..3 {
                    for dc in 0..3 {
                        pairs.insert((3 * r + dr, 3 * s + dc));
                    }
                }
            }
        }
    }

    // Build CSR arrays directly from sorted pairs — no COO intermediate
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::with_capacity(pairs.len());
    let values         = vec![0.0f64; pairs.len()];

    for &(i, j) in &pairs {
        row_offsets[i + 1] += 1;
        col_indices.push(j);
    }
    // Prefix sum to get row offsets
    for i in 0..n {
        row_offsets[i + 1] += row_offsets[i];
    }

    let csr = CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values)
        .expect("Invalid CSR pattern in build_csr_pattern");

    // entry_map: (i,j) -> flat index in values Vec
    // BTreeSet iteration order is sorted — same order as CSR values array
    let entry_map: HashMap<(usize, usize), usize> = pairs
        .into_iter()
        .enumerate()
        .map(|(idx, pair)| (pair, idx))
        .collect();

    (csr, entry_map)
}

pub struct Assembler {
    /// Connectivity table: 4 global node indices per tetrahedron.
    connectivity: Vec<[usize; 4]>,
    /// Reference geometry cached per element.
    geometry:     Vec<ElementGeometry>,
    /// Pre-built sparse pattern — fixed for the lifetime of the assembler.
    /// Reused at every sparse assembly call to avoid repeated COO->CSR conversion.
    csr_pattern:  CsrMatrix<f64>,
    /// Maps (i, j) DOF pair to its flat index in the CSR values array.
    /// Built once at Assembler::new alongside csr_pattern.
    entry_map:    HashMap<(usize, usize), usize>,
    /// Element coloring for parallel assembly.
    /// colors[c] = list of element indices with color c.
    /// Elements of the same color share no nodes — safe to assemble in parallel.
    colors:       Option<Vec<Vec<usize>>>
}

impl Assembler {

    /// Construit l'assembleur et calcule les donnees geometriques de reference.
    /// Volumes et gradients (b,c,d) sont calcules une seule fois ici.
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
                let volume = tetra_volume(p0, p1, p2, p3);
                let (bv, cv, dv) = tetra_bcd(p0, p1, p2, p3);
                ElementGeometry { volume, b: bv, c: cv, d: dv }
            })
            .collect();

        let n = 3 * mesh.nodes.len();
        let (csr_pattern, entry_map) = build_csr_pattern(&connectivity, n);
        //let colors = build_element_colors(&connectivity);
        Assembler { connectivity, geometry, csr_pattern, entry_map, colors : None }
    }

    

    /// Deplacement du noeud i depuis le vecteur u global.
    fn node_displacement(u: &DVector<f64>, i: usize) -> Vector3<f64> {
        Vector3::new(u[3 * i], u[3 * i + 1], u[3 * i + 2])
    }

    /// Assemble la matrice de rigidite tangente K = sum_e Bt * C * B * V.
    /// u : deplacements courants — utilise pour calculer F par element.
    /// Pour u = zeros, F = I et on retombe sur le cas lineaire.
    pub fn assemble_tangent<B: BMatrix>(
        &self,
        mesh: &Mesh,
        material: &dyn MaterialLaw,
        u: &DVector<f64>,
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

            let f_grad = compute_deformation_gradient(
                &u0, &u1, &u2, &u3,
                &geo.b, &geo.c, &geo.d,
            );

            let b_mat = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
            let c_mat = material.tangent_stiffness(&f_grad);
            let ke = b_mat.transpose() * c_mat * b_mat * geo.volume;

            let global_indices = [a, b, c, d];
            for r in 0..4 {
                for s in 0..4 {
                    let mut block = k.fixed_view_mut::<3, 3>(
                        3 * global_indices[r],
                        3 * global_indices[s],
                    );
                    block += ke.fixed_view::<3, 3>(3 * r, 3 * s);
                }
            }
        }
        k
    }


    /// Assembles the sparse tangent stiffness matrix K in parallel using rayon.
    /// Uses atomic f64 additions via fetch_update — no unsafe, no coloring required.
    /// Each entry in the CSR values array is updated atomically, preventing races
    /// on shared nodes between threads.
    ///
    /// Prefer this over assemble_tangent_sparse for large meshes.
    pub fn assemble_tangent_sparse_parallel<B: BMatrix>(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
    ) -> CsrMatrix<f64> {
        use std::sync::atomic::{AtomicU64, Ordering};

        let mut k = self.csr_pattern.clone();
        let n_vals = k.values().len();

        // Atomic values array — one AtomicU64 per CSR entry
        // f64 bits are stored as u64 to enable atomic addition
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

                let f_grad = compute_deformation_gradient(
                    &u0, &u1, &u2, &u3,
                    &geo.b, &geo.c, &geo.d,
                );

                let b_mat = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
                let c_mat = material.tangent_stiffness(&f_grad);
                let ke = b_mat.transpose() * c_mat * b_mat * geo.volume;

                let global_indices = [a, b, c, d];
                for r in 0..4 {
                    for s in 0..4 {
                        let block = ke.fixed_view::<3, 3>(3 * r, 3 * s);
                        for dr in 0..3 {
                            for dc in 0..3 {
                                let i = 3 * global_indices[r] + dr;
                                let j = 3 * global_indices[s] + dc;
                                let idx = self.entry_map[&(i, j)];
                                let v = block[(dr, dc)];
                                // Atomic f64 addition: read-modify-write, no race possible
                                atomic_values[idx].fetch_update(
                                    Ordering::Relaxed,
                                    Ordering::Relaxed,
                                    |old| Some((f64::from_bits(old) + v).to_bits()),
                                ).unwrap();
                            }
                        }
                    }
                }
            });

        // Copy atomic values into CSR matrix
        let values = k.values_mut();
        for (i, av) in atomic_values.iter().enumerate() {
            values[i] = f64::from_bits(av.load(Ordering::Relaxed));
        }
        k
    }

    /// Assembles the sparse tangent stiffness matrix K = sum_e Bt * C * B * V.
    /// Uses the pre-built CSR pattern from Assembler::new — avoids COO->CSR
    /// conversion and sorting overhead. Values are written directly into the
    /// CSR values array via entry_map lookup.
    pub fn assemble_tangent_sparse<B: BMatrix>(
        &self,
        mesh:     &Mesh,
        material: &dyn MaterialLaw,
        u:        &DVector<f64>,
    ) -> CsrMatrix<f64> {
        let mut k = self.csr_pattern.clone();
        let values = k.values_mut();

        // Zero out values from previous assembly
        for v in values.iter_mut() { *v = 0.0; }

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

            let b_mat = B::compute(&geo.b, &geo.c, &geo.d, geo.volume, &f_grad);
            let c_mat = material.tangent_stiffness(&f_grad);
            let ke = b_mat.transpose() * c_mat * b_mat * geo.volume;

            let global_indices = [a, b, c, d];
            for r in 0..4 {
                for s in 0..4 {
                    let block = ke.fixed_view::<3, 3>(3 * r, 3 * s);
                    for dr in 0..3 {
                        for dc in 0..3 {
                            let i = 3 * global_indices[r] + dr;
                            let j = 3 * global_indices[s] + dc;
                            let idx = self.entry_map[&(i, j)];
                            values[idx] += block[(dr, dc)];
                        }
                    }
                }
            }
        }
        k
    }

    /// Assemble le vecteur des forces internes f_int.
    /// Pour chaque element : f_int_i = V * P^T * gradNi
    /// avec P = F * S (PK1) et S = pk2_stress(F) (PK2).
    pub fn assemble_internal_forces(
        &self,
        mesh: &Mesh,
        material: &dyn MaterialLaw,
        u: &DVector<f64>,
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

            let f_grad = compute_deformation_gradient(
                &u0, &u1, &u2, &u3,
                &geo.b, &geo.c, &geo.d,
            );

            // PK2 -> PK1 : P = F * S
            let s = material.pk2_stress(&f_grad);
            let p = f_grad * s;

            // f_int_i = V * P^T * gradNi
            let global_indices = [a, b, c, d];
            for i in 0..4 {
                let grad_n = Vector3::new(geo.b[i], geo.c[i], geo.d[i]);
                let f_node = geo.volume * p.transpose() * grad_n;
                let global_i = global_indices[i];
                f_int[3 * global_i]     += f_node[0];
                f_int[3 * global_i + 1] += f_node[1];
                f_int[3 * global_i + 2] += f_node[2];
            }
        }
        f_int
    }


    
    /// Assemble le vecteur de masse concentree (lumped mass).
    /// Chaque noeud recoit 1/4 de la masse de chaque element connecte.
    /// Stocke comme DVector diagonal (3 DDL par noeud).
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

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::SaintVenantKirchhoff;
    use crate::mesh::Mesh;

    fn unit_mesh() -> Mesh {
        Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0)
    }

    fn svk() -> SaintVenantKirchhoff {
        SaintVenantKirchhoff { youngs_modulus: 1e6, poisson_ratio: 0.3, density: 1000.0 }
    }

    #[test]
    fn test_volume_calculation() {
        let v = tetra_volume(
            &Vector3::new(0.0, 0.0, 0.0),
            &Vector3::new(1.0, 0.0, 0.0),
            &Vector3::new(0.0, 1.0, 0.0),
            &Vector3::new(0.0, 0.0, 1.0),
        );
        assert_eq!(v, 1.0 / 6.0);
    }

    /// K assemblee avec u=0 doit etre symetrique.
    #[test]
    fn test_tangent_symmetric_at_zero() {
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let u = DVector::zeros(3 * mesh.nodes.len());
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u);
        let diff = (&k - k.transpose()).abs().max();
        assert!(diff < 1e-10, "K doit etre symetrique, diff = {:.2e}", diff);
    }

    /// f_int doit etre nul quand u=0 (pas de deformation -> pas de forces internes).
    #[test]
    fn test_internal_forces_zero_at_rest() {
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let u = DVector::zeros(3 * mesh.nodes.len());
        let f_int = assembler.assemble_internal_forces(&mesh, &mat, &u);
        assert!(f_int.norm() < 1e-10, "f_int doit etre nul au repos, norm = {:.2e}", f_int.norm());
    }

    /// Pour de petits deplacements, f_int doit etre proche de K*u (linearisation).
    #[test]
    fn test_internal_forces_consistent_with_tangent() {
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let n = 3 * mesh.nodes.len();
        let u_zero = DVector::zeros(n);
        let k = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);

        let mut u_small = DVector::zeros(n);
        u_small[3] = 1e-6;
        u_small[4] = 1e-6;

        let f_int = assembler.assemble_internal_forces(&mesh, &mat, &u_small);
        let f_lin = &k * &u_small;

        let error = (&f_int - &f_lin).norm() / f_lin.norm();
        assert!(error < 1e-4, "f_int doit lineariser vers K*u pour petits u, erreur = {:.2e}", error);
    }

    /// Masse totale = density * volume.
    #[test]
    fn test_mass_assembly() {
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let mass = assembler.assemble_mass(&mesh, &mat);
        let total = mass.sum() / 3.0;
        let expected = 1000.0 * 1.0;
        assert!((total - expected).abs() / expected < 1e-10);
    }

    /// assemble_tangent_sparse = assemble_tangent
    #[test]
    fn test_assemble_method_comparaison(){
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let n = 3 * mesh.nodes.len();
        let u_zero = DVector::zeros(n);
        
        let k_dense = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let k_sparse = assembler.assemble_tangent_sparse::<LinearBMatrix>(&mesh, &mat, &u_zero);

        let k_dense_from_sparse = DMatrix::from(&k_sparse);
        let diff = (k_dense - k_dense_from_sparse).abs().max();
        assert!(diff < 1e-9, "diff = {:.2e}", diff);
    }

    /// assemble_tangent_sparse = assemble_tangent_sparse_parallel
    #[test]
    fn test_assemble_parrallel_method_comparaison(){
        let mesh = unit_mesh();
        let mat = svk();
        let assembler = Assembler::new(&mesh);
        let n = 3 * mesh.nodes.len();
        let u_zero = DVector::zeros(n);
        
        let k_par = assembler.assemble_tangent_sparse_parallel::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let k_seq = assembler.assemble_tangent_sparse::<LinearBMatrix>(&mesh, &mat, &u_zero);


        let k_seq_dense = DMatrix::from(&k_seq);
        let k_par_dense = DMatrix::from(&k_par);

        let diff = (k_seq_dense - k_par_dense).abs().max();
        assert!(diff < 1e-9, "diff = {:.2e}", diff);
    }

    ///  assemble_tangent_sparse slower than assemble_tangent_sparse_parallel ?
    #[test]
    fn test_assemble_parallel_speedup() {
        let mesh = Mesh::generate(30, 30, 30, 1.0, 1.0, 1.0);
        let mat  = svk();
        let assembler = Assembler::new(&mesh);
        let u_zero = DVector::zeros(3 * mesh.nodes.len());

        let t0 = std::time::Instant::now();
        let _ = assembler.assemble_tangent_sparse::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let t_seq = t0.elapsed();

        let t0 = std::time::Instant::now();
        let _ = assembler.assemble_tangent_sparse_parallel::<LinearBMatrix>(&mesh, &mat, &u_zero);
        let t_par = t0.elapsed();

        println!("sequential: {:.2?}  parallel: {:.2?}  speedup: {:.2}x",
            t_seq, t_par, t_seq.as_secs_f64() / t_par.as_secs_f64());
    }

    

}