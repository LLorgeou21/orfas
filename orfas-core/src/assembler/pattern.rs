// UTF-8
// assembler/pattern.rs — CSR sparsity pattern and element coloring.

use nalgebra_sparse::CsrMatrix;
use std::collections::{HashMap, HashSet};

/// Pre-computes the CSR sparsity pattern from mesh connectivity.
/// Works for any element size — connectivity is a slice of node index lists.
/// Called once in Assembler::new; result stored for the simulation lifetime.
pub fn build_csr_pattern(
    connectivity: &[Vec<usize>],
    n: usize,
) -> (CsrMatrix<f64>, HashMap<(usize, usize), usize>) {
    let mut pairs: std::collections::BTreeSet<(usize, usize)> = Default::default();

    // Derive ndof from the matrix size n and the max node index in connectivity.
    // n = ndof * n_nodes, and node indices are 0..n_nodes-1.
    // We infer ndof by finding max node index + 1 and dividing n by it.
    // This avoids threading ndof through the call site while remaining correct
    // for any DofType (Vec3Dof ndof=3, Vec6Dof ndof=6).
    let max_node = connectivity.iter()
        .flat_map(|nodes| nodes.iter().copied())
        .max()
        .unwrap_or(0);
    let n_nodes = max_node + 1;
    let ndof    = if n_nodes > 0 { n / n_nodes } else { 3 };

    for nodes in connectivity {
        for &r in nodes {
            for &s in nodes {
                for dr in 0..ndof {
                    for dc in 0..ndof {
                        pairs.insert((ndof * r + dr, ndof * s + dc));
                    }
                }
            }
        }
    }

    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::with_capacity(pairs.len());
    let values          = vec![0.0f64; pairs.len()];

    for &(i, j) in &pairs {
        row_offsets[i + 1] += 1;
        col_indices.push(j);
    }
    for i in 0..n {
        row_offsets[i + 1] += row_offsets[i];
    }

    let csr = CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values)
        .expect("Invalid CSR pattern in build_csr_pattern");

    let entry_map: HashMap<(usize, usize), usize> = pairs
        .into_iter()
        .enumerate()
        .map(|(idx, pair)| (pair, idx))
        .collect();

    (csr, entry_map)
}

/// Greedy element coloring for parallel assembly.
/// Elements of the same color share no nodes — safe for concurrent writes.
/// Works for any element size.
pub fn build_element_colors(connectivity: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n_elems = connectivity.len();
    let mut node_to_elems: HashMap<usize, Vec<usize>> = HashMap::new();

    for (elem_idx, nodes) in connectivity.iter().enumerate() {
        for &node in nodes {
            node_to_elems.entry(node).or_insert_with(Vec::new).push(elem_idx);
        }
    }

    let mut elem_color = vec![usize::MAX; n_elems];
    for elem_idx in 0..n_elems {
        let mut used_colors: HashSet<usize> = HashSet::new();
        for &node in &connectivity[elem_idx] {
            if let Some(neighbors) = node_to_elems.get(&node) {
                for &neighbor in neighbors {
                    if neighbor != elem_idx && elem_color[neighbor] != usize::MAX {
                        used_colors.insert(elem_color[neighbor]);
                    }
                }
            }
        }
        let mut color = 0;
        while used_colors.contains(&color) { color += 1; }
        elem_color[elem_idx] = color;
    }

    let n_colors = *elem_color.iter().max().unwrap_or(&0) + 1;
    let mut colors: Vec<Vec<usize>> = vec![Vec::new(); n_colors];
    for (elem_idx, &color) in elem_color.iter().enumerate() {
        colors[color].push(elem_idx);
    }
    colors
}