// UTF-8
// assembler/pattern.rs — CSR sparsity pattern and element coloring.

use nalgebra_sparse::CsrMatrix;
use std::collections::{HashMap, HashSet};

/// Pre-computes the CSR sparsity pattern from mesh connectivity.
pub fn build_csr_pattern(
    connectivity: &[[usize; 4]],
    n: usize,
) -> (CsrMatrix<f64>, HashMap<(usize, usize), usize>) {
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
/// Elements of the same color share no nodes.
pub fn build_element_colors(connectivity: &[[usize; 4]]) -> Vec<Vec<usize>> {
    let n_elems = connectivity.len();
    let mut node_to_elems: HashMap<usize, Vec<usize>> = HashMap::new();
    for (elem_idx, &[a, b, c, d]) in connectivity.iter().enumerate() {
        for &node in &[a, b, c, d] {
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