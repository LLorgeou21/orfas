// UTF-8
// orfas-core/src/element/subdivision.rs
// Temporary utility to upgrade a Tet4Mesh to a Tet10Mesh by inserting
// midside nodes on each edge. To be replaced by a proper VTK Tet10
// loader in v0.9.0.

use std::collections::HashMap;
use nalgebra::Vector3;
use crate::mesh::{Mesh, Node, Tet4Mesh, Tet10Mesh};

/// Convert a Tet4Mesh to a Tet10Mesh by subdividing each edge.
/// For each unique edge (a, b), a midside node is inserted at (p_a + p_b) / 2.
/// Edges are deduplicated via a HashMap keyed on (min(a,b), max(a,b)).
///
/// Node ordering follows basix convention:
///   [0,1,2,3] — corner nodes (same as Tet4, reoriented if needed)
///   [4] — midside of edge (2,3)
///   [5] — midside of edge (1,3)
///   [6] — midside of edge (1,2)
///   [7] — midside of edge (0,3)
///   [8] — midside of edge (0,2)
///   [9] — midside of edge (0,1)
///
/// Elements with negative orientation (signed_vol < 0) have n0 and n1
/// swapped to ensure a positive jacobian determinant in Tet10::precompute.
pub fn tet4_to_tet10(tet4: &Tet4Mesh) -> Tet10Mesh {
    let mut nodes: Vec<Node> = tet4.nodes.iter().map(|n| Node {
        position: n.position,
        velocity: Vector3::zeros(),
        mass:     n.mass,
    }).collect();

    let mut edge_map: HashMap<(usize, usize), usize> = HashMap::new();
    let mut elements: Vec<[usize; 10]>               = Vec::with_capacity(tet4.elements.len());

    for &[n0, n1, n2, n3] in &tet4.elements {
        let p0 = tet4.nodes[n0].position;
        let p1 = tet4.nodes[n1].position;
        let p2 = tet4.nodes[n2].position;
        let p3 = tet4.nodes[n3].position;

        // Ensure positive orientation for Tet10 jacobian.
        // Swap n0 and n1 if element has negative signed volume.
        let signed_vol = (p1 - p0).dot(&(p2 - p0).cross(&(p3 - p0)));
        let (a, b, c, d) = if signed_vol > 0.0 {
            (n0, n1, n2, n3)
        } else {
            (n1, n0, n2, n3)
        };

        let m_ab = get_or_insert_midnode(&mut nodes, &mut edge_map, a, b);
        let m_bc = get_or_insert_midnode(&mut nodes, &mut edge_map, b, c);
        let m_ac = get_or_insert_midnode(&mut nodes, &mut edge_map, a, c);
        let m_ad = get_or_insert_midnode(&mut nodes, &mut edge_map, a, d);
        let m_bd = get_or_insert_midnode(&mut nodes, &mut edge_map, b, d);
        let m_cd = get_or_insert_midnode(&mut nodes, &mut edge_map, c, d);

        // Basix ordering: [a,b,c,d, mid(c,d), mid(b,d), mid(b,c), mid(a,d), mid(a,c), mid(a,b)]
        elements.push([a, b, c, d, m_cd, m_bd, m_bc, m_ad, m_ac, m_ab]);
    }

    Mesh { nodes, elements }
}

/// Get the index of the midside node on edge (a, b), inserting it if needed.
/// Key is always (min, max) to avoid duplicate entries for the same edge.
fn get_or_insert_midnode(
    nodes:    &mut Vec<Node>,
    edge_map: &mut HashMap<(usize, usize), usize>,
    a: usize,
    b: usize,
) -> usize {
    let key = if a < b { (a, b) } else { (b, a) };
    if let Some(&idx) = edge_map.get(&key) {
        return idx;
    }
    let mid  = (nodes[a].position + nodes[b].position) * 0.5;
    let mass = (nodes[a].mass + nodes[b].mass) * 0.5;
    let idx  = nodes.len();
    nodes.push(Node { position: mid, velocity: Vector3::zeros(), mass });
    edge_map.insert(key, idx);
    idx
}