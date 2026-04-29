// UTF-8
// orfas-viewer/src/render.rs — node picking and depth sorting.

use orfas_core::element::FemMesh;
use crate::state::Camera;

/// Returns the index of the mesh node closest to the given screen position.
/// Delegates projection to Camera::project.
pub fn screen_to_node(
    pos:    &egui::Pos2,
    mesh:   &impl FemMesh,
    camera: &Camera,
    center: egui::Pos2,
    scale:  f32,
) -> usize {
    mesh.nodes().iter().enumerate()
        .map(|(i, node)| {
            let p  = camera.project(&node.position, center, scale);
            let dx = p.x - pos.x;
            let dy = p.y - pos.y;
            (i, dx * dx + dy * dy)
        })
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Returns node indices sorted back-to-front (painter's algorithm).
/// Nodes with larger depth (further from camera) are drawn first.
pub fn depth_sorted_nodes(mesh: &impl FemMesh, camera: &Camera) -> Vec<usize> {
    let forward = camera.view_direction();
    let mut indices: Vec<usize> = (0..mesh.n_nodes()).collect();
    indices.sort_by(|&a, &b| {
        let pa = mesh.nodes()[a].position.cast::<f32>() - camera.target;
        let pb = mesh.nodes()[b].position.cast::<f32>() - camera.target;
        let da = pa.dot(&forward);
        let db = pb.dot(&forward);
        da.partial_cmp(&db).unwrap()
    });
    indices
}