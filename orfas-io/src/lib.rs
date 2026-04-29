// UTF-8
// orfas-io/src/lib.rs — VTK mesh reader for Tet4 meshes.
// Supports ASCII VTK unstructured grid format with tetrahedral elements (cell type 10).
// Tet10 VTK loading is planned for v0.9.0.

use orfas_core::mesh::{Node, Tet4Mesh};
use nalgebra::Vector3;

/// Errors that can occur when reading a VTK file.
#[derive(Debug)]
pub enum IoError {
    /// File could not be read from disk.
    UnreadableFile,
    /// File does not conform to the expected VTK ASCII format.
    InvalidFormat,
    /// Cell type is not a linear tetrahedron (type 10). Tet10 support planned for v0.9.0.
    UnsupportedCellType,
    /// A node index in the cell list exceeds the number of nodes.
    IndexOutOfBounds,
    /// Cell count header does not match the actual number of cells parsed.
    InvalidCellCount,
}

/// Internal parser state machine for VTK ASCII format.
enum ParseState {
    Header,
    Points,
    Cells,
    CellTypes,
}

/// Read a VTK ASCII unstructured grid file and return a Tet4Mesh.
///
/// Expected format:
/// - Line 0: `# vtk DataFile Version ...`
/// - Line 2: `ASCII`
/// - Line 3: `DATASET UNSTRUCTURED_GRID`
/// - `POINTS n float` section with one `x y z` per line
/// - `CELLS n m` section with one `4 i0 i1 i2 i3` per line (Tet4 only)
///
/// # Errors
/// Returns `IoError::UnsupportedCellType` if any cell has a count != 4.
/// Tet10 meshes (10 nodes per cell) are not yet supported — use v0.9.0.
pub fn read_vtk(path: &str) -> Result<Tet4Mesh, IoError> {
    let mut nodes    = Vec::new();
    let mut elements = Vec::new();

    match std::fs::read_to_string(path) {
        Ok(lines) => {
            let mut parse_state = ParseState::Header;

            if !lines.lines().nth(0).ok_or(IoError::InvalidFormat)?.contains("# vtk") {
                return Err(IoError::InvalidFormat);
            }
            if !lines.lines().nth(2).ok_or(IoError::InvalidFormat)?.contains("ASCII") {
                return Err(IoError::InvalidFormat);
            }
            if !lines.lines().nth(3).ok_or(IoError::InvalidFormat)?.contains("DATASET UNSTRUCTURED_GRID") {
                return Err(IoError::InvalidFormat);
            }

            for line in lines.lines() {
                if line.contains("CELLS")      { parse_state = ParseState::Cells;     continue; }
                if line.contains("POINTS")     { parse_state = ParseState::Points;    continue; }
                if line.contains("CELL_TYPES") { parse_state = ParseState::CellTypes; continue; }

                match parse_state {
                    ParseState::CellTypes => {}

                    ParseState::Points => {
                        let mut coords = Vec::new();
                        for p in line.split_whitespace() {
                            coords.push(p.parse::<f64>().map_err(|_| IoError::InvalidFormat)?);
                        }
                        if coords.len() != 3 { return Err(IoError::InvalidFormat); }
                        nodes.push(Node {
                            position: Vector3::new(coords[0], coords[1], coords[2]),
                            velocity: Vector3::zeros(),
                            mass:     1.0,
                        });
                    }

                    ParseState::Cells => {
                        let mut idx = Vec::new();
                        for p in line.split_whitespace() {
                            idx.push(p.parse::<usize>().map_err(|_| IoError::InvalidFormat)?);
                        }
                        // VTK Tet4: first value is node count = 4, followed by 4 indices
                        if idx[0] != 4 { return Err(IoError::UnsupportedCellType); }
                        if idx.len() != 5 { return Err(IoError::InvalidFormat); }
                        elements.push([idx[1], idx[2], idx[3], idx[4]]);
                    }

                    ParseState::Header => {}
                }
            }
        }
        Err(_) => return Err(IoError::UnreadableFile),
    }

    Ok(Tet4Mesh { nodes, elements })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_vtk_valid() {
        // Load the reference cube mesh (8 nodes, 6 tetrahedra)
        // Same geometry as Tet4Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0)
        let mesh = read_vtk("tests/data/test_mesh.vtk").expect("Failed to load test mesh");
        assert_eq!(mesh.nodes.len(),    8, "Expected 8 nodes");
        assert_eq!(mesh.elements.len(), 6, "Expected 6 tetrahedra");
    }

    #[test]
    fn test_read_vtk_node_positions() {
        // Node 0 should be at (0,0,0) and node 7 at (1,1,1)
        let mesh = read_vtk("tests/data/test_mesh.vtk").expect("Failed to load test mesh");
        let p0 = mesh.nodes[0].position;
        let p7 = mesh.nodes[7].position;
        assert!(p0.norm() < 1e-10,             "Node 0 should be at origin");
        assert!((p7 - Vector3::new(1.0, 1.0, 1.0)).norm() < 1e-10, "Node 7 should be at (1,1,1)");
    }

    #[test]
    fn test_read_vtk_file_not_found() {
        let result = read_vtk("tests/data/does_not_exist.vtk");
        assert!(matches!(result, Err(IoError::UnreadableFile)));
    }

    #[test]
    fn test_read_vtk_invalid_format() {
        let result = read_vtk("tests/data/invalid_format.vtk");
        assert!(matches!(result, Err(IoError::InvalidFormat)));
    }
}