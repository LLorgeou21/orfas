use orfas_core::mesh::{Mesh, Tetrahedron,Node};
use nalgebra::Vector3;

#[derive(Debug)]
pub enum IoError {
    UnreadableFile,
    InvalidFormat,
    UnsupportedCellType,
    IndexOutOfBounds,
    InvalidCellCount,
}

enum ParseState {
    Header,
    Points,
    Cells,
    CellTypes
}


pub fn read_vtk(path : &str)->Result<Mesh,IoError>{
    let mut nodes = Vec::new();
    let mut elements = Vec::new();
    match std::fs::read_to_string(path) {
        Ok(lines) => {
            let mut parsestate = ParseState::Header;
            if !lines.lines().nth(0).ok_or(IoError::InvalidFormat)?.contains("# vtk"){return Err(IoError::InvalidFormat)}
            if !lines.lines().nth(2).ok_or(IoError::InvalidFormat)?.contains("ASCII"){return Err(IoError::InvalidFormat)}
            if !lines.lines().nth(3).ok_or(IoError::InvalidFormat)?.contains("DATASET UNSTRUCTURED_GRID"){return Err(IoError::InvalidFormat)}
            for line in lines.lines() {
                if line.contains("CELLS"){parsestate = ParseState::Cells;continue;}
                if line.contains("POINTS"){parsestate = ParseState::Points;continue;}
                if line.contains("CELL_TYPES"){parsestate = ParseState::CellTypes;continue;}
                match parsestate{
                    ParseState::CellTypes =>{},
                    ParseState::Points =>{
                        let mut coords = Vec::new();
                        for p in line.split(" "){
                            coords.push(p.parse::<f64>().map_err(|_| IoError::InvalidFormat)?);
                        }
                        if coords.len()!=3{return Err(IoError::InvalidFormat);}
                        let position = Vector3::new(coords[0], coords[1], coords[2]);
                        nodes.push(Node { position, velocity : Vector3::zeros(), mass : 1. });
                    },
                    ParseState::Cells =>{
                        let mut indices_nodes = Vec::new();
                        for p in line.split(" "){
                            indices_nodes.push(p.parse::<usize>().map_err(|_| IoError::InvalidFormat)?);
                        }
                        if indices_nodes[0]!=4{return Err(IoError::UnsupportedCellType);}
                        if indices_nodes.len()!=5{return Err(IoError::InvalidFormat);}
                        elements.push(Tetrahedron { indices: [indices_nodes[1], indices_nodes[2], indices_nodes[3], indices_nodes[4]] });
                    },
                    ParseState::Header=>{}
                }
            }
        },  
        Err(_)=>{return Err(IoError::UnreadableFile);}
    }
    Ok(Mesh { nodes, elements })
}


#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn test_read_vtk_valid() {
        // Load the reference cube mesh (8 nodes, 6 tetrahedra)
        // This is the same geometry as Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0)
        let mesh = read_vtk("tests/data/test_mesh.vtk").expect("Failed to load test mesh");
        assert_eq!(mesh.nodes.len(), 8, "Expected 8 nodes");
        assert_eq!(mesh.elements.len(), 6, "Expected 6 tetrahedra");
    }
 
    #[test]
    fn test_read_vtk_node_positions() {
        // Verify that node positions are parsed correctly
        // Node 0 should be at (0, 0, 0) and node 7 at (1, 1, 1)
        let mesh = read_vtk("tests/data/test_mesh.vtk").expect("Failed to load test mesh");
        let p0 = mesh.nodes[0].position;
        let p7 = mesh.nodes[7].position;
        assert!((p0.x - 0.0).abs() < 1e-10 && (p0.y - 0.0).abs() < 1e-10 && (p0.z - 0.0).abs() < 1e-10, "Node 0 should be at origin");
        assert!((p7.x - 1.0).abs() < 1e-10 && (p7.y - 1.0).abs() < 1e-10 && (p7.z - 1.0).abs() < 1e-10, "Node 7 should be at (1,1,1)");
    }
 
    #[test]
    fn test_read_vtk_file_not_found() {
        // A non-existent file must return IoError::UnreadableFile
        let result = read_vtk("tests/data/does_not_exist.vtk");
        assert!(matches!(result, Err(IoError::UnreadableFile)));
    }
 
    #[test]
    fn test_read_vtk_invalid_format() {
        // A file that does not start with "# vtk" must return IoError::InvalidFormat
        let result = read_vtk("tests/data/invalid_format.vtk");
        assert!(matches!(result, Err(IoError::InvalidFormat)));
    }
}
