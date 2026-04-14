use nalgebra::Vector3;



/// A node in the FEM mesh.
/// Stores the 3D position in world space and whether
/// the node is fixed (zero displacement boundary condition).
pub struct Node {
    pub position  : Vector3<f64>,
    pub fixed : bool
}

/// A tetrahedral element in the FEM mesh.
/// Defined by 4 indices pointing into the mesh's node list.
/// Positions are not stored directly — this avoids duplication
/// and ensures shared nodes between adjacent elements stay consistent.
pub struct Tetrahedron {
    pub indices: [usize; 4],
}

/// A FEM mesh composed of nodes and tetrahedral elements.
/// Nodes store 3D positions, elements store indices into the node list.
/// Together they define the geometry and topology of the simulated object.
pub struct Mesh {
    pub nodes: Vec<Node>,
    pub elements: Vec<Tetrahedron>,
}

impl Mesh {

    /// Fn to generate a mesh for FEM
    /// Used in early version for debugging and to speed the dev
    /// Define the mesh by giving the numbers of nodes in each dimensions and the distance between them in each dimension
    /// No point is fixed by default
    /// Uses an alternating 5-tetrahedra decomposition per cube to avoid degenerate elements
    pub fn generate(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Mesh {
        let mut nodes = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    let position = Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    nodes.push(Node { position, fixed: false });
                }
            }
        }
        let mut tetrahedron = Vec::new();
        for i in 0..nx-1 {
            for j in 0..ny-1 {
                for k in 0..nz-1 {
                    let p0 = i       + j     * nx + k     * nx * ny;
                    let p1 = (i+1)   + j     * nx + k     * nx * ny;
                    let p2 = i       + (j+1) * nx + k     * nx * ny;
                    let p3 = (i+1)   + (j+1) * nx + k     * nx * ny;
                    let p4 = i       + j     * nx + (k+1) * nx * ny;
                    let p5 = (i+1)   + j     * nx + (k+1) * nx * ny;
                    let p6 = i       + (j+1) * nx + (k+1) * nx * ny;
                    let p7 = (i+1)   + (j+1) * nx + (k+1) * nx * ny;

                    if (i + j + k) % 2 == 0 {
                        tetrahedron.push(Tetrahedron { indices: [p0, p1, p3, p7] });
                        tetrahedron.push(Tetrahedron { indices: [p0, p1, p5, p7] });
                        tetrahedron.push(Tetrahedron { indices: [p0, p2, p3, p7] });
                        tetrahedron.push(Tetrahedron { indices: [p0, p2, p6, p7] });
                        tetrahedron.push(Tetrahedron { indices: [p0, p4, p5, p7] });
                        tetrahedron.push(Tetrahedron { indices: [p0, p4, p6, p7] });
                    } else {
                        tetrahedron.push(Tetrahedron { indices: [p1, p0, p4, p6] });
                        tetrahedron.push(Tetrahedron { indices: [p1, p0, p2, p6] });
                        tetrahedron.push(Tetrahedron { indices: [p1, p5, p4, p6] });
                        tetrahedron.push(Tetrahedron { indices: [p1, p5, p7, p6] });
                        tetrahedron.push(Tetrahedron { indices: [p1, p3, p2, p6] });
                        tetrahedron.push(Tetrahedron { indices: [p1, p3, p7, p6] });
                    }
                }
            }
        }
        Mesh { nodes, elements: tetrahedron }
    }

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_mesh() {
        let mesh = Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements.len(), 6);
    }
    #[test]
    fn test_generate_mesh2() {
        let mesh = Mesh::generate(3, 3, 3, 1.0, 1.0, 1.0);
        assert_eq!(mesh.nodes.len(), 27);
        assert_eq!(mesh.elements.len(), 48);
    }
}

