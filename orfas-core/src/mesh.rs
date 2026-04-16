use nalgebra::Vector3;



/// A node in the FEM mesh.
/// Stores the 3D position in world space the velocity
/// and the mass
pub struct Node {
    pub position  : Vector3<f64>,
    pub velocity  : Vector3<f64>,
    pub mass      : f64
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
    /// Node index formula: i + j*nx + k*nx*ny
    /// where i in 0..nx (x-axis), j in 0..ny (y-axis), k in 0..nz (z-axis)
    /// Loop order matches the formula: k outermost, j middle, i innermost
    pub fn generate(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Mesh {
        let mut nodes = Vec::new();
        // Loop order matches index formula: i + j*nx + k*nx*ny
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let position = Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    nodes.push(Node { position, velocity : Vector3::zeros(), mass : 1. });
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

    #[test]
    fn test_node_positions() {
        // Verify that node index formula i + j*nx + k*nx*ny matches actual positions
        let nx = 3; let ny = 3; let nz = 3;
        let mesh = Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let idx = i + j * nx + k * nx * ny;
                    let pos = &mesh.nodes[idx].position;
                    assert_eq!(pos.x, i as f64, "bad x at idx={}", idx);
                    assert_eq!(pos.y, j as f64, "bad y at idx={}", idx);
                    assert_eq!(pos.z, k as f64, "bad z at idx={}", idx);
                }
            }
        }
    }
}