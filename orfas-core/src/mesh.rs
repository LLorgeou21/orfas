// UTF-8
// orfas-core/src/mesh.rs

use nalgebra::Vector3;
use crate::element::traits::FemMesh;

// ---------------------------------------------------------------------------
// Node
// ---------------------------------------------------------------------------

/// A node in the FEM mesh.
/// Stores the 3D reference position, velocity, and mass.
pub struct Node {
    pub position: Vector3<f64>,
    pub velocity: Vector3<f64>,
    pub mass:     f64,
}

// ---------------------------------------------------------------------------
// Mesh<N>
// ---------------------------------------------------------------------------

/// A FEM mesh composed of nodes and elements.
/// Generic over N: the number of nodes per element.
/// Each element stores exactly N indices into the node list.
/// Use type aliases Tet4Mesh / Tet10Mesh for clarity at call sites.
pub struct Mesh<const N: usize> {
    pub nodes:    Vec<Node>,
    pub elements: Vec<[usize; N]>,
}

/// Linear tetrahedral mesh (4 nodes per element).
pub type Tet4Mesh = Mesh<4>;

/// Quadratic tetrahedral mesh (10 nodes per element).
pub type Tet10Mesh = Mesh<10>;

/// Linear hexahedral mesh (8 nodes per element).
pub type Hex8Mesh = Mesh<8>;

pub type Beam2Mesh = Mesh<2>;

// ---------------------------------------------------------------------------
// FemMesh implementations
// ---------------------------------------------------------------------------

impl FemMesh for Mesh<4> {
    fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// Returns connectivity as Vec<Vec<usize>> for Assembler consumption.
    fn connectivity(&self) -> Vec<Vec<usize>> {
        self.elements.iter().map(|e| e.to_vec()).collect()
    }
}

impl FemMesh for Mesh<10> {
    fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    fn connectivity(&self) -> Vec<Vec<usize>> {
        self.elements.iter().map(|e| e.to_vec()).collect()
    }
}

impl FemMesh for Mesh<8> {
    fn nodes(&self) -> &[Node] {
        &self.nodes
    }
 
    /// Returns connectivity as Vec<Vec<usize>> for Assembler consumption.
    fn connectivity(&self) -> Vec<Vec<usize>> {
        self.elements.iter().map(|e| e.to_vec()).collect()
    }
}

impl FemMesh for Mesh<2> {
    fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    fn connectivity(&self) -> Vec<Vec<usize>> {
        self.elements.iter().map(|e| e.to_vec()).collect()
    }
}

// ---------------------------------------------------------------------------
// Mesh<4> constructors
// ---------------------------------------------------------------------------

impl Mesh<4> {
    /// Generate a structured Tet4 mesh.
    /// Node index formula: i + j*nx + k*nx*ny
    /// where i in 0..nx (x-axis), j in 0..ny (y-axis), k in 0..nz (z-axis).
    /// Loop order matches the formula: k outermost, j middle, i innermost.
    /// Each cube is split into 6 tetrahedra using alternating decomposition
    /// based on (i+j+k) % 2 to avoid degenerate elements.
    pub fn generate(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Tet4Mesh {
        let mut nodes = Vec::new();
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let position = Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    nodes.push(Node { position, velocity: Vector3::zeros(), mass: 1.0 });
                }
            }
        }

        let mut elements = Vec::new();
        for i in 0..nx - 1 {
            for j in 0..ny - 1 {
                for k in 0..nz - 1 {
                    let p0 = i       + j       * nx + k       * nx * ny;
                    let p1 = (i + 1) + j       * nx + k       * nx * ny;
                    let p2 = i       + (j + 1) * nx + k       * nx * ny;
                    let p3 = (i + 1) + (j + 1) * nx + k       * nx * ny;
                    let p4 = i       + j       * nx + (k + 1) * nx * ny;
                    let p5 = (i + 1) + j       * nx + (k + 1) * nx * ny;
                    let p6 = i       + (j + 1) * nx + (k + 1) * nx * ny;
                    let p7 = (i + 1) + (j + 1) * nx + (k + 1) * nx * ny;

                    if (i + j + k) % 2 == 0 {
                        elements.push([p0, p1, p3, p7]);
                        elements.push([p0, p1, p5, p7]);
                        elements.push([p0, p2, p3, p7]);
                        elements.push([p0, p2, p6, p7]);
                        elements.push([p0, p4, p5, p7]);
                        elements.push([p0, p4, p6, p7]);
                    } else {
                        elements.push([p1, p0, p4, p6]);
                        elements.push([p1, p0, p2, p6]);
                        elements.push([p1, p5, p4, p6]);
                        elements.push([p1, p5, p7, p6]);
                        elements.push([p1, p3, p2, p6]);
                        elements.push([p1, p3, p7, p6]);
                    }
                }
            }
        }

        Mesh { nodes, elements }
    }
}

impl Mesh<8> {
    /// Generate a structured Hex8 mesh filling a box of size dx*dy*dz.
    /// Produces (nx-1)*(ny-1)*(nz-1) hexahedral elements.
    /// Node index formula: i + j*nx + k*nx*ny
    /// where i in 0..nx, j in 0..ny, k in 0..nz.
    ///
    /// Node ordering per element follows VTK Hex8 convention:
    ///   bottom face (k) CCW: p0,p1,p2,p3 then top face (k+1) CCW: p4,p5,p6,p7.
    /// This matches Hex8::NODE_RST and ensures positive Jacobian determinants.
    pub fn generate(nx: usize, ny: usize, nz: usize, dx: f64, dy: f64, dz: f64) -> Hex8Mesh {
        // Build node list in the same i+j*nx+k*nx*ny order as Tet4Mesh
        let mut nodes = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let position = Vector3::new(i as f64 * dx, j as f64 * dy, k as f64 * dz);
                    nodes.push(Node { position, velocity: Vector3::zeros(), mass: 1.0 });
                }
            }
        }
 
        // Build element connectivity: one hex per (i,j,k) cell
        let mut elements = Vec::with_capacity((nx - 1) * (ny - 1) * (nz - 1));
        for k in 0..nz - 1 {
            for j in 0..ny - 1 {
                for i in 0..nx - 1 {
                    // 8 corners of the cell, bottom face first (CCW from -z side)
                    let p0 = i       + j       * nx + k       * nx * ny;
                    let p1 = (i + 1) + j       * nx + k       * nx * ny;
                    let p2 = (i + 1) + (j + 1) * nx + k       * nx * ny;
                    let p3 = i       + (j + 1) * nx + k       * nx * ny;
                    let p4 = i       + j       * nx + (k + 1) * nx * ny;
                    let p5 = (i + 1) + j       * nx + (k + 1) * nx * ny;
                    let p6 = (i + 1) + (j + 1) * nx + (k + 1) * nx * ny;
                    let p7 = i       + (j + 1) * nx + (k + 1) * nx * ny;
 
                    elements.push([p0, p1, p2, p3, p4, p5, p6, p7]);
                }
            }
        }
 
        Mesh { nodes, elements }
    }
}


// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_mesh() {
        let mesh = Tet4Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements.len(), 6);
    }

    #[test]
    fn test_generate_mesh2() {
        let mesh = Tet4Mesh::generate(3, 3, 3, 1.0, 1.0, 1.0);
        assert_eq!(mesh.nodes.len(), 27);
        assert_eq!(mesh.elements.len(), 48);
    }

    #[test]
    fn test_node_positions() {
        // Verify that node index formula i + j*nx + k*nx*ny matches actual positions
        let nx = 3; let ny = 3; let nz = 3;
        let mesh = Tet4Mesh::generate(nx, ny, nz, 1.0, 1.0, 1.0);
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

    #[test]
    fn test_fem_mesh_trait() {
        let mesh = Tet4Mesh::generate(2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(mesh.n_nodes(), 8);
        assert_eq!(mesh.n_elements(), 6);
        assert_eq!(mesh.connectivity()[0].len(), 4);
    }

    #[test]
    fn test_generate_hex8_mesh() {
        let mesh = Mesh::<8>::generate(2, 2, 2, 1.0, 1.0, 1.0);
        // 2x2x2 nodes, 1x1x1 element
        assert_eq!(mesh.nodes.len(), 8);
        assert_eq!(mesh.elements.len(), 1);
    }
 
    #[test]
    fn test_generate_hex8_mesh_3x3x3() {
        let mesh = Mesh::<8>::generate(3, 3, 3, 1.0, 1.0, 1.0);
        // 3x3x3 nodes, 2x2x2 = 8 elements
        assert_eq!(mesh.nodes.len(), 27);
        assert_eq!(mesh.elements.len(), 8);
    }
 
    #[test]
    fn test_hex8_node_positions() {
        let nx = 3; let ny = 3; let nz = 3;
        let mesh = Mesh::<8>::generate(nx, ny, nz, 1.0, 1.0, 1.0);
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
 
    #[test]
    fn test_hex8_fem_mesh_trait() {
        let mesh = Mesh::<8>::generate(2, 2, 2, 1.0, 1.0, 1.0);
        assert_eq!(mesh.n_nodes(), 8);
        assert_eq!(mesh.n_elements(), 1);
        assert_eq!(mesh.connectivity()[0].len(), 8);
    }
}