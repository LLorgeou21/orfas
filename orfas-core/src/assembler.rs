use nalgebra::{DMatrix, Matrix3, Matrix6, SMatrix, Vector3};
use crate::material::MaterialLaw;
use crate::mesh::Mesh;
type Matrix6x12 = SMatrix<f64, 6, 12>;

/// Trait defining how the strain-displacement matrix B is computed.
/// Linear elasticity uses a constant B matrix.
/// Non-linear laws (v0.4) will override this with a deformation-dependent B.
pub trait BMatrix {
    fn compute(
        p0: &Vector3<f64>, p1: &Vector3<f64>,
        p2: &Vector3<f64>, p3: &Vector3<f64>,
        deformation: Option<&Matrix6<f64>>
    ) -> Matrix6x12;
}

/// Linear implementation of BMatrix — B is constant, deformation is ignored.
pub struct LinearBMatrix;

impl BMatrix for LinearBMatrix {
    fn compute(
        p0: &Vector3<f64>, p1: &Vector3<f64>,
        p2: &Vector3<f64>, p3: &Vector3<f64>,
        _deformation: Option<&Matrix6<f64>>
    ) -> Matrix6x12 {
        let volume = tetra_volume(p0, p1, p2, p3);
        let (b, c, d) = tetra_bcd(p0, p1, p2, p3);
        tetra_b_matrix(&b, &c, &d, volume)
    }
}

/// Compute the volume of a tetrahedron using the formula V = |det([p1-p0, p2-p0, p3-p0])| / 6
fn tetra_volume(p0: &Vector3<f64>, p1: &Vector3<f64>, p2: &Vector3<f64>, p3: &Vector3<f64>) -> f64 {
    Matrix3::from_columns(&[p1-p0, p2-p0, p3-p0]).determinant().abs() / 6.0
}

/// Computes the shape function derivatives (b, c, d) for a linear tetrahedron.
/// These coefficients are the partial derivatives of the 4 shape functions
/// with respect to x (b), y (c), and z (d).
/// They are constant throughout the element — a key property of linear tetrahedra.
/// The 4th coefficient is derived from the partition of unity property:
/// the sum of all shape functions equals 1, so the sum of their derivatives equals 0,
/// giving b4 = -(b1+b2+b3), and similarly for c and d.
fn tetra_bcd(p0: &Vector3<f64>, p1: &Vector3<f64>, p2: &Vector3<f64>, p3: &Vector3<f64>) -> ([f64; 4], [f64; 4], [f64; 4]) {
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

fn tetra_b_matrix(b: &[f64; 4], c: &[f64; 4], d: &[f64; 4], volume: f64) -> Matrix6x12 {
    Matrix6x12::from_row_slice(&[
        b[0],0.0,0.0,b[1],0.0,0.0,b[2],0.0,0.0,b[3],0.0,0.0,
        0.0,c[0],0.0,0.0,c[1],0.0,0.0,c[2],0.0,0.0,c[3],0.0,
        0.0,0.0,d[0],0.0,0.0,d[1],0.0,0.0,d[2],0.0,0.0,d[3],
        c[0],b[0],0.0,c[1],b[1],0.0,c[2],b[2],0.0,c[3],b[3],0.0,
        0.0,d[0],c[0],0.0,d[1],c[1],0.0,d[2],c[2],0.0,d[3],c[3],
        d[0],0.0,b[0],d[1],0.0,b[1],d[2],0.0,b[2],d[3],0.0,b[3],
    ]).scale(1.0 / (6.0 * volume))
}

/// Computes the local stiffness matrix Ke = Bt * C * B * V for a tetrahedron.
/// Generic over the B matrix computation strategy via the BMatrix trait.
pub fn tetra_stiffness_matrix<B: BMatrix>(
    p0: &Vector3<f64>, p1: &Vector3<f64>,
    p2: &Vector3<f64>, p3: &Vector3<f64>,
    material: &dyn MaterialLaw,
    deformation: Option<&Matrix6<f64>>
) -> SMatrix<f64, 12, 12> {
    let volume = tetra_volume(p0, p1, p2, p3);
    if volume < 1e-10 {
        return SMatrix::zeros();
    }
    let b_mat = B::compute(p0, p1, p2, p3, deformation);
    let c_mat = material.stiffness_matrix(deformation);
    b_mat.transpose() * c_mat * b_mat * volume
}

pub struct Assembler {
    /// Connectivity table: for each tetrahedron, the 4 global node indices
    connectivity: Vec<[usize; 4]>,
}


impl Assembler {

    /// Copy the indices of the tetrahedron to build the connectivity matrix
    pub fn new(mesh: &Mesh) -> Assembler {
    let connectivity = mesh.elements.iter().map(|tetra| tetra.indices).collect();
    Assembler { connectivity }
    }

    // Assemble the matrix 
    pub fn assemble<B: BMatrix>(&self,mesh: &Mesh,material: &dyn MaterialLaw,) -> DMatrix<f64> {
        let n = 3 * mesh.nodes.len();
        let mut k : DMatrix<f64> = DMatrix::zeros(n, n);
        for (i,_) in mesh.elements.iter().enumerate(){
            let [a,b,c,d] = self.connectivity[i];
            let p0 = &mesh.nodes[a].position;
            let p1 = &mesh.nodes[b].position;
            let p2 = &mesh.nodes[c].position;
            let p3 = &mesh.nodes[d].position;
            let ke =tetra_stiffness_matrix::<B>(p0, p1, p2, p3, material, None);
            let global_indices = [a, b, c, d];
            for r in 0..4 {
                for s in 0..4 {
                    let global_r = global_indices[r];
                    let global_s = global_indices[s];
                    let mut block = k.fixed_view_mut::<3, 3>(3 * global_r, 3 * global_s);
                    block += ke.fixed_view::<3, 3>(3 * r, 3 * s);
                }
            }
        }
        k
        }
        

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_volume_calculation() {
        let v = tetra_volume(
            &Vector3::new(0.0, 0.0, 0.0),
            &Vector3::new(1.0, 0.0, 0.0),
            &Vector3::new(0.0, 1.0, 0.0),
            &Vector3::new(0.0, 0.0, 1.0)
        );
        assert_eq!(v, 1.0 / 6.0);
    }
}