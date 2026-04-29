// UTF-8
// orfas-core/src/element/traits.rs

use nalgebra::{DMatrix, DVector, Matrix3, Vector3};
use nalgebra::Vector6;
use crate::mesh::Node;

// ---------------------------------------------------------------------------
// DofType trait
// ---------------------------------------------------------------------------

/// Describes the degree-of-freedom layout for a class of finite elements.
///
/// This trait is the ORFAS equivalent of SOFA's DataTypes template:
/// it encodes the number of DOFs per node and the associated coordinate
/// and derivative types at compile time, with zero runtime cost.
///
/// Two concrete implementations are provided:
/// - [`Vec3Dof`] — 3 DOFs per node (translations only). Used by all
///   volumetric elements: Tet4, Tet10, Hex8.
/// - [`Vec6Dof`] — 6 DOFs per node (3 translations + 3 rotations).
///   Used by structural elements: Beam2, Shell.
///
/// Downstream types (`MechanicalState`, `Assembler`, `BoundaryConditions`,
/// `ImplicitEulerIntegrator`) are generic over `D: DofType` and use
/// `D::N_DOF` wherever the number of DOFs per node is needed, replacing
/// all hard-coded `3`s in the previous architecture.
pub trait DofType: Send + Sync + 'static {
    /// Number of degrees of freedom per node.
    /// Used by the assembler to compute global matrix indices:
    /// node i occupies rows/cols [i*N_DOF .. i*N_DOF + N_DOF].
    const N_DOF: usize;

    /// Coordinate type for one node: the state vector entry.
    /// Vec3Dof -> Vector3<f64> (position in 3D).
    /// Vec6Dof -> Vector6<f64> (position + rotation angles).
    type Coord: Send + Sync + Clone;

    /// Derivative type for one node: velocity / angular velocity.
    /// Same dimension as Coord for all current implementations.
    type Deriv: Send + Sync + Clone;
}

// ---------------------------------------------------------------------------
// Concrete DofType implementations
// ---------------------------------------------------------------------------

/// 3 DOFs per node — translations only.
/// Used by all volumetric elements (Tet4, Tet10, Hex8).
/// Coord and Deriv are both Vector3<f64>.
pub struct Vec3Dof;

impl DofType for Vec3Dof {
    const N_DOF: usize = 3;
    type Coord = Vector3<f64>;
    type Deriv = Vector3<f64>;
}

/// 6 DOFs per node — 3 translations + 3 rotations.
/// Used by structural elements: Beam2, Shell (Mindlin-Reissner / DKT).
/// Coord and Deriv are both Vector6<f64>.
///
/// This type covers small-rotation formulations where rotational DOFs
/// are additive (q += dq), which is valid for beams and shells under
/// the assumptions of structural mechanics.
///
/// Large-rotation rigid body dynamics (v0.13.x) will introduce a
/// separate `RigidDof` type using SO(3) / quaternion kinematics
/// (7 DOFs per node: 3 translations + 1 unit quaternion).
/// `Vec6Dof` will NOT be modified — the two types coexist, following
/// the same pattern as SOFA's Vec6d and Rigid3d DataTypes.
pub struct Vec6Dof;

impl DofType for Vec6Dof {
    const N_DOF: usize = 6;
    type Coord = Vector6<f64>;
    type Deriv = Vector6<f64>;
}

// ---------------------------------------------------------------------------
// FemMesh trait
// ---------------------------------------------------------------------------

/// Abstraction over any finite element mesh.
/// Allows Assembler<E> to work with Mesh<4>, Mesh<10>, and any future
/// Mesh<N> without depending on the const generic N directly.
/// Implement this for each Mesh<N> variant used in practice.
pub trait FemMesh: Sync {
    /// Slice of all nodes in the mesh.
    fn nodes(&self) -> &[Node];

    /// Connectivity: each entry is a slice of node indices for one element.
    fn connectivity(&self) -> Vec<Vec<usize>>;

    /// Total number of nodes.
    fn n_nodes(&self) -> usize {
        self.nodes().len()
    }

    /// Total number of elements.
    fn n_elements(&self) -> usize {
        self.connectivity().len()
    }
}

// ---------------------------------------------------------------------------
// ElementGeometry marker trait
// ---------------------------------------------------------------------------

/// Marker trait for precomputed per-element geometric data.
/// Each FiniteElement implementation defines its own Geometry type
/// storing exactly what it needs — no more.
/// Tet4: analytical gradients (b, c, d) + volume.
/// Tet10: Gauss-point gradients + jacobian determinants.
pub trait ElementGeometry: Send + Sync {}

// ---------------------------------------------------------------------------
// FiniteElement trait
// ---------------------------------------------------------------------------

/// Core abstraction for a finite element type.
/// All element-specific logic lives here: shape functions,
/// gradients, integration points, B-matrix construction,
/// and volume/jacobian queries.
/// The Assembler is generic over E: FiniteElement and knows
/// nothing about the internals of each element.
pub trait FiniteElement: Send + Sync + 'static {
    /// Precomputed geometric data for this element type.
    type Geometry: ElementGeometry;

    /// DOF layout for this element type.
    /// Volumetric elements (Tet4, Tet10, Hex8) use Vec3Dof.
    /// Structural elements (Beam2, Shell) use Vec6Dof.
    /// All index arithmetic in the assembler derives from Dof::N_DOF
    /// — no hard-coded 3s anywhere in the pipeline.
    type Dof: DofType;

    /// Number of nodes per element (4 for Tet4, 10 for Tet10, 8 for Hex8).
    const N_NODES: usize;

    /// Number of DOFs per element: N_NODES * Dof::N_DOF.
    /// Replaces the former hard-coded N_NODES * 3.
    const N_DOFS: usize = Self::N_NODES * <Self::Dof as DofType>::N_DOF;

    /// Precompute and cache all geometric data for one element.
    /// Called once in Assembler::new — result stored for the
    /// lifetime of the assembler.
    ///
    /// # Arguments
    /// * `nodes` - reference positions of the element nodes, length N_NODES
    fn precompute(nodes: &[Vector3<f64>]) -> Self::Geometry;

    /// Evaluate shape functions at a reference-space point xi.
    /// Returns a vector of length N_NODES.
    ///
    /// # Arguments
    /// * `xi` - coordinates in reference space (barycentric or isoparametric)
    fn shape_functions(xi: &Vector3<f64>) -> DVector<f64>;

    /// Evaluate shape function gradients with respect to physical
    /// coordinates at a given Gauss point, using precomputed geometry.
    /// Returns a matrix of shape (N_NODES x 3): each row is grad(N_i).
    ///
    /// # Arguments
    /// * `geometry`    - precomputed data for this element
    /// * `gauss_index` - index of the Gauss point being evaluated
    fn shape_gradients(geometry: &Self::Geometry, gauss_index: usize) -> DMatrix<f64>;

    /// Returns the list of (xi, weight) integration points in
    /// reference space. Used by the assembler to loop over Gauss points.
    fn integration_points() -> Vec<(Vector3<f64>, f64)>;

    /// Build the strain-displacement matrix B from shape function
    /// gradients evaluated at one Gauss point.
    ///
    /// Dimensions depend on the element DofType:
    /// - Vec3Dof (volumetric): 6 x (N_NODES * 3), standard Voigt B-matrix.
    /// - Vec6Dof (structural): element-specific dimensions encoding
    ///   membrane, bending, and shear contributions.
    ///
    /// Voigt ordering (for volumetric): [11, 22, 33, 12, 23, 13],
    /// consistent with hooke_voigt and nh_tangent_voigt throughout
    /// the codebase. Do not add extra shear factors in B.
    ///
    /// # Arguments
    /// * `grad_n` - shape function gradients (N_NODES x 3)
    /// * `f`      - deformation gradient at this Gauss point
    fn b_matrix(grad_n: &DMatrix<f64>, f: &Matrix3<f64>) -> DMatrix<f64>;

    /// Total volume of the element.
    /// For Tet4: geo.volume directly.
    /// For Tet10: sum of det_j * weight over all Gauss points.
    fn element_volume(geometry: &Self::Geometry) -> f64;

    /// Jacobian determinant at Gauss point gauss_index.
    /// Used to scale the integrand: integral ~ sum_g f(g) * det_J(g) * w(g).
    /// For Tet4: constant, equals geo.volume (single Gauss point, weight = 1).
    /// For Tet10: det_j stored per Gauss point in Tet10Geometry.
    fn gauss_det_j(geometry: &Self::Geometry, gauss_index: usize) -> f64;

    /// Returns the element stiffness matrix directly, bypassing the Gauss-point
    /// loop in the assembler. Used by analytical elements (Beam2, Shell) whose
    /// stiffness is not derived from a B-matrix formulation.
    ///
    /// Returns `None` for all volumetric elements (Tet4, Tet10, Hex8) — they
    /// use the standard B^T C B integration path. The default implementation
    /// returns `None` so existing elements require no change.
    ///
    /// Returns `Some(k_e)` for structural elements that override this method.
    /// The returned matrix must be in the **global** coordinate system and sized
    /// (N_NODES * N_DOF) x (N_NODES * N_DOF) — e.g. 12x12 for Beam2.
    /// It is scattered directly into K by the assembler, skipping Gauss points.
    ///
    /// # Arguments
    /// * `nodes`    - reference positions of the element nodes in global coords
    /// * `material` - material law providing E, nu, density
    fn element_stiffness(
        nodes:    &[Vector3<f64>],
        material: &dyn crate::material::MaterialLaw,
    ) -> Option<DMatrix<f64>> {
        let _ = (nodes, material);
        None
    }

    /// Compute the deformation gradient F at a Gauss point.
    /// F = I + sum_i u_i (x) grad(N_i)
    /// Default implementation — valid for all volumetric element types.
    /// Structural elements (Beam2, Shell) override this if needed.
    ///
    /// # Arguments
    /// * `grad_n`        - shape function gradients (N_NODES x 3)
    /// * `displacements` - displacement vectors for each node, length N_NODES
    fn deformation_gradient(
        grad_n: &DMatrix<f64>,
        displacements: &[Vector3<f64>],
    ) -> Matrix3<f64> {
        let mut f = Matrix3::identity();
        for (i, u) in displacements.iter().enumerate() {
            let grad_ni = grad_n.row(i);
            for row in 0..3 {
                for col in 0..3 {
                    f[(row, col)] += u[row] * grad_ni[col];
                }
            }
        }
        f
    }
}