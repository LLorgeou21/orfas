// UTF-8
// orfas-core/src/element/mod.rs

pub mod traits;
pub mod tet4;
pub mod tet10;
pub mod subdivision;
pub mod hex8;
pub mod beam2;
pub mod shell4;

pub use traits::{ElementGeometry, FemMesh, FiniteElement};
pub use tet4::{Tet4, Tet4Geometry};
pub use tet10::{Tet10, Tet10Geometry};
pub use subdivision::tet4_to_tet10;
pub use hex8::Hex8;
pub use beam2::Beam2;
pub use shell4::Shell4;

#[cfg(test)]
mod tests;