// UTF-8
// material/anisotropic/mod.rs — anisotropic material models.

pub mod holzapfel_ogden;
pub mod no_anisotropy;

pub use holzapfel_ogden::HolzapfelOgden;
pub use no_anisotropy::NoAnisotropy;