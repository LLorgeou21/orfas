// UTF-8
// material/isochoric/mod.rs — re-exports all isochoric material laws.
// To add a new isochoric law: create a new file here and add a line below.

mod neo_hookean;
mod mooney_rivlin;
mod ogden;

pub use neo_hookean::NeoHookeanIso;
pub use mooney_rivlin::MooneyRivlinIso;
pub use ogden::OgdenIso;