// UTF-8
// material/context.rs — evaluation context passed to material laws.

use nalgebra::Vector3;
use super::internal_variables::{ElementInternalVars, InternalVariables};
use super::fiber_fields::FiberField;

// ─── MaterialContext ──────────────────────────────────────────────────────────

/// Per-element evaluation context passed to material laws at each quadrature point.
///
/// Two variants for internal variables:
/// - `iv_ref` — read-only borrow, used by `pk2_stress` and `assemble_internal_forces`
/// - `iv`     — mutable borrow, used by `update_state` and `update_internal_variables`
pub struct MaterialContext<'a> {
    /// Time step size in seconds. Zero for static and quasi-static problems.
    pub dt: f64,
    /// Fiber directions for this element in the reference configuration.
    pub fiber_dirs: &'a [Vector3<f64>],
    /// Read-only internal variables — used by `pk2_stress` to read sum_Q.
    /// None for purely elastic materials.
    pub iv_ref: Option<&'a ElementInternalVars>,
    /// Mutable internal variables — used by `update_state` to write Q_alpha.
    /// None for purely elastic materials and during Newton iterations.
    pub iv: Option<&'a mut ElementInternalVars>,
}

impl<'a> MaterialContext<'a> {
    /// Full constructor — for viscoelastic anisotropic materials.
    pub fn new(
        fiber_field: &'a FiberField,
        elem_idx:    usize,
        dt:          f64,
        iv:          Option<&'a mut ElementInternalVars>,
    ) -> Self {
        MaterialContext {
            dt,
            fiber_dirs: fiber_field.directions_for(elem_idx),
            iv_ref:     None,
            iv,
        }
    }

    /// Constructor for anisotropic elastic materials — no internal variables, no dt.
    pub fn from_fiber_field(fiber_field: &'a FiberField, elem_idx: usize) -> Self {
        MaterialContext {
            dt:         0.0,
            fiber_dirs: fiber_field.directions_for(elem_idx),
            iv_ref:     None,
            iv:         None,
        }
    }

    /// Constructor for direct fiber directions — useful in tests.
    pub fn from_fiber_dirs(fiber_dirs: &'a [Vector3<f64>]) -> Self {
        MaterialContext {
            dt:         0.0,
            fiber_dirs,
            iv_ref:     None,
            iv:         None,
        }
    }

    /// Constructor for isotropic elastic materials with a time step — no fibers, no iv.
    pub fn from_dt(dt: f64) -> Self {
        MaterialContext {
            dt,
            fiber_dirs: &[],
            iv_ref:     None,
            iv:         None,
        }
    }
}

impl Default for MaterialContext<'_> {
    fn default() -> Self {
        MaterialContext {
            dt:         0.0,
            fiber_dirs: &[],
            iv_ref:     None,
            iv:         None,
        }
    }
}

// ─── SimulationContext ────────────────────────────────────────────────────────

pub struct SimulationContext {
    pub fiber_field: FiberField,
    pub dt:          f64,
    pub iv:          Option<InternalVariables>,
}

impl SimulationContext {
    pub fn new(fiber_field: FiberField, dt: f64) -> Self {
        SimulationContext { fiber_field, dt, iv: None }
    }

    pub fn static_elastic(fiber_field: FiberField) -> Self {
        SimulationContext { fiber_field, dt: 0.0, iv: None }
    }

    pub fn isotropic_static(n_elements: usize) -> Self {
        SimulationContext {
            fiber_field: FiberField::empty(n_elements),
            dt:          0.0,
            iv:          None,
        }
    }

    pub fn isotropic_dynamic(n_elements: usize, dt: f64) -> Self {
        SimulationContext {
            fiber_field: FiberField::empty(n_elements),
            dt,
            iv: None,
        }
    }

    pub fn viscoelastic(fiber_field: FiberField, dt: f64, iv: InternalVariables) -> Self {
        SimulationContext { fiber_field, dt, iv: Some(iv) }
    }

    pub fn isotropic_viscoelastic(n_elements: usize, dt: f64, iv: InternalVariables) -> Self {
        SimulationContext {
            fiber_field: FiberField::empty(n_elements),
            dt,
            iv: Some(iv),
        }
    }

    /// Read-only context — passes iv_ref for pk2_stress.
    /// Used by assemble_internal_forces and assemble_tangent.
    pub fn material_context_for(&self, elem_idx: usize) -> MaterialContext {
        MaterialContext {
            dt:         self.dt,
            fiber_dirs: self.fiber_field.directions_for(elem_idx),
            iv_ref:     self.iv.as_ref().map(|iv| iv.get(elem_idx)),
            iv:         None,
        }
    }

    /// Mutable context — passes iv for update_state.
    /// Used by update_internal_variables only.
    pub fn material_context_for_mut(&mut self, elem_idx: usize) -> MaterialContext {
        MaterialContext {
            dt:         self.dt,
            fiber_dirs: self.fiber_field.directions_for(elem_idx),
            iv_ref:     None,
            iv:         self.iv.as_mut().map(|iv| iv.get_mut(elem_idx)),
        }
    }
}