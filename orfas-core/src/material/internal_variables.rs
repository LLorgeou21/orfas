// UTF-8
// material/internal_variables.rs — per-element internal variables for history-dependent materials.
//
// Layout of ElementInternalVars.data for ViscoelasticMaterial<I, A, V>:
//
//   [0..6*m_iso]                              Q_iso[0..m_iso]     (Prony tensors, isochoric)
//   [6*m_iso..6*(m_iso+m_aniso)]              Q_aniso[0..m_aniso] (Prony tensors, anisotropic)
//   [6*(m_iso+m_aniso)..+6]                   S_iso_prev          (previous isochoric PK2)
//   [6*(m_iso+m_aniso)+6..+6]                 S_aniso_prev        (previous anisotropic PK2)
//   [6*(m_iso+m_aniso)+12..+6]                sum_Q               (precomputed sum of Q_alpha)
//
// Total size per element: 6*(m_iso + m_aniso + 3)
//
// Access pattern:
//   pk2_stress    — reads  sum_Q only (O(1), called every Newton iteration)
//   update_state  — reads/writes Q_iso, Q_aniso, S_iso_prev, S_aniso_prev, sum_Q
//                   (called once per time step after Newton convergence)

use nalgebra::DVector;

// ─── ElementInternalVars ─────────────────────────────────────────────────────

/// Per-element internal variable storage for viscoelastic materials.
///
/// Stores Prony series tensors Q_alpha (isochoric and anisotropic contributions),
/// previous-step elastic stresses, and the precomputed sum of all Q_alpha tensors.
/// All tensors are stored in Voigt notation [11,22,33,12,23,13] — 6 components each.
///
/// `sum_Q` is the only field read by `pk2_stress` at every Newton iteration.
/// All other fields are written by `update_state` once per time step.
pub struct ElementInternalVars {
    /// Flat storage — layout described in module header.
    pub data:    DVector<f64>,
    /// Number of isochoric Prony processes.
    pub m_iso:   usize,
    /// Number of anisotropic Prony processes.
    pub m_aniso: usize,
}

impl ElementInternalVars {
    /// Creates zero-initialized internal variables for given process counts.
    /// Total data size: 6*(m_iso + m_aniso + 3).
    pub fn zeros(m_iso: usize, m_aniso: usize) -> Self {
        ElementInternalVars {
            data:    DVector::zeros(6 * (m_iso + m_aniso + 3)),
            m_iso,
            m_aniso,
        }
    }

    // ── Offset helpers ────────────────────────────────────────────────────────

    fn offset_q_iso(&self, alpha: usize) -> usize {
        6 * alpha
    }

    fn offset_q_aniso(&self, alpha: usize) -> usize {
        6 * self.m_iso + 6 * alpha
    }

    fn offset_s_iso_prev(&self) -> usize {
        6 * (self.m_iso + self.m_aniso)
    }

    fn offset_s_aniso_prev(&self) -> usize {
        6 * (self.m_iso + self.m_aniso) + 6
    }

    fn offset_sum_q(&self) -> usize {
        6 * (self.m_iso + self.m_aniso) + 12
    }

    // ── Q_iso accessors ───────────────────────────────────────────────────────

    /// Returns Q_iso[alpha] as a slice of 6 Voigt components.
    pub fn q_iso(&self, alpha: usize) -> &[f64] {
        let start = self.offset_q_iso(alpha);
        &self.data.as_slice()[start..start + 6]
    }

    /// Returns a mutable slice for Q_iso[alpha].
    pub fn q_iso_mut(&mut self, alpha: usize) -> &mut [f64] {
        let start = self.offset_q_iso(alpha);
        &mut self.data.as_mut_slice()[start..start + 6]
    }

    // ── Q_aniso accessors ─────────────────────────────────────────────────────

    /// Returns Q_aniso[alpha] as a slice of 6 Voigt components.
    pub fn q_aniso(&self, alpha: usize) -> &[f64] {
        let start = self.offset_q_aniso(alpha);
        &self.data.as_slice()[start..start + 6]
    }

    /// Returns a mutable slice for Q_aniso[alpha].
    pub fn q_aniso_mut(&mut self, alpha: usize) -> &mut [f64] {
        let start = self.offset_q_aniso(alpha);
        &mut self.data.as_mut_slice()[start..start + 6]
    }

    // ── S_prev accessors ──────────────────────────────────────────────────────

    /// Returns S_iso_prev as a slice of 6 Voigt components.
    pub fn s_iso_prev(&self) -> &[f64] {
        let start = self.offset_s_iso_prev();
        &self.data.as_slice()[start..start + 6]
    }

    /// Returns S_aniso_prev as a slice of 6 Voigt components.
    pub fn s_aniso_prev(&self) -> &[f64] {
        let start = self.offset_s_aniso_prev();
        &self.data.as_slice()[start..start + 6]
    }

    /// Writes S_iso_prev from a Voigt slice.
    pub fn set_s_iso_prev(&mut self, s: &[f64]) {
        let start = self.offset_s_iso_prev();
        self.data.as_mut_slice()[start..start + 6].copy_from_slice(s);
    }

    /// Writes S_aniso_prev from a Voigt slice.
    pub fn set_s_aniso_prev(&mut self, s: &[f64]) {
        let start = self.offset_s_aniso_prev();
        self.data.as_mut_slice()[start..start + 6].copy_from_slice(s);
    }

    // ── sum_Q accessors ───────────────────────────────────────────────────────

    /// Returns the precomputed sum of all Q_alpha tensors as a Voigt slice.
    /// Read by `pk2_stress` at every Newton iteration — O(1).
    pub fn sum_q(&self) -> &[f64] {
        let start = self.offset_sum_q();
        &self.data.as_slice()[start..start + 6]
    }

    /// Writes the precomputed sum_Q from a Voigt slice.
    /// Called by `update_state` once per time step after Newton convergence.
    pub fn set_sum_q(&mut self, s: &[f64]) {
        let start = self.offset_sum_q();
        self.data.as_mut_slice()[start..start + 6].copy_from_slice(s);
    }
}

// ─── InternalVariables ───────────────────────────────────────────────────────

/// Global storage of per-element internal variables for a simulation.
///
/// One `ElementInternalVars` per mesh element — all initialized with the same
/// m_iso and m_aniso counts matching the `ViscoelasticMaterial` parameters.
pub struct InternalVariables {
    data: Vec<ElementInternalVars>,
}

impl InternalVariables {
    /// Creates zero-initialized internal variables for all elements.
    pub fn new(n_elements: usize, m_iso: usize, m_aniso: usize) -> Self {
        InternalVariables {
            data: (0..n_elements)
                .map(|_| ElementInternalVars::zeros(m_iso, m_aniso))
                .collect(),
        }
    }

    /// Returns a reference to the internal variables of element `elem_idx`.
    pub fn get(&self, elem_idx: usize) -> &ElementInternalVars {
        &self.data[elem_idx]
    }

    /// Returns a mutable reference to the internal variables of element `elem_idx`.
    pub fn get_mut(&mut self, elem_idx: usize) -> &mut ElementInternalVars {
        &mut self.data[elem_idx]
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.data.len()
    }
}