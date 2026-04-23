// UTF-8
// material/traits.rs — core hyperelastic traits.

use nalgebra::{Matrix3, Matrix6};

// ─── MaterialLaw ──────────────────────────────────────────────────────────────

/// Defines a hyperelastic material law in the Lagrangian frame.
///
/// All methods take the deformation gradient F (3x3) as input.
/// Each implementation derives what it needs internally (E, C, J, ...).
///
/// Conventions:
/// - `pk2_stress` returns S (2nd Piola-Kirchhoff, symmetric, reference frame).
/// - The assembler computes P = F*S (1st Piola-Kirchhoff) internally.
/// - `tangent_stiffness` returns dS/dE in Voigt notation (6x6).
/// - Voigt ordering: [11, 22, 33, 12, 23, 13].
/// - Engineering shear convention: no extra Voigt factors in C — they live in B.
///
/// Implementors:
/// - `SaintVenantKirchhoff`        — direct impl, E-based, no iso/vol split
/// - `CompressibleMaterial<I, V>`  — isochoric/volumetric split, F-based
pub trait MaterialLaw: Sync {
    /// Mass density (kg/m^3) in the reference configuration.
    fn density(&self) -> f64;

    /// Strain energy density W(F).
    fn strain_energy(&self, f: &Matrix3<f64>) -> f64;

    /// 2nd Piola-Kirchhoff stress tensor S = dW/dE (3x3, symmetric).
    fn pk2_stress(&self, f: &Matrix3<f64>) -> Matrix3<f64>;

    /// Material tangent stiffness C = dS/dE in Voigt notation (6x6).
    /// Constant for SVK. Depends on F for all isochoric/volumetric materials.
    fn tangent_stiffness(&self, f: &Matrix3<f64>) -> Matrix6<f64>;
}

// ─── IsochoricPart ────────────────────────────────────────────────────────────

/// Isochoric (volume-preserving) part of a decoupled hyperelastic energy.
///
/// Works with the isochoric right Cauchy-Green tensor:
///   C_bar = J^{-2/3} * C    where C = F^T F,  J = det(F)
///
/// Isochoric invariants:
///   I1_bar = tr(C_bar)  = J^{-2/3} * tr(C)
///   I2_bar = 0.5*(tr(C_bar)^2 - tr(C_bar^2))
///
/// All methods take F as input — the implementor derives C, J, C_bar internally.
/// The full PK2 stress and tangent are assembled by `CompressibleMaterial`.
pub trait IsochoricPart: Sync {
    /// Isochoric strain energy density W_iso(F).
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64;

    /// Isochoric PK2 stress contribution S_iso (3x3, symmetric).
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64>;

    /// Isochoric tangent stiffness contribution in Voigt notation (6x6).
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64>;
}

// ─── VolumetricPart ───────────────────────────────────────────────────────────

/// Volumetric part of a decoupled hyperelastic energy W_vol(J).
///
/// Provides the three scalar derivatives needed to assemble the volumetric
/// PK2 stress and material tangent in `CompressibleMaterial`:
///
///   S_vol = dW_vol/dJ * J * C^{-1}
///
///   C_vol = J*(J*d2W + dW) * (C^{-1} tensor C^{-1})
///         - 2*J*dW         * (C^{-1} odot C^{-1})
///
/// where J = det(F) > 0.
pub trait VolumetricPart: Sync {
    /// Volumetric strain energy density W_vol(J).
    fn strain_energy_vol(&self, j: f64) -> f64;

    /// First derivative dW_vol/dJ.
    fn dw_vol(&self, j: f64) -> f64;

    /// Second derivative d^2W_vol/dJ^2.
    fn d2w_vol(&self, j: f64) -> f64;
}