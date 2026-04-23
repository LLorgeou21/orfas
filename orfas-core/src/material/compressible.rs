// UTF-8
// material/compressible.rs — generic compressible hyperelastic material.

use nalgebra::{Matrix3, Matrix6};
use super::traits::{MaterialLaw, IsochoricPart, VolumetricPart};
use super::helpers::cinv_tangent_voigt;

// ─── CompressibleMaterial ─────────────────────────────────────────────────────

/// Generic compressible hyperelastic material with isochoric/volumetric split.
///
/// Strain energy:
///   W(F) = W_iso(F) + W_vol(J)
///
/// where:
///   J     = det(F)          (volume ratio)
///
/// PK2 stress (Holzapfel 2000, eq. 6.88):
///   S = S_iso + S_vol
///
///   S_vol = dW_vol/dJ * J * C^{-1}
///
/// Material tangent (Holzapfel 2000, eq. 6.168):
///   C_tangent = C_iso + C_vol
///
///   C_vol = J*(J*d2W + dW) * (C^{-1} tensor C^{-1})
///         - 2*J*dW         * (C^{-1} odot C^{-1})
///
/// where dW = dW_vol/dJ, d2W = d^2W_vol/dJ^2.
///
/// Type parameters:
///   I — isochoric part, implements `IsochoricPart`
///   V — volumetric part, implements `VolumetricPart`
///
/// Breaking change from v0.6:
///   Old: `NeoHookean { youngs_modulus, poisson_ratio, density }`
///   New: `CompressibleMaterial { iso: NeoHookeanIso { mu }, vol: VolumetricLnJ { kappa }, density }`
///   See CHANGELOG for migration guide.
pub struct CompressibleMaterial<I: IsochoricPart, V: VolumetricPart> {
    /// Isochoric (deviatoric) part of the strain energy.
    pub iso: I,
    /// Volumetric part of the strain energy.
    pub vol: V,
    /// Mass density (kg/m^3) in the reference configuration.
    pub density: f64,
}

impl<I: IsochoricPart, V: VolumetricPart> MaterialLaw for CompressibleMaterial<I, V> {

    fn density(&self) -> f64 {
        self.density
    }

    /// W = W_iso(F) + W_vol(J)
    fn strain_energy(&self, f: &Matrix3<f64>) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        self.iso.strain_energy_iso(f) + self.vol.strain_energy_vol(j)
    }

    /// S = S_iso + S_vol
    ///
    /// S_vol = dW_vol/dJ * dJ/dC * 2 = dW_vol/dJ * J * C^{-1}
    ///
    /// Note: the factor J (not J/2) comes from the chain rule
    ///   dJ/dC = J/2 * C^{-1}  and  S = 2 * dW/dC,
    /// which gives S_vol = 2 * dW/dJ * J/2 * C^{-1} = dW/dJ * J * C^{-1}.
    fn pk2_stress(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };

        // S_vol = dW/dJ * J * C^{-1}
        let s_vol = self.vol.dw_vol(j) * j * c_inv;
        let s_iso = self.iso.pk2_stress_iso(f);

        s_iso + s_vol
    }

    /// C_tangent = C_iso + C_vol
    ///
    /// C_vol = J*(J*d2W + dW) * (C^{-1} tensor C^{-1})
    ///       - 2*J*dW         * (C^{-1} odot C^{-1})
    fn tangent_stiffness(&self, f: &Matrix3<f64>) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        let dw  = self.vol.dw_vol(j);
        let d2w = self.vol.d2w_vol(j);

        // Volumetric tangent coefficients (Holzapfel eq. 6.168)
        let a_vol =  j * (j * d2w + dw); // tensor coefficient
        let b_vol = -j * dw;       // odot coefficient

        let c_vol = cinv_tangent_voigt(&c_inv, a_vol, b_vol);
        let c_iso = self.iso.tangent_iso(f);

        c_iso + c_vol
    }
}