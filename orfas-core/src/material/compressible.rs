// UTF-8
// material/compressible.rs — generic compressible hyperelastic material.

use nalgebra::{Matrix3, Matrix6};
use super::traits::{MaterialLaw, IsochoricPart, VolumetricPart, AnisotropicPart};
use crate::material::MaterialContext;
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
    fn strain_energy(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        self.iso.strain_energy_iso(f) + self.vol.strain_energy_vol(j)
    }

    /// S = S_iso + S_vol
    ///
    /// S_vol = dW_vol/dJ * dJ/dC * 2 = dW_vol/dJ * J * C^{-1}
    fn pk2_stress(&self, f: &Matrix3<f64>, _ctx: &mut MaterialContext) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };

        let s_vol = self.vol.dw_vol(j) * j * c_inv;
        let s_iso = self.iso.pk2_stress_iso(f);

        s_iso + s_vol
    }

    /// C_tangent = C_iso + C_vol
    fn tangent_stiffness(&self, f: &Matrix3<f64>, _ctx: &MaterialContext) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        let dw  = self.vol.dw_vol(j);
        let d2w = self.vol.d2w_vol(j);

        let a_vol =  j * (j * d2w + dw);
        let b_vol = -j * dw;

        let c_vol = cinv_tangent_voigt(&c_inv, a_vol, b_vol);
        let c_iso = self.iso.tangent_iso(f);

        c_iso + c_vol
    }
}

// ─── CompressibleAnisotropicMaterial ─────────────────────────────────────────

/// Generic compressible hyperelastic material with isochoric/anisotropic/volumetric split.
///
/// Strain energy:
///   W(F, fiber_dirs) = W_iso(F) + W_aniso(F, fiber_dirs) + W_vol(J)
///
/// PK2 stress:
///   S = S_iso + S_aniso + S_vol
///
/// Material tangent:
///   C_tangent = C_iso + C_aniso + C_vol
///
/// This is the natural extension of `CompressibleMaterial` for fiber-reinforced
/// biological tissues (arteries, myocardium, ligaments).
///
/// Type parameters:
///   I — isochoric ground matrix part, implements `IsochoricPart`
///   A — anisotropic fiber part,       implements `AnisotropicPart`
///   V — volumetric penalty part,      implements `VolumetricPart`
///
/// Usage example (HGO arterial model):
///   CompressibleAnisotropicMaterial {
///       iso:     NeoHookeanIso { mu },
///       aniso:   HolzapfelOgden { k1, k2 },
///       vol:     VolumetricLnJ { kappa },
///       density: rho,
///   }
pub struct CompressibleAnisotropicMaterial<I: IsochoricPart,A: AnisotropicPart,V: VolumetricPart,>{
    /// Isochoric ground matrix contribution (e.g. NeoHookeanIso).
    pub iso: I,
    /// Anisotropic fiber contribution (e.g. HolzapfelOgden).
    pub aniso: A,
    /// Volumetric penalty contribution (e.g. VolumetricLnJ).
    pub vol: V,
    /// Mass density (kg/m^3) in the reference configuration.
    pub density: f64,
}

impl<I: IsochoricPart, A: AnisotropicPart, V: VolumetricPart> MaterialLaw
    for CompressibleAnisotropicMaterial<I, A, V>
{
    fn density(&self) -> f64 {
        self.density
    }

    /// W = W_iso(F) + W_aniso(F, fiber_dirs) + W_vol(J)
    fn strain_energy(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        self.iso.strain_energy_iso(f)
            + self.aniso.strain_energy_aniso(f, &ctx)
            + self.vol.strain_energy_vol(j)
    }

    /// S = S_iso + S_aniso + S_vol
    ///
    /// S_vol = dW_vol/dJ * J * C^{-1}
    /// S_iso and S_aniso are computed by their respective parts.
    fn pk2_stress(&self, f: &Matrix3<f64>, ctx: &mut MaterialContext) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };

        let s_vol   = self.vol.dw_vol(j) * j * c_inv;
        let s_iso   = self.iso.pk2_stress_iso(f);
        let s_aniso = self.aniso.pk2_stress_aniso(f, ctx);

        s_iso + s_aniso + s_vol
    }

    /// C_tangent = C_iso + C_aniso + C_vol
    fn tangent_stiffness(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        let dw  = self.vol.dw_vol(j);
        let d2w = self.vol.d2w_vol(j);

        let a_vol =  j * (j * d2w + dw);
        let b_vol = -j * dw;

        let c_vol   = cinv_tangent_voigt(&c_inv, a_vol, b_vol);
        let c_iso   = self.iso.tangent_iso(f);
        let c_aniso = self.aniso.tangent_aniso(f, &ctx);

        c_iso + c_aniso + c_vol
    }
}