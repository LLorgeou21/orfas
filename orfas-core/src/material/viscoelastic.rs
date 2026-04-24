// UTF-8
// material/viscoelastic.rs — viscoelastic material wrapper using Prony series.
//
// Reference: Holzapfel & Gasser (2001), Box 1.
//
// Wraps any CompressibleAnisotropicMaterial<I, A, V> and adds viscoelastic
// dissipation via a Prony series on the isochoric and anisotropic contributions.
//
// The volumetric part is purely elastic — no dissipation on volume changes.
//
// Algorithmic update at each time step (Box 1):
//   delta_alpha_a = beta_alpha_a * exp(-dt / 2*tau_alpha_a)
//   H_alpha_n     = exp(-dt/2*tau) * [exp(-dt/2*tau)*Q_n - beta*S_prev]
//   Q_alpha_n+1   = H_alpha_n + delta_alpha_a * S_n+1
//   S_n+1         = S_iso + S_aniso + S_vol + sum_alpha Q_alpha_n+1
//   C_n+1         = C_vol + (1 + delta_a_iso)*C_iso + (1 + delta_a_aniso)*C_aniso

use nalgebra::{Matrix3, Matrix6};
use crate::material::traits::{MaterialLaw, IsochoricPart, VolumetricPart, AnisotropicPart};
use crate::material::context::MaterialContext;
use crate::material::helpers::{cinv_tangent_voigt, mat3_to_voigt, voigt_to_mat3};

// ─── ViscoelasticMaterial ─────────────────────────────────────────────────────

/// Viscoelastic material with Prony series on isochoric and anisotropic contributions.
///
/// Strain energy:
///   W = W_iso(F) + W_aniso(F, fiber_dirs) + W_vol(J)
///
/// Total PK2 stress (algorithmic, Box 1):
///   S_n+1 = S_iso_inf + S_aniso_inf + S_vol_inf + sum_alpha Q_alpha_n+1
///
/// Algorithmic tangent (Box 1):
///   C_n+1 = C_vol_inf + (1 + delta_iso)*C_iso_inf + (1 + delta_aniso)*C_aniso_inf
///
/// Internal variables stored in `ctx.iv` (ElementInternalVars):
///   Q_iso[alpha]   — isochoric Prony tensors (Voigt, 6 components each)
///   Q_aniso[alpha] — anisotropic Prony tensors (Voigt, 6 components each)
///   S_iso_prev     — previous step isochoric PK2 (Voigt)
///   S_aniso_prev   — previous step anisotropic PK2 (Voigt)
///
/// Type parameters:
///   I — isochoric part,    implements `IsochoricPart`
///   A — anisotropic part,  implements `AnisotropicPart` (use `NoAnisotropy` if none)
///   V — volumetric part,   implements `VolumetricPart`
pub struct ViscoelasticMaterial<
    I: IsochoricPart,
    A: AnisotropicPart,
    V: VolumetricPart,
> {
    /// Isochoric ground matrix (e.g. NeoHookeanIso).
    pub iso:        I,
    /// Anisotropic fiber part (e.g. HolzapfelOgden or NoAnisotropy).
    pub aniso:      A,
    /// Volumetric penalty (e.g. VolumetricLnJ).
    pub vol:        V,
    /// Mass density (kg/m^3) in the reference configuration.
    pub density:    f64,
    /// Relaxation times for isochoric Prony processes (seconds).
    pub tau_iso:    Vec<f64>,
    /// Free-energy factors for isochoric Prony processes (dimensionless).
    pub beta_iso:   Vec<f64>,
    /// Relaxation times for anisotropic Prony processes (seconds).
    pub tau_aniso:  Vec<f64>,
    /// Free-energy factors for anisotropic Prony processes (dimensionless).
    pub beta_aniso: Vec<f64>,
}

impl<I: IsochoricPart, A: AnisotropicPart, V: VolumetricPart>
    ViscoelasticMaterial<I, A, V>
{
    /// Number of isochoric Prony processes.
    pub fn m_iso(&self) -> usize { self.tau_iso.len() }

    /// Number of anisotropic Prony processes.
    pub fn m_aniso(&self) -> usize { self.tau_aniso.len() }

    /// Computes delta_alpha_a = beta_alpha_a * exp(-dt / 2*tau_alpha_a).
    /// Used both in the stress update and the tangent scaling.
    fn delta(beta: f64, tau: f64, dt: f64) -> f64 {
        beta * (-dt / (2.0 * tau)).exp()
    }

    /// Computes the total delta_a = sum_alpha delta_alpha_a for tangent scaling.
    fn delta_total(betas: &[f64], taus: &[f64], dt: f64) -> f64 {
        betas.iter().zip(taus.iter())
            .map(|(&b, &t)| Self::delta(b, t, dt))
            .sum()
    }

    /// Updates Q_alpha for one contribution (iso or aniso) and returns the sum
    /// of all Q_alpha_n+1 as a Matrix3.
    ///
    /// Algorithm (Box 1):
    ///   H_alpha_n   = exp(-dt/tau)*Q_alpha_n - delta_alpha*S_prev_n
    ///   Q_alpha_n+1 = H_alpha_n + delta_alpha * S_n+1
    ///   sum Q       = sum_alpha Q_alpha_n+1
    ///
    /// `q_get` and `q_set` are closures to read/write Q_alpha from ElementInternalVars.
    fn update_q(
        s_cur:  &Matrix3<f64>,   // S_inf at current step n+1
        s_prev: &[f64],          // S_inf at previous step n (Voigt)
        taus:   &[f64],
        betas:  &[f64],
        dt:     f64,
        q_slices: &mut Vec<[f64; 6]>, // Q_alpha stored as Voigt
    ) -> Matrix3<f64> {
        let s_cur_v  = mat3_to_voigt(s_cur);
        let mut sum_q = Matrix3::zeros();

        for alpha in 0..taus.len() {
            let tau   = taus[alpha];
            let beta  = betas[alpha];
            let decay = (-dt / tau).exp();
            let delta = Self::delta(beta, tau, dt);

            // H_alpha_n = exp(-dt/tau)*Q_alpha_n - delta*S_prev_n
            // Q_alpha_n+1 = H_alpha_n + delta*S_cur_n+1
            let q_new: [f64; 6] = std::array::from_fn(|i| {
                decay * q_slices[alpha][i]
                    - delta * s_prev[i]
                    + delta * s_cur_v[i]
            });

            q_slices[alpha] = q_new;
            sum_q += voigt_to_mat3(&q_new);
        }
        sum_q
    }
}

impl<I: IsochoricPart, A: AnisotropicPart, V: VolumetricPart> MaterialLaw
    for ViscoelasticMaterial<I, A, V>
{
    fn density(&self) -> f64 { self.density }

    /// W = W_iso + W_aniso + W_vol (elastic equilibrium energy — no dissipation term).
    fn strain_energy(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        self.iso.strain_energy_iso(f)
            + self.aniso.strain_energy_aniso(f, ctx)
            + self.vol.strain_energy_vol(j)
    }

    /// S_n+1 = S_iso_inf + S_aniso_inf + S_vol_inf + sum_Q_n
    ///
    /// Reads sum_Q from ctx.iv — O(1), does NOT update iv.
    /// If ctx.iv is None, returns elastic stress only (fallback).
    fn pk2_stress(&self, f: &Matrix3<f64>, ctx: &mut MaterialContext) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };

        let s_iso   = self.iso.pk2_stress_iso(f);
        let s_aniso = self.aniso.pk2_stress_aniso(f, ctx);
        let s_vol   = self.vol.dw_vol(j) * j * c_inv;
        let s_eq    = s_iso + s_aniso + s_vol;

        // Read sum_Q from iv — O(1)
        let sum_q = match ctx.iv_ref.as_ref() {
            Some(iv) => voigt_to_mat3(iv.sum_q()),
            None     => Matrix3::zeros(),
        };

        s_eq + sum_q
    }

    /// C_n+1 = C_vol_inf + (1 + delta_iso)*C_iso_inf + (1 + delta_aniso)*C_aniso_inf
    ///
    /// Read-only — does not update internal variables.
    /// If ctx.iv is None, returns the elastic tangent.
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
        let c_aniso = self.aniso.tangent_aniso(f, ctx);

        // Elastic fallback — no internal variables
        if ctx.dt == 0.0 {
            return c_vol + c_iso + c_aniso;
        }

        let dt = ctx.dt;

        // Algorithmic scaling (Box 1): C_n+1 = C_vol + (1+delta_iso)*C_iso + (1+delta_aniso)*C_aniso
        let delta_iso   = Self::delta_total(&self.beta_iso,   &self.tau_iso,   dt);
        let delta_aniso = Self::delta_total(&self.beta_aniso, &self.tau_aniso, dt);

        c_vol + (1.0 + delta_iso) * c_iso + (1.0 + delta_aniso) * c_aniso
    }


    fn update_state(&self, f: &Matrix3<f64>, ctx: &mut MaterialContext) {
        if ctx.dt == 0.0 { return; }
        let dt = ctx.dt;
        let j = f.determinant();
        if j <= 0.0 { return; }

        // Calculate stresses before borrowing iv
        let s_iso   = self.iso.pk2_stress_iso(f);
        let s_aniso = self.aniso.pk2_stress_aniso(f, ctx);

        let iv = match ctx.iv.as_mut() {
            Some(iv) => iv,
            None     => return,
        };

        let m_iso   = self.m_iso();
        let m_aniso = self.m_aniso();

        let mut q_iso_slices:   Vec<[f64; 6]> = (0..m_iso).map(|a| iv.q_iso(a).try_into().unwrap()).collect();
        let mut q_aniso_slices: Vec<[f64; 6]> = (0..m_aniso).map(|a| iv.q_aniso(a).try_into().unwrap()).collect();
        let s_iso_prev:   [f64; 6] = iv.s_iso_prev().try_into().unwrap();
        let s_aniso_prev: [f64; 6] = iv.s_aniso_prev().try_into().unwrap();

        // Reuse update_q
        let sum_q_iso   = Self::update_q(&s_iso,   &s_iso_prev,   &self.tau_iso,   &self.beta_iso,   dt, &mut q_iso_slices);
        let sum_q_aniso = Self::update_q(&s_aniso, &s_aniso_prev, &self.tau_aniso, &self.beta_aniso, dt, &mut q_aniso_slices);

        // Write back
        for alpha in 0..m_iso   { iv.q_iso_mut(alpha).copy_from_slice(&q_iso_slices[alpha]); }
        for alpha in 0..m_aniso { iv.q_aniso_mut(alpha).copy_from_slice(&q_aniso_slices[alpha]); }
        iv.set_s_iso_prev(&mat3_to_voigt(&s_iso));
        iv.set_s_aniso_prev(&mat3_to_voigt(&s_aniso));

        // sum_Q = sum_Q_iso + sum_Q_aniso — precomputed for pk2_stress
        let sum_q_total = mat3_to_voigt(&(sum_q_iso + sum_q_aniso));
        iv.set_sum_q(&sum_q_total);
    }
}