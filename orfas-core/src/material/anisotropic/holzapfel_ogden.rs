// UTF-8
// material/anisotropic/holzapfel_ogden.rs — Holzapfel-Gasser-Ogden anisotropic fiber model.
//
// Reference: Cheng & Zhang (2018), eq. (37c), (49), (56).
//
// This struct implements only the anisotropic isochoric fiber contribution.
// The isotropic ground matrix (NeoHookeanIso) and volumetric penalty (VolumetricLnJ)
// are handled separately by CompressibleAnisotropicMaterial<I, A, V>.
//
// Fiber directions are passed at call time via fiber_dirs: &[Vector3<f64>].
// Each entry is a unit vector a0i in the reference configuration.
// Fibers only contribute under tension: if I_bar_i <= 1.0, contribution is zero.

use nalgebra::{Matrix3, Matrix6, Vector3};
use crate::material::traits::AnisotropicPart;
use crate::material::MaterialContext;
use crate::material::helpers::{
    a_tensor_a_voigt,
    a_tensor_b_plus_b_tensor_a,
    cinv_tangent_voigt,
};

// ─── HolzapfelOgden ───────────────────────────────────────────────────────────

/// Anisotropic isochoric fiber contribution from the Holzapfel-Gasser-Ogden model.
///
/// Strain energy (Cheng & Zhang eq. 37c):
///   W_aniso = k1/(2*k2) * sum_i [ exp(k2*(I_bar_i - 1)^2) - 1 ]
///
/// where:
///   I_bar_i = J^{-2/3} * (a0i . C . a0i)   (modified pseudo-invariant, eq. 43)
///   a0i     = fiber direction unit vector in reference configuration
///
/// Parameters:
///   k1 > 0  stress-like parameter (Pa)
///   k2 > 0  dimensionless exponential parameter
#[derive(serde::Serialize, serde::Deserialize)]
pub struct HolzapfelOgden {
    /// Stress-like fiber parameter (Pa). Must be > 0.
    pub k1: f64,
    /// Dimensionless exponential fiber parameter. Must be > 0.
    pub k2: f64,
}

impl HolzapfelOgden {
    /// Compute the modified pseudo-invariant I_bar_i = J^{-2/3} * (a0i . C . a0i).
    /// Returns (I_i, I_bar_i) — both are needed in the tangent computation.
    fn invariants(j: f64, c: &Matrix3<f64>, a0i: &Vector3<f64>) -> (f64, f64) {
        // I_i = a0i . C . a0i  (unmodified, eq. 41-42)
        let i_i = a0i.dot(&(c * a0i));
        // I_bar_i = J^{-2/3} * I_i  (modified, eq. 43)
        let i_bar_i = j.powf(-2.0 / 3.0) * i_i;
        (i_i, i_bar_i)
    }

    /// First derivative dPsi/dI_bar_i (inline below eq. 49, Cheng & Zhang).
    /// Returns 0 if I_bar_i <= 1 (fibers only resist tension).
    pub fn dw(k1: f64, k2: f64, i_bar: f64) -> f64 {
        if i_bar <= 1.0 {
            return 0.0;
        }
        k1 * (i_bar - 1.0) * ((k2 * (i_bar - 1.0).powi(2)).exp())
    }

    /// Second derivative d2Psi/dI_bar_i^2 (derived from dw).
    /// Returns 0 if I_bar_i <= 1 (fibers only resist tension).
    pub fn d2w(k1: f64, k2: f64, i_bar: f64) -> f64 {
        if i_bar <= 1.0 {
            return 0.0;
        }
        let exp = (k2 * (i_bar - 1.0).powi(2)).exp();
        k1 * exp * (1.0 + 2.0 * k2 * (i_bar - 1.0).powi(2))
    }
}

impl AnisotropicPart for HolzapfelOgden {

    /// Anisotropic isochoric strain energy density (Cheng & Zhang eq. 37c):
    ///   W_aniso = k1/(2*k2) * sum_i [ exp(k2*(I_bar_i - 1)^2) - 1 ]
    /// Fiber contribution is zero if I_bar_i <= 1 (compression).
    fn strain_energy_aniso(&self, f: &Matrix3<f64>, ctx :&MaterialContext) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }

        let c = f.transpose() * f;
        let mut w = 0.0;

        for a0i in ctx.fiber_dirs {
            let (_, i_bar) = Self::invariants(j, &c, a0i);
            if i_bar <= 1.0 { continue; }
            w += (self.k1 / (2.0 * self.k2))
                * ((self.k2 * (i_bar - 1.0).powi(2)).exp() - 1.0);
        }
        w
    }

    /// Anisotropic isochoric PK2 stress (Cheng & Zhang eq. 49):
    ///   S_aniso = 2*J^{-2/3} * sum_i dW * (A0i - 1/3 * I_i * C^{-1})
    ///
    /// Note: the 1/3 term uses the unmodified invariant I_i, not I_bar_i.
    fn pk2_stress_aniso(&self, f: &Matrix3<f64>,ctx: &mut MaterialContext) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };

        let j_23 = j.powf(-2.0 / 3.0);
        let mut s = Matrix3::zeros();

        for a0i in ctx.fiber_dirs {
            let (i_i, i_bar) = Self::invariants(j, &c, a0i);
            let dw = Self::dw(self.k1, self.k2, i_bar);
            if dw == 0.0 { continue; }

            // Structural tensor A0i = a0i outer a0i
            let a0i_mat = a0i * a0i.transpose();

            // S_aniso += 2*J^{-2/3} * dW * (A0i - 1/3 * I_i * C^{-1})
            s += 2.0 * j_23 * dw * (a0i_mat - (1.0 / 3.0) * i_i * c_inv);
        }
        s
    }

    /// Anisotropic isochoric tangent stiffness (Cheng & Zhang eq. 56).
    ///
    /// For each fiber family i, four terms are summed and multiplied by J^{-4/3}:
    ///
    /// Term 1: 4 * d2W * (A0i tensor A0i)
    /// Term 2: -4/3 * (I_bar_i*d2W + dW) * (C_bar^{-1} tensor A0i + A0i tensor C_bar^{-1})
    /// Term 3: +4/9 * (I_bar_i^2*d2W + I_bar_i*dW) * (C_bar^{-1} tensor C_bar^{-1})
    /// Term 4: +4/3 * I_bar_i * dW * (C_bar^{-1} odot C_bar^{-1})
    ///
    /// where C_bar^{-1} = J^{2/3} * C^{-1}
    fn tangent_aniso(&self, f: &Matrix3<f64>, ctx: &MaterialContext) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        // C_bar^{-1} = J^{2/3} * C^{-1}
        let c_bar_inv = j.powf(2.0 / 3.0) * c_inv;
        let j_43      = j.powf(-4.0 / 3.0);

        let mut c_aniso = Matrix6::zeros();

        for a0i in ctx.fiber_dirs {
            let (_, i_bar) = Self::invariants(j, &c, a0i);
            let dw  = Self::dw(self.k1, self.k2, i_bar);
            let d2w = Self::d2w(self.k1, self.k2, i_bar);
            if dw == 0.0 && d2w == 0.0 { continue; }

            // Structural tensor A0i = a0i outer a0i
            let a0i_mat = a0i * a0i.transpose();

            // Term 1: 4 * d2W * (A0i tensor A0i)
            let term1 = 4.0 * d2w * a_tensor_a_voigt(&a0i_mat);

            // Term 2: -4/3 * (I_bar*d2W + dW) * (C_bar^{-1} tensor A0i + A0i tensor C_bar^{-1})
            let coeff2 = -4.0 / 3.0 * (i_bar * d2w + dw);
            let term2  = coeff2 * a_tensor_b_plus_b_tensor_a(&c_bar_inv, &a0i_mat);

            // Term 3: +4/9 * (I_bar^2*d2W + I_bar*dW) * (C_bar^{-1} tensor C_bar^{-1})
            let coeff3 = 4.0 / 9.0 * (i_bar.powi(2) * d2w + i_bar * dw);
            let term3  = cinv_tangent_voigt(&c_bar_inv, coeff3, 0.0);

            // Term 4: +4/3 * I_bar * dW * (C_bar^{-1} odot C_bar^{-1})
            // Coefficient halved to 2/3 because cinv_tangent_voigt b-term computes
            // 2*(A odot A) — see helpers.rs convention (b coeff = 2*odot coeff).
            let coeff4 = 2.0 / 3.0 * i_bar * dw;
            let term4  = cinv_tangent_voigt(&c_bar_inv, 0.0, coeff4);

            // Sum all terms and multiply by J^{-4/3}
            c_aniso += j_43 * (term1 + term2 + term3 + term4);
        }
        c_aniso
    }
}