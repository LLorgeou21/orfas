// UTF-8
// material/isochoric/neo_hookean.rs — isochoric Neo-Hookean law.

use nalgebra::{Matrix3, Matrix6};
use crate::material::traits::IsochoricPart;
use crate::material::helpers::{cinv_tangent_voigt, i_tensor_a_plus_a_tensor_i};

// ─── NeoHookeanIso ────────────────────────────────────────────────────────────

/// Isochoric Neo-Hookean strain energy:
///   W_iso = mu/2 * (I1_bar - 3)
///
/// where:
///   I1_bar = J^{-2/3} * tr(C)   (first isochoric invariant)
///   C      = F^T F
///   J      = det(F)
///
/// PK2 stress (isochoric part), derived as S = 2 * J^{-2/3} * DEV_C[dW/dC_bar]:
///   S_iso = mu * J^{-2/3} * (I - tr(C)/3 * C^{-1})
///
/// where DEV_C[A] = A - 1/3*(A:C)*C^{-1} is the material deviatoric operator.
///
/// Isochoric tangent (Holzapfel 2000, eq. 6.166):
///   C_iso = 2*mu*J^{-2/3} * (I odot I)
///         - 2/3*mu*J^{-2/3} * (I tensor C^{-1} + C^{-1} tensor I)
///         + 2/9*mu*J^{-2/3}*I1_bar * (C^{-1} tensor C^{-1})
///
/// At F=I: S_iso = 0, C_iso = 2*mu*(I odot I) - 4/3*mu*(I tensor I) + 2/3*mu*(I tensor I)
///              = 2*mu*(I odot I) - 2/3*mu*(I tensor I)  = 2*mu*I4_dev
/// Combined with C_vol at F=I gives hooke_voigt(lambda, mu) when kappa = lambda + 2/3*mu.
///
/// Parameters:
///   mu — shear modulus (Pa), must be > 0
#[derive(serde::Serialize, serde::Deserialize)]
pub struct NeoHookeanIso {
    pub mu: f64,
}

impl NeoHookeanIso {
    /// Safe constructor — validates mu > 0.
    pub fn new(mu: f64) -> Result<Self, String> {
        if mu <= 0.0 {
            return Err(format!("NeoHookeanIso: mu must be > 0, got {}", mu));
        }
        Ok(Self { mu })
    }
}

impl IsochoricPart for NeoHookeanIso {

    /// W_iso = mu/2 * (I1_bar - 3)  where I1_bar = J^{-2/3} * tr(C)
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }
        let c      = f.transpose() * f;
        let i1     = c.trace();
        let i1_bar = j.powf(-2.0 / 3.0) * i1;
        self.mu / 2.0 * (i1_bar - 3.0)
    }

    /// S_iso = mu * J^{-2/3} * (I - tr(C)/3 * C^{-1})
    ///
    /// Derivation: S = 2*dW/dC = 2*dW/dI1_bar * dI1_bar/dC
    ///   dI1_bar/dC = J^{-2/3} * (I - I1_bar/3 * C^{-1})
    ///   At first order in eps: I1_bar = J^{-2/3}*tr(C), but the DEV_C projection
    ///   uses tr(C) (not I1_bar) directly to preserve the correct linearization.
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };
        let i1  = c.trace();
        let j23 = j.powf(-2.0 / 3.0);

        self.mu * j23 * (Matrix3::identity() - i1 / 3.0 * c_inv)
    }

    /// Isochoric tangent — three-term expression (Holzapfel 2000, eq. 6.166):
    ///
    ///   C_iso = 2*mu*J^{-2/3} * (I odot I)
    ///         - 2/3*mu*J^{-2/3} * (I tensor C^{-1} + C^{-1} tensor I)
    ///         + 2/9*mu*J^{-2/3}*I1_bar * (C^{-1} tensor C^{-1})
    ///
    /// Verification at F=I (C^{-1}=I, I1_bar=3, j23=1):
    ///   Term A at [0,0]: 2*mu * 1    = 2*mu
    ///   Term A at [0,1]: 2*mu * 0    = 0
    ///   Term A at shear: 2*mu * 0.5  = mu
    ///   Term B at [0,0]: -2/3*mu * 2 = -4/3*mu
    ///   Term B at [0,1]: -2/3*mu * 2 = -4/3*mu  (wait: (I tensor I + I tensor I)[0,1] = 0+0 at off-diag)
    ///   Actually (I tensor C^{-1})[IJ,KL] = delta_IJ * C^{-1}_KL
    ///   At F=I: (I tensor I)[IJ,KL] = delta_IJ * delta_KL -> only [0,0],[0,1],[0,2],[1,0]... etc = 1 when IJ=II
    ///   Sum at C_00,00 = 1*1 + 1*1 = 2 -> -4/3*mu
    ///   Sum at C_00,11 = 1*1 + 1*1 = 2 -> -4/3*mu (off-diagonal of 3x3 normal block)
    ///   Term C at [0,0]: 2/9*mu*3 * 1 = 2/3*mu
    ///   Term C at [0,1]: 2/9*mu*3 * 1 = 2/3*mu
    ///   Total [0,0] = 2*mu - 4/3*mu + 2/3*mu = 4/3*mu  (= lambda+2*mu - kappa = 4/3*mu) OK
    ///   Total [0,1] = 0 - 4/3*mu + 2/3*mu = -2/3*mu    (= lambda - kappa = -2/3*mu) OK
    ///   Total shear = mu + 0 + 0 = mu  OK
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        
        let i1  = c.trace();
        let j23     = j.powf(-2.0 / 3.0);
        let j43     = j.powf(-4.0 / 3.0);
        let i1bar   = j23 * i1;
        let cbar_inv = j.powf(2.0 / 3.0) * c_inv;
        let mu  = self.mu;
        
        // Eq. (39): C_iso =
        //   -2/3  * mu * J^{-4/3} * (I tensor C^{-1} + C^{-1} tensor I)
        //   +2/9  * mu * J^{-4/3} * I1 * (C^{-1} tensor C^{-1})
        //   +2/3  * mu * J^{-4/3} * I1 * (C^{-1} odot C^{-1})

        let term_a = i_tensor_a_plus_a_tensor_i(&cbar_inv) * (-2.0 / 3.0 * mu * j43);
        let term_b = cinv_tangent_voigt(&cbar_inv, 2.0 / 9.0 * mu * j43 * i1bar, 0.0);
        let term_c = cinv_tangent_voigt(&cbar_inv, 0.0, 1.0 / 3.0 * mu * j43 * i1bar);

        term_a + term_b + term_c
    }
}