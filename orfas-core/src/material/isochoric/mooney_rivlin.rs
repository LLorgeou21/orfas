// UTF-8
// material/isochoric/mooney_rivlin.rs — isochoric Mooney-Rivlin law.

use nalgebra::{Matrix3, Matrix6};
use crate::material::traits::IsochoricPart;
use crate::material::helpers::{
    cinv_tangent_voigt,
    i_tensor_a_plus_a_tensor_i,
    a_tensor_b_plus_b_tensor_a,
};

// ─── MooneyRivlinIso ──────────────────────────────────────────────────────────

/// Isochoric Mooney-Rivlin strain energy:
///   W_iso = c1*(I1_bar - 3) + c2*(I2_bar - 3)
///
/// where:
///   I1_bar = J^{-2/3} * tr(C)
///   I2_bar = J^{-4/3} * 0.5*(tr(C)^2 - tr(C^2))
///   C      = F^T F,   J = det(F)
///
/// Reduces to NeoHookeanIso when c2 = 0 (with mu = 2*c1).
///
/// PK2 stress (isochoric part), Bonet & Wood 2nd ed. eq. 6.29:
///   S_iso = 2*J^{-2/3} * c1 * (I - tr(C)/3 * C^{-1})
///         + 2*J^{-4/3} * c2 * (tr(C)*I - C - 2/3*I2*C^{-1})
///
/// where I2 = 0.5*(tr(C)^2 - tr(C^2)) is the second invariant (not barred).
///
/// Isochoric tangent: linearization of S_iso with respect to E,
/// composed of the c1 contribution (identical to NeoHookeanIso * 2) and
/// the c2 contribution (Holzapfel 2000, Section 6.6).
///
/// Parameters:
///   c1 — first Mooney-Rivlin constant (Pa), must be > 0
///   c2 — second Mooney-Rivlin constant (Pa), can be negative but c1+c2 > 0
#[derive(serde::Serialize, serde::Deserialize)]
pub struct MooneyRivlinIso {
    pub c1: f64,
    pub c2: f64,
}

impl MooneyRivlinIso {
    /// Safe constructor.
    /// Validates c1 > 0 and c1 + c2 > 0 (stability condition at small strains).
    pub fn new(c1: f64, c2: f64) -> Result<Self, String> {
        if c1 <= 0.0 {
            return Err(format!("MooneyRivlinIso: c1 must be > 0, got {}", c1));
        }
        if c1 + c2 <= 0.0 {
            return Err(format!(
                "MooneyRivlinIso: c1 + c2 must be > 0 for stability, got {}",
                c1 + c2
            ));
        }
        Ok(Self { c1, c2 })
    }
}

impl IsochoricPart for MooneyRivlinIso {

    /// W_iso = c1*(I1_bar - 3) + c2*(I2_bar - 3)
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }

        let c      = f.transpose() * f;
        let i1     = c.trace();
        let i2     = 0.5 * (i1 * i1 - (c * c).trace());
        let i1_bar = j.powf(-2.0 / 3.0) * i1;
        let i2_bar = j.powf(-4.0 / 3.0) * i2;

        self.c1 * (i1_bar - 3.0) + self.c2 * (i2_bar - 3.0)
    }

    /// S_iso = 2*c1 * s1 + 2*c2 * s2
    ///
    /// s1 = J^{-2/3} * (I - tr(C)/3 * C^{-1})                  <- c1 term, same as NH with mu=2*c1
    /// s2 = J^{-4/3} * (tr(C)*I - C - 2/3*I2 * C^{-1})         <- c2 term
    ///
    /// This is the DEV_C projection of dW/dC (Bonet & Wood 2008, eq. 6.29).
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix3::zeros(),
        };
        let i1  = c.trace();
        let i2  = 0.5 * (i1 * i1 - (c * c).trace());
        let j23 = j.powf(-2.0 / 3.0);
        let j43 = j.powf(-4.0 / 3.0);
        let id  = Matrix3::identity();

        let s1 = j23 * (id - i1 / 3.0 * c_inv);
        let s2 = j43 * (i1 * id - c - 2.0 / 3.0 * i2 * c_inv);

        2.0 * (self.c1 * s1 + self.c2 * s2)
    }

    /// Isochoric tangent for Mooney-Rivlin.
    ///
    /// C_iso = C_c1 + C_c2
    ///
    /// C_c1 is the NeoHookeanIso tangent with mu -> 2*c1 (three terms).
    ///
    /// C_c2 (Holzapfel 2000, Section 6.6 / Bonet & Wood 2008, eq. 6.32):
    ///   C_c2 = 4*c2*J^{-4/3} * [
    ///       I1_bar * (I odot I)
    ///     - (C odot I)_sym
    ///     - 2/3*I1_bar * (I tensor C^{-1} + C^{-1} tensor I)
    ///     + 2/3 * (C tensor C^{-1} + C^{-1} tensor C)
    ///     + 4/9*I2_bar * (C^{-1} tensor C^{-1})
    ///     - 2/3*I2_bar * (C^{-1} odot C^{-1})       <- note: this term is missing in some refs
    ///   ]
    ///
    /// Note: I1_bar and I2_bar here are the *barred* (isochoric) invariants.
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let c     = f.transpose() * f;
        let c_inv = match c.try_inverse() {
            Some(inv) => inv,
            None      => return Matrix6::zeros(),
        };

        let j23      = j.powf(-2.0 / 3.0);
        let j43      = j.powf(-4.0 / 3.0);
        let i1       = c.trace();
        let i2       = 0.5 * (i1 * i1 - (c * c).trace());
        let i1bar    = j23 * i1;
        let i2bar    = j43 * i2;
        let cbar     = j23 * c;
        let cbar_inv = j.powf(2.0 / 3.0) * c_inv;

        let c1 = 2.*self.c1;
        let c2 = 2.*self.c2;

        // Term 1: 2*J^{-4/3}*c2 * (I tensor I - I4)
        // I tensor I : delta_IJ * delta_KL
        // I4 = I odot I : 0.5*(delta_IK*delta_JL + delta_IL*delta_JK)  <- avec 0.5 car c'est la def standard
        // Notre cinv_tangent_voigt sans 0.5 -> pour avoir I4 il faut coefficient 0.5
        let mut term1 = Matrix6::zeros();
        let idx = [(0usize, 0usize), (1,1), (2,2), (0,1), (1,2), (0,2)];
        for i in 0..6 {
            for jj in 0..6 {
                let (p, q) = idx[i];
                let (r, s) = idx[jj];
                let d_ij = if p == q { 1.0 } else { 0.0 };
                let d_kl = if r == s { 1.0 } else { 0.0 };
                // I tensor I
                let i_x_i = d_ij * d_kl;
                // I4 = I odot I avec le vrai 0.5
                let i4 = 0.5 * (
                    (if p==r && q==s { 1.0 } else { 0.0 })
                + (if p==s && q==r { 1.0 } else { 0.0 })
                );
                term1[(i, jj)] = i_x_i - i4;
            }
        }
        term1 *= 2.0 * j43 * c2;

        // Term 2: -2/3*J^{-4/3}*(c1 + 2*c2*I1bar) * (Cbar^{-1} tensor I + I tensor Cbar^{-1})
        let term2 = i_tensor_a_plus_a_tensor_i(&cbar_inv)
            * (-2.0 / 3.0 * j43 * (c1 + 2.0 * c2 * i1bar));

        // Term 3: 4/3*J^{-4/3}*c2 * (Cbar^{-1} tensor Cbar + Cbar tensor Cbar^{-1})
        let term3 = a_tensor_b_plus_b_tensor_a(&cbar_inv, &cbar)
            * (4.0 / 3.0 * j43 * c2);

        // Term 4: 2/9*J^{-4/3}*(c1*I1bar + 4*c2*I2bar) * (Cbar^{-1} tensor Cbar^{-1})
        let term4 = cinv_tangent_voigt(&cbar_inv,
            2.0 / 9.0 * j43 * (c1 * i1bar + 4.0 * c2 * i2bar),
            0.0);

        // Term 5: 2/3*J^{-4/3}*(c1*I1bar + 2*c2*I2bar) * (Cbar^{-1} odot Cbar^{-1})
        // Notre convention sans 0.5 -> coefficient / 2
        let term5 = cinv_tangent_voigt(&cbar_inv,
            0.0,
            1.0 / 3.0 * j43 * (c1 * i1bar + 2.0 * c2 * i2bar));

        term1 + term2 + term3 + term4 + term5
    }
}