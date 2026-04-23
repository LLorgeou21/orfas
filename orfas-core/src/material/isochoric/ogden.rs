// UTF-8
// material/isochoric/ogden.rs — isochoric Ogden law.

use nalgebra::{Matrix3, Matrix6, SymmetricEigen};
use crate::material::traits::IsochoricPart;

// ─── OgdenIso ─────────────────────────────────────────────────────────────────

/// Isochoric Ogden strain energy (N-term sum):
///   W_iso = sum_i (mu_i/alpha_i) * (lam1_bar^alpha_i + lam2_bar^alpha_i + lam3_bar^alpha_i - 3)
///
/// where lam_k_bar = J^{-1/3} * lam_k are the isochoric principal stretches
/// and lam_k = sqrt(eigenvalue_k(C)) are the principal stretches from C = F^T F.
///
/// Reduces to NeoHookeanIso when N=1, alpha_1=2, mu_1=2*mu.
///
/// Reference: Connolly, Mackenzie, Gorash (2019),
/// "Isotropic hyperelasticity in principal stretches: explicit elasticity
/// tensors and numerical implementation", Computational Mechanics 64:1273-1288.
/// Equations (11), (12), (22), (23), (35).
///
/// Parameters:
///   mu    — shear-like moduli (Pa), Vec of length N
///   alpha — exponent parameters (dimensionless), Vec of length N
///   Stability condition: sum_i mu_i * alpha_i > 0
pub struct OgdenIso {
    pub mu:    Vec<f64>,
    pub alpha: Vec<f64>,
}

impl OgdenIso {
    /// Safe constructor.
    /// Validates mu.len() == alpha.len(), N >= 1, and sum(mu_i * alpha_i) > 0.
    pub fn new(mu: Vec<f64>, alpha: Vec<f64>) -> Result<Self, String> {
        if mu.len() != alpha.len() {
            return Err(format!(
                "OgdenIso: mu and alpha must have the same length, got {} and {}",
                mu.len(), alpha.len()
            ));
        }
        if mu.is_empty() {
            return Err("OgdenIso: must have at least one term (N >= 1)".to_string());
        }
        let stability: f64 = mu.iter().zip(alpha.iter()).map(|(m, a)| m * a).sum();
        if stability <= 0.0 {
            return Err(format!(
                "OgdenIso: sum(mu_i * alpha_i) must be > 0 for stability, got {}",
                stability
            ));
        }
        Ok(Self { mu, alpha })
    }

    /// Compute isochoric principal stretches lam_bar_a = J^{-1/3} * lam_a
    /// and principal directions N_a from the spectral decomposition of C = F^T F.
    ///
    /// Returns (eigenvalues_squared, eigenvectors, lam_bar)
    /// where eigenvalues_squared[a] = lam_a^2 and eigenvectors[:,a] = N_a.
    fn spectral_decomposition(f: &Matrix3<f64>)
        -> (nalgebra::Vector3<f64>, Matrix3<f64>, nalgebra::Vector3<f64>)
    {
        let j   = f.determinant();
        let c   = f.transpose() * f;
        let eig = SymmetricEigen::new(c);

        // eigenvalues are lam_a^2, eigenvectors[:,a] = N_a
        let lam2 = eig.eigenvalues;
        let n    = eig.eigenvectors;

        // isochoric stretches: lam_bar_a = J^{-1/3} * sqrt(lam_a^2)
        let j13 = j.powf(-1.0 / 3.0);
        let lam_bar = nalgebra::Vector3::new(
            j13 * lam2[0].sqrt(),
            j13 * lam2[1].sqrt(),
            j13 * lam2[2].sqrt(),
        );

        (lam2, n, lam_bar)
    }

    /// First derivative of W_iso with respect to isochoric principal stretch lam_bar_a:
    ///   dW/d(lam_bar_a) = sum_i mu_i * lam_bar_a^(alpha_i - 1)
    fn dw_dlam(&self, lam_bar: f64) -> f64 {
        self.mu.iter().zip(self.alpha.iter())
            .map(|(m, a)| m * lam_bar.powf(a - 1.0))
            .sum()
    }

    /// Second derivative of W_iso with respect to lam_bar_a (diagonal only, a == b):
    ///   d^2W/d(lam_bar_a)^2 = sum_i mu_i * (alpha_i - 1) * lam_bar_a^(alpha_i - 2)
    fn d2w_dlam2(&self, lam_bar: f64) -> f64 {
        self.mu.iter().zip(self.alpha.iter())
            .map(|(m, a)| m * (a - 1.0) * lam_bar.powf(a - 2.0))
            .sum()
    }

    /// Stress coefficient beta_a (Connolly eq. 12):
    ///   beta_a = lam_bar_a * dW/d(lam_bar_a) - 1/3 * sum_b lam_bar_b * dW/d(lam_bar_b)
    fn beta(&self, lam_bar: &nalgebra::Vector3<f64>) -> nalgebra::Vector3<f64> {
        let dw: nalgebra::Vector3<f64> = nalgebra::Vector3::new(
            self.dw_dlam(lam_bar[0]),
            self.dw_dlam(lam_bar[1]),
            self.dw_dlam(lam_bar[2]),
        );
        // sum_b lam_bar_b * dW/d(lam_bar_b)
        let trace_term: f64 = (0..3).map(|b| lam_bar[b] * dw[b]).sum();

        nalgebra::Vector3::new(
            lam_bar[0] * dw[0] - trace_term / 3.0,
            lam_bar[1] * dw[1] - trace_term / 3.0,
            lam_bar[2] * dw[2] - trace_term / 3.0,
        )
    }

    /// Elasticity coefficient gamma_ab (Connolly eq. 35).
    ///
    /// For Ogden, d^2W/(d_lam_a * d_lam_b) = 0 for a != b (separable energy),
    /// so the cross terms vanish and gamma_ab simplifies significantly.
    ///
    /// gamma_ab = (d^2W/d_lam_a^2 * lam_bar_a^2 + dW/d_lam_a * lam_bar_a) * delta_ab
    ///          + 1/9 * sum_{c} (d^2W/d_lam_c^2 * lam_bar_c^2 + dW/d_lam_c * lam_bar_c)
    ///          - 1/3 * (d^2W/d_lam_a^2 * lam_bar_a^2 + dW/d_lam_a * lam_bar_a) * delta_ab ... (see eq 35)
    fn gamma(&self, lam_bar: &nalgebra::Vector3<f64>) -> nalgebra::Matrix3<f64> {
        // Helper: lam_bar_a * d(lam_bar_a * dW/d_lam_bar_a) / d_lam_bar_a
        //       = d^2W/d_lam_a^2 * lam_bar_a^2 + dW/d_lam_a * lam_bar_a  (for a==b)
        //       = 0                                                          (for a!=b, Ogden)
        let h: nalgebra::Vector3<f64> = nalgebra::Vector3::new(
            self.d2w_dlam2(lam_bar[0]) * lam_bar[0] * lam_bar[0]
                + self.dw_dlam(lam_bar[0]) * lam_bar[0],
            self.d2w_dlam2(lam_bar[1]) * lam_bar[1] * lam_bar[1]
                + self.dw_dlam(lam_bar[1]) * lam_bar[1],
            self.d2w_dlam2(lam_bar[2]) * lam_bar[2] * lam_bar[2]
                + self.dw_dlam(lam_bar[2]) * lam_bar[2],
        );

        let sum_h: f64 = h.iter().sum();

        let mut gam = nalgebra::Matrix3::zeros();
        for a in 0..3 {
            for b in 0..3 {
                // eq. 35 with cross terms = 0 (Ogden separability)
                let diag_ab = if a == b { h[a] } else { 0.0 };
                let row_a   = h[a];
                let col_b   = h[b];
                gam[(a, b)] = diag_ab
                    + sum_h / 9.0
                    - row_a / 3.0
                    - col_b / 3.0;
            }
        }
        gam
    }
}

impl IsochoricPart for OgdenIso {

    /// W_iso = sum_i (mu_i/alpha_i) * (lam1_bar^alpha_i + lam2_bar^alpha_i + lam3_bar^alpha_i - 3)
    fn strain_energy_iso(&self, f: &Matrix3<f64>) -> f64 {
        let j = f.determinant();
        if j <= 0.0 { return f64::INFINITY; }

        let (_, _, lam_bar) = Self::spectral_decomposition(f);

        self.mu.iter().zip(self.alpha.iter()).map(|(m, a)| {
            (m / a) * (lam_bar[0].powf(*a) + lam_bar[1].powf(*a) + lam_bar[2].powf(*a) - 3.0)
        }).sum()
    }

    /// S_iso = sum_a beta_a * lam_a^{-2} * (N_a tensor N_a)
    ///
    /// Connolly eq. (11). beta_a defined in eq. (12).
    fn pk2_stress_iso(&self, f: &Matrix3<f64>) -> Matrix3<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix3::zeros(); }

        let (lam2, n, lam_bar) = Self::spectral_decomposition(f);
        let beta = self.beta(&lam_bar);

        let mut s = Matrix3::zeros();
        for a in 0..3 {
            let na = n.column(a);
            let na_outer = na * na.transpose(); // N_a tensor N_a
            s += beta[a] / lam2[a] * na_outer;
        }
        s
    }

    /// C_iso — material elasticity tensor in reference configuration.
    ///
    /// Connolly eq. (22), two contributions:
    ///
    /// Term 1 (diagonal a==b and off-diagonal a!=b normal part):
    ///   sum_{a,b} (gamma_ab * lam_a^{-2} * lam_b^{-2} - 2*delta_ab*beta_a*lam_a^{-4})
    ///             * (N_a tensor N_a tensor N_b tensor N_b)
    ///
    /// Term 2 (off-diagonal a!=b shear part):
    ///   sum_{a!=b} (beta_b*lam_b^{-2} - beta_a*lam_a^{-2}) / (lam_b^2 - lam_a^2)
    ///              * [(N_a tensor N_b) tensor (N_a tensor N_b + N_b tensor N_a)]
    ///
    /// When lam_a^2 ~ lam_b^2, L'Hopital's rule is applied (eq. 25):
    ///   limit = lam_b^{-4} * (gamma_bb/2 - beta_b) - lam_a^{-2}*lam_b^{-2}*gamma_ab/2
    ///
    /// Numerical tolerance for equal/similar eigenvalues: 1e-6 (Connolly Sect. 3.3).
    fn tangent_iso(&self, f: &Matrix3<f64>) -> Matrix6<f64> {
        let j = f.determinant();
        if j <= 0.0 { return Matrix6::zeros(); }

        let (lam2, n, lam_bar) = Self::spectral_decomposition(f);
        let beta  = self.beta(&lam_bar);
        let gamma = self.gamma(&lam_bar);

        // Voigt index map: 0->(0,0), 1->(1,1), 2->(2,2), 3->(0,1), 4->(1,3... wait:
        // Voigt ordering: [11,22,33,12,13,23] -> indices [(0,0),(1,1),(2,2),(0,1),(0,2),(1,2)]
        let voigt = [(0usize,0usize),(1,1),(2,2),(0,1),(1,2),(0,2)];

        // Build symmetric dyadic products sym(N_a tensor N_b) in Voigt (eq. 37)
        // sym(N_a tensor N_b)_IJ = 0.5*(N_a_I*N_b_J + N_a_J*N_b_I)
        // For a==b: just N_a tensor N_a (already symmetric)
        // For a!=b: three unique pairs (0,1), (0,2), (1,2)
        let sym_nn = |a: usize, b: usize| -> nalgebra::Vector6<f64> {
            let na = n.column(a);
            let nb = n.column(b);
            let mut v = nalgebra::Vector6::zeros();
            for (i, &(p, q)) in voigt.iter().enumerate() {
                v[i] = 0.5 * (na[p]*nb[q] + na[q]*nb[p]);
            }
            v
        };

        // Precompute sym(N_a tensor N_b) for all pairs
        let nn: Vec<Vec<nalgebra::Vector6<f64>>> = (0..3).map(|a| {
            (0..3).map(|b| sym_nn(a, b)).collect()
        }).collect();

        let mut c_iso = Matrix6::zeros();

        // ── Term 1: sum_{a,b} coeff_ab * (nn_a tensor nn_b) ─────────────────
        // coeff_ab = gamma_ab * lam_a^{-2} * lam_b^{-2} - 2*delta_ab*beta_a*lam_a^{-4}
        for a in 0..3 {
            for b in 0..3 {
                let coeff = gamma[(a,b)] / (lam2[a] * lam2[b])
                    - if a == b { 2.0 * beta[a] / (lam2[a] * lam2[a]) } else { 0.0 };

                // outer product nn_a (x) nn_b in 6x6 Voigt
                for i in 0..6 {
                    for jj in 0..6 {
                        c_iso[(i, jj)] += coeff * nn[a][a][i] * nn[b][b][jj];
                    }
                }
            }
        }

        // ── Term 2: sum_{a!=b} off-diagonal shear contribution ───────────────
        // coeff = (beta_b*lam_b^{-2} - beta_a*lam_a^{-2}) / (lam_b^2 - lam_a^2)
        // with L'Hopital when |lam_a^2 - lam_b^2| < tol
        let tol = 1e-6_f64;

        for a in 0..3 {
            for b in 0..3 {
                if a == b { continue; }

                let dlam2 = lam2[b] - lam2[a];
                let coeff = if dlam2.abs() < tol {
                    // L'Hopital's rule (Connolly eq. 25)
                    lam2[b].powi(-2) * (0.5 * gamma[(b,b)] - beta[b])
                    - 0.5 / (lam2[a] * lam2[b]) * gamma[(a,b)]
                } else {
                    (beta[b] / lam2[b] - beta[a] / lam2[a]) / dlam2
                };

                // [(N_a tensor N_b) tensor (N_a tensor N_b + N_b tensor N_a)]
                // = sym(N_a tensor N_b) tensor sym(N_a tensor N_b) * 2
                // (by symmetry argument from Connolly eq. 38)
                let sab = &nn[a][b]; // sym(N_a tensor N_b)
                for i in 0..6 {
                    for jj in 0..6 {
                        c_iso[(i, jj)] += coeff * sab[i] * sab[jj] * 2.0;
                    }
                }
            }
        }

        c_iso
    }
}