// UTF-8
// material/volumetric.rs — volumetric strain energy implementations.

use super::traits::VolumetricPart;

// ─── VolumetricLnJ ────────────────────────────────────────────────────────────

/// Logarithmic volumetric strain energy: W_vol = kappa/2 * (ln J)^2
///
/// This is the volumetric term used in the classic compressible Neo-Hookean
/// formulation (Ciarlet). Default volumetric part in ORFAS.
///
/// Properties:
/// - W_vol(J=1) = 0              (stress-free reference)
/// - dW/dJ(J=1) = 0              (no volumetric stress at rest)
/// - d^2W/dJ^2(J=1) = kappa      (bulk modulus at rest)
/// - W_vol -> +inf as J -> 0     (correct penalization of full compression)
///
/// Parameters:
///   kappa — bulk modulus (Pa), must be > 0
pub struct VolumetricLnJ {
    pub kappa: f64,
}

impl VolumetricLnJ {
    /// Safe constructor — validates kappa > 0.
    pub fn new(kappa: f64) -> Result<Self, String> {
        if kappa <= 0.0 {
            return Err(format!("VolumetricLnJ: kappa must be > 0, got {}", kappa));
        }
        Ok(Self { kappa })
    }
}

impl VolumetricPart for VolumetricLnJ {
    /// W_vol = kappa/2 * (ln J)^2
    fn strain_energy_vol(&self, j: f64) -> f64 {
        let ln_j = j.ln();
        self.kappa / 2.0 * ln_j * ln_j
    }

    /// dW/dJ = kappa * ln(J) / J
    fn dw_vol(&self, j: f64) -> f64 {
        self.kappa * j.ln() / j
    }

    /// d^2W/dJ^2 = kappa * (1 - ln J) / J^2
    fn d2w_vol(&self, j: f64) -> f64 {
        self.kappa * (1.0 - j.ln()) / (j * j)
    }
}

// ─── VolumetricQuad ──────────────────────────────────────────────────────────

/// Quadratic volumetric strain energy: W_vol = kappa/2 * (J - 1)^2
///
/// Alternative formulation (Penn, Ogden). Simpler derivatives but does not
/// penalize J -> 0 as strongly as VolumetricLnJ.
///
/// Properties:
/// - W_vol(J=1) = 0              (stress-free reference)
/// - dW/dJ(J=1) = 0              (no volumetric stress at rest)
/// - d^2W/dJ^2 = kappa           (constant bulk modulus, independent of J)
///
/// Parameters:
///   kappa — bulk modulus (Pa), must be > 0
pub struct VolumetricQuad {
    pub kappa: f64,
}

impl VolumetricQuad {
    /// Safe constructor — validates kappa > 0.
    pub fn new(kappa: f64) -> Result<Self, String> {
        if kappa <= 0.0 {
            return Err(format!("VolumetricQuad: kappa must be > 0, got {}", kappa));
        }
        Ok(Self { kappa })
    }
}

impl VolumetricPart for VolumetricQuad {
    /// W_vol = kappa/2 * (J - 1)^2
    fn strain_energy_vol(&self, j: f64) -> f64 {
        self.kappa / 2.0 * (j - 1.0) * (j - 1.0)
    }

    /// dW/dJ = kappa * (J - 1)
    fn dw_vol(&self, j: f64) -> f64 {
        self.kappa * (j - 1.0)
    }

    /// d^2W/dJ^2 = kappa
    fn d2w_vol(&self, _j: f64) -> f64 {
        self.kappa
    }
}