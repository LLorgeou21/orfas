use nalgebra::Vector3;

/// Stores fiber directions for each element in the mesh.
/// Each element has a list of fiber family directions (unit vectors in reference configuration).
/// For isotropic materials, use `FiberField::empty()` which provides empty slices.
pub struct FiberField {
    /// directions[element_idx][family_idx] = unit vector in reference configuration
    pub directions: Vec<Vec<Vector3<f64>>>,
}

impl FiberField {
    /// Creates a FiberField directly from a precomputed directions array.
    /// The caller is responsible for ensuring all vectors are normalized.
    /// `directions[i]` is the list of fiber family directions for element i.
    pub fn new(directions: Vec<Vec<Vector3<f64>>>) -> FiberField {
        FiberField { directions }
    }

    /// Creates a FiberField where every element shares the same fiber directions.
    /// Each direction in `dirs` is normalized at construction time.
    /// Panics if any direction has norm < 1e-10 (zero or near-zero vector).
    pub fn uniform(n_elements: usize, dirs: &[Vector3<f64>]) -> FiberField {
        let normalized: Vec<Vector3<f64>> = dirs
            .iter()
            .map(|d| {
                let norm = d.norm();
                assert!(
                    norm > 1e-10,
                    "FiberField::uniform — fiber direction has near-zero norm: {}",
                    norm
                );
                d / norm
            })
            .collect();

        FiberField {
            directions: vec![normalized; n_elements],
        }
    }

    /// Creates a FiberField with no fiber directions for each element.
    /// Used for isotropic materials — `directions_for` returns an empty slice.
    pub fn empty(n_elements: usize) -> FiberField {
        FiberField {
            directions: vec![vec![]; n_elements],
        }
    }

    /// Returns the fiber directions for a given element as a slice.
    /// Returns an empty slice for isotropic elements (created via `empty`).
    /// Panics if `elem_idx` is out of bounds.
    pub fn directions_for(&self, elem_idx: usize) -> &[Vector3<f64>] {
        &self.directions[elem_idx]
    }

    /// Creates a FiberField with one fiber family per element arranged in a helix.
    /// The helix is defined by a longitudinal `axis` (normalized) and a winding `angle_deg`
    /// measured from `axis` toward a reference perpendicular direction `up`.
    /// Useful for single-family fiber models (tendons, ligaments).
    /// Panics if `axis` or `up` has near-zero norm, or if `axis` and `up` are parallel.
    pub fn helix(
        n_elements: usize,
        axis: Vector3<f64>,
        up: Vector3<f64>,
        angle_deg: f64,
    ) -> FiberField {
        let axis_norm = axis.norm();
        assert!(
            axis_norm > 1e-10,
            "FiberField::helix — axis has near-zero norm"
        );
        let up_norm = up.norm();
        assert!(
            up_norm > 1e-10,
            "FiberField::helix — up has near-zero norm"
        );

        let axis = axis / axis_norm;
        let up = up / up_norm;

        // Build orthonormal basis: axis + two perpendicular directions
        let perp = {
            let candidate = up - axis * axis.dot(&up);
            let norm = candidate.norm();
            assert!(
                norm > 1e-10,
                "FiberField::helix — axis and up are parallel, cannot build helix basis"
            );
            candidate / norm
        };

        let angle_rad = angle_deg.to_radians();
        // Fiber direction: rotate axis toward perp by angle_deg
        let fiber = axis * angle_rad.cos() + perp * angle_rad.sin();

        FiberField {
            directions: vec![vec![fiber]; n_elements],
        }
    }

    /// Creates a FiberField with two symmetric fiber families arranged in a double helix.
    /// Family 1: +angle_deg from axis, Family 2: -angle_deg from axis.
    /// This is the classical HGO configuration for arterial walls and myocardium.
    /// Panics under the same conditions as `helix`.
    pub fn helix_two_families(
        n_elements: usize,
        axis: Vector3<f64>,
        up: Vector3<f64>,
        angle_deg: f64,
    ) -> FiberField {
        let axis_norm = axis.norm();
        assert!(
            axis_norm > 1e-10,
            "FiberField::helix_two_families — axis has near-zero norm"
        );
        let up_norm = up.norm();
        assert!(
            up_norm > 1e-10,
            "FiberField::helix_two_families — up has near-zero norm"
        );

        let axis = axis / axis_norm;
        let up = up / up_norm;

        let perp = {
            let candidate = up - axis * axis.dot(&up);
            let norm = candidate.norm();
            assert!(
                norm > 1e-10,
                "FiberField::helix_two_families — axis and up are parallel"
            );
            candidate / norm
        };

        let angle_rad = angle_deg.to_radians();
        let fiber_pos = axis * angle_rad.cos() + perp * angle_rad.sin();
        let fiber_neg = axis * angle_rad.cos() - perp * angle_rad.sin();

        FiberField {
            directions: vec![vec![fiber_pos, fiber_neg]; n_elements],
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_helix_two_families_angle_zero_aligns_with_axis() {
        let axis: nalgebra::Matrix<f64, nalgebra::Const<3>, nalgebra::Const<1>, nalgebra::ArrayStorage<f64, 3, 1>> = Vector3::new(1.0, 0.0, 0.0);
        let up = Vector3::new(0.0, 1.0, 0.0);
        let ff = FiberField::helix_two_families(1, axis, up, 0.0);
        let dirs = ff.directions_for(0);
        assert!((dirs[0] - axis).norm() < 1e-10);
        assert!((dirs[1] - axis).norm() < 1e-10);
    }

    #[test]
    fn test_helix_two_families_angle_90_perpendicular_to_axis() {
        let axis = Vector3::new(1.0, 0.0, 0.0);
        let up = Vector3::new(0.0, 1.0, 0.0);
        let ff = FiberField::helix_two_families(1, axis, up, 90.0);
        let dirs = ff.directions_for(0);
        assert!(dirs[0].dot(&axis).abs() < 1e-10);
        assert!(dirs[1].dot(&axis).abs() < 1e-10);
        // The two families must be symmetric (opposite y components)
        assert!((dirs[0] + dirs[1]).norm() < 1e-10);
    }
}