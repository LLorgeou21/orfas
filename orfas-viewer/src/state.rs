use nalgebra::{DMatrix, DVector, Vector3};
use orfas_core::{
    boundary::{BoundaryConditionResult, Constraint, Load},
    material::{CompressibleMaterial, MaterialLaw, 
           SaintVenantKirchhoff, NeoHookeanIso, MooneyRivlinIso, OgdenIso,
           VolumetricLnJ, lame},
    mechanical_state::MechanicalState,
    mesh::Mesh,
};

// ─── Enums ────────────────────────────────────────────────────────────────────

#[derive(PartialEq)]
pub enum MaterialChoice {
    SaintVenantKirchhoff,
    NeoHookean,
    MooneyRivlin,
    Ogden,
}

#[derive(PartialEq)]
pub enum SolverChoice {
    Direct,
    Newton,
    NewtonCachedK,
    /// Sparse Newton-Raphson with conjugate gradient solver — sequential assembly.
    /// Preferred for large meshes (> 1000 nodes).
    NewtonSparse,
    /// Sparse Newton-Raphson with conjugate gradient solver — parallel assembly.
    /// Recommended for very large meshes (> 8000 nodes) where parallelism pays off.
    NewtonSparseParallel,
}

#[derive(PartialEq)]
pub enum BoundaryChoice {
    Penalty,
    Elimination,
}

#[derive(PartialEq)]
pub enum SimulationMode {
    Static,
    Dynamic,
}

// ─── Camera ───────────────────────────────────────────────────────────────────

/// Orbit camera with explicit target point.
///
/// Rotation is expressed as yaw (around world Y) and pitch (around local X).
/// The camera always looks at `target`. Pan translates `target` in the screen
/// plane. Zoom moves the camera along the view direction.
///
/// Coordinate convention:
///   - yaw   : rotation around Y axis (left/right)
///   - pitch : rotation around local X axis (up/down), clamped to avoid gimbal
///   - distance : distance from target to eye, > 0
///   - target : world-space point the camera orbits around
pub struct Camera {
    pub yaw:      f32,
    pub pitch:    f32,
    pub distance: f32,
    pub target:   Vector3<f32>,
}

impl Camera {
    pub fn default() -> Self {
        Camera {
            yaw:      0.5,
            pitch:    0.3,
            distance: 5.0,
            target:   Vector3::zeros(),
        }
    }

    /// Focus the camera on a mesh by computing its bounding box and setting
    /// target to its center and distance to fit it in view.
    pub fn focus_on_mesh(&mut self, nodes: &[orfas_core::mesh::Node]) {
        if nodes.is_empty() { return; }

        let mut min = nodes[0].position;
        let mut max = nodes[0].position;
        for n in nodes {
            min = min.zip_map(&n.position, f64::min);
            max = max.zip_map(&n.position, f64::max);
        }
        let center = ((min + max) / 2.0).cast::<f32>();
        let extent = (max - min).norm() as f32;

        self.target   = center;
        self.distance = (extent * 1.5).max(0.5);
    }

    /// Returns the eye position in world space.
    pub fn eye(&self) -> Vector3<f32> {
        let dir = self.view_direction();
        self.target - dir * self.distance
    }

    /// Returns the unit vector pointing from eye to target (view direction).
    pub fn view_direction(&self) -> Vector3<f32> {
        let cos_yaw   = self.yaw.cos();
        let sin_yaw   = self.yaw.sin();
        let cos_pitch = self.pitch.cos();
        let sin_pitch = self.pitch.sin();
        Vector3::new(
            sin_yaw * cos_pitch,
            -sin_pitch,
            cos_yaw * cos_pitch,
        )
    }

    /// Right vector in world space (perpendicular to view and world-up).
    pub fn right(&self) -> Vector3<f32> {
        let forward = self.view_direction();
        let world_up = Vector3::new(0.0f32, 1.0, 0.0);
        forward.cross(&world_up).normalize()
    }

    /// Up vector in world space (perpendicular to view and right).
    pub fn up(&self) -> Vector3<f32> {
        let forward = self.view_direction();
        let right = self.right();
        right.cross(&forward).normalize()
    }

    /// Projects a world-space point onto the screen using perspective division.
    /// `center` is the screen center, `scale` is the pixel-per-unit factor.
    pub fn project(&self, point: &Vector3<f64>, center: egui::Pos2, scale: f32) -> egui::Pos2 {
        let p = point.cast::<f32>() - self.target;

        let right   = self.right();
        let up      = self.up();
        let forward = self.view_direction();

        let rx =  p.dot(&right);
        let ry =  p.dot(&up);
        let rz =  p.dot(&forward);

        let depth = self.distance - rz;
        let depth = depth.max(1e-3);

        let sx = center.x + rx / depth * scale;
        let sy = center.y - ry / depth * scale;
        egui::Pos2::new(sx, sy)
    }

    /// Pan: translate target in the screen plane by `delta` pixels.
    /// `scale` is the same pixel-per-unit factor used in project().
    pub fn pan(&mut self, delta: egui::Vec2, scale: f32) {
        let units_per_pixel = self.distance / scale;
        let right = self.right();
        let up    = self.up();
        self.target -= right * delta.x * units_per_pixel;
        self.target += up   * delta.y * units_per_pixel;
    }
}

// ─── AppState ─────────────────────────────────────────────────────────────────

pub struct AppState {
    // Mesh generation
    pub nx: usize, pub ny: usize, pub nz: usize,
    pub dx: f64,   pub dy: f64,   pub dz: f64,
    // Material
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    pub density:        f64,
    pub c1: f64,
    pub c2: f64,
    pub ogden_mu:    Vec<f64>,
    pub ogden_alpha: Vec<f64>,
    // Dynamic parameters
    pub alpha: f64,
    pub beta:  f64,
    pub dt:    f64,
    // Choices
    pub material:        MaterialChoice,
    pub solver:          SolverChoice,
    pub boundary:        BoundaryChoice,
    pub simulation_mode: SimulationMode,
    // Simulation data
    pub mesh:              Option<Mesh>,
    pub displacements:     Option<DVector<f64>>,
    pub initial_positions: Option<DVector<f64>>,
    pub mechanical_state:  Option<MechanicalState>,
    pub mass:              Option<DVector<f64>>,
    pub c_cached:          Option<DMatrix<f64>>,
    pub f_cached:          Option<DVector<f64>>,
    pub bc_result_cached:  Option<BoundaryConditionResult>,
    pub n_dofs_total:      usize,
    pub running:           bool,
    // Scene
    pub constraint:    Constraint,
    pub loads:         Vec<Load>,
    pub selected_load: Option<usize>,
    pub selected_node: Option<usize>,
    // View
    pub camera:            Camera,
    pub deformation_scale: f64,
    pub mesh_path:         Option<String>,
}

impl AppState {
    pub fn new() -> Self {
        AppState {
            nx: 2, ny: 2, nz: 2,
            dx: 1., dy: 1., dz: 1.,
            youngs_modulus: 1e6,
            poisson_ratio:  0.3,
            density:        1000.0,
            c1: 500.0,
            c2: 50.0,
            ogden_mu:    vec![1000.0],
            ogden_alpha: vec![2.0],
            alpha: 0.1,
            beta:  0.01,
            dt:    0.01,
            material:        MaterialChoice::SaintVenantKirchhoff,
            solver:          SolverChoice::Direct,
            boundary:        BoundaryChoice::Penalty,
            simulation_mode: SimulationMode::Static,
            mesh:              None,
            displacements:     None,
            initial_positions: None,
            mechanical_state:  None,
            mass:              None,
            c_cached:          None,
            f_cached:          None,
            bc_result_cached:  None,
            n_dofs_total:      0,
            running:           false,
            constraint:    Constraint { list: Vec::new() },
            loads:         Vec::new(),
            selected_load: None,
            selected_node: None,
            camera:            Camera::default(),
            deformation_scale: 30.,
            mesh_path:         None,
        }
    }
}

// ─── Material factory ─────────────────────────────────────────────────────────

/// Builds the material from the current AppState choices.
pub fn make_material(state: &AppState) -> Box<dyn MaterialLaw> {
    match state.material {
        MaterialChoice::SaintVenantKirchhoff => Box::new(SaintVenantKirchhoff {
            youngs_modulus: state.youngs_modulus,
            poisson_ratio:  state.poisson_ratio,
            density:        state.density,
        }),
        MaterialChoice::NeoHookean => {
            let (lambda, mu) = lame(state.youngs_modulus, state.poisson_ratio);
            let kappa = lambda + 2.0 / 3.0 * mu;
            Box::new(CompressibleMaterial {
                iso:     NeoHookeanIso   { mu },
                vol:     VolumetricLnJ   { kappa },
                density: state.density,
            })
        },
        MaterialChoice::MooneyRivlin => {
            let (lambda, mu) = lame(state.youngs_modulus, state.poisson_ratio);
            let kappa = lambda + 2.0 / 3.0 * mu;
            // c1 + c2 = mu/2 for small strain equivalence
            // split 90/10 by default
            let c1 = 0.45 * mu;
            let c2 = 0.05 * mu;
            Box::new(CompressibleMaterial {
                iso:     MooneyRivlinIso { c1, c2 },
                vol:     VolumetricLnJ   { kappa },
                density: state.density,
            })
        },
        MaterialChoice::Ogden => {
            let (lambda, mu) = lame(state.youngs_modulus, state.poisson_ratio);
            let kappa = lambda + 2.0 / 3.0 * mu;
            Box::new(CompressibleMaterial {
                iso:     OgdenIso::new(vec![mu], vec![2.0])
                            .unwrap_or_else(|_| OgdenIso::new(vec![1000.0], vec![2.0]).unwrap()),
                vol:     VolumetricLnJ { kappa },
                density: state.density,
            })
        },
    }
}