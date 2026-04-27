// UTF-8
// viewer/src/state.rs — application state and material factory.
//
// Reorganized in v0.7.2:
//   - MaterialConfig replaces the flat MaterialChoice + per-param fields
//   - MaterialSource separates Manual and TissuePreset modes
//   - make_material reads from MaterialConfig exclusively
//   - Camera unchanged

use nalgebra::{DMatrix, DVector, Vector3};
use orfas_core::{
    boundary::{BoundaryConditionResult, Constraint, Load},
    MaterialLaw,
    CompressibleMaterial, CompressibleAnisotropicMaterial,
    SaintVenantKirchhoff, NeoHookeanIso, MooneyRivlinIso, OgdenIso,
    VolumetricLnJ, HolzapfelOgden, lame,
};
use orfas_tissues::presets::all_presets;
use orfas_core::mechanical_state::MechanicalState;
use orfas_core::mesh::Mesh;

// ─── Material model selector ──────────────────────────────────────────────────

/// Which mathematical model is used for the material.
///
/// Used in both Manual and TissuePreset modes — a preset carries
/// a model type that determines which parameter fields are shown.
#[derive(PartialEq, Clone)]
pub enum MaterialModel {
    SaintVenantKirchhoff,
    NeoHookean,
    MooneyRivlin,
    Ogden,
    HolzapfelOgden,
}

impl MaterialModel {
    pub fn label(&self) -> &'static str {
        match self {
            MaterialModel::SaintVenantKirchhoff => "Saint Venant-Kirchhoff",
            MaterialModel::NeoHookean           => "Neo-Hookean",
            MaterialModel::MooneyRivlin         => "Mooney-Rivlin",
            MaterialModel::Ogden                => "Ogden",
            MaterialModel::HolzapfelOgden       => "Holzapfel-Ogden",
        }
    }
}

// ─── Material source ──────────────────────────────────────────────────────────

/// Whether the material comes from manual parameter entry or a tissue preset.
#[derive(PartialEq, Clone)]
pub enum MaterialSource {
    /// User enters parameters manually.
    Manual,
    /// Parameters loaded from a tissue preset (may be further edited).
    TissuePreset,
}

// ─── Material parameters ──────────────────────────────────────────────────────

/// All material parameters in one place.
///
/// Only the fields relevant to the selected `MaterialModel` are used
/// by `make_material`. The rest are preserved so switching models
/// does not reset unrelated fields.
#[derive(Clone)]
pub struct MaterialParams {
    // SVK
    pub youngs_modulus: f64,
    pub poisson_ratio:  f64,
    // Shared
    pub density: f64,
    // NeoHookean / MooneyRivlin / Ogden / HGO — derived from E/nu or set directly
    pub mu:    f64,
    pub kappa: f64,
    // MooneyRivlin
    pub c1: f64,
    pub c2: f64,
    // Ogden
    pub ogden_mu:    Vec<f64>,
    pub ogden_alpha: Vec<f64>,
    // HolzapfelOgden
    pub k1:              f64,
    pub k2:              f64,
    pub fiber_angle_deg: f64,
}

impl Default for MaterialParams {
    fn default() -> Self {
        MaterialParams {
            youngs_modulus:  1e6,
            poisson_ratio:   0.3,
            density:         1000.0,
            mu:              384615.0,
            kappa:           833333.0,
            c1:              173077.0,  
            c2:              19231.0,    
            ogden_mu:        vec![384615.0],
            ogden_alpha:     vec![2.0],
            k1:              2363.2,
            k2:              0.8393,
            fiber_angle_deg: 49.98,
        }
    }
}

// ─── Material configuration ───────────────────────────────────────────────────

/// Full material configuration: source, model, parameters, and preset tracking.
pub struct MaterialConfig {
    /// Manual or TissuePreset.
    pub source: MaterialSource,
    /// Selected mathematical model.
    pub model: MaterialModel,
    /// Current parameter values (edited by user or loaded from preset).
    pub params: MaterialParams,
    /// Index of the selected preset in `all_presets()` (used when source = TissuePreset).
    pub preset_index: usize,
    /// True if any parameter has been modified after loading a preset.
    pub preset_modified: bool,
}

impl Default for MaterialConfig {
    fn default() -> Self {
        MaterialConfig {
            source:           MaterialSource::Manual,
            model:            MaterialModel::SaintVenantKirchhoff,
            params:           MaterialParams::default(),
            preset_index:     0,
            preset_modified:  false,
        }
    }
}

impl MaterialConfig {
    /// Load a tissue preset by index into the current config.
    ///
    /// Fills `params` from the preset's nominal parameter values,
    /// sets `model` to match the preset's material type,
    /// and resets `preset_modified` to false.
    pub fn load_preset(&mut self, index: usize) {
        let presets = all_presets();
        if index >= presets.len() { return; }

        let preset = &presets[index];
        let meta   = preset.metadata();

        self.preset_index    = index;
        self.preset_modified = false;
        self.source          = MaterialSource::TissuePreset;

        // Map model string to MaterialModel
        self.model = match meta.model {
            "Neo-Hookean"        => MaterialModel::NeoHookean,
            "Mooney-Rivlin"      => MaterialModel::MooneyRivlin,
            "Holzapfel-Ogden"    => MaterialModel::HolzapfelOgden,
            _                    => MaterialModel::SaintVenantKirchhoff,
        };

        // Fill params from preset nominal values
        let ci = &meta.confidence_intervals;
        let p  = &mut self.params;

        p.density = preset.material().density();

        match self.model {
            MaterialModel::NeoHookean => {
                p.mu    = ci.get("mu")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.mu);
                p.kappa = ci.get("kappa").map(|c| (c.min + c.max) / 2.0).unwrap_or(p.kappa);
                // Use nominal center — actual nominal stored in material()
                // Re-extract from the built material via a test deformation is not
                // practical here; use midpoint of CI as a reasonable default.
                // Presets with exact nominal values should override below.
            }
            MaterialModel::MooneyRivlin => {
                p.c1    = ci.get("c1")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.c1);
                p.c2    = ci.get("c2")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.c2);
                p.kappa = ci.get("kappa").map(|c| (c.min + c.max) / 2.0).unwrap_or(p.kappa);
            }
            MaterialModel::HolzapfelOgden => {
                p.mu    = ci.get("mu")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.mu);
                p.k1    = ci.get("k1")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.k1);
                p.k2    = ci.get("k2")   .map(|c| (c.min + c.max) / 2.0).unwrap_or(p.k2);
                p.kappa = ci.get("kappa").map(|c| (c.min + c.max) / 2.0).unwrap_or(p.kappa);
            }
            _ => {}
        }
    }

    /// Returns true if any current parameter is outside its confidence interval
    /// for the active preset.
    ///
    /// Always false when source = Manual.
    pub fn is_outside_confidence(&self) -> bool {
        if self.source != MaterialSource::TissuePreset { return false; }
        let presets = all_presets();
        if self.preset_index >= presets.len() { return false; }
        let meta = presets[self.preset_index].metadata();
        let ci   = &meta.confidence_intervals;
        let p    = &self.params;

        let check = |key: &str, val: f64| -> bool {
            ci.get(key).map_or(false, |c| !c.contains(val))
        };

        match self.model {
            MaterialModel::NeoHookean     => check("mu", p.mu) || check("kappa", p.kappa),
            MaterialModel::MooneyRivlin   => check("c1", p.c1) || check("c2", p.c2) || check("kappa", p.kappa),
            MaterialModel::HolzapfelOgden => check("mu", p.mu) || check("k1", p.k1) || check("k2", p.k2),
            _                             => false,
        }
    }
}

// ─── Other enums ──────────────────────────────────────────────────────────────

#[derive(PartialEq)]
pub enum SolverChoice {
    Direct,
    Newton,
    NewtonCachedK,
    /// Sparse Newton-Raphson with conjugate gradient — sequential assembly.
    NewtonSparse,
    /// Sparse Newton-Raphson with conjugate gradient — parallel assembly.
    NewtonSparseParallel,
}

impl SolverChoice {
    pub fn label(&self) -> &'static str {
        match self {
            SolverChoice::Direct               => "Direct LU",
            SolverChoice::Newton               => "Newton-Raphson",
            SolverChoice::NewtonCachedK        => "Newton (cached K)",
            SolverChoice::NewtonSparse         => "Newton (sparse CG)",
            SolverChoice::NewtonSparseParallel => "Newton (sparse CG parallel)",
        }
    }
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
pub struct Camera {
    pub yaw:      f32,
    pub pitch:    f32,
    pub distance: f32,
    pub target:   Vector3<f32>,
}

impl Camera {
    pub fn default() -> Self {
        Camera { yaw: 0.5, pitch: 0.3, distance: 5.0, target: Vector3::zeros() }
    }

    /// Focus the camera on a mesh bounding box.
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

    pub fn eye(&self) -> Vector3<f32> {
        self.target - self.view_direction() * self.distance
    }

    pub fn view_direction(&self) -> Vector3<f32> {
        let (sy, cy) = (self.yaw.sin(),   self.yaw.cos());
        let (sp, cp) = (self.pitch.sin(), self.pitch.cos());
        Vector3::new(sy * cp, -sp, cy * cp)
    }

    pub fn right(&self) -> Vector3<f32> {
        self.view_direction().cross(&Vector3::new(0.0, 1.0, 0.0)).normalize()
    }

    pub fn up(&self) -> Vector3<f32> {
        self.right().cross(&self.view_direction()).normalize()
    }

    pub fn project(&self, point: &Vector3<f64>, center: egui::Pos2, scale: f32) -> egui::Pos2 {
        let p       = point.cast::<f32>() - self.target;
        let rx      =  p.dot(&self.right());
        let ry      =  p.dot(&self.up());
        let rz      =  p.dot(&self.view_direction());
        let depth   = (self.distance - rz).max(1e-3);
        egui::Pos2::new(
            center.x + rx / depth * scale,
            center.y - ry / depth * scale,
        )
    }

    pub fn pan(&mut self, delta: egui::Vec2, scale: f32) {
        let upp = self.distance / scale;
        self.target -= self.right() * delta.x * upp;
        self.target += self.up()    * delta.y * upp;
    }
}

// ─── AppState ─────────────────────────────────────────────────────────────────

pub struct AppState {
    // Mesh generation
    pub nx: usize, pub ny: usize, pub nz: usize,
    pub dx: f64,   pub dy: f64,   pub dz: f64,

    // Material (unified config)
    pub material: MaterialConfig,

    // Solver / simulation
    pub solver:          SolverChoice,
    pub boundary:        BoundaryChoice,
    pub simulation_mode: SimulationMode,

    // Dynamic parameters
    pub alpha: f64,
    pub beta:  f64,
    pub dt:    f64,

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

    // Json reading
    pub last_error: Option<String>,
    pub loaded_preset_name: Option<String>,
}

impl AppState {
    pub fn new() -> Self {
        AppState {
            nx: 2, ny: 2, nz: 2,
            dx: 1., dy: 1., dz: 1.,
            material:        MaterialConfig::default(),
            solver:          SolverChoice::Direct,
            boundary:        BoundaryChoice::Penalty,
            simulation_mode: SimulationMode::Static,
            alpha: 0.1,
            beta:  0.01,
            dt:    0.01,
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
            last_error: None,
            loaded_preset_name: None,
        }
    }
}

// ─── Material factory ─────────────────────────────────────────────────────────

/// Build the active material law from the current MaterialConfig.
///
/// Reads exclusively from `state.material` — no other AppState fields
/// are used. This function is called by both static and dynamic simulation
/// paths and must remain pure (no side effects on state).
pub fn make_material(state: &AppState) -> Box<dyn MaterialLaw> {
    let p = &state.material.params;
    match state.material.model {
        MaterialModel::SaintVenantKirchhoff => Box::new(SaintVenantKirchhoff {
            youngs_modulus: p.youngs_modulus,
            poisson_ratio:  p.poisson_ratio,
            density:        p.density,
        }),
        MaterialModel::NeoHookean => Box::new(CompressibleMaterial {
            iso:     NeoHookeanIso { mu: p.mu },
            vol:     VolumetricLnJ { kappa: p.kappa },
            density: p.density,
        }),
        MaterialModel::MooneyRivlin => Box::new(CompressibleMaterial {
            iso:     MooneyRivlinIso { c1: p.c1, c2: p.c2 },
            vol:     VolumetricLnJ   { kappa: p.kappa },
            density: p.density,
        }),
        MaterialModel::Ogden => Box::new(CompressibleMaterial {
            iso:     OgdenIso::new(p.ogden_mu.clone(), p.ogden_alpha.clone())
                         .unwrap_or_else(|_| OgdenIso::new(vec![1000.0], vec![2.0]).unwrap()),
            vol:     VolumetricLnJ { kappa: p.kappa },
            density: p.density,
        }),
        MaterialModel::HolzapfelOgden => Box::new(CompressibleAnisotropicMaterial {
            iso:     NeoHookeanIso  { mu: p.mu },
            aniso:   HolzapfelOgden { k1: p.k1, k2: p.k2 },
            vol:     VolumetricLnJ  { kappa: p.kappa },
            density: p.density,
        }),
    }
}