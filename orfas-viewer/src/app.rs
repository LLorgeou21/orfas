// UTF-8
// viewer/src/app.rs — main application and UI panels.
//
// UI structure (v0.7.2):
//   Left panel (scrollable):
//     - Section: Mesh
//     - Section: Material  (tabs: Manual | Tissue Preset)
//     - Section: Solver & Simulation
//     - Section: Boundary Conditions
//     - Section: Node Inspector
//   Central panel:
//     - 3D viewport

use nalgebra::{DVector, Vector3};
use orfas_core::{
    assembler::{Assembler},
    boundary::{Constraint, FixedNode, Load},
    integrator::{ImplicitEulerIntegrator, IntegratorMethod},
    solver::DirectSolver,
    element::Tet4,
};
use orfas_io::read_vtk;
use orfas_tissues::presets::all_presets;

use crate::render::{depth_sorted_nodes, screen_to_node};
use crate::simulation::{init_simulation_dynamic, run_simulation_static};
use crate::state::{
    AppState, BoundaryChoice, MaterialModel, MaterialSource,
    SimulationMode, SolverChoice, make_material,
};

// ─── App struct ───────────────────────────────────────────────────────────────

pub struct MyEguiApp {
    pub state: AppState,
}

impl MyEguiApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self { state: AppState::new() }
    }
}

// ─── eframe::App ──────────────────────────────────────────────────────────────

impl eframe::App for MyEguiApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.dynamic_step(ui);

        egui::Panel::left("left_panel").show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                self.ui_section_mesh(ui);
                ui.add_space(4.0);
                self.ui_section_material(ui);
                ui.add_space(4.0);
                self.ui_section_solver(ui);
                ui.add_space(4.0);
                self.ui_section_boundary(ui);
                ui.add_space(4.0);
                self.ui_section_node_inspector(ui);
            });
        });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            self.ui_viewport(ui);
        });
    }
}

// ─── Dynamic step ─────────────────────────────────────────────────────────────

impl MyEguiApp {
    fn dynamic_step(&mut self, ui: &egui::Ui) {
        if !self.state.running { return; }

        let material = make_material(&self.state);
        if let (
            Some(ref mut mech),
            Some(ref mass),
            Some(ref c),
            Some(ref f),
            Some(ref bc_result),
            Some(ref init),
        ) = (
            self.state.mechanical_state.as_mut(),
            self.state.mass.as_ref(),
            self.state.c_cached.as_ref(),
            self.state.f_cached.as_ref(),
            self.state.bc_result_cached.as_ref(),
            self.state.initial_positions.as_ref(),
        ) {
            if let Some(ref mesh) = self.state.mesh {
                let assembler = Assembler::<Tet4>::new(mesh);
                let integrator = ImplicitEulerIntegrator::default();
                match integrator.step(
                    mech, mass, c, f, self.state.dt,
                    &assembler, mesh, material.as_ref(), bc_result, &DirectSolver,
                ) {
                    Ok(_) => {
                        let disp_reduced = &mech.position - *init;
                        let disp_full    = bc_result.reconstruct_ref(
                            &DVector::from_vec(disp_reduced.iter().copied().collect()),
                            self.state.n_dofs_total,
                        );
                        self.state.displacements = Some(disp_full);
                    }
                    Err(e) => {
                        println!("Integration error: {:?}", e);
                        self.state.running = false;
                    }
                }
            }
        }
        ui.ctx().request_repaint();
    }
}

// ─── Section: Mesh ────────────────────────────────────────────────────────────

impl MyEguiApp {
    /// Mesh section — generation parameters and file loading.
    fn ui_section_mesh(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Mesh", |ui| {
            ui.group(|ui| {
                ui.label("Generated mesh");
                egui::Grid::new("mesh_grid")
                    .num_columns(2)
                    .spacing([8.0, 4.0])
                    .show(ui, |ui| {
                        ui.label("nx"); ui.add(egui::DragValue::new(&mut self.state.nx).speed(1).range(1..=100)); ui.end_row();
                        ui.label("ny"); ui.add(egui::DragValue::new(&mut self.state.ny).speed(1).range(1..=100)); ui.end_row();
                        ui.label("nz"); ui.add(egui::DragValue::new(&mut self.state.nz).speed(1).range(1..=100)); ui.end_row();
                        ui.label("dx (m)"); ui.add(egui::DragValue::new(&mut self.state.dx).speed(0.01)); ui.end_row();
                        ui.label("dy (m)"); ui.add(egui::DragValue::new(&mut self.state.dy).speed(0.01)); ui.end_row();
                        ui.label("dz (m)"); ui.add(egui::DragValue::new(&mut self.state.dz).speed(0.01)); ui.end_row();
                    });
            });

            ui.group(|ui| {
                ui.label("Load from file (.vtk)");
                if ui.button("Browse...").clicked() {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("VTK", &["vtk"])
                        .pick_file()
                    {
                        self.state.mesh_path = Some(path.display().to_string());
                    }
                }
                if let Some(ref path) = self.state.mesh_path.clone() {
                    ui.label(egui::RichText::new(path).small());
                    if ui.button("Load mesh").clicked() {
                        match read_vtk(path) {
                            Ok(mesh) => {
                                self.state.constraint    = Constraint { list: Vec::new() };
                                self.state.loads         = Vec::new();
                                self.state.selected_load = None;
                                self.state.camera.focus_on_mesh(&mesh.nodes);
                                self.state.mesh          = Some(mesh);
                            }
                            Err(e) => println!("VTK load error: {:?}", e),
                        }
                    }
                }
            });
        });
    }
}

// ─── Section: Material ────────────────────────────────────────────────────────

impl MyEguiApp {
    /// Material section — Manual tab and Tissue Preset tab.
    fn ui_section_material(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Material", |ui| {
            // ── Source tabs ───────────────────────────────────────────────────
            ui.horizontal(|ui| {
                if ui.selectable_label(
                    self.state.material.source == MaterialSource::Manual,
                    "Manual",
                ).clicked() {
                    self.state.material.source = MaterialSource::Manual;
                }
                if ui.selectable_label(
                    self.state.material.source == MaterialSource::TissuePreset,
                    "Tissue preset",
                ).clicked() {
                    self.state.material.source = MaterialSource::TissuePreset;
                }
            });

            ui.separator();

            if let Some(ref err) = self.state.last_error {
                ui.colored_label(egui::Color32::RED, err);
            }

            match self.state.material.source {
                MaterialSource::Manual       => self.ui_material_manual(ui),
                MaterialSource::TissuePreset => self.ui_material_preset(ui),
            }
        });
    }

    /// Manual material tab — model selector + parameter sliders.
    fn ui_material_manual(&mut self, ui: &mut egui::Ui) {
        // Model selector
        egui::ComboBox::from_label("Model")
            .selected_text(self.state.material.model.label())
            .show_ui(ui, |ui| {
                for model in [
                    MaterialModel::SaintVenantKirchhoff,
                    MaterialModel::NeoHookean,
                    MaterialModel::MooneyRivlin,
                    MaterialModel::Ogden,
                    MaterialModel::HolzapfelOgden,
                ] {
                    let label = model.label();
                    ui.selectable_value(&mut self.state.material.model, model, label);
                }
            });

        ui.add_space(4.0);

        // Parameter fields — only relevant to current model
        egui::Grid::new("mat_params_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                self.ui_material_params(ui);
            });
    }

    /// Tissue preset tab — preset selector, metadata tooltip, parameter sliders.
    fn ui_material_preset(&mut self, ui: &mut egui::Ui) {
        let presets = all_presets();

        // Preset selector
        let current_name = presets
            .get(self.state.material.preset_index)
            .map(|p| p.metadata().name)
            .unwrap_or("—");

        egui::ComboBox::from_label("Tissue")
            .selected_text(current_name)
            .show_ui(ui, |ui| {
                for (i, preset) in presets.iter().enumerate() {
                    let meta = preset.metadata();
                    let resp = ui.selectable_value(
                        &mut self.state.material.preset_index,
                        i,
                        meta.name,
                    );
                    // Tooltip with reference on hover
                    resp.on_hover_ui(|ui| {
                        ui.label(egui::RichText::new(meta.reference).small());
                        if let Some(doi) = meta.doi {
                            ui.label(egui::RichText::new(format!("DOI: {}", doi)).small());
                        }
                        ui.label(egui::RichText::new(format!("Protocol: {}", meta.protocol)).small());
                        if !meta.notes.is_empty() {
                            ui.label(egui::RichText::new(meta.notes).small().italics());
                        }
                    });
                }
            });

        // Load preset button
        if ui.button("Load preset").clicked() {
            let idx = self.state.material.preset_index;
            self.state.loaded_preset_name = None;
            self.state.material.load_preset(idx);
        }

        if ui.button("Load from JSON...").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("JSON", &["json"])
                .pick_file()
            {
                match orfas_tissues::load_preset_from_file(&path) {
                    Ok(loaded) => {
                        self.state.last_error = None;
                        self.state.loaded_preset_name = Some(loaded.metadata.name.clone());
                        let p = &mut self.state.material.params;
                        let params = &loaded.metadata.parameters;
                        p.density = *params.get("density").unwrap_or(&p.density);
                        self.state.material.model = match loaded.metadata.model.as_str() {
                            "neo_hookean"             => {
                                p.mu    = *params.get("mu")   .unwrap_or(&p.mu);
                                p.kappa = *params.get("kappa").unwrap_or(&p.kappa);
                                MaterialModel::NeoHookean
                            }
                            "mooney_rivlin"           => {
                                p.c1    = *params.get("c1")   .unwrap_or(&p.c1);
                                p.c2    = *params.get("c2")   .unwrap_or(&p.c2);
                                p.kappa = *params.get("kappa").unwrap_or(&p.kappa);
                                MaterialModel::MooneyRivlin
                            }
                            "holzapfel_ogden"         => {
                                p.mu    = *params.get("mu")   .unwrap_or(&p.mu);
                                p.k1    = *params.get("k1")   .unwrap_or(&p.k1);
                                p.k2    = *params.get("k2")   .unwrap_or(&p.k2);
                                p.kappa = *params.get("kappa").unwrap_or(&p.kappa);
                                MaterialModel::HolzapfelOgden
                            }
                            "saint_venant_kirchhoff"  => {
                                p.youngs_modulus = *params.get("youngs_modulus").unwrap_or(&p.youngs_modulus);
                                p.poisson_ratio  = *params.get("poisson_ratio") .unwrap_or(&p.poisson_ratio);
                                MaterialModel::SaintVenantKirchhoff
                            }
                            other => {
                                self.state.last_error = Some(format!("Unknown model: '{}'", other));
                                return;
                            }
                        };
                        self.state.material.source       = MaterialSource::TissuePreset;
                        self.state.material.preset_modified = false;
                    }
                    Err(e) => {
                        self.state.last_error = Some(format!("JSON load error: {}", e));
                    }
                }
            }
        }

        if let Some(ref name) = self.state.loaded_preset_name {
            ui.colored_label(egui::Color32::GREEN, format!("✓ {}", name));
        }
        ui.add_space(4.0);

        // Modified badge — shown if params are outside CI
        if self.state.material.is_outside_confidence() {
            ui.horizontal(|ui| {
                ui.label(
                    egui::RichText::new("⚠ modified")
                        .small()
                        .color(egui::Color32::from_rgb(220, 180, 60)),
                );
            });
        }

        // Metadata display (model name)
        if let Some(preset) = presets.get(self.state.material.preset_index) {
            let meta = preset.metadata();
            ui.label(
                egui::RichText::new(format!("Model: {}", meta.model))
                    .small()
                    .color(egui::Color32::GRAY),
            );
            // Full reference on hover over the label
            ui.label(egui::RichText::new(meta.reference).small())
                .on_hover_ui(|ui| {
                    if let Some(doi) = meta.doi {
                        ui.label(format!("DOI: {}", doi));
                    }
                    ui.label(format!("Protocol: {}", meta.protocol));
                });
        }

        ui.add_space(4.0);

        // Editable parameters (same fields as manual, but pre-filled)
        egui::Grid::new("preset_params_grid")
            .num_columns(2)
            .spacing([8.0, 4.0])
            .show(ui, |ui| {
                self.ui_material_params(ui);
            });
    }

    /// Shared parameter grid — rendered inside both Manual and Preset tabs.
    ///
    /// Shows only the fields relevant to the current model.
    /// Must be called inside an `egui::Grid`.
    fn ui_material_params(&mut self, ui: &mut egui::Ui) {
        let p = &mut self.state.material.params;

        // Density — always shown
        ui.label("Density (kg/m³)");
        ui.add(egui::DragValue::new(&mut p.density).speed(1.0));
        ui.end_row();

        match self.state.material.model {
            MaterialModel::SaintVenantKirchhoff => {
                ui.label("Young's modulus (Pa)");
                ui.add(egui::DragValue::new(&mut p.youngs_modulus).speed(100.0));
                ui.end_row();
                ui.label("Poisson ratio");
                ui.add(egui::DragValue::new(&mut p.poisson_ratio).speed(0.01).range(0.0..=0.499));
                ui.end_row();
            }
            MaterialModel::NeoHookean => {
                ui.label("mu (Pa)");
                ui.add(egui::DragValue::new(&mut p.mu).speed(10.0));
                ui.end_row();
                ui.label("kappa (Pa)");
                ui.add(egui::DragValue::new(&mut p.kappa).speed(100.0));
                ui.end_row();
            }
            MaterialModel::MooneyRivlin => {
                ui.label("c1 (Pa)");
                ui.add(egui::DragValue::new(&mut p.c1).speed(10.0));
                ui.end_row();
                ui.label("c2 (Pa)");
                ui.add(egui::DragValue::new(&mut p.c2).speed(10.0));
                ui.end_row();
                ui.label("kappa (Pa)");
                ui.add(egui::DragValue::new(&mut p.kappa).speed(100.0));
                ui.end_row();
            }
            MaterialModel::Ogden => {
                ui.label("mu[0] (Pa)");
                if let Some(v) = p.ogden_mu.first_mut() {
                    ui.add(egui::DragValue::new(v).speed(10.0));
                }
                ui.end_row();
                ui.label("alpha[0]");
                if let Some(v) = p.ogden_alpha.first_mut() {
                    ui.add(egui::DragValue::new(v).speed(0.1));
                }
                ui.end_row();
                ui.label("kappa (Pa)");
                ui.add(egui::DragValue::new(&mut p.kappa).speed(100.0));
                ui.end_row();
            }
            MaterialModel::HolzapfelOgden => {
                ui.label("mu (Pa)");
                ui.add(egui::DragValue::new(&mut p.mu).speed(1.0));
                ui.end_row();
                ui.label("k1 (Pa)");
                ui.add(egui::DragValue::new(&mut p.k1).speed(10.0));
                ui.end_row();
                ui.label("k2");
                ui.add(egui::DragValue::new(&mut p.k2).speed(0.001));
                ui.end_row();
                ui.label("kappa (Pa)");
                ui.add(egui::DragValue::new(&mut p.kappa).speed(100.0));
                ui.end_row();
                ui.label("Fiber angle (deg)");
                ui.add(egui::DragValue::new(&mut p.fiber_angle_deg).speed(0.5).range(0.0..=90.0));
                ui.end_row();
            }
        }

        // Deformation scale — always shown at bottom of params
        ui.label("Deformation scale");
        ui.add(egui::DragValue::new(&mut self.state.deformation_scale).speed(0.5));
        ui.end_row();
    }
}

// ─── Section: Solver & Simulation ─────────────────────────────────────────────

impl MyEguiApp {
    /// Solver and simulation section — mode, solver choice, dynamic params, run buttons.
    fn ui_section_solver(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Solver & Simulation", |ui| {
            egui::Grid::new("solver_grid")
                .num_columns(2)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    // Simulation mode
                    ui.label("Mode");
                    egui::ComboBox::from_id_source("mode_combo")
                        .selected_text(match self.state.simulation_mode {
                            SimulationMode::Static  => "Static",
                            SimulationMode::Dynamic => "Dynamic",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.state.simulation_mode, SimulationMode::Static,  "Static");
                            ui.selectable_value(&mut self.state.simulation_mode, SimulationMode::Dynamic, "Dynamic");
                        });
                    ui.end_row();

                    // Solver
                    ui.label("Solver");
                    egui::ComboBox::from_id_source("solver_combo")
                        .selected_text(self.state.solver.label())
                        .show_ui(ui, |ui| {
                            for s in [
                                SolverChoice::Direct,
                                SolverChoice::Newton,
                                SolverChoice::NewtonCachedK,
                                SolverChoice::NewtonSparse,
                                SolverChoice::NewtonSparseParallel,
                            ] {
                                let label = s.label();
                                ui.selectable_value(&mut self.state.solver, s, label);
                            }
                        });
                    ui.end_row();

                    // Boundary method
                    ui.label("Boundary");
                    egui::ComboBox::from_id_source("boundary_combo")
                        .selected_text(match self.state.boundary {
                            BoundaryChoice::Penalty     => "Penalty",
                            BoundaryChoice::Elimination => "Elimination",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Penalty,     "Penalty");
                            ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Elimination, "Elimination");
                        });
                    ui.end_row();

                    // Dynamic-only params
                    if self.state.simulation_mode == SimulationMode::Dynamic {
                        ui.label("dt (s)");
                        ui.add(egui::DragValue::new(&mut self.state.dt).speed(0.001).range(1e-6..=1.0));
                        ui.end_row();
                        ui.label("Rayleigh alpha");
                        ui.add(egui::DragValue::new(&mut self.state.alpha).speed(0.001));
                        ui.end_row();
                        ui.label("Rayleigh beta");
                        ui.add(egui::DragValue::new(&mut self.state.beta).speed(0.0001));
                        ui.end_row();
                    }
                });

            ui.add_space(4.0);

            // Run buttons
            ui.horizontal(|ui| {
                match self.state.simulation_mode {
                    SimulationMode::Static => {
                        if ui.button("▶ Run static").clicked() {
                            run_simulation_static(&mut self.state);
                        }
                    }
                    SimulationMode::Dynamic => {
                        if ui.button("⚙ Init dynamic").clicked() {
                            init_simulation_dynamic(&mut self.state);
                        }
                        let run_label = if self.state.running { "⏹ Stop" } else { "▶ Run" };
                        if ui.button(run_label).clicked() {
                            self.state.running = !self.state.running;
                        }
                    }
                }
                if ui.button("↺ Reset").clicked() {
                    self.state.displacements     = None;
                    self.state.mechanical_state  = None;
                    self.state.mass              = None;
                    self.state.c_cached          = None;
                    self.state.f_cached          = None;
                    self.state.bc_result_cached  = None;
                    self.state.initial_positions = None;
                    self.state.running           = false;
                    self.state.last_error        = None;
                    self.state.loaded_preset_name = None;
                }
            });
        });
    }
}

// ─── Section: Boundary Conditions ─────────────────────────────────────────────

impl MyEguiApp {
    /// Boundary conditions section — constraints and loads.
    fn ui_section_boundary(&mut self, ui: &mut egui::Ui) {
        ui.collapsing("Boundary Conditions", |ui| {
            // Constraints
            ui.group(|ui| {
                ui.label(format!("Fixed nodes: {}", self.state.constraint.list.len()));
                if !self.state.constraint.list.is_empty() {
                    if ui.small_button("Clear all").clicked() {
                        self.state.constraint.list.clear();
                    }
                }
            });

            ui.add_space(2.0);

            // Loads
            ui.group(|ui| {
                ui.horizontal(|ui| {
                    ui.label("Loads");
                    if ui.small_button("+ Add load").clicked() {
                        self.state.loads.push(Load {
                            list:  Vec::new(),
                            force: Vector3::zeros(),
                        });
                    }
                });

                let mut to_remove: Option<usize> = None;
                for (i, load) in self.state.loads.iter_mut().enumerate() {
                    let selected = self.state.selected_load == Some(i);
                    if ui.selectable_label(selected, format!("Load {}", i)).clicked() {
                        self.state.selected_load = Some(i);
                    }
                    if selected {
                        egui::Grid::new(format!("load_grid_{}", i))
                            .num_columns(2)
                            .spacing([8.0, 2.0])
                            .show(ui, |ui| {
                                ui.label("fx (N)"); ui.add(egui::DragValue::new(&mut load.force.x).speed(1.0)); ui.end_row();
                                ui.label("fy (N)"); ui.add(egui::DragValue::new(&mut load.force.y).speed(1.0)); ui.end_row();
                                ui.label("fz (N)"); ui.add(egui::DragValue::new(&mut load.force.z).speed(1.0)); ui.end_row();
                            });
                        ui.label(format!("{} nodes assigned", load.list.len()));
                        if ui.small_button("Remove load").clicked() {
                            to_remove = Some(i);
                        }
                    }
                }
                if let Some(i) = to_remove {
                    self.state.loads.remove(i);
                    if self.state.selected_load == Some(i) {
                        self.state.selected_load = None;
                    }
                }
            });
        });
    }
}

// ─── Section: Node Inspector ──────────────────────────────────────────────────

impl MyEguiApp {
    /// Node inspector section — shown when a node is selected in the viewport.
    fn ui_section_node_inspector(&mut self, ui: &mut egui::Ui) {
        let idx = match self.state.selected_node {
            Some(i) => i,
            None    => return,
        };
        if self.state.mesh.is_none() { return; }

        ui.collapsing("Node Inspector", |ui| {
            let pos = self.state.mesh.as_ref().unwrap().nodes[idx].position;
            ui.label(format!("Node {}", idx));
            ui.label(format!("x = {:.4}  y = {:.4}  z = {:.4}", pos.x, pos.y, pos.z));

            if let Some(ref displacements) = self.state.displacements {
                let ux = displacements[3 * idx];
                let uy = displacements[3 * idx + 1];
                let uz = displacements[3 * idx + 2];
                let amp = (ux*ux + uy*uy + uz*uz).sqrt();
                ui.label(format!("u = ({:.4e}, {:.4e}, {:.4e})  |u| = {:.4e}", ux, uy, uz, amp));
            }

            ui.add_space(4.0);

            if ui.button("Add to constraint").clicked() {
                if !self.state.constraint.list.iter().any(|n| n.indice == idx) {
                    self.state.constraint.list.push(FixedNode::all(idx));
                }
            }
            if self.state.constraint.list.iter().any(|n| n.indice == idx) {
                ui.label(egui::RichText::new("In constraint").color(egui::Color32::RED));
                if ui.small_button("Remove from constraint").clicked() {
                    self.state.constraint.list.retain(|n| n.indice != idx);
                }
            }

            if let Some(load_idx) = self.state.selected_load {
                if load_idx < self.state.loads.len() {
                    if ui.button(format!("Add to Load {}", load_idx)).clicked() {
                        if !self.state.loads[load_idx].list.contains(&idx) {
                            self.state.loads[load_idx].list.push(idx);
                        }
                    }
                    if self.state.loads[load_idx].list.contains(&idx) {
                        ui.label(egui::RichText::new(format!("In Load {}", load_idx))
                            .color(egui::Color32::YELLOW));
                        if ui.small_button("Remove from load").clicked() {
                            self.state.loads[load_idx].list.retain(|&n| n != idx);
                        }
                    }
                }
            }
        });
    }
}

// ─── Viewport ─────────────────────────────────────────────────────────────────

impl MyEguiApp {
    fn ui_viewport(&mut self, ui: &mut egui::Ui) {
        let (response, painter) = ui.allocate_painter(
            ui.available_size(),
            egui::Sense::drag() | egui::Sense::click(),
        );
        let rect   = response.rect;
        let center = rect.center();
        let scale  = 100.0f32;

        // Camera controls
        if response.dragged() {
            let delta     = response.drag_delta();
            let modifiers = ui.input(|i| i.modifiers);
            if modifiers.shift {
                self.state.camera.pan(delta, scale);
            } else {
                self.state.camera.yaw   += delta.x * 0.01;
                self.state.camera.pitch += delta.y * 0.01;
                self.state.camera.pitch  = self.state.camera.pitch.clamp(
                    -std::f32::consts::FRAC_PI_2 + 0.05,
                     std::f32::consts::FRAC_PI_2 - 0.05,
                );
            }
        }
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            self.state.camera.distance *= (1.0 - scroll * 0.001).max(0.01);
            self.state.camera.distance  = self.state.camera.distance.max(0.01);
        }

        // Node picking
        if response.clicked() {
            if let Some(pos) = response.interact_pointer_pos() {
                if let Some(ref mesh) = self.state.mesh {
                    self.state.selected_node =
                        Some(screen_to_node(&pos, mesh, &self.state.camera, center, scale));
                }
            }
        }

        let mesh = match self.state.mesh.as_ref() {
            Some(m) => m,
            None    => return,
        };

        let fixed_indices: std::collections::HashSet<usize> =
            self.state.constraint.list.iter().map(|n| n.indice).collect();
        let loaded_indices: std::collections::HashSet<usize> =
            self.state.loads.iter().flat_map(|l| l.list.iter().copied()).collect();

        // Draw mesh
        if let Some(ref displacements) = self.state.displacements {
            let scale_def    = self.state.deformation_scale;
            let deformed_pos = |i: usize| -> Vector3<f64> {
                mesh.nodes[i].position + Vector3::new(
                    displacements[3*i]     * scale_def,
                    displacements[3*i + 1] * scale_def,
                    displacements[3*i + 2] * scale_def,
                )
            };

            // Edges
            for tetra in &mesh.elements {
                let pts: Vec<egui::Pos2> = tetra.iter()
                    .map(|&i| self.state.camera.project(&deformed_pos(i), center, scale))
                    .collect();
                for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                    painter.line_segment([pts[a], pts[b]],
                        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                }
            }

            // Displacement colormap
            let amplitudes: Vec<f64> = (0..mesh.nodes.len()).map(|i| {
                let ux = displacements[3*i];
                let uy = displacements[3*i + 1];
                let uz = displacements[3*i + 2];
                (ux*ux + uy*uy + uz*uz).sqrt()
            }).collect();
            let max_amp = amplitudes.iter().cloned().fold(0.0f64, f64::max);

            let sorted = depth_sorted_nodes(mesh, &self.state.camera);
            for i in sorted {
                let p      = self.state.camera.project(&deformed_pos(i), center, scale);
                let t      = if max_amp > 1e-15 { (amplitudes[i] / max_amp) as f32 } else { 0.0 };
                let color  = egui::Color32::from_rgb((t * 255.0) as u8, 0, ((1.0 - t) * 255.0) as u8);
                let radius = if self.state.selected_node == Some(i) { 6.0 } else { 3.0 };
                painter.circle_filled(p, radius, color);
            }

        } else {
            // Reference mesh
            for tetra in &mesh.elements {
                let pts: Vec<egui::Pos2> = tetra.iter()
                    .map(|&i| self.state.camera.project(&mesh.nodes[i].position, center, scale))
                    .collect();
                for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                    painter.line_segment([pts[a], pts[b]],
                        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                }
            }

            let sorted = depth_sorted_nodes(mesh, &self.state.camera);
            for i in sorted {
                let p = self.state.camera.project(&mesh.nodes[i].position, center, scale);
                if fixed_indices.contains(&i)  { painter.circle_filled(p, 5.0, egui::Color32::RED);    }
                if loaded_indices.contains(&i) { painter.circle_filled(p, 5.0, egui::Color32::YELLOW); }
                let color = if self.state.selected_node == Some(i) {
                    egui::Color32::WHITE
                } else {
                    egui::Color32::GREEN
                };
                painter.circle_filled(p, 3.0, color);
            }
        }

        // Status bar
        let status = format!(
            "{} nodes  |  {} elements  |  {}{}",
            mesh.nodes.len(),
            mesh.elements.len(),
            self.state.material.model.label(),
            if self.state.material.is_outside_confidence() { "  ⚠" } else { "" },
        );
        painter.text(
            egui::Pos2::new(rect.left() + 8.0, rect.bottom() - 18.0),
            egui::Align2::LEFT_CENTER,
            status,
            egui::FontId::proportional(11.0),
            egui::Color32::GRAY,
        );
    }
}