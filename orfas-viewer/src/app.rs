use nalgebra::{DVector, Vector3};
use orfas_core::{
    assembler::{Assembler, LinearBMatrix},
    boundary::{Constraint, FixedNode, Load},
    integrator::{ImplicitEulerIntegrator, IntegratorMethod},
    solver::DirectSolver,
};
use orfas_io::read_vtk;

use crate::render::{depth_sorted_nodes, screen_to_node};
use crate::simulation::{init_simulation_dynamic, run_simulation_static};
use crate::state::{
    AppState, BoundaryChoice, MaterialChoice, SimulationMode, SolverChoice, make_material,
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

        egui::Panel::left("params").show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                self.ui_params(ui);
                self.ui_mesh_loading(ui);
                self.ui_simulation_buttons(ui);
                self.ui_loads(ui);
                self.ui_node_inspector(ui);
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
                let assembler  = Assembler::new(mesh);
                let integrator = ImplicitEulerIntegrator::default();
                match integrator.step::<LinearBMatrix>(
                    mech, mass, c, f, self.state.dt,
                    &assembler, mesh, material.as_ref(), bc_result, &DirectSolver,
                ) {
                    Ok(_) => {
                        let disp_reduced = &mech.position - *init;
                        let disp_full = bc_result.reconstruct_ref(
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

// ─── Left panel ───────────────────────────────────────────────────────────────

impl MyEguiApp {
    fn ui_params(&mut self, ui: &mut egui::Ui) {
        ui.group(|ui| {
            egui::ComboBox::from_label("Material")
                .selected_text(match self.state.material {
                    MaterialChoice::SaintVenantKirchhoff => "Saint Venant-Kirchhoff",
                    MaterialChoice::NeoHookean           => "Neo-Hookean",
                    MaterialChoice::MooneyRivlin         => "Mooney-Rivlin",
                    MaterialChoice::Ogden                => "Ogden",
                    MaterialChoice::HolzapfelOgden       => "Holzapfel-Ogden"
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.state.material,
                        MaterialChoice::SaintVenantKirchhoff,
                        "Saint Venant-Kirchhoff",
                    );
                    ui.selectable_value(
                        &mut self.state.material,
                        MaterialChoice::NeoHookean,
                        "Neo-Hookean",
                    );
                    ui.selectable_value(
                        &mut self.state.material,
                        MaterialChoice::MooneyRivlin,
                        "Mooney-Rivlin",
                    );
                    ui.selectable_value(
                        &mut self.state.material,
                        MaterialChoice::Ogden,
                        "Ogden",
                    );
                    ui.selectable_value(
                        &mut self.state.material,
                        MaterialChoice::HolzapfelOgden,
                        "Holzapfel-Ogden",
                    );
                });

            egui::ComboBox::from_label("Solver")
                .selected_text(match self.state.solver {
                    SolverChoice::Direct        => "Direct LU",
                    SolverChoice::Newton        => "Newton-Raphson",
                    SolverChoice::NewtonCachedK => "Newton (cached K)",
                    SolverChoice::NewtonSparse  => "Newton (sparse CG)",
                    SolverChoice::NewtonSparseParallel => "Newton (sparse CG parallel)",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.state.solver, SolverChoice::Direct,        "Direct LU");
                    ui.selectable_value(&mut self.state.solver, SolverChoice::Newton,        "Newton-Raphson");
                    ui.selectable_value(&mut self.state.solver, SolverChoice::NewtonCachedK, "Newton (cached K)");
                    ui.selectable_value(&mut self.state.solver, SolverChoice::NewtonSparse,  "Newton (sparse CG)");
                    ui.selectable_value(&mut self.state.solver, SolverChoice::NewtonSparseParallel, "Newton (sparse CG parallel)");
                });

            egui::ComboBox::from_label("Boundary")
                .selected_text(match self.state.boundary {
                    BoundaryChoice::Penalty     => "Penalty",
                    BoundaryChoice::Elimination => "Elimination",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Penalty,     "Penalty");
                    ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Elimination, "Elimination");
                });

            egui::ComboBox::from_label("Mode")
                .selected_text(match self.state.simulation_mode {
                    SimulationMode::Static  => "Static",
                    SimulationMode::Dynamic => "Dynamic",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.state.simulation_mode, SimulationMode::Static,  "Static");
                    ui.selectable_value(&mut self.state.simulation_mode, SimulationMode::Dynamic, "Dynamic");
                });
        });

        ui.group(|ui| {
            ui.vertical(|ui| {
                ui.label("Young's modulus (Pa)");
                ui.add(egui::DragValue::new(&mut self.state.youngs_modulus).speed(100.0));
                ui.label("Poisson ratio");
                ui.add(egui::DragValue::new(&mut self.state.poisson_ratio).speed(0.01));
                ui.label("Density (kg/m3)");
                ui.add(egui::DragValue::new(&mut self.state.density).speed(1.0));
                ui.label("Deformation scale");
                ui.add(egui::DragValue::new(&mut self.state.deformation_scale).speed(1.0));
            });
        });

        if self.state.simulation_mode == SimulationMode::Dynamic {
            ui.group(|ui| {
                ui.vertical(|ui| {
                    ui.label("dt (s)");
                    ui.add(egui::DragValue::new(&mut self.state.dt).speed(0.001));
                    ui.label("Rayleigh alpha");
                    ui.add(egui::DragValue::new(&mut self.state.alpha).speed(0.01));
                    ui.label("Rayleigh beta");
                    ui.add(egui::DragValue::new(&mut self.state.beta).speed(0.001));
                });
            });
        }

        ui.group(|ui| {
            ui.vertical(|ui| {
                ui.label("Number of Nodes");
                ui.label("x"); ui.add(egui::DragValue::new(&mut self.state.nx).speed(1));
                ui.label("y"); ui.add(egui::DragValue::new(&mut self.state.ny).speed(1));
                ui.label("z"); ui.add(egui::DragValue::new(&mut self.state.nz).speed(1));
            });
        });

        ui.group(|ui| {
            ui.vertical(|ui| {
                ui.label("Distance between Nodes");
                ui.label("x"); ui.add(egui::DragValue::new(&mut self.state.dx).speed(2));
                ui.label("y"); ui.add(egui::DragValue::new(&mut self.state.dy).speed(2));
                ui.label("z"); ui.add(egui::DragValue::new(&mut self.state.dz).speed(2));
            });
        });
    }

    fn ui_mesh_loading(&mut self, ui: &mut egui::Ui) {
        if ui.button("Browse...").clicked() {
            if let Some(path) = rfd::FileDialog::new()
                .add_filter("VTK", &["vtk"])
                .pick_file()
            {
                self.state.mesh_path = Some(path.to_string_lossy().to_string());
                self.state.mesh = match read_vtk(&path.to_string_lossy()) {
                    Ok(mesh) => {
                        println!("Loaded: {} nodes, {} elements", mesh.nodes.len(), mesh.elements.len());
                        // Auto-focus camera on the loaded mesh
                        self.state.camera.focus_on_mesh(&mesh.nodes);
                        self.state.displacements     = None;
                        self.state.initial_positions = None;
                        self.state.mechanical_state  = None;
                        self.state.bc_result_cached  = None;
                        self.state.n_dofs_total      = 0;
                        self.state.running           = false;
                        self.state.selected_node     = None;
                        self.state.constraint        = Constraint { list: Vec::new() };
                        self.state.loads             = Vec::new();
                        self.state.selected_load     = None;
                        Some(mesh)
                    }
                    Err(e) => { println!("Error: {:?}", e); None }
                };
            }
        }
    }

    fn ui_simulation_buttons(&mut self, ui: &mut egui::Ui) {
        match self.state.simulation_mode {
            SimulationMode::Static => {
                if ui.button("Simulate").clicked() {
                    run_simulation_static(&mut self.state);
                }
            }
            SimulationMode::Dynamic => {
                if ui.button("Initialize").clicked() {
                    init_simulation_dynamic(&mut self.state);
                }
                ui.horizontal(|ui| {
                    let label = if self.state.running { "Pause" } else { "Play" };
                    if ui.button(label).clicked() {
                        self.state.running = !self.state.running;
                    }
                    if ui.button("Reset").clicked() {
                        init_simulation_dynamic(&mut self.state);
                        self.state.running = false;
                    }
                });
            }
        }
    }

    fn ui_loads(&mut self, ui: &mut egui::Ui) {
        ui.separator();
        ui.label("Loads");

        if ui.button("Add Load").clicked() {
            self.state.loads.push(Load { list: Vec::new(), force: Vector3::zeros() });
            self.state.selected_load = Some(self.state.loads.len() - 1);
        }

        for (i, load) in self.state.loads.iter_mut().enumerate() {
            let selected = self.state.selected_load == Some(i);
            if ui.selectable_label(selected, format!("Load {}", i)).clicked() {
                self.state.selected_load = Some(i);
            }
            if selected {
                ui.horizontal(|ui| { ui.label("fx"); ui.add(egui::DragValue::new(&mut load.force.x).speed(1.0)); });
                ui.horizontal(|ui| { ui.label("fy"); ui.add(egui::DragValue::new(&mut load.force.y).speed(1.0)); });
                ui.horizontal(|ui| { ui.label("fz"); ui.add(egui::DragValue::new(&mut load.force.z).speed(1.0)); });
                ui.label(format!("{} nodes", load.list.len()));
            }
        }
    }

    fn ui_node_inspector(&mut self, ui: &mut egui::Ui) {
        let idx = match self.state.selected_node {
            Some(i) => i,
            None    => return,
        };
        if self.state.mesh.is_none() { return; }

        ui.separator();
        ui.label(format!("Node {}", idx));

        let pos = self.state.mesh.as_ref().unwrap().nodes[idx].position;
        ui.label(format!("x={:.3}  y={:.3}  z={:.3}", pos.x, pos.y, pos.z));

        if ui.button("Add to constraint").clicked() {
            if !self.state.constraint.list.iter().any(|n| n.indice == idx) {
                self.state.constraint.list.push(FixedNode::all(idx));
            }
        }
        if self.state.constraint.list.iter().any(|n| n.indice == idx) {
            ui.label("In constraint");
            if ui.button("Remove from constraint").clicked() {
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
                    ui.label(format!("In Load {}", load_idx));
                    if ui.button("Remove from load").clicked() {
                        self.state.loads[load_idx].list.retain(|&n| n != idx);
                    }
                }
            }
        }
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

        // ── Camera controls ────────────────────────────────────────────────
        if response.dragged() {
            let delta = response.drag_delta();
            let modifiers = ui.input(|i| i.modifiers);

            if modifiers.shift {
                // Shift + drag = pan
                self.state.camera.pan(delta, scale);
            } else {
                // Drag = orbit
                self.state.camera.yaw   += delta.x * 0.01;
                self.state.camera.pitch += delta.y * 0.01;
                // Clamp pitch to avoid gimbal flip
                self.state.camera.pitch = self.state.camera.pitch
                    .clamp(-std::f32::consts::FRAC_PI_2 + 0.05, std::f32::consts::FRAC_PI_2 - 0.05);
            }
        }

        // Scroll = zoom
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll != 0.0 {
            self.state.camera.distance *= (1.0 - scroll * 0.001).max(0.01);
            self.state.camera.distance  = self.state.camera.distance.max(0.01);
        }

        // ── Node picking ───────────────────────────────────────────────────
        if response.clicked() {
            if let Some(pos) = response.interact_pointer_pos() {
                if let Some(ref mesh) = self.state.mesh {
                    let i = screen_to_node(&pos, mesh, &self.state.camera, center, scale);
                    self.state.selected_node = Some(i);
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

        // ── Draw mesh ──────────────────────────────────────────────────────
        if let Some(ref displacements) = self.state.displacements {
            let scale_def = self.state.deformation_scale;
            let deformed_pos = |i: usize| -> Vector3<f64> {
                mesh.nodes[i].position + Vector3::new(
                    displacements[3*i]     * scale_def,
                    displacements[3*i + 1] * scale_def,
                    displacements[3*i + 2] * scale_def,
                )
            };

            // Edges
            for tetra in &mesh.elements {
                let pts: Vec<egui::Pos2> = tetra.indices.iter()
                    .map(|&i| self.state.camera.project(&deformed_pos(i), center, scale))
                    .collect();
                for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                    painter.line_segment(
                        [pts[a], pts[b]],
                        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY),
                    );
                }
            }

            // Displacement amplitude colormap
            let amplitudes: Vec<f64> = (0..mesh.nodes.len()).map(|i| {
                let ux = displacements[3*i];
                let uy = displacements[3*i + 1];
                let uz = displacements[3*i + 2];
                (ux*ux + uy*uy + uz*uz).sqrt()
            }).collect();
            let max_amp = amplitudes.iter().cloned().fold(0.0f64, f64::max);

            // Draw nodes back-to-front (depth sort)
            let sorted = depth_sorted_nodes(mesh, &self.state.camera);
            for i in sorted {
                let p = self.state.camera.project(&deformed_pos(i), center, scale);
                let t = if max_amp > 1e-15 { (amplitudes[i] / max_amp) as f32 } else { 0.0 };
                let color  = egui::Color32::from_rgb((t * 255.0) as u8, 0, ((1.0 - t) * 255.0) as u8);
                let radius = if self.state.selected_node == Some(i) { 6.0 } else { 3.0 };
                painter.circle_filled(p, radius, color);
            }

        } else {
            // Reference mesh (no simulation yet)
            for tetra in &mesh.elements {
                let pts: Vec<egui::Pos2> = tetra.indices.iter()
                    .map(|&i| self.state.camera.project(&mesh.nodes[i].position, center, scale))
                    .collect();
                for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                    painter.line_segment(
                        [pts[a], pts[b]],
                        egui::Stroke::new(1.0, egui::Color32::DARK_GRAY),
                    );
                }
            }

            // Draw nodes back-to-front
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
    }
}