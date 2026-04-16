use orfas_core::{
    mesh::Mesh,
    material::SaintVenantKirchhoff,
    assembler::{Assembler, LinearBMatrix},
    boundary::{BoundaryConditionResult, BoundaryConditions, Constraint, EliminationMethod, FixedNode, Load, PenaltyMethod},
    solver::{DirectSolver, NewtonRaphson, NonlinearSolver, Solver},
    damping::{DampingModel, RayleighDamping},
    integrator::{ImplicitEulerIntegrator, IntegratorMethod},
    mechanical_state::MechanicalState,
};
use nalgebra::{DMatrix, DVector, Vector3};
use orfas_io::read_vtk;

fn main() {
    let native_options = eframe::NativeOptions::default();
    let _ = eframe::run_native("ORFAS", native_options, Box::new(|cc| Ok(Box::new(MyEguiApp::new(cc)))));
}

// ─── Enums ────────────────────────────────────────────────────────────────────

#[derive(PartialEq)]
enum MaterialChoice {
    SaintVenantKirchhoff,
}

#[derive(PartialEq)]
enum SolverChoice {
    Direct,
    Newton,
}

#[derive(PartialEq)]
enum BoundaryChoice {
    Penalty,
    Elimination,
}

#[derive(PartialEq)]
enum SimulationMode {
    Static,
    Dynamic,
}

// ─── State ────────────────────────────────────────────────────────────────────

struct Camera {
    yaw:      f32,
    pitch:    f32,
    distance: f32,
}

struct AppState {
    // Mesh generation
    nx: usize, ny: usize, nz: usize,
    dx: f64,   dy: f64,   dz: f64,
    // Material
    youngs_modulus: f64,
    poisson_ratio:  f64,
    density:        f64,
    // Dynamic parameters
    alpha: f64,
    beta:  f64,
    dt:    f64,
    // Choices
    material:        MaterialChoice,
    solver:          SolverChoice,
    boundary:        BoundaryChoice,
    simulation_mode: SimulationMode,
    // Simulation data
    mesh:              Option<Mesh>,
    displacements:     Option<DVector<f64>>,
    initial_positions: Option<DVector<f64>>,
    mechanical_state:  Option<MechanicalState>,
    mass:              Option<DVector<f64>>,
    c_cached:          Option<DMatrix<f64>>,
    f_cached:          Option<DVector<f64>>,
    // Boundary condition result cached for the dynamic loop
    bc_result_cached:  Option<BoundaryConditionResult>,
    n_dofs_total:      usize,
    running:           bool,
    // Scene
    constraint:    Constraint,
    loads:         Vec<Load>,
    selected_load: Option<usize>,
    selected_node: Option<usize>,
    // View
    camera:           Camera,
    deformation_scale: f64,
    mesh_path:        Option<String>,
    mesh_center:      Option<Vector3<f64>>,
}

struct MyEguiApp {
    state: AppState,
}

impl MyEguiApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            state: AppState {
                nx: 2, ny: 2, nz: 2,
                dx: 1., dy: 1., dz: 1.,
                youngs_modulus: 1e6,
                poisson_ratio:  0.3,
                density:        1000.0,
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
                camera: Camera { yaw: 0.5, pitch: 0.3, distance: 5.0 },
                deformation_scale: 30.,
                mesh_path:    None,
                mesh_center:  None,
            }
        }
    }
}

// ─── App ──────────────────────────────────────────────────────────────────────

impl eframe::App for MyEguiApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {

        // ── Dynamic step — advance one step per frame if running ───────────
        if self.state.running {
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
                    let assembler = Assembler::new(mesh);
                    let integrator = ImplicitEulerIntegrator::default();
                    match integrator.step::<LinearBMatrix>(
                        mech, mass, c, f, self.state.dt,
                        &assembler, mesh, &material, bc_result, &DirectSolver,
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

        egui::Panel::left("params").show_inside(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {

            // ── Simulation parameters ──────────────────────────────────────
            ui.group(|ui| {
                egui::ComboBox::from_label("Material")
                    .selected_text(match self.state.material {
                        MaterialChoice::SaintVenantKirchhoff => "Saint Venant-Kirchhoff",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut self.state.material,
                            MaterialChoice::SaintVenantKirchhoff,
                            "Saint Venant-Kirchhoff",
                        );
                    });
                egui::ComboBox::from_label("Solver")
                    .selected_text(match self.state.solver {
                        SolverChoice::Direct => "Direct LU",
                        SolverChoice::Newton => "Newton-Raphson",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.state.solver, SolverChoice::Direct, "Direct LU");
                        ui.selectable_value(&mut self.state.solver, SolverChoice::Newton, "Newton-Raphson");
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

            // Dynamic parameters — only shown in dynamic mode
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

            // ── Mesh loading ───────────────────────────────────────────────
            if ui.button("Browse...").clicked() {
                if let Some(path) = rfd::FileDialog::new()
                    .add_filter("VTK", &["vtk"])
                    .pick_file()
                {
                    self.state.mesh_path = Some(path.to_string_lossy().to_string());
                    self.state.mesh = match read_vtk(&path.to_string_lossy()) {
                        Ok(mesh) => {
                            println!("Loaded: {} nodes, {} elements", mesh.nodes.len(), mesh.elements.len());
                            let center = mesh.nodes.iter().fold(Vector3::zeros(), |acc, n| acc + n.position)
                                / mesh.nodes.len() as f64;
                            self.state.mesh_center      = Some(center);
                            self.state.displacements    = None;
                            self.state.initial_positions = None;
                            self.state.mechanical_state = None;
                            self.state.bc_result_cached = None;
                            self.state.n_dofs_total     = 0;
                            self.state.running          = false;
                            self.state.selected_node    = None;
                            self.state.constraint       = Constraint { list: Vec::new() };
                            self.state.loads            = Vec::new();
                            self.state.selected_load    = None;
                            Some(mesh)
                        }
                        Err(e) => { println!("Error: {:?}", e); None }
                    };
                }
            }

            // ── Simulation buttons ─────────────────────────────────────────
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

            // ── Loads ──────────────────────────────────────────────────────
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

            // ── Node inspector ─────────────────────────────────────────────
            if let Some(idx) = self.state.selected_node {
                if let Some(ref mesh) = self.state.mesh {
                    ui.separator();
                    ui.label(format!("Node {}", idx));
                    let pos = mesh.nodes[idx].position;
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
        });});

        // ── 3D viewport ────────────────────────────────────────────────────
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let (response, painter) = ui.allocate_painter(
                ui.available_size(),
                egui::Sense::drag() | egui::Sense::click(),
            );
            let rect   = response.rect;
            let center = rect.center();
            let scale  = 100.0f32;

            if response.dragged() {
                let delta = response.drag_delta();
                self.state.camera.yaw   += delta.x * 0.01;
                self.state.camera.pitch += delta.y * 0.01;
            }
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            self.state.camera.distance -= scroll * 0.01;
            self.state.camera.distance  = self.state.camera.distance.max(0.1);

            if response.clicked() {
                if let Some(pos) = response.interact_pointer_pos() {
                    if let Some(ref mesh) = self.state.mesh {
                        let mesh_center = self.state.mesh_center.unwrap();
                        let i = screen_to_node(&pos, mesh, &self.state.camera, center, scale, &mesh_center);
                        self.state.selected_node = Some(i);
                    }
                }
            }

            if let Some(mesh) = &self.state.mesh {
                let center_mesh = self.state.mesh_center.unwrap();

                let fixed_indices: std::collections::HashSet<usize> =
                    self.state.constraint.list.iter().map(|n| n.indice).collect();
                let loaded_indices: std::collections::HashSet<usize> =
                    self.state.loads.iter().flat_map(|l| l.list.iter().copied()).collect();

                if let Some(displacements) = &self.state.displacements {
                    let scale_def = self.state.deformation_scale;
                    let deformed_pos = |i: usize| -> Vector3<f64> {
                        mesh.nodes[i].position + Vector3::new(
                            displacements[3*i]     * scale_def,
                            displacements[3*i + 1] * scale_def,
                            displacements[3*i + 2] * scale_def,
                        ) - center_mesh
                    };

                    for tetra in &mesh.elements {
                        let pts: Vec<egui::Pos2> = tetra.indices.iter()
                            .map(|&i| project(&deformed_pos(i), &self.state.camera, center, scale))
                            .collect();
                        for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }

                    let amplitudes: Vec<f64> = (0..mesh.nodes.len()).map(|i| {
                        let ux = displacements[3*i];
                        let uy = displacements[3*i + 1];
                        let uz = displacements[3*i + 2];
                        (ux*ux + uy*uy + uz*uz).sqrt()
                    }).collect();
                    let max_amp = amplitudes.iter().cloned().fold(0.0f64, f64::max);

                    for i in 0..mesh.nodes.len() {
                        let p = project(&deformed_pos(i), &self.state.camera, center, scale);
                        let t = if max_amp > 1e-15 { (amplitudes[i] / max_amp) as f32 } else { 0.0 };
                        let color  = egui::Color32::from_rgb((t * 255.0) as u8, 0, ((1.0 - t) * 255.0) as u8);
                        let radius = if self.state.selected_node == Some(i) { 6.0 } else { 3.0 };
                        painter.circle_filled(p, radius, color);
                    }

                } else {
                    for tetra in &mesh.elements {
                        let pts: Vec<egui::Pos2> = tetra.indices.iter()
                            .map(|&i| project(&(mesh.nodes[i].position - center_mesh), &self.state.camera, center, scale))
                            .collect();
                        for (a, b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }
                    for (i, node) in mesh.nodes.iter().enumerate() {
                        let pos = node.position - center_mesh;
                        let p   = project(&pos, &self.state.camera, center, scale);
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
        });
    }
}

// ─── Helpers ──────────────────────────────────────────────────────────────────

fn make_material(state: &AppState) -> SaintVenantKirchhoff {
    SaintVenantKirchhoff {
        youngs_modulus: state.youngs_modulus,
        poisson_ratio:  state.poisson_ratio,
        density:        state.density,
    }
}

fn project(point: &Vector3<f64>, camera: &Camera, center: egui::Pos2, scale: f32) -> egui::Pos2 {
    let cos_yaw   = camera.yaw.cos();
    let sin_yaw   = camera.yaw.sin();
    let rx = point.x as f32 * cos_yaw - point.z as f32 * sin_yaw;
    let ry = point.y as f32;
    let rz = point.x as f32 * sin_yaw + point.z as f32 * cos_yaw;

    let cos_pitch = camera.pitch.cos();
    let sin_pitch = camera.pitch.sin();
    let ry2 = ry * cos_pitch - rz * sin_pitch;
    let rz2 = ry * sin_pitch + rz * cos_pitch;

    let d  = camera.distance + rz2;
    let sx = center.x + rx / d * scale;
    let sy = center.y - ry2 / d * scale;
    egui::Pos2::new(sx, sy)
}

fn screen_to_node(
    pos: &egui::Pos2, mesh: &Mesh, camera: &Camera,
    center: egui::Pos2, scale: f32, mesh_center: &Vector3<f64>,
) -> usize {
    mesh.nodes.iter().enumerate()
        .map(|(i, node)| {
            let p  = project(&(node.position - mesh_center), camera, center, scale);
            let dx = p.x - pos.x;
            let dy = p.y - pos.y;
            (i, dx * dx + dy * dy)
        })
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0)
}

// ─── Simulation ───────────────────────────────────────────────────────────────

/// Builds mesh, material, assembler, K, mass, C, bc_result.
/// Retourne (mesh, mass_red, c_red, bc_result).
fn build_simulation(state: &mut AppState)
    -> Option<(Mesh, DVector<f64>, DMatrix<f64>, BoundaryConditionResult)>
{
    let mesh = match state.mesh.take() {
        Some(m) => m,
        None => {
            let m = Mesh::generate(state.nx, state.ny, state.nz, state.dx, state.dy, state.dz);
            state.constraint  = Constraint { list: Vec::new() };
            state.loads       = Vec::new();
            state.selected_load = None;
            m
        }
    };

    let material = make_material(state);
    let center = mesh.nodes.iter().fold(Vector3::zeros(), |acc, n| acc + n.position)
        / mesh.nodes.len() as f64;
    state.mesh_center = Some(center);

    let assembler = Assembler::new(&mesh);
    let u_zero    = DVector::zeros(3 * mesh.nodes.len());
    let k_full    = assembler.assemble_tangent::<LinearBMatrix>(&mesh, &material, &u_zero);
    let mass      = assembler.assemble_mass(&mesh, &material);
    let c_full    = RayleighDamping { alpha: state.alpha, beta: state.beta }.compute(&mass, &k_full);

    let bc = match state.boundary {
        BoundaryChoice::Penalty     => BoundaryConditions::new(state.constraint.clone(), state.loads.clone(), Box::new(PenaltyMethod)),
        BoundaryChoice::Elimination => BoundaryConditions::new(state.constraint.clone(), state.loads.clone(), Box::new(EliminationMethod)),
    };
    let bc_result = bc.apply(&k_full, mesh.nodes.len());

    // Reduce mass and c to the free DOF space
    let (mass_red, c_red) = match &bc_result.free_dofs {
        Some(free) => (
            DVector::from_iterator(free.len(), free.iter().map(|&i| mass[i])),
            DMatrix::from_fn(free.len(), free.len(), |r, s| c_full[(free[r], free[s])]),
        ),
        None => (mass, c_full),
    };

    Some((mesh, mass_red, c_red, bc_result))
}

fn run_simulation_static(state: &mut AppState) {
    if let Some((mesh, _mass, _c, bc_result)) = build_simulation(state) {
        let material  = make_material(state);
        let assembler = Assembler::new(&mesh);

        let result = match state.solver {
            SolverChoice::Direct => {
                DirectSolver.solve(&bc_result.k, &bc_result.f)
                    .map(|u_red| bc_result.reconstruct(u_red))
            }
            SolverChoice::Newton => {
                NewtonRaphson::default()
                    .solve::<LinearBMatrix>(&assembler, &mesh, &material, &bc_result, &DirectSolver)
                    .map(|u_red| bc_result.reconstruct(u_red))
            }
        };

        match result {
            Ok(u) => {
                state.displacements = Some(u);
                state.mesh = Some(mesh);
            }
            Err(e) => println!("Simulation error: {:?}", e),
        }
    }
}

fn init_simulation_dynamic(state: &mut AppState) {
    if let Some((mesh, mass_red, c_red, bc_result)) = build_simulation(state) {
        let n_red        = mass_red.len();
        let n_dofs_total = 3 * mesh.nodes.len();
        let mech         = MechanicalState::new(n_red / 3);
        let init_pos     = mech.position.clone();
        let f_ext        = bc_result.f.clone();

        state.initial_positions = Some(init_pos);
        state.n_dofs_total      = n_dofs_total;
        state.mechanical_state  = Some(mech);
        state.mass              = Some(mass_red);
        state.c_cached          = Some(c_red);
        state.f_cached          = Some(f_ext);
        state.bc_result_cached  = Some(bc_result);
        state.displacements     = None;
        state.mesh              = Some(mesh);
    }
}