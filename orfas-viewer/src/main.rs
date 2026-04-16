use orfas_core::{mesh::Mesh, material::LinearElastic,
    boundary::{BoundaryConditions, PenaltyMethod, EliminationMethod, FixedNode, Constraint, Load},
    assembler::{Assembler, LinearBMatrix},
    solver::{DirectSolver, Solver},
    damping::{RayleighDamping, DampingModel},
    integrator::{ImplicitEulerIntegrator, IntegratorMethod},
    mechanical_state::MechanicalState};
use nalgebra::Vector3;
use orfas_io::read_vtk;


fn main() {
    let native_options = eframe::NativeOptions::default();
    let _ = eframe::run_native("MyApp", native_options, Box::new(|cc| Ok(Box::new(MyEguiApp::new(cc)))));
}


#[derive(PartialEq)]
enum MaterialChoice {
    LinearElastic,
}

#[derive(PartialEq)]
enum SolverChoice {
    DirectCholesky,
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

struct Camera {
    yaw: f32,
    pitch: f32,
    distance: f32,
}

struct AppState {
    nx: usize, ny: usize, nz: usize,
    dx: f64, dy: f64, dz: f64,
    youngs_modulus: f64,
    poisson_ratio: f64,
    density: f64,
    alpha: f64,
    beta: f64,
    dt: f64,
    material: MaterialChoice,
    solver: SolverChoice,
    boundary: BoundaryChoice,
    simulation_mode: SimulationMode,
    mesh: Option<Mesh>,
    displacements: Option<nalgebra::DVector<f64>>,
    initial_positions: Option<nalgebra::DVector<f64>>,
    mechanical_state: Option<MechanicalState>,
    mass: Option<nalgebra::DVector<f64>>,
    k_cached: Option<nalgebra::DMatrix<f64>>,
    c_cached: Option<nalgebra::DMatrix<f64>>,
    f_cached: Option<nalgebra::DVector<f64>>,
    free_dofs: Option<Vec<usize>>,
    n_dofs_total: usize,
    running: bool,
    constraint: Constraint,
    loads: Vec<Load>,
    selected_load: Option<usize>,
    selected_node: Option<usize>,
    camera: Camera,
    deformation_scale: f64,
    mesh_path: Option<String>,
    mesh_center: Option<Vector3<f64>>,
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
                poisson_ratio: 0.3,
                density: 1000.0,
                alpha: 0.1,
                beta: 0.01,
                dt: 0.01,
                material: MaterialChoice::LinearElastic,
                solver: SolverChoice::DirectCholesky,
                boundary: BoundaryChoice::Penalty,
                simulation_mode: SimulationMode::Static,
                mesh: None,
                displacements: None,
                initial_positions: None,
                mechanical_state: None,
                mass: None,
                k_cached: None,
                c_cached: None,
                f_cached: None,
                free_dofs: None,
                n_dofs_total: 0,
                running: false,
                constraint: Constraint { list: Vec::new() },
                loads: Vec::new(),
                selected_load: None,
                selected_node: None,
                camera: Camera { yaw: 0.5, pitch: 0.3, distance: 5.0 },
                deformation_scale: 30.,
                mesh_path: None,
                mesh_center: None,
            }
        }
    }
}

impl eframe::App for MyEguiApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {

        // Dynamic step — advance one step per frame if running
        if self.state.running {
            if let (Some(ref mut mech), Some(ref mass), Some(ref k), Some(ref c), Some(ref f), Some(ref init)) = (
                self.state.mechanical_state.as_mut(),
                self.state.mass.as_ref(),
                self.state.k_cached.as_ref(),
                self.state.c_cached.as_ref(),
                self.state.f_cached.as_ref(),
                self.state.initial_positions.as_ref(),
            ) {
                let integrator = ImplicitEulerIntegrator;
                let solver = DirectSolver {};
                match integrator.step(mech, mass, k, c, f, self.state.dt, &solver) {
                    Ok(_) => {
                        // Compute relative displacement from initial positions
                        let disp_reduced = &mech.position - *init;

                        // Reconstruct full displacement vector for rendering
                        let mut disp_full = nalgebra::DVector::zeros(self.state.n_dofs_total);
                        match &self.state.free_dofs {
                            Some(dofs) => {
                                for (j, &dof) in dofs.iter().enumerate() {
                                    disp_full[dof] = disp_reduced[j];
                                }
                            }
                            None => disp_full = disp_reduced,
                        }
                        self.state.displacements = Some(disp_full);
                    }
                    Err(e) => {
                        println!("Integration error: {:?}", e);
                        self.state.running = false;
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
                        MaterialChoice::LinearElastic => "Linear Elastic",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.state.material, MaterialChoice::LinearElastic, "Linear Elastic");
                    });
                egui::ComboBox::from_label("Solver")
                    .selected_text(match self.state.solver {
                        SolverChoice::DirectCholesky => "Direct Cholesky",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.state.solver, SolverChoice::DirectCholesky, "Direct Cholesky");
                    });
                egui::ComboBox::from_label("Boundary")
                    .selected_text(match self.state.boundary {
                        BoundaryChoice::Penalty => "Penalty",
                        BoundaryChoice::Elimination => "Elimination",
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Penalty, "Penalty");
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
                    .pick_file() {
                    self.state.mesh_path = Some(path.to_string_lossy().to_string());
                    let path_str = path.to_string_lossy();
                    self.state.mesh = match read_vtk(&path_str) {
                        Ok(mesh) => {
                            println!("Loaded: {} nodes, {} elements", mesh.nodes.len(), mesh.elements.len());
                            let center = mesh.nodes.iter().fold(
                                Vector3::zeros(), |acc, n| acc + n.position
                            ) / mesh.nodes.len() as f64;
                            self.state.mesh_center = Some(center);
                            self.state.displacements = None;
                            self.state.initial_positions = None;
                            self.state.mechanical_state = None;
                            self.state.free_dofs = None;
                            self.state.n_dofs_total = 0;
                            self.state.running = false;
                            self.state.selected_node = None;
                            self.state.constraint = Constraint { list: Vec::new() };
                            self.state.loads = Vec::new();
                            self.state.selected_load = None;
                            Some(mesh)
                        }
                        Err(e) => { println!("Error: {:?}", e); None }
                    };
                }
            }

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
                    ui.horizontal(|ui| {
                        ui.label("fx"); ui.add(egui::DragValue::new(&mut load.force.x).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("fy"); ui.add(egui::DragValue::new(&mut load.force.y).speed(1.0));
                    });
                    ui.horizontal(|ui| {
                        ui.label("fz"); ui.add(egui::DragValue::new(&mut load.force.z).speed(1.0));
                    });
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
                        let already = self.state.constraint.list.iter().any(|n| n.indice == idx);
                        if !already {
                            self.state.constraint.list.push(FixedNode::all(idx));
                        }
                    }
                    let in_constraint = self.state.constraint.list.iter().any(|n| n.indice == idx);
                    if in_constraint {
                        ui.label("In constraint");
                        if ui.button("Remove from constraint").clicked() {
                            self.state.constraint.list.retain(|n| n.indice != idx);
                        }
                    }

                    if let Some(load_idx) = self.state.selected_load {
                        if load_idx < self.state.loads.len() {
                            if ui.button(format!("Add to Load {}", load_idx)).clicked() {
                                let already = self.state.loads[load_idx].list.contains(&idx);
                                if !already {
                                    self.state.loads[load_idx].list.push(idx);
                                }
                            }
                            let in_load = self.state.loads[load_idx].list.contains(&idx);
                            if in_load {
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

        egui::CentralPanel::default().show_inside(ui, |ui| {
            let (response, painter) = ui.allocate_painter(
                ui.available_size(),
                egui::Sense::drag() | egui::Sense::click()
            );
            let rect = response.rect;
            let center = rect.center();
            let scale = 100.0f32;

            if response.dragged() {
                let delta = response.drag_delta();
                self.state.camera.yaw   += delta.x * 0.01;
                self.state.camera.pitch += delta.y * 0.01;
            }
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            self.state.camera.distance -= scroll * 0.01;
            self.state.camera.distance = self.state.camera.distance.max(0.1);

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
                        let edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];
                        for (a, b) in edges {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }

                    let amplitudes: Vec<f64> = (0..mesh.nodes.len()).map(|i| {
                        let ux = displacements[3*i];
                        let uy = displacements[3*i + 1];
                        let uz = displacements[3*i + 2];
                        (ux*ux + uy*uy + uz*uz).sqrt()
                    }).collect();
                    let max_amplitude = amplitudes.iter().cloned().fold(0.0f64, f64::max);

                    for i in 0..mesh.nodes.len() {
                        let p = project(&deformed_pos(i), &self.state.camera, center, scale);
                        let t = if max_amplitude > 1e-15 { (amplitudes[i] / max_amplitude) as f32 } else { 0.0 };
                        let color = egui::Color32::from_rgb((t * 255.0) as u8, 0, ((1.0 - t) * 255.0) as u8);
                        let radius = if self.state.selected_node == Some(i) { 6.0 } else { 3.0 };
                        painter.circle_filled(p, radius, color);
                    }

                } else {
                    for tetra in &mesh.elements {
                        let pts: Vec<egui::Pos2> = tetra.indices.iter()
                            .map(|&i| project(&(mesh.nodes[i].position - center_mesh), &self.state.camera, center, scale))
                            .collect();
                        let edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];
                        for (a, b) in edges {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }
                    for (i, node) in mesh.nodes.iter().enumerate() {
                        let pos = node.position - center_mesh;
                        let p = project(&pos, &self.state.camera, center, scale);
                        if fixed_indices.contains(&i) {
                            painter.circle_filled(p, 5.0, egui::Color32::RED);
                        }
                        if loaded_indices.contains(&i) {
                            painter.circle_filled(p, 5.0, egui::Color32::YELLOW);
                        }
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


fn project(point: &Vector3<f64>, camera: &Camera, center: egui::Pos2, scale: f32) -> egui::Pos2 {
    let cos_yaw = camera.yaw.cos();
    let sin_yaw = camera.yaw.sin();
    let rx = point.x as f32 * cos_yaw - point.z as f32 * sin_yaw;
    let ry = point.y as f32;
    let rz = point.x as f32 * sin_yaw + point.z as f32 * cos_yaw;

    let cos_pitch = camera.pitch.cos();
    let sin_pitch = camera.pitch.sin();
    let ry2 = ry * cos_pitch - rz * sin_pitch;
    let rz2 = ry * sin_pitch + rz * cos_pitch;

    let d = camera.distance + rz2;
    let sx = center.x + rx / d * scale;
    let sy = center.y - ry2 / d * scale;
    egui::Pos2::new(sx, sy)
}

fn screen_to_node(pos: &egui::Pos2, mesh: &Mesh, camera: &Camera, center: egui::Pos2, scale: f32, mesh_center: &Vector3<f64>) -> usize {
    mesh.nodes.iter()
        .enumerate()
        .map(|(i, node)| {
            let pos_centered = node.position - mesh_center;
            let projected = project(&pos_centered, camera, center, scale);
            let dx = projected.x - pos.x;
            let dy = projected.y - pos.y;
            (i, dx * dx + dy * dy)
        })
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0)
}

fn build_mesh_and_matrices(state: &mut AppState) -> Option<(Mesh, nalgebra::DVector<f64>, nalgebra::DMatrix<f64>, nalgebra::DMatrix<f64>, nalgebra::DVector<f64>, Option<Vec<usize>>)> {
    let mesh = match state.mesh.take() {
        Some(mesh) => mesh,
        None => {
            let mesh = Mesh::generate(state.nx, state.ny, state.nz, state.dx, state.dy, state.dz);
            state.constraint = Constraint { list: Vec::new() };
            state.loads = Vec::new();
            state.selected_load = None;
            mesh
        }
    };

    let material = LinearElastic {
        youngs_modulus: state.youngs_modulus,
        poisson_ratio: state.poisson_ratio,
        density: state.density,
    };

    let center_mesh = mesh.nodes.iter().fold(
        Vector3::zeros(), |acc, n| acc + n.position
    ) / mesh.nodes.len() as f64;
    state.mesh_center = Some(center_mesh);

    let assembler = Assembler::new(&mesh);
    let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
    let mass = assembler.assemble_mass(&mesh, &material);

    let damping = RayleighDamping { alpha: state.alpha, beta: state.beta };
    let c = damping.compute(&mass, &k);

    let bc = match state.boundary {
        BoundaryChoice::Penalty     => BoundaryConditions::new(state.constraint.clone(), state.loads.clone(), Box::new(PenaltyMethod)),
        BoundaryChoice::Elimination => BoundaryConditions::new(state.constraint.clone(), state.loads.clone(), Box::new(EliminationMethod)),
    };
    let bc_result = bc.apply(&k, mesh.nodes.len());

    // Reduce mass and c using the same free_dofs as k when using elimination method
    let free_dofs = bc_result.free_dofs.clone();
    let (mass_reduced, c_reduced) = match &bc_result.free_dofs {
        Some(free_dofs) => {
            let m = nalgebra::DVector::from_iterator(
                free_dofs.len(),
                free_dofs.iter().map(|&i| mass[i])
            );
            let c = nalgebra::DMatrix::from_fn(free_dofs.len(), free_dofs.len(), |r, s| {
                c[(free_dofs[r], free_dofs[s])]
            });
            (m, c)
        }
        None => (mass, c),
    };

    Some((mesh, mass_reduced, bc_result.k, c_reduced, bc_result.f, free_dofs))
}

fn run_simulation_static(state: &mut AppState) {
    if let Some((mesh, _mass, k, _c, f, _free_dofs)) = build_mesh_and_matrices(state) {
        let solver = DirectSolver {};
        match solver.solve(&k, &f) {
            Ok(displacements) => {
                state.displacements = Some(displacements);
                state.mesh = Some(mesh);
            }
            Err(e) => println!("Simulation error: {:?}", e),
        }
    }
}

fn init_simulation_dynamic(state: &mut AppState) {
    if let Some((mesh, mass, k, c, f, free_dofs)) = build_mesh_and_matrices(state) {
        let n_dofs_total = 3 * mesh.nodes.len();
        let n = mass.len();
        let mut mech = MechanicalState::new(n / 3);



        let initial_positions = mech.position.clone();

        state.initial_positions = Some(initial_positions);
        state.free_dofs = free_dofs;
        state.n_dofs_total = n_dofs_total;
        state.mechanical_state = Some(mech);
        state.mass     = Some(mass);
        state.k_cached = Some(k);
        state.c_cached = Some(c);
        state.f_cached = Some(f);
        state.displacements = None;
        state.mesh = Some(mesh);
    }
}