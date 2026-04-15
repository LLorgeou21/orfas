use orfas_core::{mesh::Mesh,material::LinearElastic,
    boundary::BoundaryConditions,boundary::PenaltyMethod,boundary::EliminationMethod,boundary::FixedNode,boundary::BoundaryConditionResult,
    assembler::Assembler,assembler::LinearBMatrix,solver::DirectSolver,solver::Solver};
use nalgebra::{DVector, Vector3};
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
    material: MaterialChoice,
    solver: SolverChoice,
    boundary: BoundaryChoice,
    mesh: Option<Mesh>,
    displacements: Option<DVector<f64>>,
    node_forces: Vec<Vector3<f64>>,
    selected_node: Option<usize>,
    camera: Camera,
    deformation_scale: f64,
    mesh_path: Option<String>,
    mesh_center: Option<Vector3<f64>>,
}

struct MyEguiApp {
    state: AppState
}

impl MyEguiApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            state: AppState {
                nx: 2, ny: 2, nz: 2,
                dx: 1., dy: 1., dz: 1.,
                youngs_modulus: 1e6,
                poisson_ratio: 0.3,
                material: MaterialChoice::LinearElastic,
                solver: SolverChoice::DirectCholesky,
                boundary: BoundaryChoice::Penalty,
                mesh: None,
                displacements: None,
                node_forces: Vec::new(),
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
        egui::Panel::left("params").show_inside(ui, |ui| {

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
            });

            ui.group(|ui| {
                ui.vertical(|ui| {
                    ui.label("Young's modulus (Pa)");
                    ui.add(egui::DragValue::new(&mut self.state.youngs_modulus).speed(100.0));
                    ui.label("Poisson ratio");
                    ui.add(egui::DragValue::new(&mut self.state.poisson_ratio).speed(0.01));
                    ui.label("Deformation scale");
                    ui.add(egui::DragValue::new(&mut self.state.deformation_scale).speed(1.0));
                });
            });

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
                            self.state.selected_node = None;
                            self.state.node_forces = vec![Vector3::zeros(); mesh.nodes.len()];
                            Some(mesh)
                        }
                        Err(e) => { println!("Error: {:?}", e); None }
                    };
                }
            }

            if ui.button("Simulate").clicked() {
                run_simulation(&mut self.state);
            }

            // ── Node inspector ─────────────────────────────────────────────
            if let Some(idx) = self.state.selected_node {
                if let Some(ref mut mesh) = self.state.mesh {
                    ui.separator();
                    ui.label(format!("Node {}", idx));
                    let pos = mesh.nodes[idx].position;
                    ui.label(format!("x={:.3}  y={:.3}  z={:.3}", pos.x, pos.y, pos.z));

                    ui.checkbox(&mut mesh.nodes[idx].fixed, "Fixed (all directions)");

                    if idx < self.state.node_forces.len() {
                        ui.label("Applied force:");
                        ui.horizontal(|ui| {
                            ui.label("fx");
                            ui.add(egui::DragValue::new(&mut self.state.node_forces[idx].x).speed(1.0));
                        });
                        ui.horizontal(|ui| {
                            ui.label("fy");
                            ui.add(egui::DragValue::new(&mut self.state.node_forces[idx].y).speed(1.0));
                        });
                        ui.horizontal(|ui| {
                            ui.label("fz");
                            ui.add(egui::DragValue::new(&mut self.state.node_forces[idx].z).speed(1.0));
                        });
                    }
                }
            }
        });

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

            // Click selects the nearest node
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

                if let Some(displacements) = &self.state.displacements {
                    let scale_def = self.state.deformation_scale;

                    let deformed_pos = |i: usize| -> Vector3<f64> {
                        mesh.nodes[i].position + Vector3::new(
                            displacements[3*i]     * scale_def,
                            displacements[3*i + 1] * scale_def,
                            displacements[3*i + 2] * scale_def,
                        ) - center_mesh
                    };

                    // Draw edges using deformed positions
                    for tetra in &mesh.elements {
                        let pts: Vec<egui::Pos2> = tetra.indices.iter()
                            .map(|&i| project(&deformed_pos(i), &self.state.camera, center, scale))
                            .collect();
                        let edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];
                        for (a, b) in edges {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }

                    // Displacement colormap
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
                    // No simulation yet
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
                        if node.fixed {
                            painter.circle_filled(p, 5.0, egui::Color32::RED);
                        }
                        // Force indicator — yellow if a force is applied
                        if i < self.state.node_forces.len() && self.state.node_forces[i].norm() > 1e-10 {
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

fn run_simulation(state: &mut AppState) {
    let mut mesh = match state.mesh.take() {
        Some(mesh) => mesh,
        None => {
            let mesh = Mesh::generate(state.nx, state.ny, state.nz, state.dx, state.dy, state.dz);
            // Initialize node forces for generated mesh
            state.node_forces = vec![Vector3::zeros(); mesh.nodes.len()];
            mesh
        }
    };

    let material = LinearElastic { youngs_modulus: state.youngs_modulus, poisson_ratio: state.poisson_ratio };

    let center_mesh = mesh.nodes.iter().fold(
        Vector3::zeros(), |acc, n| acc + n.position
    ) / mesh.nodes.len() as f64;
    state.mesh_center = Some(center_mesh);

    // Build fixed nodes from mesh
    let fixed_nodes: Vec<FixedNode> = mesh.nodes.iter()
        .enumerate()
        .filter(|(_, node)| node.fixed)
        .map(|(i, _)| FixedNode::all(i))
        .collect();

    let bc = match state.boundary {
        BoundaryChoice::Penalty     => BoundaryConditions::new(fixed_nodes, Box::new(PenaltyMethod)),
        BoundaryChoice::Elimination => BoundaryConditions::new(fixed_nodes, Box::new(EliminationMethod)),
    };

    // Build force vector from node_forces
    let n = 3 * mesh.nodes.len();
    let mut f = DVector::zeros(n);
    for (i, force) in state.node_forces.iter().enumerate() {
        if i < mesh.nodes.len() {
            f[3 * i]     = force.x;
            f[3 * i + 1] = force.y;
            f[3 * i + 2] = force.z;
        }
    }

    let assembler = Assembler::new(&mesh);
    let k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
    let bc_result: BoundaryConditionResult = bc.apply(&k, &f);

    let solver = match state.solver {
        SolverChoice::DirectCholesky => DirectSolver {},
    };

    match solver.solve(&bc_result.k, &bc_result.f) {
        Ok(displacements) => {
            let u = bc_result.reconstruct(displacements);
            state.mesh = Some(mesh);
            state.displacements = Some(u);
        }
        Err(e) => {
            println!("Simulation error: {:?}", e);
        }
    }
}