use orfas_core::{mesh::Mesh,material::LinearElastic,
    boundary::BoundaryConditions,boundary::PenaltyMethod,
    assembler::Assembler,assembler::LinearBMatrix,solver::DirectSolver,solver::Solver};
use nalgebra::DVector;
use nalgebra::Vector3;



fn main() {





    let native_options = eframe::NativeOptions::default();
    let _ = eframe::run_native("MyApp", native_options, Box::new(|cc| Ok(Box::new(MyEguiApp::new(cc)))));
}





#[derive(PartialEq)]
enum MaterialChoice {
    LinearElastic,
    // NeoHookean,  // v0.4
}

#[derive(PartialEq)]
enum SolverChoice {
    DirectCholesky,
    // ConjugateGradient,  // v0.8
}

#[derive(PartialEq)]
enum BoundaryChoice {
    Penalty,
    // Elimination,  // v0.2
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

    camera: Camera,
    mesh_center: Option<nalgebra::Vector3<f64>>,
}


struct MyEguiApp{
    state : AppState
}


impl MyEguiApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self {
            state : AppState { 
                nx:2, 
                ny:2, 
                nz:2, 
                dx:1., 
                dy:1., 
                dz:1., 
                youngs_modulus:1e6, 
                poisson_ratio:0.3, 
                material:MaterialChoice::LinearElastic, 
                solver:SolverChoice::DirectCholesky, 
                boundary:BoundaryChoice::Penalty, 
                mesh: None, 
                displacements: None, 
                camera: Camera { yaw: 0.5, pitch: 0.3, distance: 5.0 },
                mesh_center : None }
        }
    }
}

impl eframe::App for MyEguiApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        egui::Panel::left("params").show_inside(ui, |ui| {
            ui.group(|ui| {
                egui::ComboBox::from_label("Material")
                    .selected_text(match self.state.material {
                        MaterialChoice::LinearElastic => "Linear Elastic",
                        // MaterialChoice::NeoHookean => "Neo-Hookean",
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
                    })
                    .show_ui(ui, |ui| {
                        ui.selectable_value(&mut self.state.boundary, BoundaryChoice::Penalty, "Penalty");
                    });
            });
            ui.group(|ui| {
                ui.vertical(|ui| {
                    ui.label("Young's modulus (Pa)");
                    ui.add(egui::DragValue::new(&mut self.state.youngs_modulus).speed(100.0));
                    ui.label("Poisson ratio");
                    ui.add(egui::DragValue::new(&mut self.state.poisson_ratio).speed(0.01));
                });
            });
            ui.group(|ui| {
                ui.vertical(|ui| {
                    ui.label("Number of Nodes");
                    ui.label("x dimension");
                    ui.add(egui::DragValue::new(&mut self.state.nx).speed(1));
                    ui.label("y dimension");
                    ui.add(egui::DragValue::new(&mut self.state.ny).speed(1));
                    ui.label("z dimension");
                    ui.add(egui::DragValue::new(&mut self.state.nz).speed(1));
                });
            });
            ui.group(|ui| {
                ui.vertical(|ui| {
                    ui.label("Distance between Nodes");
                    ui.label("x dimension");
                    ui.add(egui::DragValue::new(&mut self.state.dx).speed(2));
                    ui.label("y dimension");
                    ui.add(egui::DragValue::new(&mut self.state.dy).speed(2));
                    ui.label("z dimension");
                    ui.add(egui::DragValue::new(&mut self.state.dz).speed(2));
                });
            });
            if ui.button("Simulate").clicked() {
                run_simulation(&mut self.state);
            }
        });
        egui::CentralPanel::default().show_inside(ui, |ui| {
            let (response, painter) = ui.allocate_painter(
                ui.available_size(),
                egui::Sense::drag()
            );
            let rect = response.rect;
            let center = rect.center();
            let scale = 100.0f32;
            if response.dragged() {
                let delta = response.drag_delta();
                self.state.camera.yaw += delta.x * 0.01;
                self.state.camera.pitch += delta.y * 0.01;
            }

            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            self.state.camera.distance -= scroll * 0.1;
            self.state.camera.distance = self.state.camera.distance.max(10.);

            if let Some(mesh) = &self.state.mesh {
                if let Some(displacements) = &self.state.displacements {
                    let center_mesh = self.state.mesh_center.unwrap();

                    for tetra in &mesh.elements {
                        let pts: Vec<egui::Pos2> = tetra.indices.iter()
                            .map(|&i| {
                                let pos = mesh.nodes[i].position - center_mesh;
                                project(&pos, &self.state.camera, center, scale)
                            })
                            .collect();
                        
                        let edges = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];
                        for (a, b) in edges {
                            painter.line_segment([pts[a], pts[b]], egui::Stroke::new(1.0, egui::Color32::DARK_GRAY));
                        }
                    }
                    let mut amplitudes = Vec::new();
                    let mut max_amplitude : f64 = 0.;
                    for (i,_) in mesh.nodes.iter().enumerate(){
                        let ux = displacements[3*i];
                        let uy = displacements[3*i + 1];
                        let uz = displacements[3*i + 2];
                        let amplitude = (ux*ux + uy*uy + uz*uz).sqrt();
                        amplitudes.push(amplitude);
                        if amplitude>max_amplitude{
                            max_amplitude = amplitude
                        }
                    }
                    for (i,node) in mesh.nodes.iter().enumerate() {
                        
                        let pos = node.position - center_mesh;
                        let p = project(&pos, &self.state.camera, center, scale);
                        let t = (amplitudes[i] / max_amplitude) as f32;
                        let color = egui::Color32::from_rgb(
                            (t * 255.0) as u8,         // rouge — max displacement
                            0,
                            ((1.0 - t) * 255.0) as u8  // bleu — no displacement
                        );
                        painter.circle_filled(p, 3.0, color);
                    }
                }
            }
        });
    }
}


fn project(point: &Vector3<f64>, camera: &Camera, center: egui::Pos2, scale: f32) -> egui::Pos2 {
    // Rotation yaw (autour de Y)
    let cos_yaw = camera.yaw.cos();
    let sin_yaw = camera.yaw.sin();
    let rx = point.x as f32 * cos_yaw - point.z as f32 * sin_yaw;
    let ry = point.y as f32;
    let rz = point.x as f32 * sin_yaw + point.z as f32 * cos_yaw;

    // Rotation pitch (autour de X)
    let cos_pitch = camera.pitch.cos();
    let sin_pitch = camera.pitch.sin();
    let ry2 = ry * cos_pitch - rz * sin_pitch;
    let rz2 = ry * sin_pitch + rz * cos_pitch;

    // Projection perspective
    let d = camera.distance + rz2;
    let sx = center.x + rx / d * scale;
    let sy = center.y - ry2 / d * scale;
    egui::Pos2::new(sx, sy)
}

fn run_simulation(state: &mut AppState) {
    let mut mesh = Mesh::generate(state.nx, state.ny, state.nz, state.dx, state.dy, state.dz);
    let material = LinearElastic { youngs_modulus: state.youngs_modulus, poisson_ratio: state.poisson_ratio };
    let center_mesh = mesh.nodes.iter().fold(
        nalgebra::Vector3::zeros(),
        |acc, n| acc + n.position
    ) / mesh.nodes.len() as f64;
    state.mesh_center = Some(center_mesh);

    let mut fixed_nodes = Vec::new();
    for i in 0..state.nx {
        for j in 0..state.ny {
            let idx = i + j * state.nx;
            mesh.nodes[idx].fixed = true;
            fixed_nodes.push(idx);
        }
    }
    
    let bc = match state.boundary {
        BoundaryChoice::Penalty => BoundaryConditions::new(fixed_nodes, Box::new(PenaltyMethod)),
    };
    
    let n = 3 * mesh.nodes.len();
    let mut f = DVector::zeros(n);
    

    let force = -100.0;
    for i in 0..state.nx {
        for j in 0..state.ny {
            let idx = i + j * state.nx + (state.nz - 1) * state.nx * state.ny;
            f[3 * idx + 2] = force;
        }
    }
    
    let assembler = Assembler::new(&mesh);
    let mut k = assembler.assemble::<LinearBMatrix>(&mesh, &material);
    bc.apply(&mut k, &mut f);
    
    let solver = match state.solver {
        SolverChoice::DirectCholesky => DirectSolver {},
    };
    
    match solver.solve(&k, &f) {
        Ok(displacements) => {
            state.mesh = Some(mesh);
            state.displacements = Some(displacements);
        }
        Err(e) => {
            println!("Simulation error: {:?}", e);
        }
    }
}









