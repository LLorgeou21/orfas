mod app;
mod render;
mod simulation;
mod state;

fn main() {
    let native_options = eframe::NativeOptions::default();
    let _ = eframe::run_native(
        "ORFAS",
        native_options,
        Box::new(|cc| Ok(Box::new(app::MyEguiApp::new(cc)))),
    );
}