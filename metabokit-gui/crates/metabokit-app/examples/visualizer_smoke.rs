//! Exercise the visualizer data service without launching a window.
//!
//! ```text
//! cargo run -p metabokit-app --example visualizer_smoke -- <results> <sample> <mz> <rt>
//! ```

use metabokit_app_lib::visualizer::{
    visualizer_feature_data, visualizer_open, visualizer_overview, VisualizerState,
};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        eprintln!("usage: visualizer_smoke <results folder> <sample> <mz> <rt>");
        std::process::exit(2);
    }
    let output = args[0].clone();
    let sample = args[1].clone();
    let mz: f32 = args[2].parse().expect("mz must be numeric");
    let rt: f32 = args[3].parse().expect("rt must be numeric");

    let started = std::time::Instant::now();
    let session = visualizer_open(output.clone()).expect("open visualizer session");
    let session = serde_json::to_value(session).unwrap();
    println!(
        "session: {} samples, matches={}",
        session["samples"].as_array().map(Vec::len).unwrap_or(0),
        session["matchesAvailable"]
    );
    println!(
        "session loaded in {:.1} ms",
        started.elapsed().as_secs_f64() * 1000.0
    );
    let started = std::time::Instant::now();
    let overview = visualizer_overview(output.clone(), sample.clone()).expect("load overview");
    let overview = serde_json::to_value(overview).unwrap();
    println!(
        "overview: {} plotted / {} total ({})",
        overview["mz"].as_array().map(Vec::len).unwrap_or(0),
        overview["total"],
        overview["source"]
    );
    println!(
        "overview loaded in {:.1} ms",
        started.elapsed().as_secs_f64() * 1000.0
    );
    let state = VisualizerState::default();
    let started = std::time::Instant::now();
    let detail = visualizer_feature_data(output.clone(), sample.clone(), mz, rt, 0.1, &state)
        .expect("load feature detail");
    let detail = serde_json::to_value(detail).unwrap();
    println!(
        "detail: {} XIC points, {} spectra, {} mirrors",
        detail["chromatogram"].as_array().map(Vec::len).unwrap_or(0),
        detail["spectra"].as_array().map(Vec::len).unwrap_or(0),
        detail["mirrors"].as_array().map(Vec::len).unwrap_or(0)
    );
    println!(
        "cold detail loaded in {:.1} ms",
        started.elapsed().as_secs_f64() * 1000.0
    );
    let started = std::time::Instant::now();
    let _ = visualizer_feature_data(output, sample, mz, rt, 0.1, &state)
        .expect("reload feature detail");
    println!(
        "warm detail loaded in {:.1} ms",
        started.elapsed().as_secs_f64() * 1000.0
    );
}
