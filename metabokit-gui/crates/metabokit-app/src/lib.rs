//! Tauri bridge.
//!
//! This layer owns no analysis logic. It exposes `metabokit-core` to the web
//! frontend as commands, forwards engine progress as window events, and keeps
//! one run at a time alive on a dedicated thread so the UI never blocks.

pub mod visualizer;

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use metabokit_core::discover::{self, DatasetScan, DetectedLibrary};
use metabokit_core::library;
use metabokit_core::params::{LibrarySource, Params, Problem};
use metabokit_core::pipeline;
use metabokit_core::progress::{Cancel, Event, Reporter};
use tauri::{AppHandle, Emitter, State};
use tauri_plugin_opener::OpenerExt;

/// Engine progress. Payload is a `metabokit_core::progress::Event`.
const EV_PROGRESS: &str = "mk://event";
/// Run finished successfully. Payload is a `RunOutcome`.
const EV_FINISHED: &str = "mk://finished";
/// Run stopped at the user's request.
const EV_CANCELLED: &str = "mk://cancelled";
/// Run failed. Payload is the error message.
const EV_FAILED: &str = "mk://failed";

/// Forwards engine events to the window.
struct UiReporter {
    app: AppHandle,
}

impl Reporter for UiReporter {
    fn emit(&self, event: Event) {
        // A failed emit means the window is gone; the run will notice the
        // cancellation flag or simply finish into the void.
        let _ = self.app.emit(EV_PROGRESS, event);
    }
}

#[derive(Default)]
struct RunState {
    cancel: Cancel,
    running: Arc<AtomicBool>,
}

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

#[tauri::command]
fn default_params() -> Params {
    Params::default()
}

#[tauri::command]
fn validate_params(params: Params) -> Vec<Problem> {
    params.validate()
}

#[tauri::command]
fn load_preset(path: String) -> Result<Params, String> {
    Params::load(path).map_err(|e| e.to_string())
}

#[tauri::command]
fn save_preset(path: String, params: Params) -> Result<(), String> {
    params.save(path).map_err(|e| e.to_string())
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct BuiltinInfo {
    name: String,
    positive: bool,
    negative: bool,
}

/// Which built-in libraries are actually installed, per polarity.
#[tauri::command]
fn builtin_libraries(params: Params) -> Vec<BuiltinInfo> {
    let dir = params.resolve_libs_dir();
    library::available_builtins(dir.as_ref())
        .into_iter()
        .map(|(name, positive, negative)| BuiltinInfo {
            name,
            positive,
            negative,
        })
        .collect()
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct SystemInfo {
    version: String,
    os: String,
    available_threads: usize,
    libs_dir: Option<String>,
}

#[tauri::command]
fn system_info(params: Params) -> SystemInfo {
    SystemInfo {
        version: metabokit_core::VERSION.to_string(),
        os: std::env::consts::OS.to_string(),
        available_threads: std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
        libs_dir: params
            .resolve_libs_dir()
            .map(|p| p.to_string_lossy().into_owned()),
    }
}

#[derive(serde::Serialize)]
#[serde(rename_all = "camelCase")]
struct FileInfo {
    path: String,
    name: String,
    bytes: u64,
    exists: bool,
}

/// Metadata for the sample list. Deliberately cheap — no file is opened.
#[tauri::command]
fn describe_files(paths: Vec<String>) -> Vec<FileInfo> {
    paths
        .into_iter()
        .map(|p| {
            let path = PathBuf::from(&p);
            let meta = std::fs::metadata(&path).ok();
            FileInfo {
                name: path
                    .file_name()
                    .map(|x| x.to_string_lossy().into_owned())
                    .unwrap_or_else(|| p.clone()),
                bytes: meta.as_ref().map(|m| m.len()).unwrap_or(0),
                exists: meta.is_some(),
                path: p,
            }
        })
        .collect()
}

/// Read a file's declared scan polarity, so the UI can confirm the automatic
/// choice before a long run starts.
#[tauri::command]
fn sniff_polarity(path: String) -> Result<Option<bool>, String> {
    metabokit_core::mzml::sniff_polarity(std::path::Path::new(&path)).map_err(|e| e.to_string())
}

/// Scan a dataset folder and derive a complete, ready-to-run configuration.
///
/// This is the app's main entry point: everything the old parameter tabs asked
/// for is inferred here, and reported back with provenance so the user can see
/// what was assumed.
#[tauri::command]
fn scan_dataset(path: String) -> Result<DatasetScan, String> {
    discover::scan(std::path::Path::new(&path)).map_err(|e| e.to_string())
}

/// Cheap startup check used before restoring the last dataset. Keeping this in
/// Rust avoids granting the webview general filesystem access.
#[tauri::command]
fn directory_exists(path: String) -> bool {
    std::path::Path::new(&path).is_dir()
}

/// Open a completed run's output folder with the platform file manager.
///
/// The opener plugin's web command is scope-restricted and cannot safely be
/// pre-scoped to arbitrary user-selected output folders. A native command can
/// validate the exact directory first and then hand it to the OS directly.
#[tauri::command]
fn open_output_directory(app: AppHandle, path: String) -> Result<(), String> {
    if path.trim().is_empty() {
        return Err("No output folder is available yet.".to_string());
    }

    let directory = std::path::Path::new(&path);
    if !directory.is_dir() {
        return Err(format!("Output folder does not exist: {path}"));
    }

    app.opener()
        .open_path(path, None::<String>)
        .map_err(|e| format!("could not open output folder: {e}"))
}

/// Re-derive the library list after the user points at a `libs/` folder or adds
/// a file by hand. Cheaper than re-walking the dataset.
#[tauri::command]
fn relink_libraries(mut params: Params) -> (Params, Vec<DetectedLibrary>) {
    let detected = discover::refresh_libraries(&mut params);
    (params, detected)
}

/// Add user library files to the configuration, skipping duplicates.
#[tauri::command]
fn add_libraries(mut params: Params, paths: Vec<String>, kind: String) -> Params {
    for path in paths {
        let path = PathBuf::from(path);
        let source = if kind == "csv" {
            LibrarySource::Csv(path)
        } else {
            LibrarySource::Msp(path)
        };
        if !params.libraries.contains(&source) {
            params.libraries.push(source);
        }
    }
    params
}

#[tauri::command]
fn is_running(state: State<'_, RunState>) -> bool {
    state.running.load(Ordering::SeqCst)
}

#[tauri::command]
fn cancel_run(state: State<'_, RunState>) {
    state.cancel.cancel();
}

#[tauri::command]
fn start_run(app: AppHandle, state: State<'_, RunState>, params: Params) -> Result<(), String> {
    // `swap` rather than load-then-store: two rapid clicks must not both win.
    if state.running.swap(true, Ordering::SeqCst) {
        return Err("A run is already in progress.".to_string());
    }
    state.cancel.reset();
    let cancel = state.cancel.clone();
    let running = state.running.clone();

    let spawned = std::thread::Builder::new()
        .name("metabokit-run".to_string())
        // The aligner recurses through group formation on large runs; the
        // default 2 MB is tight on Windows.
        .stack_size(16 * 1024 * 1024)
        .spawn(move || {
            let reporter = UiReporter { app: app.clone() };
            let result = pipeline::run(&params, &reporter, &cancel);
            running.store(false, Ordering::SeqCst);
            match result {
                Ok(outcome) => {
                    let _ = app.emit(EV_FINISHED, outcome);
                }
                Err(e) if e.is_cancelled() => {
                    let _ = app.emit(EV_CANCELLED, ());
                }
                Err(e) => {
                    let _ = app.emit(EV_FAILED, e.to_string());
                }
            }
        });

    if let Err(e) = spawned {
        state.running.store(false, Ordering::SeqCst);
        return Err(format!("could not start run: {e}"));
    }
    Ok(())
}

// ---------------------------------------------------------------------------

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .manage(RunState::default())
        .manage(visualizer::VisualizerState::default())
        .invoke_handler(tauri::generate_handler![
            scan_dataset,
            directory_exists,
            open_output_directory,
            relink_libraries,
            add_libraries,
            default_params,
            validate_params,
            load_preset,
            save_preset,
            builtin_libraries,
            system_info,
            describe_files,
            sniff_polarity,
            is_running,
            start_run,
            cancel_run,
            visualizer::visualizer_open,
            visualizer::visualizer_overview,
            visualizer::visualizer_feature,
        ])
        .run(tauri::generate_context!())
        .expect("failed to start MetaboKit DDA");
}
