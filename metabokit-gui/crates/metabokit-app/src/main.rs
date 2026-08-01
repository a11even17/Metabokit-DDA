// Suppress the console window that Windows would otherwise attach to a GUI
// binary in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

fn main() {
    metabokit_app_lib::run();
}
