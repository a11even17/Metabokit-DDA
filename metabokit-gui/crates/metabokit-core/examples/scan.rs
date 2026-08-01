//! Print what dataset discovery makes of a folder, without launching the GUI.
//!
//! ```text
//! cargo run -p metabokit-core --example scan -- /path/to/dataset
//! ```

fn main() {
    let Some(dir) = std::env::args().nth(1) else {
        eprintln!("usage: scan <dataset folder>");
        std::process::exit(2);
    };

    match metabokit_core::discover::scan(std::path::Path::new(&dir)) {
        Ok(scan) => {
            println!("root       {}", scan.root);
            println!("polarity   {}", scan.polarity);
            println!("output     {}", scan.output_dir);
            println!("libs dir   {}", scan.libs_dir.as_deref().unwrap_or("—"));
            println!(
                "settings   {}",
                scan.imported_settings.as_deref().unwrap_or("defaults")
            );
            println!("ready      {}", scan.ready);
            if let Some(run) = &scan.previous_run {
                println!(
                    "previous   {} samples, {} features, {} identified, {} compounds",
                    run.samples,
                    run.summary.features,
                    run.summary.identified,
                    run.summary.compounds
                );
            }

            println!("\nsamples ({})", scan.samples.len());
            for s in &scan.samples {
                println!(
                    "  {:<70} {:>8.1} MB  {}",
                    s.name,
                    s.bytes as f64 / 1e6,
                    s.subfolder.as_deref().unwrap_or("")
                );
            }

            println!("\nlibraries ({})", scan.libraries.len());
            for l in &scan.libraries {
                println!("  [{}] {} — {}", l.kind, l.label, l.detail);
            }

            println!("\nnotes");
            for n in &scan.notes {
                println!("  {:?} {:<10} {}", n.level, n.topic, n.message);
            }
        }
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    }
}
