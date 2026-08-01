//! End-to-end smoke test against real mzML.
//!
//! Parses a sample, synthesises a small `.msp` from its own strongest MS2
//! spectra so that matches are guaranteed, then runs the full pipeline over the
//! whole dataset. It proves the plumbing — parse, detect, score, align, report,
//! gap fill — not the biology: the "identifications" are the data matched
//! against itself.
//!
//! ```text
//! cargo run --release -p metabokit-core --example smoke -- <dataset folder>
//! ```

use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

use metabokit_core::progress::{Cancel, Event, Level, Reporter};
use metabokit_core::{discover, mzml, params::LibrarySource, pipeline};

struct Console;

impl Reporter for Console {
    fn emit(&self, event: Event) {
        match event {
            Event::Stage { label, .. } => println!("== {label}"),
            Event::Log { level, message } => match level {
                Level::Warn => println!("   ! {message}"),
                Level::Error => println!("   E {message}"),
                Level::Info => println!("   · {message}"),
            },
            Event::Metric { key, value } => println!("   · {key} = {value}"),
            _ => {}
        }
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let Some(dir) = std::env::args().nth(1) else {
        eprintln!("usage: smoke <dataset folder>");
        std::process::exit(2);
    };
    let root = PathBuf::from(&dir);
    let cancel = Cancel::new();

    let mut scan = discover::scan(&root)?;
    if scan.samples.is_empty() {
        eprintln!("no samples found");
        std::process::exit(1);
    }

    // ---- parse one sample and report what came out ------------------------
    let first = scan.params.mzml_files[0].clone();
    println!("parsing {}", first.file_name().unwrap().to_string_lossy());
    let t = Instant::now();
    let data = mzml::parse(&first, &cancel)?;
    let parse_secs = t.elapsed().as_secs_f64();

    let rts = data.ms1.rts();
    println!(
        "  MS1 {} scans, {} points | MS2 {} scans | RT {:.2}–{:.2} min",
        data.ms1.len(),
        data.ms1.point_count(),
        data.ms2.len(),
        rts.first().copied().unwrap_or(0.0),
        rts.last().copied().unwrap_or(0.0),
    );
    println!(
        "  {:.1} MB resident, {:.1} s ({:.0} MB/s of mzML)",
        (data.ms1.heap_bytes() + data.ms2.heap_bytes()) as f64 / 1e6,
        parse_secs,
        std::fs::metadata(&first)?.len() as f64 / 1e6 / parse_secs,
    );
    if data.rt_converted_from_seconds {
        println!("  (scan times were in seconds and were converted)");
    }

    // ---- synthesise a library from the data's own strongest spectra --------
    let mut ms2 = data.ms2;
    ms2.sort_by_precursor();
    let mut ranked: Vec<(usize, f32)> = (0..ms2.len())
        .map(|i| {
            let (_, ints) = ms2.scan(i);
            (i, ints.iter().copied().fold(0.0f32, f32::max))
        })
        .filter(|(_, v)| *v > 0.0)
        .collect();
    ranked.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    ranked.truncate(40);

    let lib_path = std::env::temp_dir().join("metabokit-smoke.msp");
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&lib_path)?);
        for (n, (i, _)) in ranked.iter().enumerate() {
            let (mzs, ints) = ms2.scan(*i);
            let mut peaks: Vec<(f32, f32)> =
                mzs.iter().copied().zip(ints.iter().copied()).collect();
            peaks.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            peaks.truncate(10);
            if peaks.len() < 3 {
                continue;
            }
            writeln!(w, "NAME: SMOKE-{n:03}")?;
            writeln!(w, "PRECURSORTYPE: [M+H]+")?;
            writeln!(w, "PRECURSORMZ: {:.4}", ms2.prec_mz(*i))?;
            writeln!(w, "INCHIKEY: SMOKETEST{n:03}AAAAAAAAAAAAA")?;
            writeln!(w, "Num Peaks: {}", peaks.len())?;
            for (mz, inten) in peaks {
                writeln!(w, "{mz:.4} {inten:.1}")?;
            }
            writeln!(w)?;
        }
    }
    println!("\nsynthetic library at {}", lib_path.display());

    // ---- run the pipeline over the whole dataset --------------------------
    scan.params.libraries = vec![LibrarySource::Msp(lib_path)];
    scan.params.output_dir = root.join("results-smoke");
    // The synthetic entries are the spectra themselves, so a low bar is fine.
    scan.params.min_peaks = 2;
    scan.params.ms2_score = 0.3;

    println!("\nrunning pipeline over {} samples\n", scan.samples.len());
    let reporter = Console;
    let outcome = pipeline::run(&scan.params, &reporter, &cancel)?;

    println!("\n=== result ===");
    println!("elapsed        {:.1} s", outcome.elapsed_seconds);
    println!("library        {} entries", outcome.library_entries);
    println!("features       {}", outcome.summary.features);
    println!("identified     {}", outcome.summary.identified);
    println!("compounds      {}", outcome.summary.compounds);
    println!("isf rows       {}", outcome.summary.isf_rows);
    for s in &outcome.samples {
        println!(
            "  {:<58} {:>6} feat {:>6} ann {:>6.1} s",
            &s.name[..s.name.len().min(58)],
            s.features,
            s.annotations,
            s.seconds
        );
    }
    println!("\noutput -> {}", outcome.output_dir);
    Ok(())
}
