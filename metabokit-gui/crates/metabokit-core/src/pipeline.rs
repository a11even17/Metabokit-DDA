//! Run orchestration.
//!
//! # Parallelism and memory
//!
//! 0.1 ran `mzml_fs.par_iter()`: every worker thread decoded, transformed and
//! scored a *different* sample, so peak resident memory was
//! `threads × sample size`. On a 16-core machine with 2 GB runs that is not a
//! workload, it is an out-of-memory error.
//!
//! Here the parallelism has moved one level down — into the m/z slices of
//! feature detection, which is where the time actually goes. Samples are
//! processed in small batches (`max_files_in_flight`, default 2) so the
//! high-water mark is set by file size rather than core count, while parsing
//! and scoring of neighbouring samples still overlap.

use std::path::{Path, PathBuf};
use std::time::Instant;

use rayon::prelude::*;

use crate::align;
use crate::error::{Error, IoContext, Result};
use crate::features;
use crate::fill;
use crate::library::{self, Library};
use crate::mzml;
use crate::params::{Params, Polarity};
use crate::progress::{Cancel, Reporter, Stage};
use crate::report::{self, ReportSummary};
use crate::scans::{
    feature_cache_name, ms1_cache_name, ms2_cache_name, write_feature_cache, write_ms1_cache,
    write_ms2_cache,
};
use crate::score::{self, Ann, SpecSet};

/// Per-sample statistics, surfaced in the UI when a run finishes.
#[derive(Clone, Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SampleStat {
    pub name: String,
    pub ms1_scans: usize,
    pub ms2_scans: usize,
    pub features: usize,
    pub annotations: usize,
    pub seconds: f64,
}

#[derive(Clone, Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunOutcome {
    pub samples: Vec<SampleStat>,
    pub summary: ReportSummary,
    pub polarity: String,
    pub library_entries: usize,
    pub elapsed_seconds: f64,
    pub output_dir: String,
    pub reports: Vec<String>,
}

struct SampleOutput {
    stem: String,
    timestamp: String,
    annotations: Vec<Ann>,
    spectra: SpecSet,
    stat: SampleStat,
}

/// Execute a full run.
pub fn run(params: &Params, reporter: &dyn Reporter, cancel: &Cancel) -> Result<RunOutcome> {
    let started = Instant::now();
    params.check()?;

    reporter.stage(Stage::Preparing);
    std::fs::create_dir_all(&params.output_dir).at(&params.output_dir)?;
    let misc = params.misc_dir();
    std::fs::create_dir_all(&misc).at(&misc)?;
    clear_stale_cache(&misc);

    let positive = resolve_polarity(params, reporter)?;
    reporter.metric(
        "polarity",
        if positive { "positive" } else { "negative" },
    );

    reporter.stage(Stage::Library);
    let lib = library::load(params, positive, reporter)?;
    let meta = library::load_metadata(params, reporter);
    cancel.check()?;

    // A private pool rather than `build_global`: a GUI can start a second run
    // in the same process, and the global pool can only be installed once.
    let threads = params.effective_threads();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .thread_name(|i| format!("metabokit-{i}"))
        .build()
        .map_err(|e| Error::param(format!("thread pool: {e}")))?;
    reporter.metric("threads", threads.to_string());
    let library_entries = lib.len();

    let outputs = pool.install(|| {
        process_samples(params, &lib, positive, reporter, cancel)
    })?;

    cancel.check()?;

    let file_stems: Vec<String> = outputs.iter().map(|o| o.stem.clone()).collect();
    let timestamps: Vec<String> = outputs.iter().map(|o| o.timestamp.clone()).collect();
    let stats: Vec<SampleStat> = outputs.iter().map(|o| o.stat.clone()).collect();

    reporter.stage(Stage::Aligning);
    // Move the annotations out sample by sample so the per-sample vectors are
    // freed as they are absorbed, rather than doubling up at the peak.
    let mut all_annotations: Vec<Ann> = Vec::new();
    let mut spectra: Vec<SpecSet> = Vec::with_capacity(outputs.len());
    for output in outputs {
        all_annotations.extend(output.annotations);
        spectra.push(output.spectra);
    }
    reporter.metric("annotations", all_annotations.len().to_string());

    let aligned = align::align(all_annotations, params);
    reporter.metric("consensusFeatures", aligned.groups.len().to_string());
    cancel.check()?;

    reporter.stage(Stage::Reporting);
    let summary = report::write_reports(
        &aligned,
        &lib,
        &meta,
        params,
        positive,
        &file_stems,
        &timestamps,
        &spectra,
        reporter,
    )?;

    // Everything below only reads the reports and the mmap'd caches, so the
    // library, annotations and spectra can go now.
    drop(aligned);
    drop(spectra);
    drop(lib);
    drop(meta);

    let mut reports = vec![
        "Report_RTseparated.csv".to_string(),
        "Report_by_ID.csv".to_string(),
        "spec_exp.txt".to_string(),
        "spec_lib.txt".to_string(),
        "spec_reduced.txt".to_string(),
    ];

    if params.gap_fill {
        reporter.stage(Stage::GapFilling);
        for table in ["RTseparated", "by_ID"] {
            cancel.check()?;
            fill::gap_fill(params, table, &file_stems, reporter, cancel)?;
            reports.push(format!("Report_{table}_fill.csv"));
        }
    }

    fill::write_run_manifest(params, positive, &file_stems)?;

    if !params.keep_cache {
        clear_stale_cache(&misc);
        reporter.info("scan caches removed (the visualizer will need a re-run)");
    }

    Ok(RunOutcome {
        samples: stats,
        summary,
        polarity: if positive { "positive" } else { "negative" }.to_string(),
        library_entries,
        elapsed_seconds: started.elapsed().as_secs_f64(),
        output_dir: params.output_dir.to_string_lossy().into_owned(),
        reports,
    })
}

/// Parse, detect and score every sample, in bounded-size batches.
fn process_samples(
    params: &Params,
    lib: &Library,
    positive: bool,
    reporter: &dyn Reporter,
    cancel: &Cancel,
) -> Result<Vec<SampleOutput>> {
    let total = params.mzml_files.len();
    let in_flight = params.effective_in_flight();
    let lib_masses = lib.masses();
    let mut outputs: Vec<SampleOutput> = Vec::with_capacity(total);
    let mut done = 0u64;

    reporter.stage(Stage::Processing);

    for batch in params.mzml_files.chunks(in_flight) {
        cancel.check()?;
        let base = outputs.len();
        let batch_results: Vec<Result<SampleOutput>> = batch
            .par_iter()
            .enumerate()
            .map(|(offset, path)| {
                process_one(
                    (base + offset) as u32,
                    path,
                    params,
                    lib,
                    lib_masses,
                    positive,
                    reporter,
                    cancel,
                )
            })
            .collect();

        for result in batch_results {
            outputs.push(result?);
            done += 1;
            reporter.progress(done, total as u64);
        }
    }
    Ok(outputs)
}

#[allow(clippy::too_many_arguments)]
fn process_one(
    index: u32,
    path: &Path,
    params: &Params,
    lib: &Library,
    lib_masses: &[f32],
    positive: bool,
    reporter: &dyn Reporter,
    cancel: &Cancel,
) -> Result<SampleOutput> {
    let started = Instant::now();
    let stem = path
        .file_stem()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| format!("sample_{index}"));
    let display = path
        .file_name()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| stem.clone());

    reporter.emit(crate::progress::Event::Sample {
        name: display.clone(),
        index: index as usize,
        total: params.mzml_files.len(),
    });

    let mut data = mzml::parse(path, cancel)?;
    if data.rt_converted_from_seconds {
        reporter.info(format!("{display}: scan times were in seconds, converted"));
    }
    if data.unknown_mode_spectra > 0 {
        reporter.warn(format!(
            "{display}: {} spectra declare neither centroid nor profile mode; assuming centroid",
            data.unknown_mode_spectra
        ));
    }
    if data.ms1.is_empty() {
        return Err(Error::mzml(path, "no MS1 scans"));
    }
    if data.ms2.is_empty() {
        reporter.warn(format!("{display}: no MS2 scans; nothing to identify"));
    }
    data.ms2.sort_by_precursor();

    let ms1_scans = data.ms1.len();
    let ms2_scans = data.ms2.len();
    reporter.info(format!(
        "{display}: {ms1_scans} MS1 + {ms2_scans} MS2 scans, {:.0} MB",
        (data.ms1.heap_bytes() + data.ms2.heap_bytes()) as f64 / 1e6
    ));

    // Persist before the in-memory copies are dropped; gap filling and the
    // visualizer both read these back through a memory map.
    let misc = params.misc_dir();
    write_ms1_cache(&misc.join(ms1_cache_name(&stem)), &data.ms1, &data.timestamp)?;
    write_ms2_cache(&misc.join(ms2_cache_name(&stem)), &data.ms2)?;

    cancel.check()?;
    let peaks = features::detect(&data.ms1, &data.ms2, lib_masses, params, cancel)?;
    reporter.info(format!("{display}: {} features", peaks.len()));
    write_feature_cache(&misc.join(feature_cache_name(&stem)), &peaks)?;

    cancel.check()?;
    let result = score::score_sample(
        index,
        &data.ms1,
        &data.ms2,
        &peaks,
        lib,
        params,
        cancel,
    )?;

    report::write_unknowns(
        &misc.join(format!("u_{stem}.txt")),
        index as usize,
        &result.spectra,
        positive,
    )?;

    let stat = SampleStat {
        name: display,
        ms1_scans,
        ms2_scans,
        features: peaks.len(),
        annotations: result.annotations.len(),
        seconds: started.elapsed().as_secs_f64(),
    };

    Ok(SampleOutput {
        stem,
        timestamp: data.timestamp,
        annotations: result.annotations,
        spectra: result.spectra,
        stat,
    })
}

/// Decide the run's polarity, sniffing the first file when set to auto.
fn resolve_polarity(params: &Params, reporter: &dyn Reporter) -> Result<bool> {
    match params.polarity {
        Polarity::Positive => Ok(true),
        Polarity::Negative => Ok(false),
        Polarity::Auto => {
            let first = params
                .mzml_files
                .first()
                .ok_or_else(|| Error::param("no mzML files"))?;
            match mzml::sniff_polarity(first)? {
                Some(p) => Ok(p),
                None => {
                    reporter.warn(
                        "no scan-polarity cvParam found; assuming positive mode. Set it \
                         explicitly if that is wrong.",
                    );
                    Ok(true)
                }
            }
        }
    }
}

/// Remove caches from a previous run so a stale sample cannot leak into this
/// one. Only files this pipeline writes are touched — the output directory may
/// legitimately hold other work.
fn clear_stale_cache(misc: &Path) {
    let Ok(entries) = std::fs::read_dir(misc) else {
        return;
    };
    for entry in entries.flatten() {
        let path: PathBuf = entry.path();
        let Some(name) = path.file_name().and_then(|x| x.to_str()) else {
            continue;
        };
        let ours = (name.starts_with("ms1_")
            || name.starts_with("ms2_")
            || name.starts_with("features_"))
            && name.ends_with(".mkc")
            || (name.starts_with("u_") && name.ends_with(".txt"))
            || name == "run.json";
        if ours {
            let _ = std::fs::remove_file(&path);
        }
    }
}
