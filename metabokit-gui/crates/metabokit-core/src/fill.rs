//! Gap filling.
//!
//! A consensus feature is rarely detected in every sample. For each blank the
//! area is re-integrated directly from the raw MS1 signal at the group's median
//! m/z and retention time, so a missing value becomes a measured zero rather
//! than a hole.
//!
//! 0.1 re-parsed the whole report once per sample, seeking back to the data
//! start each time — N passes over a file that can hold hundreds of thousands
//! of rows. Because the MS1 caches are memory-mapped here, all samples can be
//! open simultaneously at no cost in resident memory, so one streaming pass
//! fills every column at once.

use std::fs::File;
use std::io::{BufWriter, Write};

use crate::error::{IoContext, Result};
use crate::params::Params;
use crate::progress::{Cancel, Reporter};
use crate::scans::{ms1_cache_name, Ms1View};

/// Fill blanks in `Report_{table}.csv`, writing `Report_{table}_fill.csv`.
pub fn gap_fill(
    params: &Params,
    table: &str,
    file_stems: &[String],
    reporter: &dyn Reporter,
    cancel: &Cancel,
) -> Result<()> {
    let src = params.output_dir.join(format!("Report_{table}.csv"));
    let dst = params.output_dir.join(format!("Report_{table}_fill.csv"));
    if !src.is_file() {
        reporter.warn(format!("{} not found; skipping gap fill", src.display()));
        return Ok(());
    }

    let misc = params.misc_dir();
    // Mapping is address space, not resident memory: the OS pages in only the
    // scans an imputation actually touches.
    let views: Vec<Option<Ms1View>> = file_stems
        .iter()
        .map(|stem| Ms1View::open(&misc.join(ms1_cache_name(stem))).ok())
        .collect();
    let missing = views.iter().filter(|v| v.is_none()).count();
    if missing > 0 {
        reporter.warn(format!(
            "{missing} sample cache(s) unavailable; those columns cannot be imputed"
        ));
    }

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(&src)?;
    let mut records = rdr.records();

    let Some(header_top) = records.next().transpose()? else {
        return Ok(());
    };
    let Some(header) = records.next().transpose()? else {
        return Ok(());
    };

    let find = |name: &str| header.iter().position(|x| x.trim() == name);
    let (Some(detected_at), Some(mz_at), Some(rt_at)) =
        (find("%detected"), find("feature_m/z"), find("Median RT"))
    else {
        reporter.warn(format!(
            "{} is missing expected columns; skipping gap fill",
            src.display()
        ));
        return Ok(());
    };
    let area_at = detected_at + 1;
    let n_files = file_stems.len();

    let mut wtr = csv::WriterBuilder::new().from_path(&dst)?;
    wtr.write_record(&header_top)?;
    wtr.write_record(&header)?;

    let mut chrom: Vec<(f32, f32)> = Vec::new();
    let mut out: Vec<String> = Vec::new();
    let mut filled = 0usize;
    let mut rows = 0u64;

    for record in records {
        let record = record?;
        rows += 1;
        if rows % 4096 == 0 {
            cancel.check()?;
        }

        let mz: Option<f32> = record.get(mz_at).and_then(|x| x.trim().parse().ok());
        let rt: Option<f32> = record.get(rt_at).and_then(|x| x.trim().parse().ok());

        out.clear();
        out.extend(record.iter().take(area_at).map(str::to_string));

        for (i, view) in views.iter().enumerate() {
            let cell = record.get(area_at + i).unwrap_or("").trim();
            // A cell holds a plain number only when exactly one feature was
            // measured in that sample. Anything else — blank, a `*`-prefixed
            // imputed area, or a comma-separated list — is re-integrated, which
            // is what 0.1 did and what keeps the column numeric.
            if let Ok(value) = cell.parse::<f32>() {
                out.push(format!("{value:.1}"));
                continue;
            }
            let area = match (view, mz, rt) {
                (Some(view), Some(mz), Some(rt)) => {
                    view.xic(
                        mz,
                        rt,
                        params.impute_width,
                        params.integration_mz,
                        &mut chrom,
                    );
                    filled += 1;
                    trapezoid(&chrom) * 30.0
                }
                _ => 0.0,
            };
            out.push(format!("{area:.1}"));
        }

        out.extend(record.iter().skip(area_at + n_files).map(str::to_string));
        wtr.write_record(&out)?;
    }

    wtr.flush()?;
    reporter.info(format!(
        "gap fill ({table}): {filled} value(s) imputed across {rows} rows"
    ));
    Ok(())
}

fn trapezoid(chrom: &[(f32, f32)]) -> f32 {
    let mut acc = 0.0f32;
    for k in 0..chrom.len().saturating_sub(1) {
        acc += (chrom[k].1 + chrom[k + 1].1) * (chrom[k + 1].0 - chrom[k].0);
    }
    acc
}

/// Write the small binary blob the visualizer reads for its axis defaults.
/// Kept as a separate file so the UI does not have to parse the reports just
/// to learn the integration settings a run used.
pub fn write_run_manifest(params: &Params, positive: bool, stems: &[String]) -> Result<()> {
    let path = params.misc_dir().join("run.json");
    let manifest = serde_json::json!({
        "version": 2,
        "polarity": if positive { "positive" } else { "negative" },
        "ms1MsPair": params.ms1_ms2_pair,
        "integrationRt": params.integration_rt,
        "integrationMz": params.integration_mz,
        "imputeWidth": params.impute_width,
        "samples": stems,
    });
    let file = File::create(&path).at(&path)?;
    let mut w = BufWriter::new(file);
    w.write_all(manifest.to_string().as_bytes()).at(&path)?;
    w.flush().at(&path)?;
    Ok(())
}
