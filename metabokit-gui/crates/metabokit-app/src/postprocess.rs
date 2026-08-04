//! Streaming quality-control filters for completed CSV reports.
//!
//! Reports can be wide and large, so preview and export both process one CSV
//! record at a time. Only the small on-screen preview is retained in memory;
//! the original report is never modified.

use std::fs::{self, File};
use std::path::{Path, PathBuf};

use csv::{ReaderBuilder, StringRecord, WriterBuilder};
use serde::{Deserialize, Serialize};

const PREVIEW_ROWS: usize = 100;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PostSession {
    reports: Vec<String>,
}

#[derive(Clone, Deserialize)]
#[serde(rename_all = "camelCase", default)]
pub struct PostOptions {
    pub minimum_detected_percent: f32,
    pub minimum_peak_shape: f32,
    pub minimum_score: f32,
    pub minimum_sn: f32,
    pub minimum_matching_peaks: f32,
    pub maximum_cv_percent: Option<f32>,
    pub identified_only: bool,
    pub msms_only: bool,
    pub remove_isf: bool,
}

impl Default for PostOptions {
    fn default() -> Self {
        Self {
            minimum_detected_percent: 0.0,
            minimum_peak_shape: 0.0,
            minimum_score: 0.0,
            minimum_sn: 0.0,
            minimum_matching_peaks: 0.0,
            maximum_cv_percent: None,
            identified_only: false,
            msms_only: false,
            remove_isf: true,
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PostRow {
    group: String,
    name: String,
    adduct: String,
    mz: Option<f32>,
    detected_percent: Option<f32>,
    peak_shape: Option<f32>,
    score: Option<f32>,
    sn: Option<f32>,
    cv_percent: Option<f32>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct FilterCount {
    criterion: String,
    rows: usize,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PostPreview {
    source: String,
    total_rows: usize,
    kept_rows: usize,
    preview_truncated: bool,
    rows: Vec<PostRow>,
    removed_by: Vec<FilterCount>,
}

#[derive(Default)]
struct Columns {
    group: Option<usize>,
    name: Option<usize>,
    adduct: Option<usize>,
    mz: Option<usize>,
    detected: Option<usize>,
    shape: Option<usize>,
    matching_peaks: Option<usize>,
    confidence: Option<usize>,
    isf: Vec<usize>,
    scores: Vec<usize>,
    sns: Vec<usize>,
    areas: Vec<usize>,
}

impl Columns {
    fn from_header(header: &StringRecord) -> Self {
        let find = |name: &str| header.iter().position(|value| value == name);
        let prefixed = |prefix: &str| {
            header
                .iter()
                .enumerate()
                .filter_map(|(i, value)| value.starts_with(prefix).then_some(i))
                .collect()
        };
        Self {
            group: find("group"),
            name: find("name"),
            adduct: find("adduct"),
            mz: find("feature_m/z"),
            detected: find("%detected"),
            shape: find("peak_shape (median)"),
            matching_peaks: find("matching_peaks (median)"),
            confidence: find("confidence level"),
            isf: [find("ISF"), find("ISF of")]
                .into_iter()
                .flatten()
                .collect(),
            scores: prefixed("SCORE_"),
            sns: prefixed("S/N_"),
            areas: prefixed("AREA_"),
        }
    }
}

struct Metrics {
    detected_percent: Option<f32>,
    peak_shape: Option<f32>,
    score: Option<f32>,
    sn: Option<f32>,
    matching_peaks: Option<f32>,
    cv_percent: Option<f32>,
    identified: bool,
    msms: bool,
    isf: bool,
}

impl Metrics {
    fn read(record: &StringRecord, columns: &Columns) -> Self {
        let number = |at: Option<usize>| at.and_then(|i| parse(record.get(i)));
        let mut detected_percent = number(columns.detected);
        if detected_percent.is_some_and(|value| value <= 1.0) {
            detected_percent = detected_percent.map(|value| value * 100.0);
        }
        Self {
            detected_percent,
            peak_shape: number(columns.shape),
            score: median_values(record, &columns.scores),
            sn: median_values(record, &columns.sns),
            matching_peaks: number(columns.matching_peaks),
            cv_percent: coefficient_of_variation(record, &columns.areas),
            identified: columns
                .name
                .and_then(|i| record.get(i))
                .is_some_and(|value| {
                    !value.trim().is_empty() && !value.eq_ignore_ascii_case("unknown")
                }),
            msms: columns
                .confidence
                .and_then(|i| record.get(i))
                .is_some_and(|value| value.to_ascii_uppercase().contains("MSMS")),
            isf: columns
                .isf
                .iter()
                .any(|&i| record.get(i).is_some_and(|value| !value.trim().is_empty())),
        }
    }

    fn rejected_by(&self, options: &PostOptions) -> Option<&'static str> {
        if options.remove_isf && self.isf {
            return Some("In-source fragments");
        }
        if options.identified_only && !self.identified {
            return Some("Unidentified");
        }
        if options.msms_only && !self.msms {
            return Some("Not MS/MS confidence");
        }
        if below(self.detected_percent, options.minimum_detected_percent) {
            return Some("Detection frequency");
        }
        if below(self.peak_shape, options.minimum_peak_shape) {
            return Some("Peak shape");
        }
        if below(self.score, options.minimum_score) {
            return Some("MS/MS score");
        }
        if below(self.sn, options.minimum_sn) {
            return Some("Signal/noise");
        }
        if below(self.matching_peaks, options.minimum_matching_peaks) {
            return Some("Matching fragments");
        }
        if let Some(maximum) = options.maximum_cv_percent.filter(|value| *value > 0.0) {
            if self.cv_percent.map_or(true, |value| value > maximum) {
                return Some("Area CV");
            }
        }
        None
    }
}

fn parse(value: Option<&str>) -> Option<f32> {
    value?.trim().trim_start_matches('*').parse().ok()
}

fn below(value: Option<f32>, minimum: f32) -> bool {
    minimum > 0.0 && value.map_or(true, |actual| actual < minimum)
}

fn median_values(record: &StringRecord, columns: &[usize]) -> Option<f32> {
    let mut values: Vec<f32> = columns
        .iter()
        .filter_map(|&i| parse(record.get(i)))
        .collect();
    if values.is_empty() {
        return None;
    }
    values.sort_unstable_by(f32::total_cmp);
    let middle = values.len() / 2;
    Some(if values.len() % 2 == 0 {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    })
}

fn coefficient_of_variation(record: &StringRecord, columns: &[usize]) -> Option<f32> {
    let values: Vec<f64> = columns
        .iter()
        .filter_map(|&i| parse(record.get(i)).map(f64::from))
        .collect();
    if values.len() < 2 {
        return None;
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    if mean.abs() <= f64::EPSILON {
        return None;
    }
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    Some((variance.sqrt() / mean.abs() * 100.0) as f32)
}

fn output_path(value: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(value);
    if !path.is_dir() {
        return Err(format!("Results folder does not exist: {value}"));
    }
    Ok(path)
}

fn report_path(output: &Path, report: &str) -> Result<PathBuf, String> {
    let name = Path::new(report);
    if name.file_name().and_then(|value| value.to_str()) != Some(report)
        || name.extension().and_then(|value| value.to_str()) != Some("csv")
    {
        return Err("Invalid report name.".to_string());
    }
    let path = output.join(name);
    if !path.is_file() {
        return Err(format!("Report does not exist: {report}"));
    }
    Ok(path)
}

fn reports(output: &Path) -> Result<Vec<String>, String> {
    let mut names: Vec<String> = fs::read_dir(output)
        .map_err(|e| format!("could not read {}: {e}", output.display()))?
        .filter_map(Result::ok)
        .filter_map(|entry| {
            let name = entry.file_name().to_string_lossy().into_owned();
            (entry.path().is_file()
                && name.starts_with("Report")
                && name.to_ascii_lowercase().ends_with(".csv"))
            .then_some(name)
        })
        .collect();
    names.sort_by_key(|name| {
        if name == "Report_RTseparated.csv" {
            (0, name.clone())
        } else if name == "Report_ID.csv" {
            (1, name.clone())
        } else {
            (2, name.clone())
        }
    });
    Ok(names)
}

#[tauri::command]
pub fn postprocess_open(output_dir: String) -> Result<PostSession, String> {
    let output = output_path(&output_dir)?;
    Ok(PostSession {
        reports: reports(&output)?,
    })
}

fn process(
    path: &Path,
    options: &PostOptions,
    mut writer: Option<&mut csv::Writer<File>>,
) -> Result<PostPreview, String> {
    let mut reader = ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(path)
        .map_err(|e| format!("could not read {}: {e}", path.display()))?;
    let mut columns = None;
    let mut total_rows = 0usize;
    let mut kept_rows = 0usize;
    let mut preview = Vec::with_capacity(PREVIEW_ROWS);
    let mut removed: Vec<(&'static str, usize)> = Vec::new();

    for result in reader.records() {
        let record = result.map_err(|e| format!("{}: {e}", path.display()))?;
        if columns.is_none() {
            if record.get(0) == Some("group") {
                columns = Some(Columns::from_header(&record));
            }
            if let Some(output) = writer.as_mut() {
                output
                    .write_record(&record)
                    .map_err(|e| format!("could not write export: {e}"))?;
            }
            continue;
        }

        total_rows += 1;
        let columns = columns.as_ref().expect("columns established above");
        let metrics = Metrics::read(&record, columns);
        if let Some(reason) = metrics.rejected_by(options) {
            if let Some((_, count)) = removed.iter_mut().find(|(name, _)| *name == reason) {
                *count += 1;
            } else {
                removed.push((reason, 1));
            }
            continue;
        }

        kept_rows += 1;
        if let Some(output) = writer.as_mut() {
            output
                .write_record(&record)
                .map_err(|e| format!("could not write export: {e}"))?;
        }
        if preview.len() < PREVIEW_ROWS {
            let field =
                |at: Option<usize>| at.and_then(|i| record.get(i)).unwrap_or("").to_string();
            preview.push(PostRow {
                group: field(columns.group),
                name: field(columns.name),
                adduct: field(columns.adduct),
                mz: columns.mz.and_then(|i| parse(record.get(i))),
                detected_percent: metrics.detected_percent,
                peak_shape: metrics.peak_shape,
                score: metrics.score,
                sn: metrics.sn,
                cv_percent: metrics.cv_percent,
            });
        }
    }
    if columns.is_none() {
        return Err("The report header could not be found.".to_string());
    }

    Ok(PostPreview {
        source: path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned(),
        total_rows,
        kept_rows,
        preview_truncated: kept_rows > preview.len(),
        rows: preview,
        removed_by: removed
            .into_iter()
            .map(|(criterion, rows)| FilterCount {
                criterion: criterion.to_string(),
                rows,
            })
            .collect(),
    })
}

#[tauri::command]
pub fn postprocess_preview(
    output_dir: String,
    report: String,
    options: PostOptions,
) -> Result<PostPreview, String> {
    let output = output_path(&output_dir)?;
    let path = report_path(&output, &report)?;
    process(&path, &options, None)
}

#[tauri::command]
pub fn postprocess_export(
    output_dir: String,
    report: String,
    destination: String,
    options: PostOptions,
) -> Result<PostPreview, String> {
    let output = output_path(&output_dir)?;
    let source = report_path(&output, &report)?;
    let destination = PathBuf::from(destination);
    if destination.as_os_str().is_empty() {
        return Err("Choose an export file first.".to_string());
    }
    if destination.exists() && source.canonicalize().ok() == destination.canonicalize().ok() {
        return Err("Choose a new file; the original report is never overwritten.".to_string());
    }
    let parent = destination
        .parent()
        .ok_or("The export location is invalid.")?;
    if !parent.is_dir() {
        return Err("The export folder does not exist.".to_string());
    }
    let file = File::create(&destination)
        .map_err(|e| format!("could not create {}: {e}", destination.display()))?;
    let mut writer = WriterBuilder::new().flexible(true).from_writer(file);
    let preview = process(&source, &options, Some(&mut writer))?;
    writer
        .flush()
        .map_err(|e| format!("could not finish export: {e}"))?;
    Ok(preview)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cv_uses_sample_standard_deviation() {
        let record = StringRecord::from(vec!["10", "12", "14"]);
        let cv = coefficient_of_variation(&record, &[0, 1, 2]).unwrap();
        assert!((cv - 16.666_666).abs() < 0.001);
    }

    #[test]
    fn a_requested_missing_metric_is_rejected() {
        assert!(below(None, 0.5));
        assert!(!below(None, 0.0));
    }

    #[test]
    fn report_filtering_is_streamed_and_keeps_the_preamble() {
        let path = std::env::temp_dir().join(format!(
            "metabokit-postprocess-{}-{}.csv",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(
            &path,
            ",,2026-01-01\ngroup,name,feature_m/z,%detected,peak_shape (median),SCORE_A,SCORE_B,S/N_A,S/N_B,AREA_A,AREA_B\n1,keep,100.0,1.0,0.9,0.8,0.6,10,12,100,110\n2,drop,200.0,0.2,0.9,0.9,0.9,20,20,100,100\n",
        )
        .unwrap();
        let options = PostOptions {
            minimum_detected_percent: 50.0,
            remove_isf: false,
            ..PostOptions::default()
        };
        let result = process(&path, &options, None).unwrap();
        assert_eq!(result.total_rows, 2);
        assert_eq!(result.kept_rows, 1);
        assert_eq!(result.rows[0].name, "keep");
        assert_eq!(result.removed_by[0].criterion, "Detection frequency");
        let destination = path.with_extension("clean.csv");
        let exported = postprocess_export(
            path.parent().unwrap().to_string_lossy().into_owned(),
            path.file_name().unwrap().to_string_lossy().into_owned(),
            destination.to_string_lossy().into_owned(),
            options.clone(),
        )
        .unwrap();
        assert_eq!(exported.kept_rows, 1);
        let clean = fs::read_to_string(&destination).unwrap();
        assert!(clean.contains("group,name,feature_m/z"));
        assert!(clean.contains("1,keep,100.0"));
        assert!(!clean.contains("2,drop,200.0"));
        assert!(postprocess_export(
            path.parent().unwrap().to_string_lossy().into_owned(),
            path.file_name().unwrap().to_string_lossy().into_owned(),
            path.to_string_lossy().into_owned(),
            options,
        )
        .is_err());
        let _ = fs::remove_file(destination);
        let _ = fs::remove_file(path);
    }
}
