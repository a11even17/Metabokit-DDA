//! Bounded-memory data service for the integrated visualizer.
//!
//! Large scan caches stay memory-mapped. IPC responses are columnar where they
//! can be large, and only one selected feature's chromatogram/spectra are ever
//! materialized. The spectrum exports are indexed by byte offset, so browsing
//! a new feature does not re-read or retain every peak in every match.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::UNIX_EPOCH;

use metabokit_core::scans::{
    feature_cache_name, ms1_cache_name, ms2_cache_name, FeatureView, Ms1View, Ms2View,
};
use serde::{Deserialize, Serialize};
use tauri::State;

const MAX_OVERVIEW_POINTS: usize = 200_000;
const MAX_SPECTRA: usize = 96;
const MAX_SPECTRUM_PEAKS: usize = 4_096;
const MAX_MIRRORS: usize = 38;

#[derive(Default)]
pub struct VisualizerState {
    match_index: Mutex<Option<MatchIndex>>,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(default, rename_all = "camelCase")]
struct RunSettings {
    ms1_ms_pair: f32,
    integration_rt: f32,
    integration_mz: f32,
    samples: Vec<String>,
}

impl Default for RunSettings {
    fn default() -> Self {
        Self {
            ms1_ms_pair: 0.01,
            integration_rt: 1.5,
            integration_mz: 0.006,
            samples: Vec::new(),
        }
    }
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizParams {
    ms1_ms_pair: f32,
    integration_rt: f32,
    integration_mz: f32,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizSample {
    id: String,
    label: String,
    cache_bytes: u64,
    has_features: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizSession {
    samples: Vec<VizSample>,
    params: VizParams,
    matches_available: bool,
}

#[tauri::command]
pub fn visualizer_open(output_dir: String) -> Result<VizSession, String> {
    let output = checked_output(&output_dir)?;
    let misc = output.join("misc");
    let settings = read_settings(&misc);
    let mut available = cache_samples(&misc)?;

    let mut samples = Vec::new();
    for wanted in &settings.samples {
        if let Some((ms1, ms2)) = available.remove(wanted) {
            samples.push(sample_info(&misc, wanted, ms1, ms2));
        }
    }
    let mut remaining: Vec<_> = available.into_iter().collect();
    remaining.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    for (sample, (ms1, ms2)) in remaining {
        samples.push(sample_info(&misc, &sample, ms1, ms2));
    }

    Ok(VizSession {
        samples,
        params: VizParams {
            ms1_ms_pair: settings.ms1_ms_pair,
            integration_rt: settings.integration_rt,
            integration_mz: settings.integration_mz,
        },
        matches_available: output.join("spec_exp.txt").is_file()
            && output.join("spec_lib.txt").is_file(),
    })
}

fn checked_output(value: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(value);
    if !path.is_dir() {
        return Err(format!("Results folder does not exist: {}", path.display()));
    }
    Ok(path)
}

fn read_settings(misc: &Path) -> RunSettings {
    File::open(misc.join("run.json"))
        .ok()
        .and_then(|file| serde_json::from_reader(file).ok())
        .unwrap_or_default()
}

fn cache_samples(misc: &Path) -> Result<HashMap<String, (u64, u64)>, String> {
    let mut out = HashMap::new();
    let entries =
        fs::read_dir(misc).map_err(|e| format!("could not read {}: {e}", misc.display()))?;
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name = name.to_string_lossy();
        let Some(sample) = name
            .strip_prefix("ms1_")
            .and_then(|value| value.strip_suffix(".mkc"))
        else {
            continue;
        };
        let ms2 = misc.join(ms2_cache_name(sample));
        if !ms2.is_file() {
            continue;
        }
        let ms1_bytes = entry.metadata().map(|m| m.len()).unwrap_or(0);
        let ms2_bytes = fs::metadata(ms2).map(|m| m.len()).unwrap_or(0);
        out.insert(sample.to_string(), (ms1_bytes, ms2_bytes));
    }
    Ok(out)
}

fn sample_info(misc: &Path, sample: &str, ms1: u64, ms2: u64) -> VizSample {
    VizSample {
        id: sample.to_string(),
        label: sample.to_string(),
        cache_bytes: ms1 + ms2,
        has_features: misc.join(feature_cache_name(sample)).is_file(),
    }
}

fn validate_sample(sample: &str) -> Result<(), String> {
    if sample.is_empty()
        || sample.contains('/')
        || sample.contains('\\')
        || sample == "."
        || sample == ".."
    {
        return Err("Invalid sample name.".to_string());
    }
    Ok(())
}

#[derive(Serialize, Default)]
#[serde(rename_all = "camelCase")]
pub struct VizOverview {
    mz: Vec<f32>,
    rt: Vec<f32>,
    half_width: Vec<f32>,
    coefficient: Vec<f32>,
    shape: Vec<f32>,
    smoothness: Vec<f32>,
    group: Vec<u32>,
    names: Vec<String>,
    has_ms2: Vec<bool>,
    total: usize,
    truncated: bool,
    source: String,
}

#[derive(Clone)]
struct ReportPoint {
    mz: f32,
    rt: f32,
    shape: f32,
    smoothness: f32,
    group: u32,
    name: String,
}

struct Point {
    mz: f32,
    rt: f32,
    half_width: f32,
    coefficient: f32,
    shape: f32,
    smoothness: f32,
    group: u32,
    name: String,
    has_ms2: bool,
}

impl VizOverview {
    fn consider(&mut self, point: Point) {
        self.total += 1;
        let slot = if self.mz.len() < MAX_OVERVIEW_POINTS {
            self.mz.len()
        } else {
            // Deterministic reservoir sampling keeps the full RT/mz range
            // represented without retaining an unbounded response.
            let mixed = splitmix64(self.total as u64) as usize % self.total;
            if mixed >= MAX_OVERVIEW_POINTS {
                self.truncated = true;
                return;
            }
            self.truncated = true;
            mixed
        };

        if slot == self.mz.len() {
            self.mz.push(point.mz);
            self.rt.push(point.rt);
            self.half_width.push(point.half_width);
            self.coefficient.push(point.coefficient);
            self.shape.push(point.shape);
            self.smoothness.push(point.smoothness);
            self.group.push(point.group);
            self.names.push(point.name);
            self.has_ms2.push(point.has_ms2);
        } else {
            self.mz[slot] = point.mz;
            self.rt[slot] = point.rt;
            self.half_width[slot] = point.half_width;
            self.coefficient[slot] = point.coefficient;
            self.shape[slot] = point.shape;
            self.smoothness[slot] = point.smoothness;
            self.group[slot] = point.group;
            self.names[slot] = point.name;
            self.has_ms2[slot] = point.has_ms2;
        }
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

#[tauri::command]
pub fn visualizer_overview(output_dir: String, sample: String) -> Result<VizOverview, String> {
    validate_sample(&sample)?;
    let output = checked_output(&output_dir)?;
    let misc = output.join("misc");
    let settings = read_settings(&misc);
    let ms2_path = misc.join(ms2_cache_name(&sample));
    let ms2 = Ms2View::open(&ms2_path).map_err(|e| e.to_string())?;
    let mut report = read_report_points(&output.join("Report_RTseparated.csv"), &sample)?;
    report.sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));

    let feature_path = misc.join(feature_cache_name(&sample));
    let mut overview = VizOverview::default();
    let precursors = ms2.precursors();
    let spectrum_rts = ms2.rts();
    if feature_path.is_file() {
        let features = FeatureView::open(&feature_path).map_err(|e| e.to_string())?;
        let mzs = features.mzs();
        let rts = features.rts();
        let half_widths = features.half_widths();
        let coefficients = features.coefficients();
        let shapes = features.shapes();
        let smoothness = features.smoothness();
        overview.source = "feature cache".to_string();
        for i in 0..features.len() {
            let mz = mzs[i];
            let rt = rts[i];
            if !nearby_ms2(precursors, spectrum_rts, mz, rt, settings.ms1_ms_pair, 0.5) {
                continue;
            }
            let matched = nearest_report(&report, mz, rt);
            overview.consider(Point {
                mz,
                rt,
                half_width: half_widths[i],
                coefficient: coefficients[i],
                shape: shapes[i],
                smoothness: smoothness[i],
                group: matched.map(|p| p.group).unwrap_or(0),
                name: matched.map(|p| p.name.clone()).unwrap_or_default(),
                has_ms2: true,
            });
        }
    } else {
        overview.source = "report fallback".to_string();
        for point in report {
            let has_ms2 = nearby_ms2(
                precursors,
                spectrum_rts,
                point.mz,
                point.rt,
                settings.ms1_ms_pair,
                0.5,
            );
            overview.consider(Point {
                mz: point.mz,
                rt: point.rt,
                half_width: 0.1,
                coefficient: 0.0,
                shape: point.shape,
                smoothness: point.smoothness,
                group: point.group,
                name: point.name,
                has_ms2,
            });
        }
    }
    Ok(overview)
}

fn nearby_ms2(precursors: &[f32], rts: &[f32], mz: f32, rt: f32, mz_tol: f32, rt_tol: f32) -> bool {
    let start = precursors.partition_point(|&value| value < mz - mz_tol);
    for i in start..precursors.len() {
        if precursors[i] >= mz + mz_tol {
            break;
        }
        if (rts[i] - rt).abs() < rt_tol {
            return true;
        }
    }
    false
}

fn nearest_report(points: &[ReportPoint], mz: f32, rt: f32) -> Option<&ReportPoint> {
    let start = points.partition_point(|point| point.mz < mz - 0.02);
    points[start..]
        .iter()
        .take_while(|point| point.mz < mz + 0.02)
        .filter(|point| (point.rt - rt).abs() < 0.05)
        .min_by(|a, b| {
            let da = (a.mz - mz).abs() * 10.0 + (a.rt - rt).abs();
            let db = (b.mz - mz).abs() * 10.0 + (b.rt - rt).abs();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
}

fn read_report_points(path: &Path, sample: &str) -> Result<Vec<ReportPoint>, String> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_path(path)
        .map_err(|e| format!("could not read {}: {e}", path.display()))?;
    let mut columns = None;
    let mut out = Vec::new();

    for record in reader.records() {
        let record = record.map_err(|e| format!("{}: {e}", path.display()))?;
        if columns.is_none() {
            if record.get(0) == Some("group") {
                let find = |name: &str| record.iter().position(|value| value == name);
                columns = Some((
                    find("group").unwrap_or(0),
                    find("feature_m/z").ok_or("report is missing feature_m/z")?,
                    find(&format!("RT_{sample}"))
                        .or_else(|| find("Median RT"))
                        .ok_or("report is missing retention time")?,
                    find("peak_shape (median)"),
                    find(&format!("S/N_{sample}")),
                    find("name").ok_or("report is missing compound names")?,
                ));
            }
            continue;
        }

        let (group_at, mz_at, rt_at, shape_at, smooth_at, name_at) = columns.unwrap();
        let Some(mz) = record.get(mz_at).and_then(|v| v.parse::<f32>().ok()) else {
            continue;
        };
        let Some(rt) = record.get(rt_at).and_then(|v| v.parse::<f32>().ok()) else {
            continue;
        };
        let name = record.get(name_at).unwrap_or("").trim();
        if name.starts_with("ISF of ") {
            continue;
        }
        out.push(ReportPoint {
            mz,
            rt,
            shape: shape_at
                .and_then(|at| record.get(at))
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0),
            smoothness: smooth_at
                .and_then(|at| record.get(at))
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0),
            group: record
                .get(group_at)
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
            name: name.to_string(),
        });
    }
    columns.ok_or_else(|| "report header was not found".to_string())?;
    Ok(out)
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizSpectrum {
    id: usize,
    precursor_mz: f32,
    rt: f32,
    ce: f32,
    peaks: Vec<[f32; 2]>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizPeak {
    mz: f32,
    intensity: f32,
    matched: bool,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizMirror {
    name: String,
    library_mz: f32,
    experimental_mz: f32,
    rt: f32,
    ce: f32,
    score: f32,
    shape: f32,
    experimental: Vec<VizPeak>,
    library: Vec<VizPeak>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VizFeatureDetail {
    chromatogram: Vec<[f32; 2]>,
    spectra: Vec<VizSpectrum>,
    spectra_truncated: bool,
    mirrors: Vec<VizMirror>,
    integration_half_width: f32,
}

#[tauri::command]
pub fn visualizer_feature(
    output_dir: String,
    sample: String,
    mz: f32,
    rt: f32,
    half_width: f32,
    state: State<'_, VisualizerState>,
) -> Result<VizFeatureDetail, String> {
    visualizer_feature_data(output_dir, sample, mz, rt, half_width, &state)
}

pub fn visualizer_feature_data(
    output_dir: String,
    sample: String,
    mz: f32,
    rt: f32,
    half_width: f32,
    state: &VisualizerState,
) -> Result<VizFeatureDetail, String> {
    validate_sample(&sample)?;
    if !mz.is_finite() || !rt.is_finite() {
        return Err("Feature coordinates are invalid.".to_string());
    }
    let output = checked_output(&output_dir)?;
    let misc = output.join("misc");
    let settings = read_settings(&misc);
    let window = 0.5f32
        .max(half_width * settings.integration_rt * 1.5)
        .min(5.0);

    let ms1 = Ms1View::open(&misc.join(ms1_cache_name(&sample))).map_err(|e| e.to_string())?;
    let ms2 = Ms2View::open(&misc.join(ms2_cache_name(&sample))).map_err(|e| e.to_string())?;

    let mut trace = Vec::new();
    ms1.xic(mz, rt, window, settings.integration_mz, &mut trace);
    let chromatogram = downsample_trace(&trace, 2_048);

    let mut candidates = Vec::new();
    let precursors = ms2.precursors();
    let spectrum_rts = ms2.rts();
    let collision_energies = ms2.ces();
    let start = precursors.partition_point(|&value| value < mz - settings.ms1_ms_pair);
    for i in start..ms2.len() {
        if precursors[i] >= mz + settings.ms1_ms_pair {
            break;
        }
        if (spectrum_rts[i] - rt).abs() < window {
            candidates.push(i);
        }
    }
    let spectra_truncated = candidates.len() > MAX_SPECTRA;
    if spectra_truncated {
        candidates.sort_unstable_by(|&a, &b| {
            let da = (spectrum_rts[a] - rt).abs()
                + (precursors[a] - mz).abs() / settings.ms1_ms_pair.max(1e-6);
            let db = (spectrum_rts[b] - rt).abs()
                + (precursors[b] - mz).abs() / settings.ms1_ms_pair.max(1e-6);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
        candidates.truncate(MAX_SPECTRA);
    }
    candidates.sort_unstable_by(|&a, &b| {
        spectrum_rts[a]
            .partial_cmp(&spectrum_rts[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let spectra = candidates
        .into_iter()
        .map(|i| {
            let (mzs, intensities) = ms2.scan(i);
            VizSpectrum {
                id: i,
                precursor_mz: precursors[i],
                rt: spectrum_rts[i],
                ce: collision_energies[i],
                peaks: compact_peaks(mzs, intensities, MAX_SPECTRUM_PEAKS),
            }
        })
        .collect();

    let mirrors = {
        let mut guard = state
            .match_index
            .lock()
            .map_err(|_| "visualizer match index is unavailable".to_string())?;
        let rebuild = guard
            .as_ref()
            .map(|index| !index.is_current(&output))
            .unwrap_or(true);
        if rebuild {
            *guard = MatchIndex::build(&output).ok();
        }
        match guard.as_ref() {
            Some(index) => index.matches(&sample, mz, rt)?,
            None => Vec::new(),
        }
    };

    Ok(VizFeatureDetail {
        chromatogram,
        spectra,
        spectra_truncated,
        mirrors,
        integration_half_width: half_width * settings.integration_rt,
    })
}

fn compact_peaks(mzs: &[f32], intensities: &[f32], max: usize) -> Vec<[f32; 2]> {
    let mut indices: Vec<usize> = (0..mzs.len().min(intensities.len())).collect();
    if indices.len() > max {
        indices.select_nth_unstable_by(max, |&a, &b| {
            intensities[b]
                .partial_cmp(&intensities[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        indices.truncate(max);
        indices.sort_unstable_by(|&a, &b| {
            mzs[a]
                .partial_cmp(&mzs[b])
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    indices
        .into_iter()
        .map(|i| [mzs[i], intensities[i]])
        .collect()
}

fn downsample_trace(points: &[(f32, f32)], max: usize) -> Vec<[f32; 2]> {
    if points.len() <= max {
        return points.iter().map(|&(x, y)| [x, y]).collect();
    }
    let buckets = (max / 2).max(1);
    let width = points.len() as f64 / buckets as f64;
    let mut out = Vec::with_capacity(max + 2);
    for bucket in 0..buckets {
        let from = (bucket as f64 * width) as usize;
        let to = (((bucket + 1) as f64 * width) as usize).min(points.len());
        if from >= to {
            continue;
        }
        let mut low = from;
        let mut high = from;
        for i in from + 1..to {
            if points[i].1 < points[low].1 {
                low = i;
            }
            if points[i].1 > points[high].1 {
                high = i;
            }
        }
        if low <= high {
            out.push([points[low].0, points[low].1]);
            if high != low {
                out.push([points[high].0, points[high].1]);
            }
        } else {
            out.push([points[high].0, points[high].1]);
            out.push([points[low].0, points[low].1]);
        }
    }
    out
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct FileStamp {
    len: u64,
    modified_ns: u128,
}

fn stamp(path: &Path) -> Option<FileStamp> {
    let meta = fs::metadata(path).ok()?;
    let modified_ns = meta
        .modified()
        .ok()
        .and_then(|time| time.duration_since(UNIX_EPOCH).ok())
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    Some(FileStamp {
        len: meta.len(),
        modified_ns,
    })
}

struct IndexedSpectrum {
    peaks_at: u64,
    peak_count: usize,
    name_id: u32,
    sample_id: u32,
    library_mz: f32,
    experimental_mz: f32,
    rt: f32,
    ce: f32,
    score: f32,
    shape: f32,
}

#[derive(Clone, Copy)]
struct PeakBlock {
    peaks_at: u64,
    peak_count: usize,
}

struct MatchIndex {
    output: PathBuf,
    exp_stamp: FileStamp,
    lib_stamp: FileStamp,
    strings: Vec<String>,
    experimental: Vec<IndexedSpectrum>,
    library: Vec<PeakBlock>,
}

impl MatchIndex {
    fn build(output: &Path) -> Result<MatchIndex, String> {
        let exp_path = output.join("spec_exp.txt");
        let lib_path = output.join("spec_lib.txt");
        let exp_stamp = stamp(&exp_path).ok_or("experimental spectra are missing")?;
        let lib_stamp = stamp(&lib_path).ok_or("library spectra are missing")?;
        let mut interner = Interner::default();
        let experimental = index_experimental(&exp_path, &mut interner)?;
        let library = index_peak_blocks(&lib_path)?;
        Ok(MatchIndex {
            output: output.to_path_buf(),
            exp_stamp,
            lib_stamp,
            strings: interner.values,
            experimental,
            library,
        })
    }

    fn is_current(&self, output: &Path) -> bool {
        self.output == output
            && stamp(&output.join("spec_exp.txt")) == Some(self.exp_stamp)
            && stamp(&output.join("spec_lib.txt")) == Some(self.lib_stamp)
    }

    fn matches(&self, sample: &str, mz: f32, rt: f32) -> Result<Vec<VizMirror>, String> {
        let mut selected: Vec<usize> = self
            .experimental
            .iter()
            .enumerate()
            .filter(|(_, entry)| {
                self.strings
                    .get(entry.sample_id as usize)
                    .is_some_and(|value| value == sample)
                    && (entry.experimental_mz - mz).abs() < 0.02
                    && (entry.rt - rt).abs() < 0.035
            })
            .map(|(i, _)| i)
            .collect();
        selected.sort_unstable_by(|&a, &b| {
            self.experimental[b]
                .score
                .partial_cmp(&self.experimental[a].score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        selected.truncate(MAX_MIRRORS);

        let mut exp_reader = BufReader::new(
            File::open(self.output.join("spec_exp.txt"))
                .map_err(|e| format!("could not open experimental spectra: {e}"))?,
        );
        let mut lib_reader = BufReader::new(
            File::open(self.output.join("spec_lib.txt"))
                .map_err(|e| format!("could not open library spectra: {e}"))?,
        );
        let mut out = Vec::with_capacity(selected.len());
        for i in selected {
            let Some(lib) = self.library.get(i) else {
                continue;
            };
            let exp = &self.experimental[i];
            out.push(VizMirror {
                name: self
                    .strings
                    .get(exp.name_id as usize)
                    .cloned()
                    .unwrap_or_default(),
                library_mz: exp.library_mz,
                experimental_mz: exp.experimental_mz,
                rt: exp.rt,
                ce: exp.ce,
                score: exp.score,
                shape: exp.shape,
                experimental: read_peaks(&mut exp_reader, exp.peaks_at, exp.peak_count)?,
                library: read_peaks(&mut lib_reader, lib.peaks_at, lib.peak_count)?,
            });
        }
        Ok(out)
    }
}

#[derive(Default)]
struct Interner {
    ids: HashMap<String, u32>,
    values: Vec<String>,
}

impl Interner {
    fn intern(&mut self, value: String) -> u32 {
        if let Some(&id) = self.ids.get(&value) {
            return id;
        }
        let id = self.values.len() as u32;
        self.values.push(value.clone());
        self.ids.insert(value, id);
        id
    }
}

#[derive(Default)]
struct PendingSpectrum {
    name: String,
    sample: String,
    library_mz: f32,
    mass_diff: f32,
    rt: f32,
    ce: f32,
    score: f32,
    shape: f32,
}

fn index_experimental(
    path: &Path,
    interner: &mut Interner,
) -> Result<Vec<IndexedSpectrum>, String> {
    let file = File::open(path).map_err(|e| format!("could not open {}: {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut offset = 0u64;
    let mut pending = PendingSpectrum::default();
    let mut out = Vec::new();

    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|e| format!("could not index {}: {e}", path.display()))?;
        if read == 0 {
            break;
        }
        offset += read as u64;
        let value = line.trim();
        if let Some(v) = value.strip_prefix("NAME: ") {
            pending = PendingSpectrum::default();
            pending.name = v.to_string();
        } else if let Some(v) = value.strip_prefix("SAMPLE: ") {
            pending.sample = v.to_string();
        } else if let Some(v) = value.strip_prefix("PRECURSORMZ: ") {
            pending.library_mz = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("MASS_DIFF(LIB-EXP): ") {
            pending.mass_diff = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("RETENTIONTIME: ") {
            pending.rt = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("COLLISIONENERGY: ") {
            pending.ce = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("SCORE: ") {
            pending.score = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("SHAPE: ") {
            pending.shape = v.parse().unwrap_or(0.0);
        } else if let Some(v) = value.strip_prefix("Num Peaks: ") {
            let peak_count = v.parse().unwrap_or(0);
            out.push(IndexedSpectrum {
                peaks_at: offset,
                peak_count,
                name_id: interner.intern(std::mem::take(&mut pending.name)),
                sample_id: interner.intern(std::mem::take(&mut pending.sample)),
                library_mz: pending.library_mz,
                experimental_mz: pending.library_mz - pending.mass_diff,
                rt: pending.rt,
                ce: pending.ce,
                score: pending.score,
                shape: pending.shape,
            });
        }
    }
    Ok(out)
}

fn index_peak_blocks(path: &Path) -> Result<Vec<PeakBlock>, String> {
    let file = File::open(path).map_err(|e| format!("could not open {}: {e}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let mut offset = 0u64;
    let mut out = Vec::new();
    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|e| format!("could not index {}: {e}", path.display()))?;
        if read == 0 {
            break;
        }
        offset += read as u64;
        if let Some(value) = line.trim().strip_prefix("Num Peaks: ") {
            out.push(PeakBlock {
                peaks_at: offset,
                peak_count: value.parse().unwrap_or(0),
            });
        }
    }
    Ok(out)
}

fn read_peaks(
    reader: &mut BufReader<File>,
    offset: u64,
    count: usize,
) -> Result<Vec<VizPeak>, String> {
    reader
        .seek(SeekFrom::Start(offset))
        .map_err(|e| format!("could not seek spectrum: {e}"))?;
    let mut out = Vec::with_capacity(count.min(MAX_SPECTRUM_PEAKS));
    let mut line = String::new();
    for _ in 0..count {
        line.clear();
        if reader
            .read_line(&mut line)
            .map_err(|e| format!("could not read spectrum: {e}"))?
            == 0
        {
            break;
        }
        let mut fields = line.split_whitespace();
        let Some(mz) = fields.next().and_then(|v| v.parse().ok()) else {
            continue;
        };
        let Some(intensity) = fields.next().and_then(|v| v.parse().ok()) else {
            continue;
        };
        let matched = fields.any(|value| value == "*");
        out.push(VizPeak {
            mz,
            intensity,
            matched,
        });
    }
    if out.len() > MAX_SPECTRUM_PEAKS {
        out.sort_unstable_by(|a, b| {
            b.intensity
                .partial_cmp(&a.intensity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        out.truncate(MAX_SPECTRUM_PEAKS);
        out.sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(std::cmp::Ordering::Equal));
    }
    Ok(out)
}
