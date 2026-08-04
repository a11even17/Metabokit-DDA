//! Dataset discovery.
//!
//! Given one folder, work out everything a run needs: which files are samples,
//! what polarity they were acquired in, which libraries are available, whether
//! a 0.1-era `param.txt` is sitting alongside them, and where results should
//! go.
//!
//! The scan reports *what it found and where it found it*, not just a filled-in
//! settings object. Auto-detection that cannot be audited is worse than no
//! auto-detection, because a wrong guess becomes a silently wrong result three
//! hours later.

use std::cmp::Ordering;
use std::collections::HashSet;
use std::fs::File;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};
use crate::library;
use crate::mzml;
use crate::params::{LibrarySource, Params, Polarity, BUILTIN_LIBRARIES};
use crate::report::ReportSummary;

/// Directory names never worth descending into.
const SKIP_DIRS: [&str; 6] = [
    "misc",
    "target",
    "node_modules",
    "results",
    "$RECYCLE.BIN",
    "System Volume Information",
];

/// How deep below the chosen folder to look for samples.
const MAX_DEPTH: usize = 3;

/// Default results folder, created inside the dataset folder.
pub const RESULTS_DIR: &str = "results";

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum NoteLevel {
    Ok,
    Warn,
    Blocked,
}

/// One auditable statement about what the scan concluded.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Note {
    pub level: NoteLevel,
    /// Short mono label, e.g. "samples", "polarity".
    pub topic: String,
    pub message: String,
}

impl Note {
    fn ok(topic: &str, message: impl Into<String>) -> Self {
        Note {
            level: NoteLevel::Ok,
            topic: topic.into(),
            message: message.into(),
        }
    }
    fn warn(topic: &str, message: impl Into<String>) -> Self {
        Note {
            level: NoteLevel::Warn,
            topic: topic.into(),
            message: message.into(),
        }
    }
    fn blocked(topic: &str, message: impl Into<String>) -> Self {
        Note {
            level: NoteLevel::Blocked,
            topic: topic.into(),
            message: message.into(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SampleEntry {
    pub name: String,
    pub path: String,
    pub bytes: u64,
    /// Path relative to the dataset root, when it sits in a subfolder.
    pub subfolder: Option<String>,
}

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DetectedLibrary {
    pub kind: String,
    pub label: String,
    pub detail: String,
}

/// Summary recovered from an existing, completed `results/` folder.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct PreviousRun {
    pub samples: usize,
    pub summary: ReportSummary,
    pub polarity: String,
    pub reports: Vec<String>,
    pub caches_present: bool,
    pub output_dir: String,
}

/// The result of scanning a dataset folder.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct DatasetScan {
    pub root: String,
    /// Fully populated and ready to run, subject to `notes`.
    pub params: Params,
    pub samples: Vec<SampleEntry>,
    pub total_bytes: u64,
    pub polarity: String,
    pub libraries: Vec<DetectedLibrary>,
    pub libs_dir: Option<String>,
    pub output_dir: String,
    /// A valid prior run found in the default output folder, if any.
    pub previous_run: Option<PreviousRun>,
    /// A 0.1 `param.txt` that was imported, if any.
    pub imported_settings: Option<String>,
    pub notes: Vec<Note>,
    /// True when the scan produced something that can actually be run.
    pub ready: bool,
}

/// Compare file names the way a person would: `PosIDA-2` before `PosIDA-10`.
fn natural_cmp(a: &str, b: &str) -> Ordering {
    let (ab, bb) = (a.as_bytes(), b.as_bytes());
    let (mut i, mut j) = (0usize, 0usize);
    while i < ab.len() && j < bb.len() {
        if ab[i].is_ascii_digit() && bb[j].is_ascii_digit() {
            let si = i;
            while i < ab.len() && ab[i].is_ascii_digit() {
                i += 1;
            }
            let sj = j;
            while j < bb.len() && bb[j].is_ascii_digit() {
                j += 1;
            }
            // Compare as numbers: longer digit run wins unless leading zeros
            // make them equal length.
            let na = a[si..i].trim_start_matches('0');
            let nb = b[sj..j].trim_start_matches('0');
            match na.len().cmp(&nb.len()).then_with(|| na.cmp(nb)) {
                Ordering::Equal => {}
                other => return other,
            }
        } else {
            match ab[i].to_ascii_lowercase().cmp(&bb[j].to_ascii_lowercase()) {
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                other => return other,
            }
        }
    }
    (ab.len() - i).cmp(&(bb.len() - j))
}

fn has_extension(path: &Path, wanted: &[&str]) -> bool {
    path.extension()
        .and_then(|x| x.to_str())
        .map(|e| wanted.iter().any(|w| e.eq_ignore_ascii_case(w)))
        .unwrap_or(false)
}

fn is_hidden(name: &str) -> bool {
    // Covers `.DS_Store`, `._resource-forks` and dot-directories in one test.
    name.starts_with('.')
}

/// Files found during the walk, bucketed by what they look like.
#[derive(Default)]
struct Walked {
    mzml: Vec<PathBuf>,
    msp: Vec<PathBuf>,
    csv: Vec<PathBuf>,
    param_txt: Option<PathBuf>,
    libs_dir: Option<PathBuf>,
}

fn walk(dir: &Path, depth: usize, out: &mut Walked) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Some(name) = path.file_name().and_then(|x| x.to_str()) else {
            continue;
        };
        if is_hidden(name) {
            continue;
        }
        let is_dir = entry.file_type().map(|t| t.is_dir()).unwrap_or(false);

        if is_dir {
            if name.eq_ignore_ascii_case("libs") {
                if out.libs_dir.is_none() {
                    out.libs_dir = Some(path.clone());
                }
                continue;
            }
            if SKIP_DIRS.iter().any(|s| name.eq_ignore_ascii_case(s)) {
                continue;
            }
            if depth < MAX_DEPTH {
                walk(&path, depth + 1, out);
            }
            continue;
        }

        if has_extension(&path, &["mzML", "mzXML"]) {
            out.mzml.push(path);
        } else if has_extension(&path, &["msp"]) {
            out.msp.push(path);
        } else if has_extension(&path, &["csv"]) {
            out.csv.push(path);
        } else if name.eq_ignore_ascii_case("param.txt") && out.param_txt.is_none() {
            out.param_txt = Some(path);
        }
    }
}

/// Does this CSV look like a `name, adduct, m/z, rt` library rather than, say,
/// a sample manifest or an exported report?
fn looks_like_library_csv(path: &Path) -> bool {
    use std::io::{BufRead, BufReader};
    let Ok(file) = std::fs::File::open(path) else {
        return false;
    };
    let mut header = String::new();
    if BufReader::new(file).read_line(&mut header).is_err() {
        return false;
    }
    let header = header.to_ascii_lowercase();
    let fields: Vec<&str> = header.split(',').map(str::trim).collect();
    if fields.len() < 3 {
        return false;
    }
    let has_name = fields.iter().any(|f| f.contains("name"));
    let has_adduct = fields
        .iter()
        .any(|f| f.contains("adduct") || f.contains("precursortype"));
    let has_mz = fields
        .iter()
        .any(|f| *f == "mz" || f.contains("m/z") || f.contains("precursormz"));
    has_adduct && (has_name || has_mz)
}

/// Scan a dataset folder and build a ready-to-run configuration.
pub fn scan(root: &Path) -> Result<DatasetScan> {
    if !root.is_dir() {
        return Err(Error::param(format!("{} is not a folder", root.display())));
    }

    let mut found = Walked::default();
    walk(root, 0, &mut found);

    let mut notes: Vec<Note> = Vec::new();

    // ---- settings ---------------------------------------------------------
    // A 0.1 `param.txt` alongside the data is the most authoritative source of
    // numeric settings we can have, so it wins over the defaults.
    let mut params;
    let mut imported_settings = None;
    match &found.param_txt {
        Some(p) => match std::fs::read_to_string(p)
            .ok()
            .and_then(|text| Params::from_legacy(&text, p).ok())
        {
            Some(imported) => {
                params = imported;
                imported_settings = Some(p.to_string_lossy().into_owned());
                notes.push(Note::ok(
                    "settings",
                    format!(
                        "imported analysis settings from {}",
                        display_relative(root, p)
                    ),
                ));
            }
            None => {
                params = Params::default();
                notes.push(Note::warn(
                    "settings",
                    format!(
                        "{} could not be parsed; using defaults",
                        display_relative(root, p)
                    ),
                ));
            }
        },
        None => params = Params::default(),
    }

    // ---- samples ----------------------------------------------------------
    // A `param.txt` may have brought its own file list (via `file_order.txt`).
    // Trust it only if those files actually exist.
    let ordered_from_settings: Vec<PathBuf> = params
        .mzml_files
        .iter()
        .filter(|p| p.is_file())
        .cloned()
        .collect();

    let mut mzml =
        if ordered_from_settings.len() >= found.mzml.len() && !ordered_from_settings.is_empty() {
            notes.push(Note::ok(
                "samples",
                "sample order taken from file_order.txt",
            ));
            ordered_from_settings
        } else {
            found.mzml.clone()
        };
    if mzml.len() == found.mzml.len() {
        mzml.sort_by(|a, b| {
            natural_cmp(
                &a.file_name()
                    .map(|x| x.to_string_lossy().into_owned())
                    .unwrap_or_default(),
                &b.file_name()
                    .map(|x| x.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            )
        });
    }

    let samples: Vec<SampleEntry> = mzml
        .iter()
        .map(|p| SampleEntry {
            name: p
                .file_name()
                .map(|x| x.to_string_lossy().into_owned())
                .unwrap_or_default(),
            bytes: std::fs::metadata(p).map(|m| m.len()).unwrap_or(0),
            subfolder: p
                .parent()
                .filter(|d| *d != root)
                .map(|d| display_relative(root, d)),
            path: p.to_string_lossy().into_owned(),
        })
        .collect();
    let total_bytes: u64 = samples.iter().map(|s| s.bytes).sum();

    if samples.is_empty() {
        notes.push(Note::blocked(
            "samples",
            "no .mzML or .mzXML files found in this folder or its subfolders",
        ));
    } else {
        notes.push(Note::ok(
            "samples",
            format!(
                "{} sample{} found, {:.0} MB total",
                samples.len(),
                if samples.len() == 1 { "" } else { "s" },
                total_bytes as f64 / 1e6
            ),
        ));
    }
    params.mzml_files = mzml;

    // ---- polarity ---------------------------------------------------------
    let polarity = match params.mzml_files.first() {
        Some(first) => match mzml::sniff_polarity(first) {
            Ok(Some(true)) => {
                params.polarity = Polarity::Positive;
                notes.push(Note::ok(
                    "polarity",
                    "positive mode, read from the first sample",
                ));
                "positive"
            }
            Ok(Some(false)) => {
                params.polarity = Polarity::Negative;
                notes.push(Note::ok(
                    "polarity",
                    "negative mode, read from the first sample",
                ));
                "negative"
            }
            _ => {
                params.polarity = Polarity::Auto;
                notes.push(Note::warn(
                    "polarity",
                    "the samples do not declare a scan polarity; set it manually if positive is wrong",
                ));
                "unknown"
            }
        },
        None => "unknown",
    };

    // ---- libraries --------------------------------------------------------
    let libs_dir = found
        .libs_dir
        .clone()
        .or_else(|| {
            let candidate = root.join("libs");
            candidate.is_dir().then_some(candidate)
        })
        .or_else(|| root.parent().map(|p| p.join("libs")).filter(|p| p.is_dir()))
        .or_else(|| params.resolve_libs_dir());

    params.libs_dir = libs_dir.clone();
    let mut libraries: Vec<DetectedLibrary> = Vec::new();
    let mut sources: Vec<LibrarySource> = Vec::new();
    let want_positive = polarity != "negative";

    if let Some(dir) = &libs_dir {
        let installed = library::available_builtins(Some(dir));
        for (name, has_pos, has_neg) in installed {
            let usable = if want_positive { has_pos } else { has_neg };
            // `Atlas_filtered` is a narrowed variant of `Atlas`; selecting both
            // would double every lipid entry.
            if usable && name != "Atlas_filtered" {
                libraries.push(DetectedLibrary {
                    kind: "builtin".into(),
                    label: name.clone(),
                    detail: format!(
                        "{} mode",
                        if want_positive {
                            "positive"
                        } else {
                            "negative"
                        }
                    ),
                });
                sources.push(LibrarySource::Builtin(name));
            }
        }
        if sources.is_empty() {
            notes.push(Note::warn(
                "libraries",
                format!(
                    "{} holds no library for {} mode",
                    display_relative(root, dir),
                    if want_positive {
                        "positive"
                    } else {
                        "negative"
                    }
                ),
            ));
        }
    }

    for path in &found.msp {
        libraries.push(DetectedLibrary {
            kind: "msp".into(),
            label: path
                .file_name()
                .map(|x| x.to_string_lossy().into_owned())
                .unwrap_or_default(),
            detail: display_relative(root, path),
        });
        sources.push(LibrarySource::Msp(path.clone()));
    }
    for path in &found.csv {
        if !looks_like_library_csv(path) {
            continue;
        }
        libraries.push(DetectedLibrary {
            kind: "csv".into(),
            label: path
                .file_name()
                .map(|x| x.to_string_lossy().into_owned())
                .unwrap_or_default(),
            detail: display_relative(root, path),
        });
        sources.push(LibrarySource::Csv(path.clone()));
    }

    if sources.is_empty() {
        notes.push(Note::blocked(
            "libraries",
            "no spectral library found — add an .msp file, or point the app at a libs/ folder",
        ));
    } else {
        notes.push(Note::ok(
            "libraries",
            format!(
                "{} librar{} ready",
                sources.len(),
                if sources.len() == 1 { "y" } else { "ies" }
            ),
        ));
    }
    params.libraries = sources;

    // ---- output -----------------------------------------------------------
    // Never write into the dataset folder itself: a re-scan would then pick up
    // our own outputs as inputs.
    let output_dir = root.join(RESULTS_DIR);
    params.output_dir = output_dir.clone();
    let previous_run = detect_previous_run(&output_dir, samples.len(), polarity);
    if previous_run.is_some() {
        notes.push(Note::ok(
            "output",
            format!("completed run found in {RESULTS_DIR}/"),
        ));
    } else if output_dir.is_dir() {
        notes.push(Note::warn(
            "output",
            format!("{RESULTS_DIR}/ exists but does not contain a complete run"),
        ));
    } else {
        notes.push(Note::ok(
            "output",
            format!("reports will be written to {RESULTS_DIR}/"),
        ));
    }

    let ready = !notes.iter().any(|n| n.level == NoteLevel::Blocked);

    Ok(DatasetScan {
        root: root.to_string_lossy().into_owned(),
        samples,
        total_bytes,
        polarity: polarity.to_string(),
        libraries,
        libs_dir: libs_dir.map(|p| p.to_string_lossy().into_owned()),
        output_dir: output_dir.to_string_lossy().into_owned(),
        previous_run,
        imported_settings,
        notes,
        ready,
        params,
    })
}

#[derive(Deserialize)]
struct SavedRunManifest {
    #[serde(default)]
    polarity: String,
    #[serde(default)]
    samples: Vec<String>,
}

/// Only call a folder completed when both primary reports are readable. The
/// manifest enriches the restored state when present, while older result sets
/// can fall back to the currently detected samples and polarity.
fn detect_previous_run(
    output_dir: &Path,
    fallback_samples: usize,
    fallback_polarity: &str,
) -> Option<PreviousRun> {
    let manifest_path = output_dir.join("misc").join("run.json");
    let manifest: Option<SavedRunManifest> = File::open(manifest_path)
        .ok()
        .and_then(|file| serde_json::from_reader(file).ok());
    let summary = summarize_reports(
        &output_dir.join("Report_RTseparated.csv"),
        &output_dir.join("Report_by_ID.csv"),
    )?;

    let report_names = [
        "Report_RTseparated.csv",
        "Report_by_ID.csv",
        "Report_RTseparated_fill.csv",
        "Report_by_ID_fill.csv",
        "spec_exp.txt",
        "spec_lib.txt",
        "spec_reduced.txt",
    ];
    let reports = report_names
        .iter()
        .filter(|name| output_dir.join(name).is_file())
        .map(|name| (*name).to_string())
        .collect();

    let misc = output_dir.join("misc");
    let mut has_ms1 = false;
    let mut has_ms2 = false;
    if let Ok(entries) = std::fs::read_dir(&misc) {
        for entry in entries.flatten() {
            let name = entry.file_name();
            let name = name.to_string_lossy();
            has_ms1 |= name.starts_with("ms1_") && name.ends_with(".mkc");
            has_ms2 |= name.starts_with("ms2_") && name.ends_with(".mkc");
        }
    }

    Some(PreviousRun {
        samples: manifest
            .as_ref()
            .map(|saved| saved.samples.len())
            .filter(|count| *count > 0)
            .unwrap_or(fallback_samples),
        summary,
        polarity: manifest
            .as_ref()
            .map(|saved| saved.polarity.as_str())
            .filter(|value| !value.is_empty())
            .unwrap_or(fallback_polarity)
            .to_string(),
        reports,
        caches_present: has_ms1 && has_ms2,
        output_dir: output_dir.to_string_lossy().into_owned(),
    })
}

fn summarize_reports(by_rt: &Path, by_id: &Path) -> Option<ReportSummary> {
    let mut feature_groups = HashSet::new();
    let mut isf_rows = 0usize;
    read_report(by_rt, |group, isf, _| {
        feature_groups.insert(group.to_string());
        if !isf.is_empty() {
            isf_rows += 1;
        }
    })?;

    let mut identified_groups = HashSet::new();
    let mut compounds = HashSet::new();
    read_report(by_id, |group, _, name| {
        if !name.is_empty() {
            identified_groups.insert(group.to_string());
            compounds.insert(name.to_string());
        }
    })?;

    Some(ReportSummary {
        features: feature_groups.len(),
        identified: identified_groups.len(),
        isf_rows,
        compounds: compounds.len(),
    })
}

/// Reports contain a timestamp preamble before their real CSV header. Locate
/// that header instead of assuming it is the first record so older output is
/// restored as reliably as newly-written output.
fn read_report(path: &Path, mut visit: impl FnMut(&str, &str, &str)) -> Option<()> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .from_path(path)
        .ok()?;
    let mut columns = None;

    for record in reader.records() {
        let record = record.ok()?;
        if columns.is_none() {
            if record.get(0) == Some("group") {
                let group = 0;
                let isf = record.iter().position(|value| value == "ISF")?;
                let name = record.iter().position(|value| value == "name")?;
                columns = Some((group, isf, name));
            }
            continue;
        }

        let (group_at, isf_at, name_at) = columns?;
        let group = record.get(group_at).unwrap_or("").trim();
        if group.is_empty() {
            continue;
        }
        visit(
            group,
            record.get(isf_at).unwrap_or("").trim(),
            record.get(name_at).unwrap_or("").trim(),
        );
    }

    columns.map(|_| ())
}

/// Path relative to the dataset root where possible, for compact display.
fn display_relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .map(|p| p.to_string_lossy().into_owned())
        .unwrap_or_else(|_| path.to_string_lossy().into_owned())
}

/// Re-run library detection after the user supplies a libs directory or an
/// individual file, without re-walking the whole dataset.
pub fn refresh_libraries(params: &mut Params) -> Vec<DetectedLibrary> {
    let want_positive = params.polarity != Polarity::Negative;
    let mut libraries = Vec::new();
    let mut builtins = Vec::new();

    if let Some(dir) = params.resolve_libs_dir() {
        for (name, has_pos, has_neg) in library::available_builtins(Some(&dir)) {
            let usable = if want_positive { has_pos } else { has_neg };
            if usable && name != "Atlas_filtered" {
                libraries.push(DetectedLibrary {
                    kind: "builtin".into(),
                    label: name.clone(),
                    detail: format!(
                        "{} mode",
                        if want_positive {
                            "positive"
                        } else {
                            "negative"
                        }
                    ),
                });
                builtins.push(LibrarySource::Builtin(name));
            }
        }
    }

    // Keep whatever user libraries are already configured.
    let user: Vec<LibrarySource> = params
        .libraries
        .iter()
        .filter(|l| !matches!(l, LibrarySource::Builtin(_)))
        .cloned()
        .collect();
    for lib in &user {
        libraries.push(DetectedLibrary {
            kind: match lib {
                LibrarySource::Msp(_) => "msp".into(),
                LibrarySource::Csv(_) => "csv".into(),
                LibrarySource::Builtin(_) => "builtin".into(),
            },
            label: lib.label(),
            detail: String::new(),
        });
    }

    params.libraries = builtins.into_iter().chain(user).collect();
    libraries
}

/// Names of every built-in library, for the advanced panel.
pub fn builtin_names() -> &'static [&'static str; 8] {
    &BUILTIN_LIBRARIES
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_order_puts_2_before_10() {
        let mut names = vec!["s-10.mzML", "s-2.mzML", "s-1.mzML"];
        names.sort_by(|a, b| natural_cmp(a, b));
        assert_eq!(names, vec!["s-1.mzML", "s-2.mzML", "s-10.mzML"]);
    }

    #[test]
    fn natural_order_is_case_insensitive_and_stable_on_text() {
        assert_eq!(natural_cmp("Abc", "abd"), Ordering::Less);
        assert_eq!(natural_cmp("run007", "run7"), Ordering::Equal);
    }

    #[test]
    fn hidden_files_are_ignored() {
        assert!(is_hidden(".DS_Store"));
        assert!(!is_hidden("sample.mzML"));
    }

    #[test]
    fn extension_match_is_case_insensitive() {
        assert!(has_extension(Path::new("a.MZML"), &["mzML", "mzXML"]));
        assert!(!has_extension(Path::new("a.raw"), &["mzML", "mzXML"]));
    }
}
