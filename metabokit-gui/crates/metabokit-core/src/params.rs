//! Run parameters.
//!
//! In 0.1 these lived in a `param.txt` TOML file read from the working
//! directory, with `unwrap()` on every lookup — a missing key was a panic. The
//! GUI owns them now, so the struct is `serde`-round-trippable to JSON for the
//! frontend and to TOML for on-disk presets, every field has a defensible
//! default, and `validate()` returns human-readable problems instead of
//! aborting.

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Built-in binary libraries shipped in the `libs/` directory next to the
/// executable. Order matters: it is the order they are concatenated in.
pub const BUILTIN_LIBRARIES: [&str; 8] = [
    "nist",
    "Atlas",
    "Atlas_filtered",
    "MSDIAL",
    "hmdb",
    "sling",
    "MassBank",
    "FiehnHILIC",
];

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Polarity {
    /// Sniff the first mzML file for a scan-polarity cvParam.
    Auto,
    Positive,
    Negative,
}

impl Default for Polarity {
    fn default() -> Self {
        Polarity::Auto
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TolUnit {
    Ppm,
    /// Absolute m/z.
    Mz,
}

impl Default for TolUnit {
    fn default() -> Self {
        TolUnit::Ppm
    }
}

impl TolUnit {
    /// Absolute tolerance at a given m/z.
    #[inline]
    pub fn absolute(self, tol: f32, mz: f32) -> f32 {
        match self {
            TolUnit::Ppm => tol * mz * 1e-6,
            TolUnit::Mz => tol,
        }
    }
}

/// One spectral library input.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", content = "value", rename_all = "lowercase")]
pub enum LibrarySource {
    /// One of [`BUILTIN_LIBRARIES`], resolved against the libs directory.
    Builtin(String),
    /// A NIST-style `.msp` text library.
    Msp(PathBuf),
    /// A `name,adduct,m/z,rt` CSV with no fragments (MS1-only matching).
    Csv(PathBuf),
}

impl LibrarySource {
    pub fn label(&self) -> String {
        match self {
            LibrarySource::Builtin(n) => n.clone(),
            LibrarySource::Msp(p) | LibrarySource::Csv(p) => p
                .file_name()
                .map(|x| x.to_string_lossy().into_owned())
                .unwrap_or_else(|| p.to_string_lossy().into_owned()),
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, rename_all = "camelCase")]
pub struct Params {
    // ---- inputs -----------------------------------------------------------
    pub mzml_files: Vec<PathBuf>,
    pub output_dir: PathBuf,
    /// Directory holding the built-in binary libraries plus `inchik.txt` and
    /// `name_formu.txt`. `None` searches next to the executable, then the
    /// working directory.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub libs_dir: Option<PathBuf>,
    pub polarity: Polarity,
    pub libraries: Vec<LibrarySource>,

    // ---- mass / retention tolerances --------------------------------------
    pub ms1_tol: f32,
    pub ms1_tol_unit: TolUnit,
    /// Absolute m/z tolerance for fragment matching.
    pub ms2_tol: f32,
    /// Precursor window half-width used to pair an MS1 feature with MS2 scans.
    pub ms1_ms2_pair: f32,
    /// m/z window for grouping the same feature across samples.
    pub mz_shift: f32,
    /// RT window for grouping the same feature across samples (minutes).
    pub rt_shift: f32,

    // ---- feature detection ------------------------------------------------
    /// Chromatographic peak width bounds in minutes, `(min, max)`.
    pub peak_width: (f32, f32),
    /// Minimum wavelet-vs-chromatogram cosine similarity for a real peak.
    pub peak_shape: f32,
    /// Minimum lag-1 autocorrelation (smoothness) for a real peak.
    pub sn_score: f32,
    /// Minimum signal-to-noise of the integrated feature.
    pub s_n_1: f32,

    // ---- spectral matching ------------------------------------------------
    pub ms2_score: f32,
    pub min_peaks: u8,
    /// Keep at most this many library fragments (and, absent an intensity
    /// cutoff, this many experimental fragments) per spectrum.
    pub match_n_fragments: usize,
    /// Absolute intensity floor for experimental fragments. Takes precedence
    /// over `match_n_fragments` on the experimental side.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub intensity_cutoff: Option<f32>,
    /// Score unmatched experimental fragments against zero (penalises
    /// co-isolation) instead of ignoring them.
    pub chimeric_spectra: bool,
    pub top_scoring_only: bool,
    pub exclude_adducts: Vec<String>,

    // ---- in-source fragments ---------------------------------------------
    /// Minimum precursor mass gap between an ISF and its parent.
    pub isf_parent_mass_diff: f32,
    /// RT window for calling two features the same elution (minutes).
    pub isf_rt_diff: f32,

    // ---- quantification ---------------------------------------------------
    /// Only report features; drop MS2 scans with no MS1 feature underneath.
    pub features_only: bool,
    /// Integration bounds as a multiple of the fitted peak half-width.
    pub integration_rt: f32,
    /// Half-width of the m/z window used to extract ion chromatograms.
    pub integration_mz: f32,
    /// Half-width (minutes) used when imputing a missing value.
    pub impute_width: f32,
    pub gap_fill: bool,

    // ---- execution --------------------------------------------------------
    /// Worker threads. `0` uses all available parallelism.
    pub threads: usize,
    /// How many mzML files may be resident in memory at once. `0` derives a
    /// value from `threads`. See `pipeline` for why this is capped separately.
    pub max_files_in_flight: usize,
    /// Keep the per-sample binary caches after the run (needed by the
    /// visualizer; roughly the size of the input).
    pub keep_cache: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            mzml_files: Vec::new(),
            output_dir: PathBuf::new(),
            libs_dir: None,
            polarity: Polarity::Auto,
            libraries: Vec::new(),

            ms1_tol: 10.0,
            ms1_tol_unit: TolUnit::Ppm,
            ms2_tol: 0.01,
            ms1_ms2_pair: 0.5,
            mz_shift: 0.01,
            rt_shift: 0.3,

            peak_width: (0.06, 1.2),
            peak_shape: 0.7,
            sn_score: 0.4,
            s_n_1: 3.0,

            ms2_score: 0.6,
            min_peaks: 2,
            match_n_fragments: 15,
            intensity_cutoff: None,
            chimeric_spectra: false,
            top_scoring_only: false,
            exclude_adducts: Vec::new(),

            isf_parent_mass_diff: 10.0,
            isf_rt_diff: 0.04,

            features_only: false,
            integration_rt: 2.0,
            integration_mz: 0.01,
            impute_width: 0.1,
            gap_fill: true,

            threads: 0,
            max_files_in_flight: 0,
            keep_cache: true,
        }
    }
}

/// A validation problem. `fatal` problems block the run.
#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct Problem {
    pub field: String,
    pub message: String,
    pub fatal: bool,
}

impl Problem {
    fn fatal(field: &str, message: impl Into<String>) -> Self {
        Problem {
            field: field.to_string(),
            message: message.into(),
            fatal: true,
        }
    }

    fn warn(field: &str, message: impl Into<String>) -> Self {
        Problem {
            field: field.to_string(),
            message: message.into(),
            fatal: false,
        }
    }
}

impl Params {
    /// Effective worker-thread count.
    pub fn effective_threads(&self) -> usize {
        if self.threads > 0 {
            self.threads
        } else {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(4)
        }
    }

    /// How many samples may be held in RAM simultaneously.
    ///
    /// The 0.1 CLI ran `mzml_fs.par_iter()` across files, so peak RSS scaled
    /// with the thread count: eight threads meant eight fully-decoded mzML
    /// files resident at once, which is where the memory blow-ups came from.
    /// The rewrite parallelises *inside* a sample instead and keeps this small,
    /// so peak memory is bounded by file size rather than by core count.
    pub fn effective_in_flight(&self) -> usize {
        let n = if self.max_files_in_flight > 0 {
            self.max_files_in_flight
        } else {
            // Two gives the parser something to overlap with scoring without
            // meaningfully raising the high-water mark.
            2
        };
        n.max(1).min(self.mzml_files.len().max(1))
    }

    /// Resolve the directory holding the built-in libraries and metadata
    /// tables. Tries the explicit setting, then next to the executable, then
    /// the executable's parent (covers `target/debug/`), then the CWD.
    pub fn resolve_libs_dir(&self) -> Option<PathBuf> {
        if let Some(dir) = &self.libs_dir {
            return dir.is_dir().then(|| dir.clone());
        }
        let mut candidates: Vec<PathBuf> = Vec::new();
        if let Ok(exe) = std::env::current_exe() {
            if let Some(dir) = exe.parent() {
                candidates.push(dir.join("libs"));
                if let Some(up) = dir.parent() {
                    candidates.push(up.join("libs"));
                    if let Some(up2) = up.parent() {
                        candidates.push(up2.join("libs"));
                    }
                }
            }
        }
        if let Ok(cwd) = std::env::current_dir() {
            candidates.push(cwd.join("libs"));
        }
        candidates.into_iter().find(|p| p.is_dir())
    }

    pub fn misc_dir(&self) -> PathBuf {
        self.output_dir.join("misc")
    }

    /// Human-readable problems with the current settings.
    pub fn validate(&self) -> Vec<Problem> {
        let mut out = Vec::new();

        if self.mzml_files.is_empty() {
            out.push(Problem::fatal(
                "mzmlFiles",
                "Select at least one mzML file.",
            ));
        }
        for f in &self.mzml_files {
            if !f.is_file() {
                out.push(Problem::fatal(
                    "mzmlFiles",
                    format!("Not found: {}", f.display()),
                ));
            }
        }
        if self.output_dir.as_os_str().is_empty() {
            out.push(Problem::fatal("outputDir", "Choose an output directory."));
        }
        if self.libraries.is_empty() {
            out.push(Problem::fatal(
                "libraries",
                "Add at least one spectral library.",
            ));
        }
        let needs_builtin = self
            .libraries
            .iter()
            .any(|l| matches!(l, LibrarySource::Builtin(_)));
        if needs_builtin && self.resolve_libs_dir().is_none() {
            out.push(Problem::fatal(
                "libsDir",
                "Built-in libraries selected but no libs/ directory was found.",
            ));
        }
        for lib in &self.libraries {
            match lib {
                LibrarySource::Builtin(name) => {
                    if !BUILTIN_LIBRARIES
                        .iter()
                        .any(|b| b.eq_ignore_ascii_case(name))
                    {
                        out.push(Problem::fatal(
                            "libraries",
                            format!("Unknown built-in library {name:?}."),
                        ));
                    }
                }
                LibrarySource::Msp(p) | LibrarySource::Csv(p) => {
                    if !p.is_file() {
                        out.push(Problem::fatal(
                            "libraries",
                            format!("Not found: {}", p.display()),
                        ));
                    }
                }
            }
        }

        if !(self.ms1_tol > 0.0) {
            out.push(Problem::fatal("ms1Tol", "MS1 tolerance must be positive."));
        }
        if !(self.ms2_tol > 0.0) {
            out.push(Problem::fatal("ms2Tol", "MS2 tolerance must be positive."));
        }
        if !(self.peak_width.0 > 0.0) || self.peak_width.1 <= self.peak_width.0 {
            out.push(Problem::fatal(
                "peakWidth",
                "Peak width must be an increasing positive range.",
            ));
        }
        if !(0.0..=1.0).contains(&self.ms2_score) {
            out.push(Problem::fatal(
                "ms2Score",
                "MS2 score must be within 0 – 1.",
            ));
        }
        if !(0.0..=1.0).contains(&self.peak_shape) {
            out.push(Problem::fatal(
                "peakShape",
                "Peak shape must be within 0 – 1.",
            ));
        }
        if self.match_n_fragments == 0 {
            out.push(Problem::fatal(
                "matchNFragments",
                "Fragment count must be at least 1.",
            ));
        }
        if self.min_peaks == 0 {
            out.push(Problem::warn(
                "minPeaks",
                "With 0 matching peaks required, every library entry in the mass window will be reported.",
            ));
        }
        if self.match_n_fragments > 255 {
            out.push(Problem::warn(
                "matchNFragments",
                "Values above 255 are clamped when spectra are written.",
            ));
        }
        if self.ms1_tol_unit == TolUnit::Ppm && self.ms1_tol > 100.0 {
            out.push(Problem::warn(
                "ms1Tol",
                "A tolerance above 100 ppm is unusually wide.",
            ));
        }
        if self.rt_shift <= 0.0 {
            out.push(Problem::fatal("rtShift", "RT shift must be positive."));
        }
        if self.integration_rt <= 0.0 {
            out.push(Problem::fatal(
                "integrationRt",
                "Integration width must be positive.",
            ));
        }
        out
    }

    pub fn check(&self) -> Result<()> {
        let problems = self.validate();
        let fatal: Vec<String> = problems
            .iter()
            .filter(|p| p.fatal)
            .map(|p| p.message.clone())
            .collect();
        if fatal.is_empty() {
            Ok(())
        } else {
            Err(Error::param(fatal.join(" ")))
        }
    }

    /// Presets are JSON, not TOML: the frontend speaks JSON natively, and TOML
    /// cannot represent this struct's field order (an array-of-tables followed
    /// by scalars) without reordering.
    pub fn to_json(&self) -> Result<String> {
        serde_json::to_string_pretty(self).map_err(|e| Error::param(e.to_string()))
    }

    pub fn from_json(text: &str) -> Result<Params> {
        serde_json::from_str(text).map_err(|e| Error::param(e.to_string()))
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        std::fs::write(path, self.to_json()?).map_err(|e| Error::io(e, path))
    }

    /// Loads a `.json` preset, or a 0.1-era `param.txt`.
    pub fn load(path: impl AsRef<Path>) -> Result<Params> {
        let path = path.as_ref();
        let text = std::fs::read_to_string(path).map_err(|e| Error::io(e, path))?;
        match Params::from_json(&text) {
            Ok(p) => Ok(p),
            Err(json_err) => Params::from_legacy(&text, path).map_err(|_| json_err),
        }
    }

    /// Import a 0.1-era `param.txt`.
    ///
    /// Kept so existing analyses stay reproducible. Unknown keys are ignored
    /// and missing keys fall back to [`Params::default`], which is the whole
    /// point: the old reader panicked on both.
    pub fn from_legacy(text: &str, origin: &Path) -> Result<Params> {
        let table: toml::Table = text
            .parse()
            .map_err(|e| Error::param(format!("{}: {e}", origin.display())))?;

        let mut p = Params::default();

        let f32_at = |t: &toml::Table, k: &str| -> Option<f32> {
            t.get(k).and_then(|v| match v {
                toml::Value::Float(x) => Some(*x as f32),
                toml::Value::Integer(x) => Some(*x as f32),
                _ => None,
            })
        };
        let bool_at = |t: &toml::Table, k: &str| t.get(k).and_then(toml::Value::as_bool);
        let usize_at = |t: &toml::Table, k: &str| {
            t.get(k)
                .and_then(toml::Value::as_integer)
                .map(|x| x as usize)
        };

        if let Some(pattern) = table.get("mzML_files").and_then(toml::Value::as_str) {
            if let Ok(paths) = glob::glob(pattern) {
                p.mzml_files = paths.filter_map(|x| x.ok()).collect();
            }
        }
        if let Some(list) = table.get("library").and_then(toml::Value::as_array) {
            for item in list.iter().filter_map(toml::Value::as_str) {
                let item = item.trim();
                if let Some(rest) = item.strip_prefix("user ") {
                    p.libraries
                        .push(LibrarySource::Msp(PathBuf::from(rest.trim())));
                } else if let Some(rest) = item.strip_prefix("csv ") {
                    p.libraries
                        .push(LibrarySource::Csv(PathBuf::from(rest.trim())));
                } else if let Some(name) = BUILTIN_LIBRARIES
                    .iter()
                    .find(|b| b.eq_ignore_ascii_case(item))
                {
                    p.libraries
                        .push(LibrarySource::Builtin((*name).to_string()));
                }
            }
        }
        if let Some(arr) = table.get("ms1tol").and_then(toml::Value::as_array) {
            if let Some(v) = arr.first() {
                p.ms1_tol = match v {
                    toml::Value::Float(x) => *x as f32,
                    toml::Value::Integer(x) => *x as f32,
                    _ => p.ms1_tol,
                };
            }
            if let Some(u) = arr.get(1).and_then(toml::Value::as_str) {
                p.ms1_tol_unit = if u.eq_ignore_ascii_case("ppm") {
                    TolUnit::Ppm
                } else {
                    TolUnit::Mz
                };
            }
        }
        if let Some(arr) = table
            .get("length_of_ion_chromatogram")
            .and_then(toml::Value::as_array)
        {
            let mut it = arr.iter().filter_map(|v| match v {
                toml::Value::Float(x) => Some(*x as f32),
                toml::Value::Integer(x) => Some(*x as f32),
                _ => None,
            });
            if let (Some(a), Some(b)) = (it.next(), it.next()) {
                p.peak_width = (a, b);
            }
        }
        if let Some(list) = table.get("exclude_adduct").and_then(toml::Value::as_array) {
            p.exclude_adducts = list
                .iter()
                .filter_map(toml::Value::as_str)
                .map(str::to_string)
                .collect();
        }

        if let Some(v) = f32_at(&table, "ms2tol") {
            p.ms2_tol = v;
        }
        if let Some(v) = f32_at(&table, "mz_shift") {
            p.mz_shift = v;
        }
        if let Some(v) = f32_at(&table, "ISF_parent_mass_diff") {
            p.isf_parent_mass_diff = v;
        }
        if let Some(v) = f32_at(&table, "MS2_score") {
            p.ms2_score = v;
        }
        if let Some(v) = f32_at(&table, "RT_shift") {
            p.rt_shift = v;
        }
        if let Some(v) = f32_at(&table, "peak_shape") {
            p.peak_shape = v;
        }
        if let Some(v) = f32_at(&table, "sn_score") {
            p.sn_score = v;
        }
        if let Some(v) = f32_at(&table, "S_N_1") {
            p.s_n_1 = v;
        }
        if let Some(v) = f32_at(&table, "MS1_MS2_pair") {
            p.ms1_ms2_pair = v;
        }
        if let Some(v) = f32_at(&table, "impute_width") {
            p.impute_width = v;
        }
        if let Some(v) = f32_at(&table, "integration_RT") {
            p.integration_rt = v;
        }
        if let Some(v) = f32_at(&table, "integration_mz") {
            p.integration_mz = v;
        }
        p.intensity_cutoff = f32_at(&table, "intensity_cutoff");
        if let Some(v) = usize_at(&table, "min_peaks") {
            p.min_peaks = v.min(255) as u8;
        }
        if let Some(v) = usize_at(&table, "num_threads") {
            p.threads = v;
        }
        if let Some(v) = usize_at(&table, "match_n_fragments") {
            p.match_n_fragments = v;
        }
        if let Some(v) = bool_at(&table, "features_only") {
            p.features_only = v;
        }
        if let Some(v) = bool_at(&table, "top_scoring_only") {
            p.top_scoring_only = v;
        }
        if let Some(v) = bool_at(&table, "chimeric_spectra") {
            p.chimeric_spectra = v;
        }

        if let Some(dir) = origin.parent() {
            if p.output_dir.as_os_str().is_empty() {
                p.output_dir = dir.to_path_buf();
            }
            // 0.1 honoured a sibling `file_order.txt`, whose first row is a
            // base directory and whose remaining rows are file names.
            let order = dir.join("file_order.txt");
            if order.is_file() {
                if let Ok(files) = read_file_order(&order) {
                    if !files.is_empty() {
                        p.mzml_files = files;
                    }
                }
            }
        }

        Ok(p)
    }
}

fn read_file_order(path: &Path) -> Result<Vec<PathBuf>> {
    let mut rdr = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path(path)?;
    let mut records = rdr.records();
    let base = match records.next() {
        Some(row) => {
            let row = row?;
            PathBuf::from(row.get(0).unwrap_or("").to_string())
        }
        None => return Ok(Vec::new()),
    };
    let mut out = Vec::new();
    for row in records {
        let row = row?;
        if let Some(name) = row.get(0) {
            if !name.is_empty() {
                out.push(base.join(name));
            }
        }
    }
    Ok(out)
}
