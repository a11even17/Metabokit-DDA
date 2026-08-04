//! Streaming mzML reader.
//!
//! # Differences from 0.1
//!
//! * **No per-tag allocation.** The old parser pushed
//!   `e.local_name().as_ref().to_vec()` onto a `Vec<Vec<u8>>` for *every*
//!   element in the document, then popped it. An mzML file has on the order of
//!   ten elements per spectrum, so a 50k-spectrum run performed roughly half a
//!   million heap allocations purely to track nesting. Two booleans and a depth
//!   counter carry the same information.
//! * **Decode into reused buffers.** base64 goes through `decode_slice` into a
//!   buffer owned by the parser instead of a `DecoderReader` per array.
//! * **Straight into columns.** Decoded points are appended to [`Ms1Set`] /
//!   [`Ms2Set`] rather than materialising a `Vec<(f32, f32)>` per scan.
//! * **Errors, not panics.** Profile-mode data, numpress compression and
//!   missing binary-array metadata all used to `panic!` or `expect` from inside
//!   a rayon worker.

use std::path::Path;

use base64::Engine;
use quick_xml::events::{BytesStart, Event};
use quick_xml::reader::Reader;

use crate::error::{Error, IoContext, Result};
use crate::progress::Cancel;
use crate::scans::{Ms1Set, Ms2Set};

/// Everything one mzML file contributes to a run.
pub struct MzmlData {
    pub ms1: Ms1Set,
    pub ms2: Ms2Set,
    /// `run/@startTimeStamp`, used as a column header in the reports.
    pub timestamp: String,
    /// `Some(true)` for positive mode. `None` when the file never declares it.
    pub polarity: Option<bool>,
    /// True when scan start times were given in seconds and converted.
    pub rt_converted_from_seconds: bool,
    /// Spectra that declared neither centroid nor profile mode.
    pub unknown_mode_spectra: u32,
}

#[derive(Clone, Copy, PartialEq)]
enum ArrayKind {
    Mz,
    Intensity,
    Other,
}

struct Decoder {
    b64: Vec<u8>,
    inflated: Vec<u8>,
    scratch: Vec<u8>,
}

impl Decoder {
    fn new() -> Self {
        Decoder {
            b64: Vec::new(),
            inflated: Vec::new(),
            scratch: Vec::new(),
        }
    }

    /// Decode one `<binary>` payload into `out`, reusing internal buffers.
    fn decode(&mut self, text: &str, zlib: bool, f64bit: bool, out: &mut Vec<f32>) -> Result<()> {
        let engine = &base64::engine::general_purpose::STANDARD;
        let raw = text.as_bytes();

        // Worst case is 3 bytes out per 4 in; the `+ 3` covers partial groups.
        self.b64.clear();
        self.b64.resize(raw.len() / 4 * 3 + 3, 0);
        let n = match engine.decode_slice(raw, &mut self.b64) {
            Ok(n) => n,
            Err(_) => {
                // Some writers wrap the payload across lines. That is the slow
                // path precisely because it is the rare one.
                self.scratch.clear();
                self.scratch
                    .extend(raw.iter().copied().filter(|b| !b.is_ascii_whitespace()));
                self.b64.clear();
                self.b64.resize(self.scratch.len() / 4 * 3 + 3, 0);
                engine
                    .decode_slice(&self.scratch, &mut self.b64)
                    .map_err(|e| Error::Decode(format!("base64: {e}")))?
            }
        };
        self.b64.truncate(n);

        let bytes: &[u8] = if zlib {
            use std::io::Read;
            self.inflated.clear();
            flate2::read::ZlibDecoder::new(self.b64.as_slice())
                .read_to_end(&mut self.inflated)
                .map_err(|e| Error::Decode(format!("zlib: {e}")))?;
            &self.inflated
        } else {
            &self.b64
        };

        out.clear();
        if f64bit {
            out.reserve(bytes.len() / 8);
            for c in bytes.chunks_exact(8) {
                let mut b = [0u8; 8];
                b.copy_from_slice(c);
                out.push(f64::from_le_bytes(b) as f32);
            }
        } else {
            out.reserve(bytes.len() / 4);
            for c in bytes.chunks_exact(4) {
                let mut b = [0u8; 4];
                b.copy_from_slice(c);
                out.push(f32::from_le_bytes(b));
            }
        }
        Ok(())
    }
}

/// Per-spectrum state, reset at every `<spectrum>`.
struct SpectrumState {
    ms_level: u8,
    rt: f32,
    precursor_mz: f32,
    collision_energy: f32,
    centroid: Option<bool>,
}

impl SpectrumState {
    fn reset(&mut self) {
        self.ms_level = 0;
        self.rt = f32::NAN;
        self.precursor_mz = f32::NAN;
        self.collision_energy = 0.0;
        self.centroid = None;
    }
}

fn attr_str(e: &BytesStart<'_>, name: &str) -> Option<String> {
    match e.try_get_attribute(name) {
        Ok(Some(a)) => std::str::from_utf8(a.value.as_ref())
            .ok()
            .map(str::to_string),
        _ => None,
    }
}

fn attr_f32(e: &BytesStart<'_>, name: &str) -> Option<f32> {
    match e.try_get_attribute(name) {
        Ok(Some(a)) => std::str::from_utf8(a.value.as_ref())
            .ok()
            .and_then(|s| s.trim().parse::<f32>().ok()),
        _ => None,
    }
}

fn attr_bytes_eq(e: &BytesStart<'_>, name: &str, want: &[u8]) -> bool {
    matches!(e.try_get_attribute(name), Ok(Some(a)) if a.value.as_ref() == want)
}

/// Parse one mzML file.
///
/// Cancellation is checked once per spectrum: often enough that stopping feels
/// immediate, rarely enough that the atomic load is free.
pub fn parse(path: &Path, cancel: &Cancel) -> Result<MzmlData> {
    let file_len = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
    // Centroided mzML expands to roughly file_len/8 f32 values across both
    // columns. Reserving up front removes the repeated doubling-and-copy of a
    // multi-hundred-megabyte vector, capped so a huge file cannot commit an
    // absurd amount before we have read a single scan.
    let est_points = ((file_len / 8) as usize).min(32 << 20);
    let est_scans = 4096;

    let mut reader = Reader::from_file(path).map_err(|e| Error::mzml(path, e.to_string()))?;
    reader.config_mut().trim_text(true);

    let mut ms1 = Ms1Set::with_capacity(est_scans, est_points);
    let mut ms2 = Ms2Set::new();

    let mut buf: Vec<u8> = Vec::with_capacity(1 << 16);
    let mut decoder = Decoder::new();

    let mut mz_l: Vec<f32> = Vec::new();
    let mut int_l: Vec<f32> = Vec::new();
    let mut keep_mz: Vec<f32> = Vec::new();
    let mut keep_int: Vec<f32> = Vec::new();

    let mut spec = SpectrumState {
        ms_level: 0,
        rt: f32::NAN,
        precursor_mz: f32::NAN,
        collision_energy: 0.0,
        centroid: None,
    };

    // Binary-array state, reset at each `<binaryDataArray>`.
    let mut zlib: Option<bool> = None;
    let mut f64bit: Option<bool> = None;
    let mut array_kind = ArrayKind::Other;

    // Nesting, tracked without allocating.
    let mut in_spectrum_list = false;
    let mut in_binary = false;

    let mut timestamp = String::new();
    let mut polarity: Option<bool> = None;
    let mut rt_seconds = false;
    let mut unknown_mode = 0u32;
    let mut spectra_seen: u64 = 0;

    loop {
        // Read the position before the call: afterwards `buf` is borrowed by
        // the returned event, and this also points at the start of the
        // offending element rather than past it.
        let pos = reader.buffer_position();
        let event = reader
            .read_event_into(&mut buf)
            .map_err(|e| Error::mzml(path, format!("at byte {pos}: {e}")))?;

        match event {
            Event::Eof => break,

            Event::Start(e) => match e.local_name().as_ref() {
                b"run" => {
                    if let Some(ts) = attr_str(&e, "startTimeStamp") {
                        timestamp = ts;
                    }
                }
                b"spectrumList" => in_spectrum_list = true,
                b"spectrum" => {
                    spec.reset();
                    spectra_seen += 1;
                    if spectra_seen % 256 == 0 {
                        cancel.check()?;
                    }
                }
                b"binaryDataArray" => {
                    zlib = None;
                    f64bit = None;
                    array_kind = ArrayKind::Other;
                }
                b"binary" => in_binary = true,
                b"cvParam" => {
                    read_cv_param(
                        &e,
                        &mut spec,
                        &mut zlib,
                        &mut f64bit,
                        &mut array_kind,
                        &mut polarity,
                        &mut rt_seconds,
                        path,
                    )?;
                }
                _ => {}
            },

            Event::Empty(e) => {
                if e.local_name().as_ref() == b"cvParam" {
                    read_cv_param(
                        &e,
                        &mut spec,
                        &mut zlib,
                        &mut f64bit,
                        &mut array_kind,
                        &mut polarity,
                        &mut rt_seconds,
                        path,
                    )?;
                }
            }

            // `array_kind == Other` covers arrays we do not consume (ion
            // mobility, baseline, …) — they are skipped without decoding.
            Event::Text(t)
                if in_binary
                    && in_spectrum_list
                    && spec.ms_level != 0
                    && array_kind != ArrayKind::Other =>
            {
                let zlib = zlib.ok_or_else(|| {
                    Error::mzml(path, "binary array without a compression cvParam")
                })?;
                let f64bit = f64bit
                    .ok_or_else(|| Error::mzml(path, "binary array without a precision cvParam"))?;
                let text = t
                    .decode()
                    .map_err(|e| Error::mzml(path, format!("binary text: {e}")))?;
                let target = if array_kind == ArrayKind::Mz {
                    &mut mz_l
                } else {
                    &mut int_l
                };
                decoder.decode(&text, zlib, f64bit, target)?;
            }

            Event::End(e) => match e.local_name().as_ref() {
                b"binary" => in_binary = false,
                b"spectrumList" => in_spectrum_list = false,
                b"spectrum" => {
                    if spec.centroid == Some(false) {
                        return Err(Error::mzml(
                            path,
                            "profile-mode spectra are not supported; centroid the file first \
                             (e.g. msconvert --filter \"peakPicking true 1-\")",
                        ));
                    }
                    if spec.ms_level != 0 && spec.centroid.is_none() {
                        unknown_mode += 1;
                    }

                    let n = mz_l.len().min(int_l.len());
                    keep_mz.clear();
                    keep_int.clear();
                    for i in 0..n {
                        // Zero-intensity points carry no information and would
                        // otherwise widen every binary search downstream.
                        if int_l[i] > 0.0 {
                            keep_mz.push(mz_l[i]);
                            keep_int.push(int_l[i]);
                        }
                    }

                    match spec.ms_level {
                        1 if spec.rt.is_finite() => {
                            ms1.push_scan(spec.rt, &keep_mz, &keep_int);
                        }
                        2 if spec.rt.is_finite() && spec.precursor_mz.is_finite() => {
                            ms2.push_scan(
                                spec.precursor_mz,
                                spec.rt,
                                spec.collision_energy,
                                &keep_mz,
                                &keep_int,
                            );
                        }
                        _ => {}
                    }

                    mz_l.clear();
                    int_l.clear();
                    spec.reset();
                }
                _ => {}
            },

            _ => {}
        }

        buf.clear();
    }

    ms1.shrink_to_fit();

    Ok(MzmlData {
        ms1,
        ms2,
        timestamp,
        polarity,
        rt_converted_from_seconds: rt_seconds,
        unknown_mode_spectra: unknown_mode,
    })
}

#[allow(clippy::too_many_arguments)]
fn read_cv_param(
    e: &BytesStart<'_>,
    spec: &mut SpectrumState,
    zlib: &mut Option<bool>,
    f64bit: &mut Option<bool>,
    array_kind: &mut ArrayKind,
    polarity: &mut Option<bool>,
    rt_seconds: &mut bool,
    path: &Path,
) -> Result<()> {
    // Borrowed, not owned: an mzML has millions of cvParams and allocating a
    // `Vec<u8>` for each accession would dominate the parse.
    let attr = match e.try_get_attribute("accession") {
        Ok(Some(a)) => a,
        _ => return Ok(()),
    };

    match attr.value.as_ref() {
        // ---- spectrum description ----
        b"MS:1000511" => {
            if let Some(v) = attr_f32(e, "value") {
                spec.ms_level = v as u8;
            }
        }
        b"MS:1000016" => {
            if let Some(v) = attr_f32(e, "value") {
                // Scan start time may be given in minutes or seconds, and in
                // either the PSI-MS or the Unit Ontology namespace —
                // ProteoWizard writes `UO:0000031` (minute) for AB SCIEX data,
                // so checking only the `MS:` terms would miss a seconds-based
                // file entirely. 0.1 assumed minutes unconditionally.
                if attr_bytes_eq(e, "unitAccession", b"MS:1000039")
                    || attr_bytes_eq(e, "unitAccession", b"UO:0000010")
                {
                    spec.rt = v / 60.0;
                    *rt_seconds = true;
                } else {
                    spec.rt = v;
                }
            }
        }
        b"MS:1000744" => {
            if let Some(v) = attr_f32(e, "value") {
                spec.precursor_mz = v;
            }
        }
        b"MS:1000045" => {
            if let Some(v) = attr_f32(e, "value") {
                spec.collision_energy = v;
            }
        }
        b"MS:1000127" => spec.centroid = Some(true),
        b"MS:1000128" => spec.centroid = Some(false),
        b"MS:1000130" => {
            if polarity.is_none() {
                *polarity = Some(true);
            }
        }
        b"MS:1000129" => {
            if polarity.is_none() {
                *polarity = Some(false);
            }
        }

        // ---- binary array description ----
        b"MS:1000523" => *f64bit = Some(true),
        b"MS:1000521" => *f64bit = Some(false),
        b"MS:1000574" => *zlib = Some(true),
        b"MS:1000576" => *zlib = Some(false),
        b"MS:1000514" => *array_kind = ArrayKind::Mz,
        b"MS:1000515" => *array_kind = ArrayKind::Intensity,

        // ---- unsupported encodings ----
        b"MS:1002312" | b"MS:1002313" | b"MS:1002314" | b"MS:1002746" | b"MS:1002747"
        | b"MS:1002748" => {
            return Err(Error::mzml(
                path,
                "numpress-compressed binary arrays are not supported; re-convert without \
                 numpress (msconvert --zlib)",
            ));
        }
        _ => {}
    }
    Ok(())
}

/// Read just far enough into a file to learn its scan polarity.
///
/// Used before libraries are loaded, since the built-in libraries are split by
/// polarity. Bounded so a malformed file cannot make this scan gigabytes.
pub fn sniff_polarity(path: &Path) -> Result<Option<bool>> {
    use std::io::{BufRead, BufReader};

    let file = std::fs::File::open(path).at(path)?;
    let mut reader = BufReader::with_capacity(1 << 16, file);
    let mut line = String::new();
    let mut scanned: usize = 0;

    loop {
        line.clear();
        let n = match reader.read_line(&mut line) {
            Ok(0) => return Ok(None),
            Ok(n) => n,
            // Binary payloads are valid UTF-8 base64, but a stray byte should
            // not be fatal for a best-effort sniff.
            Err(_) => return Ok(None),
        };
        scanned += n;
        if line.contains("MS:1000130") {
            return Ok(Some(true));
        }
        if line.contains("MS:1000129") {
            return Ok(Some(false));
        }
        // Polarity is declared in the first spectrum; if we are 64 MB in
        // without finding it, the file does not declare it.
        if scanned > 64 << 20 {
            return Ok(None);
        }
    }
}
