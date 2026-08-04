//! Spectral library loading.
//!
//! # Why this shape
//!
//! 0.1 represented a library as `Vec<Ent>` where every entry owned four
//! `String`s and a `Vec<(f32, f32)>`. A million-entry library is therefore five
//! million heap allocations, ~120 bytes of `String`/`Vec` headers per entry
//! before a single byte of payload, and a sort that swaps 100+ byte structs.
//!
//! Here the entries are columns, all text lives in one arena addressed by
//! `(offset, len)`, adducts are interned to `u16` (there are a few dozen
//! distinct values across millions of entries), and fragments live in a single
//! pair of flat arrays. Entries are referred to downstream by `u32` index,
//! which also removes the `<'a>` lifetime that 0.1 threaded through every type
//! from `Lib` to `Ann` to the aligner.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use regex::Regex;

use crate::error::{Error, IoContext, Result};
use crate::params::{LibrarySource, Params};
use crate::progress::Reporter;

/// A slice of the string arena.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StrRef {
    start: u32,
    len: u32,
}

impl StrRef {
    pub const EMPTY: StrRef = StrRef { start: 0, len: 0 };

    pub fn is_empty(self) -> bool {
        self.len == 0
    }
}

/// Append-only string storage. Interning is by value only where it pays
/// (adducts); names are near-unique so deduplicating them would cost more in
/// hashing than it saves.
#[derive(Default, Debug)]
pub struct StrArena {
    buf: String,
}

impl StrArena {
    pub fn push(&mut self, s: &str) -> StrRef {
        if s.is_empty() {
            return StrRef::EMPTY;
        }
        let start = self.buf.len() as u32;
        self.buf.push_str(s);
        StrRef {
            start,
            len: s.len() as u32,
        }
    }

    #[inline]
    pub fn get(&self, r: StrRef) -> &str {
        // Every `StrRef` was produced by `push`, so both ends sit on char
        // boundaries and this is an O(1) bounds check.
        let a = r.start as usize;
        let b = a + r.len as usize;
        self.buf.get(a..b).unwrap_or("")
    }

    pub fn heap_bytes(&self) -> usize {
        self.buf.capacity()
    }
}

/// Per-source load statistics, surfaced in the UI.
#[derive(Clone, Debug, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct LoadedSource {
    pub label: String,
    pub entries: usize,
    pub skipped: usize,
}

/// One library entry, staged before the columns are built.
struct Staged {
    mmass: f32,
    charge: i8,
    rt: f32,
    adduct: u16,
    name: StrRef,
    inchikey: StrRef,
    formula: StrRef,
    frag_start: u32,
    frag_len: u32,
}

/// The loaded library, in column form and sorted by precursor mass.
pub struct Library {
    mmass: Vec<f32>,
    charge: Vec<i8>,
    rt: Vec<f32>,
    adduct: Vec<u16>,
    name: Vec<StrRef>,
    inchikey: Vec<StrRef>,
    formula: Vec<StrRef>,
    frag_off: Vec<u32>,
    frag_mz: Vec<f32>,
    frag_i: Vec<f32>,
    strings: StrArena,
    adduct_names: Vec<String>,
    pub sources: Vec<LoadedSource>,
}

impl Library {
    #[inline]
    pub fn len(&self) -> usize {
        self.mmass.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.mmass.is_empty()
    }

    /// Precursor masses, ascending. Everything downstream binary-searches this.
    #[inline]
    pub fn masses(&self) -> &[f32] {
        &self.mmass
    }

    #[inline]
    pub fn mass(&self, i: usize) -> f32 {
        self.mmass[i]
    }

    #[inline]
    pub fn charge(&self, i: usize) -> i8 {
        self.charge[i]
    }

    /// `None` when the entry carries no reference retention time.
    #[inline]
    pub fn rt(&self, i: usize) -> Option<f32> {
        let v = self.rt[i];
        if v.is_nan() {
            None
        } else {
            Some(v)
        }
    }

    #[inline]
    pub fn name(&self, i: usize) -> &str {
        self.strings.get(self.name[i])
    }

    #[inline]
    pub fn inchikey(&self, i: usize) -> &str {
        self.strings.get(self.inchikey[i])
    }

    #[inline]
    pub fn formula(&self, i: usize) -> &str {
        self.strings.get(self.formula[i])
    }

    #[inline]
    pub fn adduct(&self, i: usize) -> &str {
        self.adduct_names
            .get(self.adduct[i] as usize)
            .map(String::as_str)
            .unwrap_or("")
    }

    /// `(mz, intensity)` columns of the entry's reference fragments.
    #[inline]
    pub fn fragments(&self, i: usize) -> (&[f32], &[f32]) {
        let a = self.frag_off[i] as usize;
        let b = self.frag_off[i + 1] as usize;
        (&self.frag_mz[a..b], &self.frag_i[a..b])
    }

    #[inline]
    pub fn fragment_count(&self, i: usize) -> usize {
        (self.frag_off[i + 1] - self.frag_off[i]) as usize
    }

    /// Index of the first entry with mass at or above `mz`.
    #[inline]
    pub fn mass_at_or_after(&self, mz: f32) -> usize {
        self.mmass.partition_point(|&x| x < mz)
    }

    pub fn heap_bytes(&self) -> usize {
        self.mmass.capacity() * 4
            + self.charge.capacity()
            + self.rt.capacity() * 4
            + self.adduct.capacity() * 2
            + (self.name.capacity() + self.inchikey.capacity() + self.formula.capacity()) * 8
            + self.frag_off.capacity() * 4
            + (self.frag_mz.capacity() + self.frag_i.capacity()) * 4
            + self.strings.heap_bytes()
    }
}

// ---------------------------------------------------------------------------
// Builder
// ---------------------------------------------------------------------------

struct Builder {
    staged: Vec<Staged>,
    frag_mz: Vec<f32>,
    frag_i: Vec<f32>,
    strings: StrArena,
    adduct_names: Vec<String>,
    adduct_ids: HashMap<String, u16>,
}

impl Builder {
    fn new() -> Self {
        Builder {
            staged: Vec::new(),
            frag_mz: Vec::new(),
            frag_i: Vec::new(),
            strings: StrArena::default(),
            adduct_names: Vec::new(),
            adduct_ids: HashMap::new(),
        }
    }

    fn adduct_id(&mut self, adduct: &str) -> u16 {
        if let Some(&id) = self.adduct_ids.get(adduct) {
            return id;
        }
        let id = self.adduct_names.len().min(u16::MAX as usize) as u16;
        self.adduct_names.push(adduct.to_string());
        self.adduct_ids.insert(adduct.to_string(), id);
        id
    }

    /// Stage one entry. `frags` are `(mz, intensity)` in the order they should
    /// be reported.
    fn push(
        &mut self,
        mmass: f32,
        charge: i8,
        rt: Option<f32>,
        adduct: &str,
        name: &str,
        inchikey: &str,
        formula: &str,
        frags: &[(f32, f32)],
    ) {
        let frag_start = self.frag_mz.len() as u32;
        for &(mz, i) in frags {
            self.frag_mz.push(mz);
            self.frag_i.push(i);
        }
        let adduct = self.adduct_id(adduct);
        self.staged.push(Staged {
            mmass,
            charge,
            rt: rt.unwrap_or(f32::NAN),
            adduct,
            name: self.strings.push(name),
            inchikey: self.strings.push(inchikey),
            formula: self.strings.push(formula),
            frag_start,
            frag_len: self.frag_mz.len() as u32 - frag_start,
        });
    }

    /// Apply the global filters, sort by mass and flatten into columns.
    fn finish(mut self, params: &Params, sources: Vec<LoadedSource>) -> Library {
        // Drop excluded adducts by id, so the comparison is an integer test
        // rather than a string compare per entry.
        let excluded: Vec<u16> = self
            .adduct_names
            .iter()
            .enumerate()
            .filter(|(_, n)| params.exclude_adducts.iter().any(|x| x == *n))
            .map(|(i, _)| i as u16)
            .collect();
        if !excluded.is_empty() {
            self.staged.retain(|s| !excluded.contains(&s.adduct));
        }

        // Sort by precursor mass. Ties break on index so a run is reproducible
        // regardless of the sort implementation.
        let mut order: Vec<u32> = (0..self.staged.len() as u32).collect();
        order.sort_unstable_by(|&a, &b| {
            let (x, y) = (self.staged[a as usize].mmass, self.staged[b as usize].mmass);
            x.partial_cmp(&y)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.cmp(&b))
        });

        let n = order.len();
        let mut lib = Library {
            mmass: Vec::with_capacity(n),
            charge: Vec::with_capacity(n),
            rt: Vec::with_capacity(n),
            adduct: Vec::with_capacity(n),
            name: Vec::with_capacity(n),
            inchikey: Vec::with_capacity(n),
            formula: Vec::with_capacity(n),
            frag_off: Vec::with_capacity(n + 1),
            frag_mz: Vec::with_capacity(self.frag_mz.len()),
            frag_i: Vec::with_capacity(self.frag_i.len()),
            strings: StrArena::default(),
            adduct_names: Vec::new(),
            sources,
        };
        lib.frag_off.push(0);

        let keep = params.match_n_fragments;
        for &oi in &order {
            let s = &self.staged[oi as usize];
            let a = s.frag_start as usize;
            let b = a + s.frag_len as usize;
            // Fragments at or above the precursor are not informative, and
            // only the top `match_n_fragments` are scored.
            let cutoff = s.mmass - 0.3;
            let mut taken = 0usize;
            for k in a..b {
                if taken >= keep {
                    break;
                }
                if self.frag_mz[k] < cutoff {
                    lib.frag_mz.push(self.frag_mz[k]);
                    lib.frag_i.push(self.frag_i[k]);
                    taken += 1;
                }
            }
            lib.mmass.push(s.mmass);
            lib.charge.push(s.charge);
            lib.rt.push(s.rt);
            lib.adduct.push(s.adduct);
            lib.name.push(s.name);
            lib.inchikey.push(s.inchikey);
            lib.formula.push(s.formula);
            lib.frag_off.push(lib.frag_mz.len() as u32);
        }

        // The arena and adduct table move across wholesale; `StrRef`s stay valid.
        lib.strings = self.strings;
        lib.adduct_names = self.adduct_names;
        lib
    }
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

/// Cursor over a mapped byte slice. Replaces 0.1's per-field `read_exact`
/// calls, which issued one syscall-shaped read per string and per float.
struct Cursor<'a> {
    b: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(b: &'a [u8]) -> Self {
        Cursor { b, at: 0 }
    }

    fn done(&self) -> bool {
        self.at >= self.b.len()
    }

    fn take(&mut self, n: usize) -> Option<&'a [u8]> {
        let end = self.at.checked_add(n)?;
        if end > self.b.len() {
            return None;
        }
        let s = &self.b[self.at..end];
        self.at = end;
        Some(s)
    }

    fn f32(&mut self) -> Option<f32> {
        let s = self.take(4)?;
        Some(f32::from_le_bytes([s[0], s[1], s[2], s[3]]))
    }

    fn i8(&mut self) -> Option<i8> {
        Some(self.take(1)?[0] as i8)
    }

    fn u8(&mut self) -> Option<u8> {
        Some(self.take(1)?[0])
    }

    fn str_u8(&mut self) -> Option<&'a str> {
        let n = self.u8()? as usize;
        let s = self.take(n)?;
        std::str::from_utf8(s).ok()
    }
}

fn map_read_only(path: &Path) -> Result<Mmap> {
    let file = File::open(path).at(path)?;
    // SAFETY: read-only mapping of a library file. As with the scan cache, an
    // external truncation mid-read can fault; libraries are static assets that
    // ship with the app, so this is the same trade-off every mmap'd data file
    // makes.
    unsafe { Mmap::map(&file) }.at(path)
}

/// Built-in binary library format:
/// `f32 mass | i8 charge | u8+bytes name | u8+bytes adduct | u8+bytes inchikey
///  | u8 n | n×f32 mz | n×f32 intensity`
fn read_builtin(
    builder: &mut Builder,
    dir: &Path,
    name: &str,
    positive: bool,
) -> Result<LoadedSource> {
    // "Atlas_filtered" is the Atlas library plus a lipid-nomenclature filter,
    // which 0.1 only ever applied in positive mode.
    let filtered = name == "Atlas_filtered";
    let base = if filtered { "Atlas" } else { name };
    let file_name = format!("{base}_{}", if positive { "pos" } else { "neg" });
    let path = dir.join(&file_name);
    if !path.is_file() {
        return Err(Error::library(
            name,
            format!("{} not found", path.display()),
        ));
    }
    let apply_filter = filtered && positive;

    let map = map_read_only(&path)?;
    let mut cur = Cursor::new(&map);
    let chain_re = Regex::new(r"\d\d?:\d\d?").expect("static regex");

    let mut entries = 0usize;
    let mut skipped = 0usize;
    let mut frags: Vec<(f32, f32)> = Vec::new();

    while !cur.done() {
        let Some(mmass) = cur.f32() else { break };
        let Some(charge) = cur.i8() else { break };
        let Some(ent_name) = cur.str_u8() else { break };
        let Some(adduct) = cur.str_u8() else { break };
        let Some(inchikey) = cur.str_u8() else { break };
        let Some(n) = cur.u8() else { break };
        let n = n as usize;

        frags.clear();
        frags.reserve(n);
        let mz_start = cur.at;
        let i_start = mz_start + n * 4;
        if cur.take(n * 8).is_none() {
            break;
        }
        for k in 0..n {
            let m = f32::from_le_bytes([
                map[mz_start + k * 4],
                map[mz_start + k * 4 + 1],
                map[mz_start + k * 4 + 2],
                map[mz_start + k * 4 + 3],
            ]);
            let v = f32::from_le_bytes([
                map[i_start + k * 4],
                map[i_start + k * 4 + 1],
                map[i_start + k * 4 + 2],
                map[i_start + k * 4 + 3],
            ]);
            frags.push((m, v));
        }

        if apply_filter && !lipid_chain_ok(&chain_re, ent_name) {
            skipped += 1;
            continue;
        }
        if apply_filter {
            // Drop fragments above the precursor-minus-water region.
            let cutoff = mmass - 18.1;
            frags.retain(|x| x.0 < cutoff);
        }

        builder.push(mmass, charge, None, adduct, ent_name, inchikey, "", &frags);
        entries += 1;
    }

    Ok(LoadedSource {
        label: name.to_string(),
        entries,
        skipped,
    })
}

/// Lipid shorthand filter carried over from 0.1: entries whose acyl-chain
/// notation falls outside plausible ranges are dropped.
fn lipid_chain_ok(re: &Regex, name: &str) -> bool {
    if ["CAR ", "VAE ", "CE ", "CL "]
        .iter()
        .any(|p| name.starts_with(p))
    {
        return true;
    }
    let chains: Vec<(u8, u8)> = re
        .find_iter(name)
        .filter_map(|m| {
            let (c, d) = m.as_str().split_once(':')?;
            Some((c.parse::<u8>().ok()?, d.parse::<u8>().ok()?))
        })
        .collect();

    let implausible = |cs: &[(u8, u8)]| cs.iter().any(|x| !(12..=26).contains(&x.0) || x.1 > 6);

    match chains.len() {
        0 => true,
        1 => {
            if name.starts_with("LP") {
                !implausible(&chains)
            } else {
                (24..=52).contains(&chains[0].0)
            }
        }
        _ => !implausible(&chains),
    }
}

/// NIST-style `.msp` text library.
fn read_msp(builder: &mut Builder, path: &Path) -> Result<LoadedSource> {
    let file = File::open(path).at(path)?;
    let mut rdr = BufReader::with_capacity(1 << 16, file);
    let bracket_re = Regex::new(r"\[(.*)\](.)").expect("static regex");
    let label = path
        .file_name()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned());

    let mut name = String::new();
    let mut formula = String::new();
    let mut inchikey = String::new();
    let mut adduct = String::new();
    let mut rt: Option<f32> = None;
    let mut mmass = 0.0f32;
    let mut charge: i8 = 1;

    let mut line = String::new();
    let mut peak_line = String::new();
    // `(mz, intensity, starred)` — a `*` marks a curated diagnostic fragment.
    let mut frags: Vec<(f32, f32, bool)> = Vec::new();
    let mut kept: Vec<(f32, f32)> = Vec::new();
    let mut entries = 0usize;

    loop {
        line.clear();
        let n = rdr.read_line(&mut line).at(path)?;
        if n == 0 {
            break;
        }
        let Some((key, value)) = line.split_once(':') else {
            continue;
        };
        let value = value.trim();
        let key = key.trim().to_ascii_uppercase();

        match key.as_str() {
            "NAME" => name = value.to_string(),
            "FORMULA" => formula = value.to_string(),
            "PRECURSORMZ" => mmass = value.parse().unwrap_or(0.0),
            "PRECURSORTYPE" | "PRECURSOR_TYPE" => {
                let (a, c) = match bracket_re.captures(value) {
                    Some(caps) => (
                        caps.get(1)
                            .map(|m| m.as_str().to_string())
                            .unwrap_or_default(),
                        caps.get(2)
                            .and_then(|m| m.as_str().parse::<i8>().ok())
                            .unwrap_or(1),
                    ),
                    None => (value.to_string(), 1),
                };
                adduct = a;
                charge = c;
                if !name.contains(&adduct) {
                    name = format!("{name} {adduct}");
                }
            }
            "RETENTIONTIME" => rt = value.parse().ok(),
            "INCHIKEY" => inchikey = value.to_string(),
            "NUM PEAKS" => {
                frags.clear();
                loop {
                    peak_line.clear();
                    if rdr.read_line(&mut peak_line).at(path)? == 0 {
                        break;
                    }
                    if peak_line.trim().is_empty() {
                        break;
                    }
                    let mut it = peak_line.split_whitespace();
                    let (Some(mz), Some(inten)) = (it.next(), it.next()) else {
                        continue;
                    };
                    let (Ok(mz), Ok(inten)) = (mz.parse::<f32>(), inten.parse::<f32>()) else {
                        continue;
                    };
                    frags.push((mz, inten, it.next() == Some("*")));
                }

                // If any fragment is starred, only starred fragments count.
                let starred = frags.iter().any(|x| x.2);
                kept.clear();
                kept.extend(frags.iter().filter(|x| !starred || x.2).map(|x| (x.0, x.1)));
                // Descending intensity, so the later `match_n_fragments`
                // truncation keeps the most informative peaks.
                kept.sort_unstable_by(|a, b| {
                    b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal)
                });

                let display_name = format!("{name} ({label})");
                builder.push(
                    mmass,
                    charge,
                    rt,
                    &adduct,
                    &display_name,
                    &inchikey,
                    &formula,
                    &kept,
                );
                entries += 1;

                name.clear();
                formula.clear();
                inchikey.clear();
                adduct.clear();
                rt = None;
                mmass = 0.0;
                charge = 1;
            }
            _ => {}
        }
    }

    Ok(LoadedSource {
        label,
        entries,
        skipped: 0,
    })
}

/// `name, adduct, m/z, rt` CSV — MS1-only entries with no fragments.
fn read_csv_library(builder: &mut Builder, path: &Path) -> Result<LoadedSource> {
    let label = path
        .file_name()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned());
    let mut rdr = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(true)
        .trim(csv::Trim::All)
        .flexible(true)
        .from_path(path)?;

    let mut entries = 0usize;
    let mut skipped = 0usize;
    for row in rdr.records() {
        let row = row?;
        let (Some(name), Some(adduct), Some(mz)) = (row.get(0), row.get(1), row.get(2)) else {
            skipped += 1;
            continue;
        };
        let Ok(mmass) = mz.parse::<f32>() else {
            skipped += 1;
            continue;
        };
        let rt = row.get(3).and_then(|x| x.parse::<f32>().ok());
        let display_name = if name.contains(adduct) {
            name.to_string()
        } else {
            format!("{name} {adduct}")
        };
        builder.push(mmass, 1, rt, adduct, &display_name, "", "", &[]);
        entries += 1;
    }

    Ok(LoadedSource {
        label,
        entries,
        skipped,
    })
}

/// Load every configured library into one mass-sorted table.
pub fn load(params: &Params, positive: bool, reporter: &dyn Reporter) -> Result<Library> {
    let mut builder = Builder::new();
    let mut sources = Vec::new();
    let libs_dir = params.resolve_libs_dir();

    for source in &params.libraries {
        let loaded = match source {
            LibrarySource::Builtin(name) => {
                let Some(dir) = libs_dir.as_ref() else {
                    return Err(Error::library(name.clone(), "no libs/ directory found"));
                };
                read_builtin(&mut builder, dir, name, positive)
            }
            LibrarySource::Msp(p) => read_msp(&mut builder, p),
            LibrarySource::Csv(p) => read_csv_library(&mut builder, p),
        };
        match loaded {
            Ok(s) => {
                reporter.info(format!(
                    "{}: {} entries{}",
                    s.label,
                    s.entries,
                    if s.skipped > 0 {
                        format!(", {} filtered out", s.skipped)
                    } else {
                        String::new()
                    }
                ));
                sources.push(s);
            }
            // One unreadable library should not abandon a run that has other
            // usable ones; the pipeline fails later if nothing loaded.
            Err(e) => reporter.warn(format!("{e}")),
        }
    }

    let lib = builder.finish(params, sources);
    if lib.is_empty() {
        return Err(Error::library(
            "libraries",
            "no entries were loaded from any source",
        ));
    }
    reporter.info(format!(
        "library: {} entries, {:.1} MB resident",
        lib.len(),
        lib.heap_bytes() as f64 / 1e6
    ));
    Ok(lib)
}

// ---------------------------------------------------------------------------
// Metadata tables
// ---------------------------------------------------------------------------

/// One row of `inchik.txt`: HMDB and LIPID MAPS annotation for an InChIKey.
#[derive(Clone, Debug, Default)]
pub struct Annotation {
    pub inchikey: String,
    pub hmdb_accession: String,
    pub hmdb_name: String,
    pub hmdb_super_class: String,
    pub hmdb_class: String,
    pub hmdb_sub_class: String,
    pub lm_id: String,
    pub lm_name: String,
    pub lm_abbrev: String,
    pub lm_core: String,
    pub lm_main_class: String,
    pub lm_sub_class: String,
    pub lm_abbrev_chains: String,
    pub formula: String,
}

/// The stereochemistry-independent first block of an InChIKey. `get` rather
/// than slicing so a malformed, non-ASCII key cannot panic mid-report.
pub fn skeleton(key: &str) -> &str {
    key.get(..14).unwrap_or(key)
}

/// Compound metadata, sorted by InChIKey for binary search.
#[derive(Default)]
pub struct Metadata {
    pub annotations: Vec<Annotation>,
    /// First 14 characters of each InChIKey (the skeleton), parallel to
    /// `annotations`, for the fallback lookup that ignores stereochemistry.
    pub skeletons: Vec<String>,
    /// `name -> formula`, sorted by name.
    pub name_formula: Vec<(String, String)>,
}

impl Metadata {
    /// Exact InChIKey matches.
    pub fn exact(&self, key: &str) -> &[Annotation] {
        let lo = self
            .annotations
            .partition_point(|a| a.inchikey.as_str() < key);
        let mut hi = lo;
        while hi < self.annotations.len() && self.annotations[hi].inchikey == key {
            hi += 1;
        }
        &self.annotations[lo..hi]
    }

    /// Skeleton (first block) matches, used when no exact key is present.
    pub fn by_skeleton(&self, key: &str) -> &[Annotation] {
        let skel = skeleton(key);
        let lo = self.skeletons.partition_point(|s| s.as_str() < skel);
        let mut hi = lo;
        while hi < self.skeletons.len() && self.skeletons[hi] == skel {
            hi += 1;
        }
        &self.annotations[lo..hi]
    }

    pub fn formula_for(&self, name: &str) -> Option<&str> {
        self.name_formula
            .binary_search_by(|probe| probe.0.as_str().cmp(name))
            .ok()
            .map(|i| self.name_formula[i].1.as_str())
    }
}

/// Load `inchik.txt` and `name_formu.txt`. Both are optional: without them the
/// reports simply carry no ontology columns, which beats 0.1's behaviour of
/// panicking at the very end of a long run.
pub fn load_metadata(params: &Params, reporter: &dyn Reporter) -> Metadata {
    let mut meta = Metadata::default();
    let Some(dir) = params.resolve_libs_dir() else {
        reporter.warn("no libs/ directory: compound ontology columns will be blank");
        return meta;
    };

    let inchik = dir.join("inchik.txt");
    if inchik.is_file() {
        match read_annotations(&inchik) {
            Ok(mut rows) => {
                rows.sort_unstable_by(|a, b| a.inchikey.cmp(&b.inchikey));
                meta.skeletons = rows
                    .iter()
                    .map(|a| skeleton(&a.inchikey).to_string())
                    .collect();
                meta.annotations = rows;
            }
            Err(e) => reporter.warn(format!("inchik.txt: {e}")),
        }
    } else {
        reporter.warn("inchik.txt not found: ontology columns will be blank");
    }

    let name_formu = dir.join("name_formu.txt");
    if name_formu.is_file() {
        match read_name_formula(&name_formu) {
            Ok(mut rows) => {
                rows.sort_unstable_by(|a, b| a.0.cmp(&b.0));
                meta.name_formula = rows;
            }
            Err(e) => reporter.warn(format!("name_formu.txt: {e}")),
        }
    }

    meta
}

fn read_annotations(path: &Path) -> Result<Vec<Annotation>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .trim(csv::Trim::All)
        .flexible(true)
        .from_path(path)?;
    let field = |row: &csv::StringRecord, i: usize| row.get(i).unwrap_or("").to_string();
    let mut out = Vec::new();
    for row in rdr.records() {
        let row = row?;
        out.push(Annotation {
            inchikey: field(&row, 0),
            hmdb_accession: field(&row, 1),
            hmdb_name: field(&row, 2),
            hmdb_super_class: field(&row, 3),
            hmdb_class: field(&row, 4),
            hmdb_sub_class: field(&row, 5),
            lm_id: field(&row, 6),
            lm_name: field(&row, 7),
            lm_abbrev: field(&row, 8),
            lm_core: field(&row, 9),
            lm_main_class: field(&row, 10),
            lm_sub_class: field(&row, 11),
            lm_abbrev_chains: field(&row, 12),
            formula: field(&row, 13),
        });
    }
    Ok(out)
}

fn read_name_formula(path: &Path) -> Result<Vec<(String, String)>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .trim(csv::Trim::All)
        .flexible(true)
        .from_path(path)?;
    let mut out = Vec::new();
    for row in rdr.records() {
        let row = row?;
        if let (Some(a), Some(b)) = (row.get(0), row.get(1)) {
            out.push((a.to_string(), b.to_string()));
        }
    }
    Ok(out)
}

/// Built-in libraries actually present in the resolved libs directory, for
/// each polarity. Lets the UI grey out what is not installed.
pub fn available_builtins(libs_dir: Option<&PathBuf>) -> Vec<(String, bool, bool)> {
    let Some(dir) = libs_dir else {
        return Vec::new();
    };
    crate::params::BUILTIN_LIBRARIES
        .iter()
        .map(|name| {
            let base = if *name == "Atlas_filtered" {
                "Atlas"
            } else {
                name
            };
            (
                (*name).to_string(),
                dir.join(format!("{base}_pos")).is_file(),
                dir.join(format!("{base}_neg")).is_file(),
            )
        })
        .collect()
}
