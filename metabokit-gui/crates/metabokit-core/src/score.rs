//! Spectral matching: pair MS2 scans with MS1 features, score them against
//! the library, and detect in-source fragments and isotopes.
//!
//! The scoring maths is a faithful port of 0.1. The structural changes are that
//! library entries are referred to by `u32` index rather than `&'a Ent` (which
//! removes the lifetime 0.1 threaded through every downstream type and lets
//! results be moved to the UI thread), and that spectra live in a columnar
//! [`SpecSet`] instead of a `Vec<Spec>` where each entry owned two `Vec`s.

use crate::error::Result;
use crate::features::Peak;
use crate::library::Library;
use crate::params::Params;
use crate::progress::Cancel;
use crate::scans::{Ms1Set, Ms2Set};

/// A library entry matched to one experimental spectrum.
#[derive(Clone, Copy, Debug)]
pub struct LibHit {
    /// Index into [`Library`].
    pub entry: u32,
    pub score: f32,
    pub matched_peaks: u8,
    /// Precursor m/z and retention time of the spectrum that produced the
    /// match, used to find it again when writing spectra.
    pub spec_prec_mz: f32,
    pub spec_rt: f32,
}

/// A related feature: either the parent an in-source fragment came from, or
/// the monoisotopic peak of an isotope.
#[derive(Clone, Copy, Debug, Default)]
pub struct Relation {
    pub mz: f32,
    pub rt: f32,
    pub score: f32,
    pub matched_peaks: u8,
}

/// One annotated feature (or one orphan spectrum) in one sample.
#[derive(Clone, Debug)]
pub struct Ann {
    /// Index of the sample within the run.
    pub file: u32,
    pub premz: f32,
    pub rt: f32,
    pub auc: f32,
    /// False when this came from an MS2 scan with no MS1 feature under it.
    pub is_feature: bool,
    pub s_n: f32,
    pub shape: f32,
    pub hits: Vec<LibHit>,
    /// Features this one appears to be an in-source fragment of.
    pub parents: Vec<Relation>,
    /// The monoisotopic peak, when this feature is an isotope of another.
    pub mono: Option<Relation>,
}

// ---------------------------------------------------------------------------
// Filtered MS2 spectra
// ---------------------------------------------------------------------------

/// Experimental spectra after filtering, in column form and sorted by
/// precursor m/z.
///
/// Two fragment lists are kept per spectrum: the *scored* list (mz-ascending,
/// truncated by the intensity cutoff or fragment count) and the *full* list
/// (intensity-descending, capped at 255) which is what gets exported.
pub struct SpecSet {
    prec_mz: Vec<f32>,
    rt: Vec<f32>,
    ce: Vec<f32>,
    scored_off: Vec<u32>,
    scored_mz: Vec<f32>,
    scored_i: Vec<f32>,
    full_off: Vec<u32>,
    full_mz: Vec<f32>,
    full_i: Vec<f32>,
}

impl Default for SpecSet {
    fn default() -> Self {
        SpecSet {
            prec_mz: Vec::new(),
            rt: Vec::new(),
            ce: Vec::new(),
            // Offset vectors must start with a sentinel zero: `scored(i)` and
            // `full(i)` both read `off[i + 1]`.
            scored_off: vec![0],
            scored_mz: Vec::new(),
            scored_i: Vec::new(),
            full_off: vec![0],
            full_mz: Vec::new(),
            full_i: Vec::new(),
        }
    }
}

impl SpecSet {
    fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.prec_mz.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.prec_mz.is_empty()
    }

    #[inline]
    pub fn prec_mz(&self, i: usize) -> f32 {
        self.prec_mz[i]
    }

    #[inline]
    pub fn rt(&self, i: usize) -> f32 {
        self.rt[i]
    }

    #[inline]
    pub fn ce(&self, i: usize) -> f32 {
        self.ce[i]
    }

    /// Fragments used for scoring: ascending m/z.
    #[inline]
    pub fn scored(&self, i: usize) -> (&[f32], &[f32]) {
        let a = self.scored_off[i] as usize;
        let b = self.scored_off[i + 1] as usize;
        (&self.scored_mz[a..b], &self.scored_i[a..b])
    }

    /// Full fragment list: descending intensity, capped at 255.
    #[inline]
    pub fn full(&self, i: usize) -> (&[f32], &[f32]) {
        let a = self.full_off[i] as usize;
        let b = self.full_off[i + 1] as usize;
        (&self.full_mz[a..b], &self.full_i[a..b])
    }

    #[inline]
    pub fn precursors(&self) -> &[f32] {
        &self.prec_mz
    }

    /// Locate the spectrum a [`LibHit`] came from.
    pub fn find(&self, prec_mz: f32, rt: f32) -> Option<usize> {
        let from = self.prec_mz.partition_point(|&x| x < prec_mz);
        (from..self.len())
            .take_while(|&i| self.prec_mz[i] == prec_mz)
            .find(|&i| self.rt[i] == rt)
    }

    pub fn heap_bytes(&self) -> usize {
        (self.prec_mz.capacity()
            + self.rt.capacity()
            + self.ce.capacity()
            + self.scored_off.capacity()
            + self.scored_mz.capacity()
            + self.scored_i.capacity()
            + self.full_off.capacity()
            + self.full_mz.capacity()
            + self.full_i.capacity())
            * 4
    }
}

/// Filter and rank the fragments of every MS2 scan.
fn filter_spectra(ms2: &Ms2Set, params: &Params) -> SpecSet {
    let mut out = SpecSet::new();
    // Reused across scans instead of a `Vec` per spectrum.
    let mut ranked: Vec<(f32, f32)> = Vec::new();

    for i in 0..ms2.len() {
        let prec = ms2.prec_mz(i);
        let (mzs, ints) = ms2.scan(i);
        // Fragments at or above the precursor are co-isolation, not fragments.
        let cutoff = prec - 0.3;
        ranked.clear();
        for k in 0..mzs.len() {
            if mzs[k] < cutoff {
                ranked.push((mzs[k], ints[k]));
            }
        }
        ranked.sort_unstable_by(|a, b| {
            b.1.partial_cmp(&a.1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        ranked.truncate(u8::MAX as usize);
        if ranked.is_empty() {
            continue;
        }

        let full_len = ranked.len();
        let scored_len = match params.intensity_cutoff {
            Some(cut) => ranked.iter().position(|x| x.1 < cut).unwrap_or(full_len),
            None => params.match_n_fragments.min(full_len),
        };
        if scored_len == 0 {
            continue;
        }

        // Full list first, while `ranked` is still intensity-ordered.
        for &(mz, inten) in &ranked {
            out.full_mz.push(mz);
            out.full_i.push(inten);
        }
        out.full_off.push(out.full_mz.len() as u32);

        ranked.truncate(scored_len);
        ranked.sort_unstable_by(|a, b| {
            a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
        });
        for &(mz, inten) in &ranked {
            out.scored_mz.push(mz);
            out.scored_i.push(inten);
        }
        out.scored_off.push(out.scored_mz.len() as u32);

        out.prec_mz.push(prec);
        out.rt.push(ms2.rt_at(i));
        out.ce.push(ms2.ce_at(i));
    }
    out
}

// ---------------------------------------------------------------------------
// Feature ↔ spectrum assignment
// ---------------------------------------------------------------------------

/// Spectra attached to a feature.
struct FeatureSpectra {
    peak: u32,
    specs: Vec<u32>,
}

/// Pair each MS1 feature with the MS2 scans acquired on it.
///
/// A scan counts if its precursor is within `ms1_ms2_pair` of the feature m/z
/// and it elutes within the feature's half-width. If nothing lands that
/// tightly, the window is relaxed to 1.5×.
fn assign_spectra(
    specs: &SpecSet,
    peaks: &[Peak],
    params: &Params,
) -> (Vec<FeatureSpectra>, Vec<u32>) {
    let mut groups: Vec<FeatureSpectra> = Vec::new();
    let mut attached: Vec<u32> = Vec::new();
    let mut candidates: Vec<(u32, f32)> = Vec::new();
    let precursors = specs.precursors();

    for (pi, peak) in peaks.iter().enumerate() {
        let from = precursors.partition_point(|&x| x < peak.mz - params.ms1_ms2_pair);
        let hi = peak.mz + params.ms1_ms2_pair;
        candidates.clear();
        for si in from..specs.len() {
            if precursors[si] >= hi {
                break;
            }
            candidates.push((si as u32, (peak.rt - specs.rt(si)).abs()));
        }
        let relaxed = peak.half_width * 1.5;
        for &(si, dt) in &candidates {
            if dt < relaxed {
                attached.push(si);
            }
        }
        let tight = candidates.iter().any(|x| x.1 < peak.half_width);
        let limit = if tight { peak.half_width } else { relaxed };
        let chosen: Vec<u32> = candidates
            .iter()
            .filter(|x| x.1 < limit)
            .map(|x| x.0)
            .collect();
        if !chosen.is_empty() {
            groups.push(FeatureSpectra {
                peak: pi as u32,
                specs: chosen,
            });
        }
    }

    let orphans = if params.features_only {
        Vec::new()
    } else {
        attached.sort_unstable();
        attached.dedup();
        (0..specs.len() as u32)
            .filter(|si| attached.binary_search(si).is_err())
            .collect()
    };

    (groups, orphans)
}

// ---------------------------------------------------------------------------
// In-source fragments
// ---------------------------------------------------------------------------

/// For each feature group, the groups it appears to be an in-source fragment of.
///
/// A lighter co-eluting feature is called an in-source fragment of a heavier
/// one when the heavier one's spectrum contains a fragment at the lighter
/// precursor's m/z with at least 10% of base-peak intensity, and the two
/// spectra share more than one of the lighter one's top fragments.
fn match_isf(groups: &[FeatureSpectra], specs: &SpecSet, peaks: &[Peak], params: &Params) -> Vec<Vec<(u32, u8)>> {
    let mut out: Vec<Vec<(u32, u8)>> = vec![Vec::new(); groups.len()];
    let tol = params.ms2_tol;

    // Descending so that `out[j]` accumulates parents from the heaviest down,
    // matching 0.1's ordering in the reports.
    for gi in (0..groups.len()).rev() {
        let parent_peak = peaks[groups[gi].peak as usize];
        for &si in &groups[gi].specs {
            let (heavy_mz, heavy_i) = specs.scored(si as usize);
            if heavy_mz.is_empty() {
                continue;
            }
            let base = heavy_i.iter().copied().fold(0.0f32, f32::max);
            let prec = specs.prec_mz(si as usize);

            // Only lighter features, separated by at least the configured gap.
            let limit = prec - params.isf_parent_mass_diff;
            let upto = groups.partition_point(|g| peaks[g.peak as usize].mz < limit);

            for gj in 0..upto {
                let light_peak = peaks[groups[gj].peak as usize];
                if (parent_peak.rt - light_peak.rt).abs() >= params.isf_rt_diff {
                    continue;
                }
                // Is the light precursor present in the heavy spectrum?
                let from = heavy_mz.partition_point(|&x| x < light_peak.mz - tol);
                let hi = light_peak.mz + tol;
                let mut survivor = 0.0f32;
                for k in from..heavy_mz.len() {
                    if heavy_mz[k] >= hi {
                        break;
                    }
                    if heavy_i[k] > survivor {
                        survivor = heavy_i[k];
                    }
                }
                if survivor <= 0.1 * base {
                    continue;
                }

                for &sj in &groups[gj].specs {
                    let (light_mz, _) = specs.full(sj as usize);
                    let light_prec = specs.prec_mz(sj as usize) + tol;
                    let mut shared = 0u8;
                    for &fmz in light_mz.iter().filter(|&&m| m < light_prec).take(10) {
                        let p = heavy_mz.partition_point(|&x| x < fmz - tol);
                        if heavy_mz.get(p).is_some_and(|&x| x < fmz + tol) {
                            shared += 1;
                        }
                    }
                    if shared > 1 {
                        // Keep the strongest evidence per parent rather than
                        // one record per contributing spectrum.
                        if out[gj].last().map(|e| e.0) == Some(gi as u32) {
                            let e = out[gj].last_mut().expect("just checked");
                            if e.1 < shared {
                                e.1 = shared;
                            }
                        } else {
                            out[gj].push((gi as u32, shared));
                        }
                    }
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Isotopes
// ---------------------------------------------------------------------------

const NEUTRON: f32 = 1.003_355;

/// Group features into isotope envelopes, returning `(isotope, monoisotopic)`
/// peak-index pairs sorted by the isotope's m/z.
fn group_isotopes(peaks: &[Peak], rt_window: f32) -> Vec<(u32, u32)> {
    let next_isotope = |from: usize, seed: &Peak| -> Option<usize> {
        let lo = seed.mz + NEUTRON - 0.003;
        let hi = seed.mz + NEUTRON + 0.003;
        (from + 1..peaks.len())
            .skip_while(|&i| peaks[i].mz < lo)
            .take_while(|&i| peaks[i].mz < hi)
            .find(|&i| {
                let c = &peaks[i];
                (c.rt - seed.rt).abs() < rt_window
                    // An isotope is never more intense than its monoisotope,
                    // and shares its chromatographic width.
                    && c.coef < seed.coef
                    && (c.half_width - seed.half_width).abs() < 1.01f32.max(seed.half_width * 0.11)
            })
    };

    let mut pairs: Vec<(u32, u32)> = Vec::new();
    for mono in 0..peaks.len() {
        let mut at = mono;
        // At most M+3, as in 0.1.
        for _ in 0..3 {
            match next_isotope(at, &peaks[at]) {
                Some(next) => {
                    pairs.push((next as u32, mono as u32));
                    at = next;
                }
                None => break,
            }
        }
    }
    pairs.sort_unstable_by(|a, b| {
        peaks[a.0 as usize]
            .mz
            .partial_cmp(&peaks[b.0 as usize].mz)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(&b.0))
    });
    pairs
}

// ---------------------------------------------------------------------------
// Library matching
// ---------------------------------------------------------------------------

/// Spectral similarity.
///
/// A square-root-weighted cosine: `Σ√(aᵢbᵢ) / √(Σaᵢ · Σbᵢ)`. The square root
/// compresses the dynamic range so a single dominant fragment cannot carry the
/// score on its own.
fn similarity(lib_i: &[f32], exp_i: &[f32]) -> f32 {
    let lib_sum: f32 = lib_i.iter().sum();
    let exp_sum: f32 = exp_i.iter().sum();
    if lib_sum <= 0.0 || exp_sum <= 0.0 {
        return 0.0;
    }
    let numerator: f32 = lib_i
        .iter()
        .zip(exp_i)
        .map(|(&a, &b)| (a * b).sqrt())
        .sum();
    numerator / (lib_sum * exp_sum).sqrt()
}

/// Scratch buffers for one call to [`match_library`].
#[derive(Default)]
struct MatchScratch {
    /// Indices of experimental fragments claimed by some library fragment.
    claimed: Vec<u32>,
    lib_i: Vec<f32>,
    exp_i: Vec<f32>,
}

/// Score every library entry whose precursor falls in `[lo, hi)` against one
/// spectrum, appending the ones that clear the score threshold.
#[allow(clippy::too_many_arguments)]
fn match_library(
    hits: &mut Vec<LibHit>,
    lib: &Library,
    (lo_idx, hi_mass): (usize, f32),
    specs: &SpecSet,
    si: usize,
    params: &Params,
    scratch: &mut MatchScratch,
) {
    let (exp_mz, exp_int) = specs.scored(si);
    let prec = specs.prec_mz(si);
    let spec_rt = specs.rt(si);
    let tol = params.ms2_tol;

    for entry in lo_idx..lib.len() {
        if lib.mass(entry) >= hi_mass {
            break;
        }
        // Entries carrying a reference RT must elute where they should.
        if let Some(ref_rt) = lib.rt(entry) {
            if (ref_rt - spec_rt).abs() >= params.rt_shift {
                continue;
            }
        }

        let (frag_mz, frag_i) = lib.fragments(entry);
        let charge = f32::from(lib.charge(entry));
        scratch.claimed.clear();
        scratch.lib_i.clear();
        scratch.exp_i.clear();

        let mut matched = 0u32;
        for (&fmz, &fi) in frag_mz.iter().zip(frag_i) {
            // Skip library fragments at or above the precursor.
            if charge * prec - fmz <= 0.1 {
                continue;
            }
            let from = exp_mz.partition_point(|&x| x < fmz - tol);
            let hi = fmz + tol;
            let mut best = 0.0f32;
            for k in from..exp_mz.len() {
                if exp_mz[k] >= hi {
                    break;
                }
                scratch.claimed.push(k as u32);
                if exp_int[k] > best {
                    best = exp_int[k];
                }
            }
            scratch.lib_i.push(fi);
            scratch.exp_i.push(best);
            if best > 0.0 {
                matched += 1;
            }
        }

        let matched = matched.min(u8::MAX as u32) as u8;
        if matched < params.min_peaks {
            continue;
        }

        if !params.chimeric_spectra {
            // Penalise experimental fragments the entry cannot explain by
            // scoring them against zero library intensity.
            scratch.claimed.sort_unstable();
            for k in 0..exp_mz.len() {
                if charge * prec - exp_mz[k] <= 0.1 {
                    break;
                }
                if scratch.claimed.binary_search(&(k as u32)).is_err() {
                    scratch.lib_i.push(0.0);
                    scratch.exp_i.push(exp_int[k]);
                }
            }
        }

        let score = similarity(&scratch.lib_i, &scratch.exp_i);
        if score <= params.ms2_score {
            continue;
        }

        let hit = LibHit {
            entry: entry as u32,
            score,
            matched_peaks: matched,
            spec_prec_mz: prec,
            spec_rt,
        };
        if params.top_scoring_only {
            if hits.is_empty() {
                hits.push(hit);
            } else if hits[0].score < score {
                hits[0] = hit;
            }
        } else {
            hits.push(hit);
        }
    }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// Extra chromatogram width, in minutes, kept either side of the integration
/// bounds so the signal-to-noise estimate has baseline to work with.
const BASELINE_WIDTH: f32 = 0.3;

pub struct SampleResult {
    pub annotations: Vec<Ann>,
    pub spectra: SpecSet,
}

/// Score one sample.
pub fn score_sample(
    file_index: u32,
    ms1: &Ms1Set,
    ms2: &Ms2Set,
    peaks: &[Peak],
    lib: &Library,
    params: &Params,
    cancel: &Cancel,
) -> Result<SampleResult> {
    let specs = filter_spectra(ms2, params);
    let (groups, orphans) = assign_spectra(&specs, peaks, params);
    let isf = match_isf(&groups, &specs, peaks, params);
    let isotopes = group_isotopes(peaks, params.isf_rt_diff);

    let mut annotations: Vec<Ann> = Vec::new();
    let mut hits: Vec<LibHit> = Vec::new();
    let mut chrom: Vec<(f32, f32)> = Vec::new();
    let mut scratch = MatchScratch::default();

    for (gi, group) in groups.iter().enumerate() {
        if gi % 512 == 0 {
            cancel.check()?;
        }
        let peak = peaks[group.peak as usize];
        let half = peak.half_width * params.integration_rt;
        ms1.xic(
            peak.mz,
            peak.rt,
            half + BASELINE_WIDTH,
            params.integration_mz,
            &mut chrom,
        );
        let lo = chrom.partition_point(|x| x.0 < peak.rt - half);
        let hi = chrom.partition_point(|x| x.0 < peak.rt + half);
        let Some(s_n) = signal_to_noise(&chrom, lo, hi, peak.rt, peak.half_width) else {
            continue;
        };
        if s_n < params.s_n_1 {
            continue;
        }

        let auc = integrate(&chrom, lo, hi);
        let tol = params.ms1_tol_unit.absolute(params.ms1_tol, peak.mz);
        let lo_idx = lib.mass_at_or_after(peak.mz - tol);

        hits.clear();
        for &si in &group.specs {
            match_library(
                &mut hits,
                lib,
                (lo_idx, peak.mz + tol),
                &specs,
                si as usize,
                params,
                &mut scratch,
            );
        }

        let parents = &isf[gi];
        if hits.is_empty() && parents.is_empty() {
            continue;
        }

        let mono = isotope_parent(&isotopes, peaks, group.peak);
        annotations.push(Ann {
            file: file_index,
            premz: peak.mz,
            rt: peak.rt,
            auc,
            is_feature: true,
            s_n,
            shape: peak.shape,
            hits: hits.clone(),
            parents: parents
                .iter()
                .map(|&(gj, matched)| {
                    let p = peaks[groups[gj as usize].peak as usize];
                    Relation {
                        mz: p.mz,
                        rt: p.rt,
                        score: 0.0,
                        matched_peaks: matched,
                    }
                })
                .collect(),
            mono,
        });
    }

    // Spectra with no feature underneath: quantify a narrow window around the
    // scan itself so they still carry an area.
    for (n, &si) in orphans.iter().enumerate() {
        if n % 512 == 0 {
            cancel.check()?;
        }
        let si = si as usize;
        let prec = specs.prec_mz(si);
        let rt = specs.rt(si);
        ms1.xic(
            prec,
            rt,
            params.impute_width + BASELINE_WIDTH,
            params.integration_mz,
            &mut chrom,
        );
        let lo = chrom.partition_point(|x| x.0 < rt - params.impute_width);
        let hi = chrom.partition_point(|x| x.0 < rt + params.impute_width);
        let auc = integrate(&chrom, lo, hi);

        let tol = params.ms1_tol_unit.absolute(params.ms1_tol, prec);
        let lo_idx = lib.mass_at_or_after(prec - tol);
        hits.clear();
        match_library(
            &mut hits,
            lib,
            (lo_idx, prec + tol),
            &specs,
            si,
            params,
            &mut scratch,
        );
        if hits.is_empty() {
            continue;
        }
        annotations.push(Ann {
            file: file_index,
            premz: prec,
            rt,
            auc,
            is_feature: false,
            s_n: 0.0,
            shape: 0.0,
            hits: hits.clone(),
            parents: Vec::new(),
            mono: None,
        });
    }

    Ok(SampleResult {
        annotations,
        spectra: specs,
    })
}

/// The monoisotopic partner of `peak`, if it is an isotope.
fn isotope_parent(isotopes: &[(u32, u32)], peaks: &[Peak], peak: u32) -> Option<Relation> {
    let mz = peaks[peak as usize].mz;
    let from = isotopes.partition_point(|&(iso, _)| peaks[iso as usize].mz < mz);
    for &(iso, mono) in &isotopes[from..] {
        if peaks[iso as usize].mz != mz {
            break;
        }
        if iso == peak {
            let m = peaks[mono as usize];
            return Some(Relation {
                mz: m.mz,
                rt: m.rt,
                score: 0.0,
                matched_peaks: 0,
            });
        }
    }
    None
}

/// Trapezoidal area under `chrom[lo..hi]`, scaled to match 0.1's units.
fn integrate(chrom: &[(f32, f32)], lo: usize, hi: usize) -> f32 {
    if hi <= lo + 1 || hi > chrom.len() {
        return 0.0;
    }
    let mut acc = 0.0f32;
    for k in lo..hi - 1 {
        acc += (chrom[k].1 + chrom[k + 1].1) * (chrom[k + 1].0 - chrom[k].0);
    }
    acc * 30.0
}

/// Apex intensity over mean baseline intensity, floored at 1.
///
/// Returns `None` when the integration window holds no points — 0.1 unwrapped
/// here and panicked on features at the very start or end of a run.
fn signal_to_noise(
    chrom: &[(f32, f32)],
    lo: usize,
    hi: usize,
    rt: f32,
    half_width: f32,
) -> Option<f32> {
    if hi <= lo || hi > chrom.len() {
        return None;
    }
    let window = &chrom[lo..hi];
    let a = window.partition_point(|x| x.0 < rt - half_width);
    let b = window.partition_point(|x| x.0 < rt + half_width);
    if b <= a {
        return None;
    }
    let apex = window[a..b].iter().map(|x| x.1).fold(f32::MIN, f32::max);
    let flank_n = (a + window.len() - b).max(1);
    let flank_sum: f32 = window[..a]
        .iter()
        .chain(&window[b..])
        .map(|x| x.1)
        .sum::<f32>();
    let baseline = (flank_sum / flank_n as f32).max(1.0);
    Some(apex / baseline)
}
