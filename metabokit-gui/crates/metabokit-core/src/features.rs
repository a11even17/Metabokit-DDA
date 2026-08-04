//! MS1 feature detection by continuous wavelet transform.
//!
//! The algorithm is unchanged from 0.1: for each narrow m/z slice that has both
//! an MS2 scan and a library candidate, build an extracted ion chromatogram,
//! convolve it with Mexican-hat wavelets across a range of scales, link local
//! maxima into ridge lines across scales, and keep ridges whose apex looks like
//! a chromatographic peak.
//!
//! What changed is how it is executed:
//!
//! * **Parallel over slices, not over files.** 0.1 ran `mzml_fs.par_iter()`, so
//!   memory scaled with core count — eight threads meant eight decoded mzML
//!   files resident. m/z slices are independent, so parallelising here gives
//!   the same utilisation with exactly one sample in RAM. `map_init` gives each
//!   worker a scratch arena it reuses for every slice it handles.
//! * **One coefficient row instead of the full matrix.** 0.1 allocated
//!   `eic_rt.len() * scales.len()` floats per slice, filled all of it, and only
//!   then extracted maxima scale by scale. Nothing crosses scales, so a single
//!   row suffices — roughly a 40× cut in the hottest buffer, which keeps it in
//!   L2 instead of streaming from RAM.
//! * **Flat ridges.** Ridge lines were `Vec<Vec<RtScC>>`, one heap allocation
//!   per ridge per slice. Only the apex and the tail are ever read, so each
//!   ridge is now a fixed 32-byte record in one flat vector.

use rayon::prelude::*;

use crate::error::Result;
use crate::params::Params;
use crate::progress::Cancel;
use crate::scans::{Ms1Set, Ms2Set};

/// A detected chromatographic feature.
#[derive(Clone, Copy, Debug)]
pub struct Peak {
    pub mz: f32,
    pub rt: f32,
    /// Fitted wavelet scale — the peak's half-width in minutes.
    pub half_width: f32,
    /// Wavelet coefficient at the apex; used as a relative intensity rank.
    pub coef: f32,
    /// Cosine similarity between the chromatogram and the fitted wavelet.
    pub shape: f32,
    /// Lag-1 autocorrelation of the chromatogram; a smoothness proxy.
    pub smooth: f32,
}

/// Width of one m/z detection slice.
const MZ_STEP: f32 = 0.007;
/// A slice spans three steps.
const SLICE_SPAN: f32 = MZ_STEP * 3.0;
/// Minimum number of scales a ridge must persist across.
const MIN_RIDGE: usize = 3;
/// How far from a sample point an MS2 scan may sit for that point to be worth
/// transforming, in minutes. Fixed at 0.25 in 0.1 and kept fixed here.
const RT_SEARCH: f32 = 0.25;

/// Per-worker reusable buffers. Allocated once per rayon thread, not once per
/// slice: a run has tens of thousands of slices.
struct Scratch {
    /// Per-scan m/z and (noise-subtracted) intensity of the slice's trace,
    /// indexed by MS1 scan index. Columnar so the intensity pass vectorises.
    eic_mz: Vec<f32>,
    eic_i: Vec<f32>,
    trace: Vec<(u32, f32, f32)>,
    msms_rt: Vec<f32>,
    noise: Vec<f32>,
    /// Wavelet sample points, midway between consecutive scans.
    eic_rt: Vec<f32>,
    /// Whether a sample point is near enough to an MS2 scan to be worth
    /// evaluating.
    active: Vec<bool>,
    /// Coefficients for the scale currently being processed.
    row: Vec<f32>,
    ridges: Vec<Ridge>,
    conv: Vec<f32>,
    chrom: Vec<(f32, f32)>,
    model: Vec<f32>,
    keep: Vec<f32>,
    out: Vec<Peak>,
}

impl Scratch {
    fn new() -> Self {
        Scratch {
            eic_mz: Vec::new(),
            eic_i: Vec::new(),
            trace: Vec::new(),
            msms_rt: Vec::new(),
            noise: Vec::new(),
            eic_rt: Vec::new(),
            active: Vec::new(),
            row: Vec::new(),
            ridges: Vec::new(),
            conv: Vec::new(),
            chrom: Vec::new(),
            model: Vec::new(),
            keep: Vec::new(),
            out: Vec::new(),
        }
    }
}

/// A ridge line across scales.
///
/// Only the apex and the tail are ever inspected, so the line is summarised
/// rather than stored.
#[derive(Clone, Copy)]
struct Ridge {
    last_rt: f32,
    last_scale: f32,
    len: u32,
    best_coef: f32,
    best_rt: f32,
    best_scale: f32,
    /// Position of the apex within the ridge.
    best_pos: u32,
}

/// Wavelet scales, in minutes. Fine near the bottom, geometric above.
fn build_scales(max_peak_width: f32) -> (Vec<f32>, Vec<f32>) {
    let mut scales: Vec<f32> = vec![0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09];
    scales.extend((0..40).map(|i| 0.1 * 1.1f32.powi(i)));
    let cut = scales.partition_point(|&x| x < max_peak_width / 2.0) + 3;
    scales.truncate(cut.min(scales.len()));
    let sqrt = scales.iter().map(|x| x.sqrt()).collect();
    (scales, sqrt)
}

/// Detect features across the whole m/z range covered by MS2 acquisition.
pub fn detect(
    ms1: &Ms1Set,
    ms2: &Ms2Set,
    lib_masses: &[f32],
    params: &Params,
    cancel: &Cancel,
) -> Result<Vec<Peak>> {
    if ms1.is_empty() || ms2.is_empty() {
        return Ok(Vec::new());
    }
    let (scales, scale_sqrt) = build_scales(params.peak_width.1);
    if scales.is_empty() {
        return Ok(Vec::new());
    }

    let precursors = ms2.precursors();
    let start = precursors[0] - 0.02;
    let end = precursors[precursors.len() - 1];
    if !(end > start) {
        return Ok(Vec::new());
    }
    let n_slices = (((end - start) / MZ_STEP).ceil() as usize).max(1);

    // Computing `mz0` as `start + i * step` rather than accumulating keeps the
    // slice boundaries exact and independent of iteration order, so a parallel
    // run reproduces a serial one bit for bit.
    let chunks: Vec<Vec<Peak>> = (0..n_slices)
        .into_par_iter()
        .map_init(Scratch::new, |scratch, i| {
            if cancel.is_cancelled() {
                return Vec::new();
            }
            let mz0 = start + i as f32 * MZ_STEP;
            let mz1 = mz0 + SLICE_SPAN;
            scratch.out.clear();
            process_slice(
                ms1,
                ms2,
                lib_masses,
                (mz0, mz1),
                &scales,
                &scale_sqrt,
                params,
                scratch,
            );
            std::mem::take(&mut scratch.out)
        })
        .collect();

    cancel.check()?;

    let mut peaks: Vec<Peak> = chunks.into_iter().flatten().collect();
    peaks.sort_unstable_by(|a, b| {
        a.mz.partial_cmp(&b.mz)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.rt.partial_cmp(&b.rt).unwrap_or(std::cmp::Ordering::Equal))
    });
    Ok(deduplicate(peaks))
}

/// Adjacent m/z slices overlap, so the same ion is detected more than once.
/// Keep the strongest of each co-eluting cluster.
fn deduplicate(peaks: Vec<Peak>) -> Vec<Peak> {
    let n = peaks.len();
    let mut keep = vec![true; n];
    let mut by_coef: Vec<u32> = (0..n as u32).collect();
    by_coef.sort_unstable_by(|&a, &b| {
        peaks[b as usize]
            .coef
            .partial_cmp(&peaks[a as usize].coef)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });

    for &idx in &by_coef {
        let pos = idx as usize;
        if !keep[pos] {
            continue;
        }
        let peak = peaks[pos];
        let lo = peak.mz - 0.009;
        let hi = peak.mz + 0.009;

        let mut j = pos;
        while j > 0 {
            j -= 1;
            if peaks[j].mz <= lo {
                break;
            }
            if (peaks[j].rt - peak.rt).abs() < peaks[j].half_width + peak.half_width {
                keep[j] = false;
            }
        }
        for (j, other) in peaks.iter().enumerate().skip(pos + 1) {
            if other.mz >= hi {
                break;
            }
            if (other.rt - peak.rt).abs() < other.half_width + peak.half_width {
                keep[j] = false;
            }
        }
    }

    peaks
        .into_iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(p, _)| p)
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn process_slice(
    ms1: &Ms1Set,
    ms2: &Ms2Set,
    lib_masses: &[f32],
    (mz0, mz1): (f32, f32),
    scales: &[f32],
    scale_sqrt: &[f32],
    params: &Params,
    s: &mut Scratch,
) {
    // Only slices with MS2 coverage are worth the wavelet transform.
    let lo = mz0 - params.ms1_ms2_pair;
    let hi = mz1 + params.ms1_ms2_pair;
    let precursors = ms2.precursors();
    let from = ms2.prec_at_or_after(lo);
    s.msms_rt.clear();
    for i in from..precursors.len() {
        if precursors[i] > hi {
            break;
        }
        s.msms_rt.push(ms2.rt_at(i));
    }
    if s.msms_rt.is_empty() {
        return;
    }

    // …and only slices a library entry could explain.
    let lib_from = lib_masses.partition_point(|&x| x < mz0 - MZ_STEP);
    match lib_masses.get(lib_from) {
        Some(&m) if m <= mz1 + MZ_STEP => {}
        _ => return,
    }

    ms1.slice_trace(mz0, mz1, &mut s.trace);
    if s.trace.len() <= 2 {
        return;
    }

    find_ridges(ms1, (mz0, mz1), scales, scale_sqrt, params, s);
}

fn find_ridges(
    ms1: &Ms1Set,
    mz_range: (f32, f32),
    scales: &[f32],
    scale_sqrt: &[f32],
    params: &Params,
    s: &mut Scratch,
) {
    let rt_all = ms1.rts();
    let n_scans = rt_all.len();

    // Scatter the slice trace into per-scan columns; scans with no signal stay
    // at zero so the chromatogram is evenly sampled.
    s.eic_mz.clear();
    s.eic_i.clear();
    s.eic_mz.resize(n_scans, 0.0);
    s.eic_i.resize(n_scans, 0.0);
    for &(idx, mz, inten) in &s.trace {
        let idx = idx as usize;
        if idx < n_scans {
            s.eic_mz[idx] = mz;
            s.eic_i[idx] = inten;
        }
    }

    // Noise floor: the 5th percentile of the observed intensities.
    s.noise.clear();
    s.noise.extend(s.trace.iter().map(|x| x.2));
    let k = s.noise.len() / 20;
    let noise = if s.noise.is_empty() {
        0.0
    } else {
        s.noise.select_nth_unstable_by(k, |a, b| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        });
        s.noise[k]
    };

    s.msms_rt
        .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Sample the transform between scans, and only where an MS2 scan was
    // acquired nearby — evaluating the wavelet where nothing was fragmented
    // cannot produce an annotation.
    s.eic_rt.clear();
    s.active.clear();
    let mut prev_active = true;
    for i in 0..n_scans.saturating_sub(1) {
        let mid = (rt_all[i] + rt_all[i + 1]) * 0.5;
        let from = s.msms_rt.partition_point(|&x| x < mid - RT_SEARCH);
        let near = s.msms_rt.get(from).is_some_and(|&x| x < mid + RT_SEARCH);
        if near || prev_active {
            prev_active = near;
            s.eic_rt.push(mid);
            s.active.push(near);
        }
    }
    if s.eic_rt.len() < 3 {
        return;
    }
    // The endpoints have no neighbours to test for a local maximum.
    let last = s.active.len() - 1;
    s.active[0] = false;
    s.active[last] = false;

    for v in &mut s.eic_i {
        *v = (*v - noise).max(0.0);
    }

    let n_pts = s.eic_rt.len();
    s.ridges.clear();

    for (scale_idx, (&scale, &sqrt_scale)) in scales.iter().zip(scale_sqrt).enumerate() {
        // --- transform at this scale ---------------------------------------
        s.row.clear();
        s.row.resize(n_pts, 0.0);
        let inv_scale = 1.0 / scale;

        for i in 0..n_pts {
            if !s.active[i] {
                continue;
            }
            let centre = s.eic_rt[i];
            let a = rt_all.partition_point(|&x| x < centre - scale);
            let b = rt_all.partition_point(|&x| x < centre + scale);
            if b < a + 2 {
                continue;
            }
            s.conv.clear();
            s.conv.reserve(b - a);
            for k in a..b {
                let t = (rt_all[k] - centre) * inv_scale;
                let t2 = t * t;
                // Mexican hat (second derivative of a Gaussian).
                s.conv.push(s.eic_i[k] * (-t2 * 0.5).exp() * (1.0 - t2));
            }
            let mut acc = 0.0f32;
            for k in 0..s.conv.len() - 1 {
                acc += (s.conv[k] + s.conv[k + 1]) * (rt_all[a + k + 1] - rt_all[a + k]);
            }
            s.row[i] = acc / sqrt_scale;
        }

        // --- local maxima, strongest first ---------------------------------
        // Each accepted maximum suppresses everything within one scale width,
        // so a broad peak yields a single ridge point per scale.
        loop {
            let mut max_i = 0usize;
            let mut max_v = f32::NEG_INFINITY;
            for (i, &v) in s.row.iter().enumerate() {
                if v >= max_v {
                    max_v = v;
                    max_i = i;
                }
            }
            if !(max_v > 0.0) {
                break;
            }
            if max_i == 0
                || max_i + 1 >= n_pts
                || s.row[max_i - 1] <= 0.0
                || s.row[max_i + 1] <= 0.0
            {
                // A maximum sitting on a zero shoulder is an edge artefact.
                s.row[max_i] = 0.0;
                continue;
            }

            let rt = s.eic_rt[max_i];
            link_ridge(&mut s.ridges, scale_idx, scales, rt, max_v, rt_all);

            let lo = rt - scale;
            let hi = rt + scale;
            let a = s.eic_rt[..max_i].partition_point(|&x| x <= lo);
            let b = max_i + s.eic_rt[max_i..].partition_point(|&x| x < hi);
            for v in &mut s.row[a..b] {
                *v = 0.0;
            }
        }
    }

    // --- evaluate the ridges ------------------------------------------------
    // Copy out first: `evaluate_ridge` needs `&mut Scratch` for its own
    // buffers, which it cannot have while borrowing `s.ridges`.
    let ridges = std::mem::take(&mut s.ridges);
    for r in &ridges {
        if let Some(peak) = evaluate_ridge(r, ms1, mz_range, params, s) {
            s.out.push(peak);
        }
    }
    s.ridges = ridges;
}

/// Attach a local maximum to an existing ridge, or start a new one.
fn link_ridge(
    ridges: &mut Vec<Ridge>,
    scale_idx: usize,
    scales: &[f32],
    rt: f32,
    coef: f32,
    rt_all: &[f32],
) {
    if scale_idx > 0 {
        let prev_scale = scales[scale_idx - 1];
        let scan_of = |t: f32| rt_all.partition_point(|&x| x < t);
        if let Some(r) = ridges.iter_mut().find(|r| {
            r.last_scale == prev_scale
                && ((r.last_rt - rt).abs() < 0.01 || scan_of(rt).abs_diff(scan_of(r.last_rt)) < 2)
        }) {
            r.last_rt = rt;
            r.last_scale = scales[scale_idx];
            r.len += 1;
            // `>=` so ties resolve to the coarser scale, matching the `max_by`
            // this replaced.
            if coef >= r.best_coef {
                r.best_coef = coef;
                r.best_rt = rt;
                r.best_scale = scales[scale_idx];
                r.best_pos = r.len - 1;
            }
            return;
        }
    }
    ridges.push(Ridge {
        last_rt: rt,
        last_scale: scales[scale_idx],
        len: 1,
        best_coef: coef,
        best_rt: rt,
        best_scale: scales[scale_idx],
        best_pos: 0,
    });
}

fn evaluate_ridge(
    r: &Ridge,
    ms1: &Ms1Set,
    mz_range: (f32, f32),
    params: &Params,
    s: &mut Scratch,
) -> Option<Peak> {
    if (r.len as usize) < MIN_RIDGE {
        return None;
    }
    let width = r.best_scale * 2.0;
    if width < params.peak_width.0 || width > params.peak_width.1 {
        return None;
    }
    // The apex must be followed by at least two coarser scales; a ridge that
    // peaks at its own tail is noise growing with scale.
    if r.best_pos as usize + MIN_RIDGE > r.len as usize {
        return None;
    }

    let rt_all = ms1.rts();
    let a = rt_all.partition_point(|&x| x < r.best_rt - r.best_scale);
    let b = rt_all.partition_point(|&x| x < r.best_rt + r.best_scale);
    if a == 0 || b >= rt_all.len() || b <= a {
        return None;
    }

    // Intensity-weighted centroid of the apex region.
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for k in a..b {
        num += s.eic_mz[k] * s.eic_i[k];
        den += s.eic_i[k];
    }
    if !(den > 0.0) {
        // 0.1 divided regardless and let a NaN m/z flow downstream, where it
        // silently passed both range tests and produced an empty chromatogram.
        return None;
    }
    let peak_mz = num / den;
    if peak_mz < mz_range.0 + 0.005 || peak_mz > mz_range.1 - 0.005 {
        return None;
    }

    // Re-extract at the refined m/z, wider than the peak so the tails are
    // available for the smoothness estimate.
    ms1.xic(
        peak_mz,
        r.best_rt,
        3.0 * r.best_scale,
        params.integration_mz,
        &mut s.chrom,
    );
    if s.chrom.len() < 3 {
        return None;
    }

    // Smoothness: lag-1 autocorrelation over the chromatogram, ignoring points
    // in a run of three consecutive zeros (detector gaps rather than shape).
    s.keep.clear();
    let n = s.chrom.len();
    for i in 0..n {
        let live = if i == 0 || i + 1 >= n {
            true
        } else {
            s.chrom[i - 1].1 > 0.0 || s.chrom[i].1 > 0.0 || s.chrom[i + 1].1 > 0.0
        };
        if live {
            s.keep.push(s.chrom[i].1);
        }
    }
    let smooth = autocorrelation(&mut s.keep)?;
    if smooth < 0.0 {
        return None;
    }

    // Shape: cosine similarity between the chromatogram and the fitted
    // wavelet, both mean-centred, integrated over retention time.
    let left = r.best_rt - r.best_scale;
    let right = r.best_rt + r.best_scale;
    let lo = s.chrom.partition_point(|x| x.0 < left);
    let hi = s.chrom.partition_point(|x| x.0 < right);
    if hi <= lo + 1 {
        return None;
    }
    s.chrom.truncate(hi);
    s.chrom.drain(..lo);
    let m = s.chrom.len();
    if m < 2 {
        return None;
    }

    s.model.clear();
    s.model.reserve(m);
    for &(rt, _) in s.chrom.iter() {
        let t = (rt - r.best_rt) / r.best_scale;
        let t2 = t * t;
        s.model.push((1.0 - t2) * (-t2 * 0.5).exp());
    }
    let model_mean = s.model.iter().sum::<f32>() / m as f32;
    for v in &mut s.model {
        *v -= model_mean;
    }
    let chrom_mean = s.chrom.iter().map(|x| x.1).sum::<f32>() / m as f32;
    for c in s.chrom.iter_mut() {
        c.1 -= chrom_mean;
    }

    let mut dot = 0.0f32;
    let mut model_mag = 0.0f32;
    let mut chrom_mag = 0.0f32;
    for k in 0..m - 1 {
        let dt = s.chrom[k + 1].0 - s.chrom[k].0;
        let (c0, c1) = (s.chrom[k].1, s.chrom[k + 1].1);
        let (m0, m1) = (s.model[k], s.model[k + 1]);
        dot += (c0 * m0 + c1 * m1) * dt;
        model_mag += (m0 * m0 + m1 * m1) * dt;
        chrom_mag += (c0 * c0 + c1 * c1) * dt;
    }
    let denom = (model_mag * chrom_mag).sqrt();
    if !(denom > 0.0) {
        return None;
    }
    let shape = dot / denom;

    (shape > params.peak_shape && smooth > params.sn_score).then_some(Peak {
        mz: peak_mz,
        rt: r.best_rt,
        half_width: r.best_scale,
        coef: r.best_coef,
        shape,
        smooth,
    })
}

/// Lag-1 autocorrelation. Mutates `data` in place (mean-centres it).
fn autocorrelation(data: &mut [f32]) -> Option<f32> {
    if data.len() < 2 {
        return None;
    }
    let mean = data.iter().sum::<f32>() / data.len() as f32;
    for d in data.iter_mut() {
        *d -= mean;
    }
    let variance: f32 = data.iter().map(|&x| x * x).sum();
    if variance == 0.0 {
        return Some(0.0);
    }
    let covariance: f32 = data.windows(2).map(|c| c[0] * c[1]).sum();
    Some(covariance / variance)
}
