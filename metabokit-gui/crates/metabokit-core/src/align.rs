//! Cross-sample alignment.
//!
//! Annotations from every sample are pooled and grouped into consensus
//! features: same m/z within `mz_shift`, co-eluting within `rt_shift`. The same
//! greedy, area-weighted procedure then runs again over in-source-fragment
//! relationships so an ISF appears once per parent rather than once per sample.
//!
//! 0.1 did this with `Vec<&'a Ann<'a>>` everywhere, which forced the library,
//! the annotations and the aligner to share one lifetime and made the whole
//! result graph impossible to move off the worker thread. Groups here are
//! `u32` indices into one flat `Vec<Ann>`: half the size, trivially `Send`,
//! and immune to the aliasing problems that shape the original code.

use crate::params::Params;
use crate::score::Ann;

/// The pooled, grouped result of a run.
pub struct Aligned {
    /// Every annotation from every sample, sorted by precursor m/z.
    pub annotations: Vec<Ann>,
    /// Consensus features: indices into `annotations`, sorted by median m/z.
    pub groups: Vec<Vec<u32>>,
    /// `(in-source fragment, parent)` pairs, ordered by the fragment's m/z.
    pub isf_pairs: Vec<(u32, u32)>,
    /// The same pairs ordered by the *parent's* m/z, for reverse lookup.
    pub isf_by_parent: Vec<(u32, u32)>,
    /// ISF relationships grouped across samples, sorted by median parent m/z.
    pub isf_groups: Vec<Vec<(u32, u32)>>,
}

/// Pool per-sample annotations and group them.
pub fn align(mut annotations: Vec<Ann>, params: &Params) -> Aligned {
    annotations.sort_unstable_by(|a, b| {
        a.premz
            .partial_cmp(&b.premz)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.rt.partial_cmp(&b.rt)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.file.cmp(&b.file))
    });

    let groups = consensus_features(&annotations, params);
    let isf_pairs = associate_isf(&annotations, params);

    let mut isf_by_parent = isf_pairs.clone();
    isf_by_parent.sort_unstable_by(|a, b| {
        annotations[a.1 as usize]
            .premz
            .partial_cmp(&annotations[b.1 as usize].premz)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.cmp(&b.1))
    });

    let isf_groups = group_isf(&annotations, &isf_pairs, params.rt_shift);

    Aligned {
        annotations,
        groups,
        isf_pairs,
        isf_by_parent,
        isf_groups,
    }
}

/// Greedily grow consensus features, strongest first.
///
/// Working from the largest peak area down means the most reliable
/// measurement, rather than acquisition order, decides where a group's
/// retention-time window starts.
fn consensus_features(annotations: &[Ann], params: &Params) -> Vec<Vec<u32>> {
    let n = annotations.len();
    // Annotations with no library hit can still join a group, but never seed
    // one — an unidentified area should not define a consensus feature.
    let mut available: Vec<bool> = annotations.iter().map(|a| !a.hits.is_empty()).collect();

    let mut by_area: Vec<u32> = (0..n as u32).collect();
    by_area.sort_unstable_by(|&a, &b| {
        annotations[b as usize]
            .auc
            .partial_cmp(&annotations[a as usize].auc)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });

    let mut groups: Vec<(f32, Vec<u32>)> = Vec::new();
    let mut medians: Vec<f32> = Vec::new();

    for &seed in &by_area {
        let pos = seed as usize;
        if !available[pos] {
            continue;
        }
        let anchor_rt = annotations[pos].rt;
        let lo = annotations[pos].premz - params.mz_shift;
        let hi = annotations[pos].premz + params.mz_shift;
        let from = annotations[..pos].partition_point(|a| a.premz < lo);
        let to = pos + 1 + annotations[pos + 1..].partition_point(|a| a.premz < hi);

        // Repeat until no remaining window in the m/z band still covers the
        // seed's retention time.
        loop {
            let mut best: Option<(usize, f32)> = None;
            for start in from..to {
                if !available[start] {
                    continue;
                }
                let start_rt = annotations[start].rt;
                if !(start_rt <= anchor_rt && anchor_rt < start_rt + params.rt_shift) {
                    continue;
                }
                let rt_end = start_rt + params.rt_shift;
                let mut total = 0.0f32;
                for k in from..to {
                    if available[k] && annotations[k].rt >= start_rt && annotations[k].rt < rt_end {
                        total += annotations[k].auc;
                    }
                }
                // `>=` so ties resolve to the later window, matching the
                // `max_by` this replaced.
                if best.map_or(true, |(_, b)| total >= b) {
                    best = Some((start, total));
                }
            }
            let Some((start, _)) = best else { break };

            let start_rt = annotations[start].rt;
            let rt_end = start_rt + params.rt_shift;
            let mut members = Vec::new();
            medians.clear();
            for k in from..to {
                if available[k] && annotations[k].rt >= start_rt && annotations[k].rt < rt_end {
                    available[k] = false;
                    members.push(k as u32);
                    medians.push(annotations[k].premz);
                }
            }
            if members.is_empty() {
                break;
            }
            groups.push((median(&mut medians), members));
        }
    }

    groups.sort_unstable_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.first().cmp(&b.1.first()))
    });
    groups.into_iter().map(|x| x.1).collect()
}

/// Resolve each annotation's recorded parent m/z and RT to an actual
/// annotation in the same sample.
///
/// The parent must itself be confidently identified — a poorly-scoring parent
/// would otherwise drag a whole ISF family into the report.
fn associate_isf(annotations: &[Ann], params: &Params) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for (i, ann) in annotations.iter().enumerate() {
        for rel in &ann.parents {
            let from = annotations.partition_point(|a| a.premz < rel.mz);
            for j in from..annotations.len() {
                let cand = &annotations[j];
                if cand.premz != rel.mz {
                    break;
                }
                if cand.rt != rel.rt || cand.file != ann.file {
                    continue;
                }
                let confident = cand.hits.iter().any(|h| {
                    h.score > params.ms2_score && h.matched_peaks >= params.min_peaks
                });
                if confident {
                    out.push((i as u32, j as u32));
                    break;
                }
            }
        }
    }
    out
}

/// Group ISF relationships across samples: same fragment m/z, same parent m/z,
/// co-eluting parents.
fn group_isf(annotations: &[Ann], pairs: &[(u32, u32)], rt_shift: f32) -> Vec<Vec<(u32, u32)>> {
    let n = pairs.len();
    let mut available = vec![true; n];
    let mut by_area: Vec<u32> = (0..n as u32).collect();
    by_area.sort_unstable_by(|&a, &b| {
        annotations[pairs[b as usize].0 as usize]
            .auc
            .partial_cmp(&annotations[pairs[a as usize].0 as usize].auc)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });

    let mut groups: Vec<(f32, Vec<(u32, u32)>)> = Vec::new();
    let mut medians: Vec<f32> = Vec::new();

    for &seed in &by_area {
        let pos = seed as usize;
        if !available[pos] {
            continue;
        }
        let (frag_i, parent_i) = pairs[pos];
        let frag_mz = annotations[frag_i as usize].premz;
        let parent_mz = annotations[parent_i as usize].premz;
        let anchor_rt = annotations[parent_i as usize].rt;

        let lo = frag_mz - 0.01;
        let hi = frag_mz + 0.01;
        let from = pairs[..pos].partition_point(|p| annotations[p.0 as usize].premz < lo);
        let to = pos
            + 1
            + pairs[pos + 1..].partition_point(|p| annotations[p.0 as usize].premz < hi);

        // Only pairs whose parent is the *same* compound count towards a group.
        let same_parent = |k: usize| {
            (annotations[pairs[k].1 as usize].premz - parent_mz).abs() < 0.01
        };

        loop {
            let mut best: Option<(usize, f32)> = None;
            for start in from..to {
                if !available[start] || !same_parent(start) {
                    continue;
                }
                let start_rt = annotations[pairs[start].1 as usize].rt;
                if !(start_rt <= anchor_rt && anchor_rt < start_rt + rt_shift) {
                    continue;
                }
                let rt_end = start_rt + rt_shift;
                let mut total = 0.0f32;
                for k in from..to {
                    if !available[k] || !same_parent(k) {
                        continue;
                    }
                    let rt = annotations[pairs[k].1 as usize].rt;
                    if rt >= start_rt && rt < rt_end {
                        total += annotations[pairs[k].1 as usize].auc;
                    }
                }
                // `>=` so ties resolve to the later window, matching the
                // `max_by` this replaced.
                if best.map_or(true, |(_, b)| total >= b) {
                    best = Some((start, total));
                }
            }
            let Some((start, _)) = best else { break };

            let start_rt = annotations[pairs[start].1 as usize].rt;
            let rt_end = start_rt + rt_shift;
            let mut members = Vec::new();
            medians.clear();
            for k in from..to {
                if !available[k] || !same_parent(k) {
                    continue;
                }
                let rt = annotations[pairs[k].1 as usize].rt;
                if rt >= start_rt && rt < rt_end {
                    available[k] = false;
                    members.push(pairs[k]);
                    medians.push(annotations[pairs[k].1 as usize].premz);
                }
            }
            if members.is_empty() {
                break;
            }
            groups.push((median(&mut medians), members));
        }
    }

    groups.sort_unstable_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.first().cmp(&b.1.first()))
    });
    groups.into_iter().map(|x| x.1).collect()
}

/// Median, computed in place with `select_nth_unstable` (no full sort).
pub fn median(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mid = values.len() / 2;
    values.select_nth_unstable_by(mid, |a, b| {
        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
    });
    if values.len() % 2 == 1 {
        values[mid]
    } else {
        let lower = values[..mid]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        (values[mid] + lower) * 0.5
    }
}

/// `[min, median, max]`, computed in one partial sort.
pub fn min_median_max(values: &mut [f32]) -> [f32; 3] {
    if values.is_empty() {
        return [0.0; 3];
    }
    let mid = values.len() / 2;
    values.select_nth_unstable_by(mid, |a, b| {
        a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
    });
    let lower = &values[..mid.max(1)];
    let min = lower.iter().copied().fold(f32::INFINITY, f32::min);
    let max = values[mid..].iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let med = if values.len() % 2 == 1 {
        values[mid]
    } else {
        let lo = values[..mid]
            .iter()
            .copied()
            .fold(f32::NEG_INFINITY, f32::max);
        (values[mid] + lo) * 0.5
    };
    [min, med, max]
}
