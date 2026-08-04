//! Report writing: the two feature tables and the three spectral libraries.
//!
//! Output files (all under the run's output directory) keep 0.1's names and
//! column layout so existing downstream scripts — including `DDAplot.r` — keep
//! working:
//!
//! * `Report_RTseparated.csv` — one row per consensus feature.
//! * `Report_by_ID.csv` — one row per (feature, compound) pair.
//! * `spec_exp.txt` / `spec_lib.txt` / `spec_reduced.txt` — matched spectra.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::align::{median, min_median_max, Aligned};
use crate::error::{IoContext, Result};
use crate::library::{skeleton, Annotation, Library, Metadata};
use crate::params::Params;
use crate::progress::Reporter;
use crate::score::{Ann, LibHit, SpecSet};

const META_COLUMNS: [&str; 29] = [
    "group",
    "ISF",
    "name",
    "ISF of",
    "isotope",
    "super class (HMDB)",
    "class (HMDB)",
    "sub class (HMDB)",
    "accession (HMDB)",
    "name (HMDB)",
    "core (LM)",
    "main_class (LM)",
    "sub_class (LM)",
    "lm_id (LM)",
    "name (LM)",
    "abbrev (LM)",
    "abbrev_chains (LM)",
    "InChIKey",
    "formula",
    "adduct",
    "feature_m/z",
    "matching_peaks (median)",
    "#lib_peaks (median)",
    "peak_shape (median)",
    "confidence level",
    "Min. RT",
    "Median RT",
    "Max. RT",
    "%detected",
];

/// Blank meta columns emitted for an ISF row, which carries no compound
/// identity of its own.
const ISF_BLANK_COLUMNS: usize = 17;

const PER_SAMPLE_ATTRS: [&str; 5] = ["AREA", "RT", "SCORE", "MP", "S/N"];

/// Number of experimental fragments exported per unmatched spectrum.
const UNKNOWN_PEAKS: usize = 20;

#[derive(Clone, Debug, Default, serde::Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ReportSummary {
    pub features: usize,
    pub identified: usize,
    pub isf_rows: usize,
    pub compounds: usize,
}

struct Ctx<'a> {
    lib: &'a Library,
    meta: &'a Metadata,
    params: &'a Params,
    ann: &'a [Ann],
    n_files: usize,
    /// Resolved at run time — `Polarity::Auto` is decided by sniffing the
    /// first mzML file, so it is not readable off `Params`.
    positive: bool,
}

/// Write both report tables and the spectral outputs.
#[allow(clippy::too_many_arguments)]
pub fn write_reports(
    aligned: &Aligned,
    lib: &Library,
    meta: &Metadata,
    params: &Params,
    positive: bool,
    file_stems: &[String],
    timestamps: &[String],
    spectra: &[SpecSet],
    reporter: &dyn Reporter,
) -> Result<ReportSummary> {
    let dir = &params.output_dir;
    let ctx = Ctx {
        lib,
        meta,
        params,
        ann: &aligned.annotations,
        n_files: file_stems.len(),
        positive,
    };

    let mut by_rt = new_writer(&dir.join("Report_RTseparated.csv"))?;
    let mut by_id = new_writer(&dir.join("Report_by_ID.csv"))?;
    write_headers(&mut by_rt, file_stems, timestamps)?;
    write_headers(&mut by_id, file_stems, timestamps)?;

    let mut exp = new_text(&dir.join("spec_exp.txt"))?;
    let mut lib_out = new_text(&dir.join("spec_lib.txt"))?;
    let mut reduced = new_text(&dir.join("spec_reduced.txt"))?;

    let mut summary = ReportSummary::default();
    let mut compound_names: Vec<&str> = Vec::new();

    for (group_no, members) in (1usize..).zip(&aligned.groups) {
        if members.is_empty() {
            continue;
        }
        summary.features += 1;

        write_spectra(
            &ctx,
            [&mut exp, &mut lib_out, &mut reduced],
            members,
            file_stems,
            spectra,
        )?;

        // Pairs in which a member of this group is the *fragment*, so the row
        // can name what it fragmented from.
        let isf_of = pairs_where_fragment(aligned, members);

        write_row(&ctx, &mut by_rt, members, group_no, &isf_of, &[], None)?;

        let mut names: Vec<&str> = members
            .iter()
            .flat_map(|&m| ctx.ann[m as usize].hits.iter())
            .map(|h| lib.name(h.entry as usize))
            .collect();
        names.sort_unstable();
        names.dedup();
        if !names.is_empty() {
            summary.identified += 1;
        }
        for name in &names {
            compound_names.push(name);
            write_row(
                &ctx,
                &mut by_id,
                members,
                group_no,
                &isf_of,
                &[],
                Some(name),
            )?;
        }

        // Rows for fragments of this group, when it is itself a parent.
        if !members
            .iter()
            .any(|&m| is_parent(aligned, &ctx.ann[m as usize], m))
        {
            continue;
        }
        let anchor = &ctx.ann[members[0] as usize];
        let lo = anchor.premz - 0.01;
        let hi = anchor.premz + 0.01;
        let from = aligned
            .isf_groups
            .partition_point(|g| ctx.ann[g[0].1 as usize].premz < lo);

        let mut families: Vec<&Vec<(u32, u32)>> = Vec::new();
        for g in &aligned.isf_groups[from..] {
            let parent = &ctx.ann[g[0].1 as usize];
            if parent.premz >= hi {
                break;
            }
            if (parent.rt - anchor.rt).abs() < params.rt_shift {
                families.push(g);
            }
        }
        // Heaviest fragment first, as in 0.1.
        families.sort_unstable_by(|a, b| {
            ctx.ann[b[0].0 as usize]
                .premz
                .partial_cmp(&ctx.ann[a[0].0 as usize].premz)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        for family in families {
            let fragments: Vec<u32> = family.iter().map(|x| x.0).collect();
            let parents: Vec<u32> = family.iter().map(|x| x.1).collect();
            write_row(&ctx, &mut by_rt, &fragments, group_no, &[], &parents, None)?;
            write_row(&ctx, &mut by_id, &fragments, group_no, &[], &parents, None)?;
            summary.isf_rows += 1;
        }
    }

    compound_names.sort_unstable();
    compound_names.dedup();
    summary.compounds = compound_names.len();

    by_rt.flush()?;
    by_id.flush()?;
    exp.flush().at(dir.join("spec_exp.txt"))?;
    lib_out.flush().at(dir.join("spec_lib.txt"))?;
    reduced.flush().at(dir.join("spec_reduced.txt"))?;

    reporter.info(format!(
        "{} consensus features, {} identified, {} in-source-fragment rows",
        summary.features, summary.identified, summary.isf_rows
    ));
    Ok(summary)
}

fn new_writer(path: &Path) -> Result<csv::Writer<File>> {
    csv::WriterBuilder::new()
        .from_path(path)
        .map_err(crate::error::Error::Csv)
}

fn new_text(path: &Path) -> Result<BufWriter<File>> {
    Ok(BufWriter::with_capacity(
        1 << 16,
        File::create(path).at(path)?,
    ))
}

fn write_headers(
    wtr: &mut csv::Writer<File>,
    file_stems: &[String],
    timestamps: &[String],
) -> Result<()> {
    // First header row: acquisition timestamps under each per-sample block.
    let mut row: Vec<&str> = vec![""; META_COLUMNS.len()];
    for _ in 0..PER_SAMPLE_ATTRS.len() {
        row.extend(timestamps.iter().map(String::as_str));
    }
    wtr.write_record(&row)?;

    // Second header row: the actual column names.
    let mut names: Vec<String> = META_COLUMNS.iter().map(|x| (*x).to_string()).collect();
    for attr in PER_SAMPLE_ATTRS {
        names.extend(file_stems.iter().map(|bn| format!("{attr}_{bn}")));
    }
    wtr.write_record(&names)?;
    Ok(())
}

/// Pairs in which one of `members` is the in-source fragment.
fn pairs_where_fragment(aligned: &Aligned, members: &[u32]) -> Vec<(u32, u32)> {
    let mut out = Vec::new();
    for &m in members {
        let mz = aligned.annotations[m as usize].premz;
        let from = aligned
            .isf_pairs
            .partition_point(|p| aligned.annotations[p.0 as usize].premz < mz);
        for &pair in &aligned.isf_pairs[from..] {
            if aligned.annotations[pair.0 as usize].premz != mz {
                break;
            }
            if pair.0 == m {
                out.push(pair);
            }
        }
    }
    out
}

fn is_parent(aligned: &Aligned, ann: &Ann, index: u32) -> bool {
    let from = aligned
        .isf_by_parent
        .partition_point(|p| aligned.annotations[p.1 as usize].premz < ann.premz);
    for &pair in &aligned.isf_by_parent[from..] {
        if aligned.annotations[pair.1 as usize].premz != ann.premz {
            break;
        }
        if pair.1 == index {
            return true;
        }
    }
    false
}

/// Hits of one annotation, optionally restricted to a single compound name
/// (keeping that compound's best-scoring hit).
///
/// 0.1 built this by cloning every `Ann` in the group once per compound name.
/// A projection avoids the copy entirely.
fn hits_for(ctx: &Ctx<'_>, ann: &Ann, name_filter: Option<&str>) -> Vec<LibHit> {
    match name_filter {
        None => ann.hits.clone(),
        Some(name) => ann
            .hits
            .iter()
            .filter(|h| ctx.lib.name(h.entry as usize) == name)
            .max_by(|a, b| {
                a.score
                    .partial_cmp(&b.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .into_iter()
            .copied()
            .collect(),
    }
}

fn join_unique(mut values: Vec<&str>, sep: &str) -> String {
    values.retain(|x| !x.trim().is_empty());
    values.sort_unstable();
    values.dedup();
    values.join(sep)
}

#[allow(clippy::too_many_arguments)]
fn write_row(
    ctx: &Ctx<'_>,
    wtr: &mut csv::Writer<File>,
    members: &[u32],
    group_no: usize,
    isf_of: &[(u32, u32)],
    parents: &[u32],
    name_filter: Option<&str>,
) -> Result<()> {
    let is_isf_row = !parents.is_empty();

    // Hits per member, after the optional compound projection.
    let member_hits: Vec<Vec<LibHit>> = members
        .iter()
        .map(|&m| hits_for(ctx, &ctx.ann[m as usize], name_filter))
        .collect();

    // A hit is "confident" if it clears both thresholds. If nothing in the
    // group is confident, fall back to reporting everything rather than an
    // empty row.
    let any_confident = member_hits
        .iter()
        .flatten()
        .any(|h| h.score > ctx.params.ms2_score && h.matched_peaks >= ctx.params.min_peaks);
    let confident = |h: &LibHit| {
        !any_confident
            || (h.score > ctx.params.ms2_score && h.matched_peaks >= ctx.params.min_peaks)
    };

    let inchikeys: Vec<&str> = {
        let mut v: Vec<&str> = member_hits
            .iter()
            .flatten()
            .filter(|h| confident(h))
            .map(|h| ctx.lib.inchikey(h.entry as usize))
            .filter(|x| !x.is_empty())
            .collect();
        v.sort_unstable();
        v.dedup();
        v
    };

    // Exact InChIKey lookup, falling back to the stereochemistry-independent
    // skeleton when the full key is unknown to the metadata tables.
    let mut annotations: Vec<&Annotation> = inchikeys
        .iter()
        .flat_map(|k| ctx.meta.exact(k).iter())
        .collect();
    if annotations.is_empty() {
        annotations = inchikeys
            .iter()
            .flat_map(|k| ctx.meta.by_skeleton(skeleton(k)).iter())
            .collect();
    }

    wtr.write_field(group_no.to_string())?;

    if is_isf_row {
        wtr.write_field("*")?;
        let parent_names: Vec<&str> = parents
            .iter()
            .flat_map(|&p| ctx.ann[p as usize].hits.iter())
            .map(|h| ctx.lib.name(h.entry as usize))
            .collect();
        wtr.write_field(format!("ISF of ({})", join_unique(parent_names, " --- ")))?;
        for _ in 0..ISF_BLANK_COLUMNS {
            wtr.write_field("")?;
        }
    } else {
        wtr.write_field("")?;
        let names: Vec<&str> = member_hits
            .iter()
            .flatten()
            .filter(|h| confident(h))
            .map(|h| ctx.lib.name(h.entry as usize))
            .collect();
        wtr.write_field(join_unique(names, " --- "))?;

        let parent_names: Vec<&str> = isf_of
            .iter()
            .flat_map(|&(_, p)| ctx.ann[p as usize].hits.iter())
            .map(|h| ctx.lib.name(h.entry as usize))
            .collect();
        wtr.write_field(join_unique(parent_names, " --- "))?;

        let isotope_tags: Vec<String> = members
            .iter()
            .filter_map(|&m| {
                let a = &ctx.ann[m as usize];
                a.mono
                    .as_ref()
                    .map(|mono| format!("M+{}", (a.premz - mono.mz).round()))
            })
            .collect();
        let mut tags: Vec<&str> = isotope_tags.iter().map(String::as_str).collect();
        tags.sort_unstable();
        tags.dedup();
        wtr.write_field(tags.join(", "))?;

        macro_rules! ontology {
            ($field:ident, $sep:expr) => {{
                let v: Vec<&str> = annotations.iter().map(|a| a.$field.as_str()).collect();
                wtr.write_field(join_unique(v, $sep))?;
            }};
        }
        ontology!(hmdb_super_class, " --- ");
        ontology!(hmdb_class, " --- ");
        ontology!(hmdb_sub_class, " --- ");
        ontology!(hmdb_accession, " --- ");
        ontology!(hmdb_name, " --- ");
        ontology!(lm_core, " --- ");
        ontology!(lm_main_class, " --- ");
        ontology!(lm_sub_class, " --- ");
        ontology!(lm_id, " --- ");
        ontology!(lm_name, " --- ");
        ontology!(lm_abbrev, " --- ");
        ontology!(lm_abbrev_chains, " --- ");

        wtr.write_field(inchikeys.join(", "))?;

        // Formula: prefer the metadata table, fall back to the library entry's
        // own formula (or a name lookup for the built-in libraries, which
        // carry none).
        if annotations.is_empty() {
            let formulas: Vec<&str> = member_hits
                .iter()
                .flatten()
                .filter(|h| confident(h))
                .map(|h| {
                    let name = ctx.lib.name(h.entry as usize);
                    ctx.meta
                        .formula_for(name)
                        .unwrap_or_else(|| ctx.lib.formula(h.entry as usize))
                })
                .collect();
            wtr.write_field(join_unique(formulas, ", "))?;
        } else {
            ontology!(formula, ", ");
        }

        let adducts: Vec<&str> = member_hits
            .iter()
            .flatten()
            .filter(|h| confident(h))
            .map(|h| ctx.lib.adduct(h.entry as usize))
            .collect();
        wtr.write_field(join_unique(adducts, ", "))?;
    }

    // --- per-sample bucketing -------------------------------------------
    let mut by_file: Vec<Vec<usize>> = vec![Vec::new(); ctx.n_files];
    for (slot, &m) in members.iter().enumerate() {
        let file = ctx.ann[m as usize].file as usize;
        if file < ctx.n_files {
            by_file[file].push(slot);
        }
    }
    for bucket in &mut by_file {
        bucket.sort_unstable_by(|&a, &b| {
            ctx.ann[members[a] as usize]
                .rt
                .partial_cmp(&ctx.ann[members[b] as usize].rt)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // For an ISF row, the per-sample score and matched-peak counts come from
    // the fragment's record of its parent, not from a library hit.
    let mut parent_rel: Vec<Option<crate::score::Relation>> = vec![None; ctx.n_files];
    for &p in parents {
        let parent = &ctx.ann[p as usize];
        let file = parent.file as usize;
        if file >= ctx.n_files {
            continue;
        }
        parent_rel[file] = by_file[file]
            .iter()
            .flat_map(|&slot| ctx.ann[members[slot] as usize].parents.iter())
            .find(|rel| rel.mz == parent.premz && rel.rt == parent.rt)
            .copied();
    }

    let mut scratch: Vec<f32> = members.iter().map(|&m| ctx.ann[m as usize].premz).collect();
    wtr.write_field(format!("{:.4}", median(&mut scratch)))?;

    // matching peaks (median)
    scratch.clear();
    if is_isf_row {
        scratch.extend(
            parent_rel
                .iter()
                .filter_map(|x| x.as_ref())
                .map(|r| f32::from(r.matched_peaks))
                .filter(|x| *x > 0.0),
        );
    } else {
        scratch.extend(
            member_hits
                .iter()
                .flatten()
                .map(|h| f32::from(h.matched_peaks))
                .filter(|x| *x > 0.0),
        );
    }
    wtr.write_field(format!("{:.1}", median(&mut scratch)))?;

    // library peak count (median)
    if is_isf_row {
        wtr.write_field("")?;
    } else {
        scratch.clear();
        scratch.extend(
            member_hits
                .iter()
                .flatten()
                .map(|h| ctx.lib.fragment_count(h.entry as usize) as f32),
        );
        if scratch.is_empty() {
            wtr.write_field("")?;
        } else {
            wtr.write_field(format!("{:.1}", median(&mut scratch)))?;
        }
    }

    scratch.clear();
    scratch.extend(members.iter().map(|&m| ctx.ann[m as usize].shape));
    wtr.write_field(format!("{:.2}", median(&mut scratch)))?;

    wtr.write_field(if any_confident { "MSMS" } else { "MS1" })?;

    scratch.clear();
    scratch.extend(members.iter().map(|&m| ctx.ann[m as usize].rt));
    for v in min_median_max(&mut scratch) {
        wtr.write_field(format!("{v:.2}"))?;
    }

    let detected = by_file.iter().filter(|x| !x.is_empty()).count();
    wtr.write_field(format!(
        "{:.2}",
        detected as f32 / ctx.n_files.max(1) as f32
    ))?;

    // --- per-sample value blocks -----------------------------------------
    let mut cell = String::new();

    // AREA — a leading `*` marks a value imputed from an orphan spectrum.
    for bucket in &by_file {
        cell.clear();
        for (n, &slot) in bucket.iter().enumerate() {
            let a = &ctx.ann[members[slot] as usize];
            if n > 0 {
                cell.push_str(", ");
            }
            cell.push_str(&format!(
                "{}{:.1}",
                if a.is_feature { "" } else { "*" },
                a.auc
            ));
        }
        wtr.write_field(&cell)?;
    }

    // RT
    for bucket in &by_file {
        cell.clear();
        for (n, &slot) in bucket.iter().enumerate() {
            if n > 0 {
                cell.push_str(", ");
            }
            cell.push_str(&format!("{:.2}", ctx.ann[members[slot] as usize].rt));
        }
        wtr.write_field(&cell)?;
    }

    // SCORE
    for (file, bucket) in by_file.iter().enumerate() {
        cell.clear();
        for (n, &slot) in bucket.iter().enumerate() {
            if n > 0 {
                cell.push_str(", ");
            }
            let value = if is_isf_row {
                parent_rel[file].map(|r| r.score).unwrap_or(0.0)
            } else {
                member_hits[slot]
                    .iter()
                    .map(|h| h.score)
                    .fold(0.0f32, f32::max)
            };
            cell.push_str(&format!("{value:.2}"));
        }
        wtr.write_field(&cell)?;
    }

    // MP (matched peaks of the best-scoring hit)
    for (file, bucket) in by_file.iter().enumerate() {
        cell.clear();
        for (n, &slot) in bucket.iter().enumerate() {
            if n > 0 {
                cell.push_str(", ");
            }
            let value = if is_isf_row {
                parent_rel[file].map(|r| r.matched_peaks).unwrap_or(0)
            } else {
                member_hits[slot]
                    .iter()
                    .max_by(|a, b| {
                        a.score
                            .partial_cmp(&b.score)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map_or(0, |h| h.matched_peaks)
            };
            cell.push_str(&value.to_string());
        }
        wtr.write_field(&cell)?;
    }

    // S/N
    for bucket in &by_file {
        cell.clear();
        for (n, &slot) in bucket.iter().enumerate() {
            if n > 0 {
                cell.push_str(", ");
            }
            cell.push_str(&format!("{:.1}", ctx.ann[members[slot] as usize].s_n));
        }
        wtr.write_field(&cell)?;
    }

    wtr.write_record(None::<&[u8]>)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Spectral output
// ---------------------------------------------------------------------------

fn write_spectra(
    ctx: &Ctx<'_>,
    [exp, lib_out, reduced]: [&mut BufWriter<File>; 3],
    members: &[u32],
    file_stems: &[String],
    spectra: &[SpecSet],
) -> Result<()> {
    let tol = ctx.params.ms2_tol;

    // One representative spectrum per compound, best score wins.
    let mut best_per_compound: Vec<(u32, LibHit)> = Vec::new();

    for &m in members {
        let ann = &ctx.ann[m as usize];
        let file = ann.file as usize;
        let Some(specs) = spectra.get(file) else {
            continue;
        };
        for hit in &ann.hits {
            let Some(si) = specs.find(hit.spec_prec_mz, hit.spec_rt) else {
                continue;
            };
            write_header_block(exp, ctx, ann, hit, specs, si, file_stems)?;
            write_experimental(exp, ctx, hit, specs, si, tol)?;
            write_header_block(lib_out, ctx, ann, hit, specs, si, file_stems)?;
            write_library(lib_out, ctx, hit, specs, si, tol)?;

            let key = ctx.lib.inchikey(hit.entry as usize);
            if key.trim().is_empty() {
                continue;
            }
            match best_per_compound
                .iter()
                .position(|(_, h)| ctx.lib.inchikey(h.entry as usize) == key)
            {
                Some(pos) if best_per_compound[pos].1.score < hit.score => {
                    best_per_compound[pos] = (m, *hit);
                }
                Some(_) => {}
                None => best_per_compound.push((m, *hit)),
            }
        }
    }

    for (m, hit) in best_per_compound {
        let ann = &ctx.ann[m as usize];
        let file = ann.file as usize;
        let Some(specs) = spectra.get(file) else {
            continue;
        };
        let Some(si) = specs.find(hit.spec_prec_mz, hit.spec_rt) else {
            continue;
        };
        write_header_block(reduced, ctx, ann, &hit, specs, si, file_stems)?;
        write_experimental(reduced, ctx, &hit, specs, si, tol)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_header_block(
    w: &mut BufWriter<File>,
    ctx: &Ctx<'_>,
    ann: &Ann,
    hit: &LibHit,
    specs: &SpecSet,
    si: usize,
    file_stems: &[String],
) -> Result<()> {
    let e = hit.entry as usize;
    writeln!(w, "NAME: {}", ctx.lib.name(e))?;
    writeln!(w, "INCHIKEY: {}", ctx.lib.inchikey(e))?;
    writeln!(
        w,
        "PRECURSORTYPE: [{}]{}",
        ctx.lib.adduct(e),
        if ctx.positive { '+' } else { '-' }
    )?;
    writeln!(w, "PRECURSORMZ: {:.4}", ctx.lib.mass(e))?;
    writeln!(
        w,
        "SAMPLE: {}",
        file_stems
            .get(ann.file as usize)
            .map(String::as_str)
            .unwrap_or("")
    )?;
    let ce = specs.ce(si);
    if ce > 0.0 {
        writeln!(w, "COLLISIONENERGY: {ce}")?;
    }
    writeln!(w, "SCORE: {:.2}", hit.score)?;
    writeln!(w, "SHAPE: {:.2}", ann.shape)?;
    writeln!(w, "MASS_DIFF(LIB-EXP): {:.4}", ctx.lib.mass(e) - ann.premz)?;
    writeln!(w, "RETENTIONTIME: {:.2}", ann.rt)?;
    Ok(())
}

/// Experimental peaks, starring those the library entry explains.
fn write_experimental(
    w: &mut BufWriter<File>,
    ctx: &Ctx<'_>,
    hit: &LibHit,
    specs: &SpecSet,
    si: usize,
    tol: f32,
) -> Result<()> {
    let (exp_mz, exp_i) = specs.scored(si);
    let (lib_mz, _) = ctx.lib.fragments(hit.entry as usize);
    let mut sorted: Vec<f32> = lib_mz.to_vec();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    writeln!(w, "Num Peaks: {}", exp_mz.len())?;
    for (&mz, &i) in exp_mz.iter().zip(exp_i) {
        let p = sorted.partition_point(|&x| x < mz - tol);
        let star = sorted.get(p).is_some_and(|&x| x < mz + tol);
        writeln!(w, "{mz:.4} {i:.2}{}", if star { " *" } else { "" })?;
    }
    writeln!(w)?;
    Ok(())
}

/// Library peaks, starring those the experiment confirms.
fn write_library(
    w: &mut BufWriter<File>,
    ctx: &Ctx<'_>,
    hit: &LibHit,
    specs: &SpecSet,
    si: usize,
    tol: f32,
) -> Result<()> {
    let (lib_mz, lib_i) = ctx.lib.fragments(hit.entry as usize);
    let (exp_mz, _) = specs.scored(si);

    writeln!(w, "Num Peaks: {}", lib_mz.len())?;
    for (&mz, &i) in lib_mz.iter().zip(lib_i) {
        let p = exp_mz.partition_point(|&x| x < mz - tol);
        let star = exp_mz.get(p).is_some_and(|&x| x < mz + tol);
        writeln!(w, "{mz:.4} {i:.2}{}", if star { " *" } else { "" })?;
    }
    writeln!(w)?;
    Ok(())
}

/// Export every filtered spectrum of one sample as unknowns, so they can be
/// re-searched against another library.
pub fn write_unknowns(
    path: &Path,
    file_index: usize,
    specs: &SpecSet,
    positive: bool,
) -> Result<()> {
    let mut w = new_text(path)?;
    for i in 0..specs.len() {
        writeln!(w, "NAME: {file_index}_{i}")?;
        writeln!(w, "PRECURSORMZ: {:.4}", specs.prec_mz(i))?;
        writeln!(
            w,
            "PRECURSORTYPE: [unknown]{}",
            if positive { "+" } else { "-" }
        )?;
        writeln!(w, "RETENTIONTIME: {:.2}", specs.rt(i))?;
        let ce = specs.ce(i);
        if ce > 0.0 {
            writeln!(w, "COLLISIONENERGY: {ce:.1}")?;
        }
        let (mzs, ints) = specs.full(i);
        let n = mzs.len().min(UNKNOWN_PEAKS);
        writeln!(w, "Num Peaks: {n}")?;
        for k in 0..n {
            writeln!(w, "{:.4} {:.2}", mzs[k], ints[k])?;
        }
        writeln!(w)?;
    }
    w.flush().at(path)?;
    Ok(())
}
