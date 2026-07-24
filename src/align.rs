use crate::Ann;
use crate::Spec;
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};

pub fn align(
    param_t: &crate::Param,
    all_dat: Vec<Vec<Ann>>,
    all_spec: &[Vec<Spec>],
    time_stamp: &[String],
) -> Result<(), Box<dyn Error>> {
    let crate::Param { mzml_fs, .. } = param_t;
    let mut all_dat: Vec<Ann> = all_dat.into_iter().flatten().collect();
    all_dat.sort_unstable_by(|x, y| x.premz.partial_cmp(&y.premz).unwrap());
    let ref_pt_rt_sep = get_ref_pt(&all_dat, param_t);
    let isf_isf = assoc_isf(&all_dat, param_t);
    let mut isf_flat: Vec<(&Ann, &Ann)> = isf_isf.clone();
    isf_flat.sort_unstable_by(|x, y| x.1.premz.partial_cmp(&y.1.premz).unwrap());
    let isf_pt = get_isf_pt(&isf_isf, param_t.rt_shift);
    let name_formu = read_name_f();
    let iddat_d = read_inchik();
    let iddat_i14: Vec<&str> = iddat_d
        .iter()
        .map(|x| &x.inchikey[..14.min(x.inchikey.len())])
        .collect();

    let mut bufw_e = BufWriter::new(File::create("spec_exp.txt")?);
    let mut bufw_r = BufWriter::new(File::create("spec_reduced.txt")?);
    let mut bufw_l = BufWriter::new(File::create("spec_lib.txt")?);
    let mut wtr_id = csv::WriterBuilder::new().from_path("Report_by_ID.csv")?;
    let mut wtr_rt = csv::WriterBuilder::new().from_path("Report_RTseparated.csv")?;
    let meta_head = [
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
    let attr = ["AREA", "RT", "SCORE", "MP", "S/N"];
    let bns: Vec<&str> = mzml_fs
        .iter()
        .map(|x| x.file_stem()?.to_str())
        .collect::<Option<_>>()
        .unwrap();

    let iter = std::iter::repeat_n("", meta_head.len())
        .chain((0..attr.len()).flat_map(|_| time_stamp.iter().map(std::string::String::as_str)));
    wtr_rt.write_record(iter.clone())?;
    wtr_id.write_record(iter)?;

    let line: Vec<String> = attr
        .into_iter()
        .flat_map(|pref| bns.iter().map(move |bn| [pref, "_", bn].concat()))
        .collect();
    let iter = meta_head
        .iter()
        .copied()
        .chain(line.iter().map(std::string::String::as_str));
    wtr_rt.write_record(iter.clone())?;
    wtr_id.write_record(iter)?;

    let print_row = |wtr_: &mut csv::Writer<File>,
                     mem: &[&Ann],
                     n0: usize,
                     p_isf: &[(&Ann, &Ann)],
                     par_: &[&Ann]|
     -> io::Result<()> {
        let mut msms_ind: Vec<bool> = mem
            .iter()
            .flat_map(|x| &x.ent_l)
            .map(|y| y.dotp > param_t.ms2_score && y.m_peaks >= param_t.min_peaks)
            .collect();
        let msms_v = msms_ind.iter().any(|x| *x);
        if !msms_v {
            msms_ind.fill(true);
        }
        let mut inchik_s: Vec<&str> = mem
            .iter()
            .flat_map(|x| &x.ent_l)
            .zip(&msms_ind)
            .filter(|(_, x)| **x)
            .map(|(x, _)| x.ent.inchik.as_str())
            .filter(|x| !x.is_empty())
            .collect();
        inchik_s.sort_unstable();
        inchik_s.dedup();
        let mut iddats: Vec<&Iddat> = inchik_s
            .iter()
            .flat_map(|x| {
                iddat_d[iddat_d.partition_point(|y| y.inchikey.as_str() < *x)..]
                    .iter()
                    .take_while(|y| y.inchikey.as_str() == *x)
            })
            .collect();
        if iddats.is_empty() {
            iddats = inchik_s
                .iter()
                .flat_map(|x| {
                    let x = &x[..14.min(x.len())];
                    let pos0 = iddat_i14.partition_point(|y| *y < x);
                    iddat_d[pos0..]
                        .iter()
                        .zip(iddat_i14[pos0..].iter().take_while(move |i14| **i14 == x))
                        .map(|y| y.0)
                })
                .collect();
        }

        wtr_.write_field(n0.to_string())?;
        if par_.is_empty() {
            wtr_.write_field("")?;
            wtr_.write_field({
                let mut nameset = mem
                    .iter()
                    .flat_map(|x| &x.ent_l)
                    .zip(&msms_ind)
                    .filter(|(_, x)| **x)
                    .map(|(x, _)| x.ent.name.as_str())
                    .collect::<Vec<_>>();
                nameset.sort_unstable();
                nameset.dedup();
                nameset.join(" --- ")
            })?;
            wtr_.write_field({
                let mut nameset = p_isf
                    .iter()
                    .flat_map(|x| &x.1.ent_l)
                    .map(|x| x.ent.name.as_str())
                    .collect::<Vec<_>>();
                nameset.sort_unstable();
                nameset.dedup();
                nameset.join(" --- ")
            })?;
            wtr_.write_field({
                let mut iso_tag: Vec<String> = mem
                    .iter()
                    .filter_map(|x| {
                        x.mono
                            .as_ref()
                            .map(|y| format!("M+{}", (x.premz - y.pmz).round()))
                    })
                    .collect::<Vec<_>>();
                iso_tag.sort_unstable();
                iso_tag.dedup();
                iso_tag.join(", ")
            })?;
            macro_rules! dedup1 {
                ($sn:ident, $sep:expr) => {{
                    let mut xx: Vec<&str> = iddats
                        .iter()
                        .map(|x| x.$sn.as_str())
                        .filter(|x| !x.is_empty())
                        .collect();
                    xx.sort_unstable();
                    xx.dedup();
                    xx.join($sep)
                }};
            }
            wtr_.write_field(dedup1!(super_hmdb, " --- "))?;
            wtr_.write_field(dedup1!(class_hmdb, " --- "))?;
            wtr_.write_field(dedup1!(sub_hmdb, " --- "))?;
            wtr_.write_field(dedup1!(acc_hmdb, " --- "))?;
            wtr_.write_field(dedup1!(name_hmdb, " --- "))?;
            wtr_.write_field(dedup1!(core_lm, " --- "))?;
            wtr_.write_field(dedup1!(main_lm, " --- "))?;
            wtr_.write_field(dedup1!(sub_lm, " --- "))?;
            wtr_.write_field(dedup1!(lm_id, " --- "))?;
            wtr_.write_field(dedup1!(name_lm, " --- "))?;
            wtr_.write_field(dedup1!(abb_lm, " --- "))?;
            wtr_.write_field(dedup1!(abb_c_lm, " --- "))?;
            wtr_.write_field(inchik_s.join(", "))?;
            wtr_.write_field({
                if iddats.is_empty() {
                    let mut f_ = mem
                        .iter()
                        .flat_map(|x| &x.ent_l)
                        .zip(&msms_ind)
                        .filter(|(_, x)| **x)
                        .map(|(x, _)| {
                            name_formu
                                .binary_search_by_key(&&x.ent.name, |y| &y.0)
                                .map_or(x.ent.formu.as_str(), |i| name_formu[i].1.as_str())
                        })
                        .filter(|x| !x.trim().is_empty())
                        .collect::<Vec<_>>();
                    f_.sort_unstable();
                    f_.dedup();
                    f_.join(", ")
                } else {
                    dedup1!(formu, ", ")
                }
            })?;
            wtr_.write_field({
                let mut add_set = mem
                    .iter()
                    .flat_map(|x| &x.ent_l)
                    .zip(&msms_ind)
                    .filter(|(_, x)| **x)
                    .map(|(x, _)| x.ent.adduct.as_str())
                    .collect::<Vec<_>>();
                add_set.sort_unstable();
                add_set.dedup();
                add_set.join(", ")
            })?;
        } else {
            wtr_.write_field("*")?;
            wtr_.write_field(format!("ISF of ({})", {
                let mut nameset = par_
                    .iter()
                    .flat_map(|x| &x.ent_l)
                    .map(|x| x.ent.name.as_str())
                    .collect::<Vec<_>>();
                nameset.sort_unstable();
                nameset.dedup();
                nameset.join(" --- ")
            }))?;
            (0..17).for_each(|_| wtr_.write_field("").unwrap());
        }
        let mut dat: Vec<Vec<&Ann>> = vec![Vec::new(); mzml_fs.len()];
        for ann_ in mem {
            dat[ann_.nn].push(ann_);
        }
        for dat_nn in &mut dat {
            dat_nn.sort_unstable_by(|x, y| x.rt.partial_cmp(&y.rt).unwrap());
        }
        let mut par_d: Vec<Option<&crate::Par>> = vec![None; mzml_fs.len()];
        for ip in par_ {
            par_d[ip.nn] = Some(
                dat[ip.nn]
                    .iter()
                    .flat_map(|x| &x.par_l)
                    .find(|x| x.pmz == ip.premz && x.prt == ip.rt)
                    .unwrap(),
            );
        }
        wtr_.write_field(format!(
            "{:.4}",
            median(mem.iter().map(|x| x.premz).collect())
        ))?;
        wtr_.write_field({
            let mpl: Vec<_> = if par_.is_empty() {
                mem.iter()
                    .flat_map(|x| &x.ent_l)
                    .map(|x| x.m_peaks)
                    .filter(|x| *x > 0)
                    .map(f32::from)
                    .collect()
            } else {
                par_d
                    .iter()
                    .filter_map(|x| x.map(|y| y.m_peaks))
                    .filter(|x| *x > 0)
                    .map(f32::from)
                    .collect()
            };
            format!("{:.1}", if mpl.is_empty() { 0. } else { median(mpl) })
        })?;
        let npeaks: Vec<_> = mem
            .iter()
            .flat_map(|x| &x.ent_l)
            .map(|x| x.ent.mz_i.len() as f32)
            .collect();
        if par_.is_empty() && !npeaks.is_empty() {
            wtr_.write_field(format!("{:.1}", median(npeaks)))?;
        } else {
            wtr_.write_field("")?;
        }
        wtr_.write_field(format!(
            "{:.2}",
            median(mem.iter().map(|x| x.shape).collect())
        ))?;
        wtr_.write_field(if msms_v { "MSMS" } else { "MS1" })?;
        for x in &quantile(mem.iter().map(|x| x.rt).collect()) {
            wtr_.write_field(format!("{x:.2}"))?;
        }
        wtr_.write_field(format!(
            "{:.2}",
            dat.iter().filter(|x| !x.is_empty()).count() as f32 / mzml_fs.len() as f32
        ))?;
        macro_rules! q_r_s_c {
            ($f:expr) => {
                dat.iter()
                    .map(|x| x.iter().map($f).collect::<Vec<_>>().join(", "))
                    .for_each(|x| wtr_.write_field(x).unwrap())
            };
        }
        q_r_s_c!(|x| format!("{}{:.1}", if x.feat { "" } else { "*" }, x.auc));
        q_r_s_c!(|x| format!("{:.2}", x.rt));
        dat.iter()
            .enumerate()
            .map(|(nn, dat_nn)| {
                dat_nn
                    .iter()
                    .map(|x| {
                        format!(
                            "{:.2}",
                            if par_.is_empty() {
                                x.ent_l
                                    .iter()
                                    .map(|x| x.dotp)
                                    .reduce(f32::max)
                                    .unwrap_or(0.)
                            } else {
                                par_d[nn].unwrap().dotp
                            }
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(", ")
            })
            .for_each(|x| wtr_.write_field(x).unwrap());
        dat.iter()
            .enumerate()
            .map(|(nn, dat_nn)| {
                dat_nn
                    .iter()
                    .map(|x| {
                        if par_.is_empty() {
                            x.ent_l
                                .iter()
                                .max_by(|x, y| x.dotp.partial_cmp(&y.dotp).unwrap())
                                .map_or(0, |x| x.m_peaks)
                        } else {
                            par_d[nn].unwrap().m_peaks
                        }
                        .to_string()
                    })
                    .collect::<Vec<_>>()
                    .join(", ")
            })
            .for_each(|x| wtr_.write_field(x).unwrap());
        q_r_s_c!(|x| format!("{:.1}", x.s_n));
        wtr_.write_record(None::<&[u8]>)?;
        Ok(())
    };
    for (n0, mem) in (1..).zip(&ref_pt_rt_sep) {
        write_spec_ent(
            [&mut bufw_e, &mut bufw_r, &mut bufw_l],
            param_t.ms2tol,
            mem,
            &bns,
            param_t.ispos,
            all_spec,
        )?;
        let p_isf = mem
            .iter()
            .flat_map(|x| {
                let pos0 = isf_isf.partition_point(|y| y.0.premz < x.premz);
                isf_isf[pos0..]
                    .iter()
                    .take_while(|y| y.0.premz == x.premz)
                    .filter(|y| y.0.rt == x.rt && y.0.nn == x.nn)
            })
            .copied()
            .collect::<Vec<(&Ann, &Ann)>>();
        print_row(&mut wtr_rt, mem, n0, &p_isf, &[])?;
        let mut ent_set = mem
            .iter()
            .flat_map(|x| &x.ent_l)
            .map(|x| &x.ent.name)
            .collect::<Vec<_>>();
        ent_set.sort_unstable();
        ent_set.dedup();
        for u_ent in ent_set {
            let mem: Vec<Ann> = mem
                .iter()
                .filter_map(|x| {
                    x.ent_l
                        .iter()
                        .filter(|x| &x.ent.name == u_ent)
                        .max_by(|x, y| x.dotp.partial_cmp(&y.dotp).unwrap())
                        .map(|y| Ann {
                            nn: x.nn,
                            premz: x.premz,
                            rt: x.rt,
                            auc: x.auc,
                            feat: x.feat,
                            s_n: x.s_n,
                            shape: x.shape,
                            ent_l: vec![y.clone()],
                            par_l: x.par_l.clone(),
                            mono: x.mono.clone(),
                        })
                })
                .collect();
            let mem: Vec<&Ann> = mem.iter().collect();
            let p_isf = mem
                .iter()
                .flat_map(|x| {
                    let pos0 = isf_isf.partition_point(|y| y.0.premz < x.premz);
                    isf_isf[pos0..]
                        .iter()
                        .take_while(|y| y.0.premz == x.premz)
                        .filter(|y| y.0.rt == x.rt && y.0.nn == x.nn)
                })
                .copied()
                .collect::<Vec<(&Ann, &Ann)>>();
            print_row(&mut wtr_id, &mem, n0, &p_isf, &[])?;
        }
        if mem.iter().any(|x| {
            let pos0 = isf_flat.partition_point(|y| y.1.premz < x.premz);
            isf_flat[pos0..]
                .iter()
                .take_while(|y| y.1.premz == x.premz)
                .any(|y| y.1.rt == x.rt && y.1.nn == x.nn)
        }) {
            let pos0 = mem[0].premz - 0.01;
            let pos0 = isf_pt.partition_point(|x| x[0].1.premz < pos0);
            let mut isf_pt_: Vec<&Vec<(&Ann, &Ann)>> = isf_pt[pos0..]
                .iter()
                .take_while(|y| y[0].1.premz < mem[0].premz + 0.01)
                .filter(|y| (y[0].1.rt - mem[0].rt).abs() < param_t.rt_shift)
                .collect();
            isf_pt_.sort_unstable_by(|y, x| x[0].0.premz.partial_cmp(&y[0].0.premz).unwrap());
            for mem_isf in isf_pt_ {
                let isf_: Vec<&Ann> = mem_isf.iter().map(|x| x.0).collect();
                let par_: Vec<&Ann> = mem_isf.iter().map(|x| x.1).collect();
                print_row(&mut wtr_rt, &isf_, n0, &[], &par_)?;
                print_row(&mut wtr_id, &isf_, n0, &[], &par_)?;
            }
        }
    }
    Ok(())
}
fn assoc_isf<'a>(all_dat: &'a [Ann], param_t: &crate::Param) -> Vec<(&'a Ann<'a>, &'a Ann<'a>)> {
    all_dat
        .iter()
        .flat_map(|x| {
            x.par_l.iter().filter_map(move |y| {
                let pos = all_dat.partition_point(|z| z.premz < y.pmz);
                all_dat[pos..]
                    .iter()
                    .take_while(|z| z.premz == y.pmz)
                    .find(|z| {
                        z.rt == y.prt
                            && z.nn == x.nn
                            && z.ent_l.iter().any(|w| {
                                !w.ent.adduct.contains("2M")
                                    && param_t.ms2_score < w.dotp
                                    && param_t.min_peaks <= w.m_peaks
                            })
                    })
                    .map(|y| (x, y))
            })
        })
        .collect()
}
fn get_isf_pt<'a>(
    all_isf: &'a [(&Ann, &Ann)],
    rt_shift: f32,
) -> Vec<Vec<(&'a Ann<'a>, &'a Ann<'a>)>> {
    let mut mem: Vec<(f32, Vec<(&Ann, &Ann)>)> = Vec::new();
    let mut inc = vec![true; all_isf.len()];
    let mut all_isf_int: Vec<(usize, (&Ann, &Ann))> = all_isf.iter().copied().enumerate().collect();
    all_isf_int.sort_unstable_by(|y, x| x.1.0.auc.partial_cmp(&y.1.0.auc).unwrap());
    let mut mem_inc1 = Vec::<((&Ann, &Ann), bool)>::new();
    for (pos, (ann_, ann_p)) in all_isf_int {
        if inc[pos] {
            let lo = ann_.premz - 0.01;
            let up = ann_.premz + 0.01;
            let p0 = all_isf[..pos].partition_point(|x| x.0.premz < lo);
            let p1 = pos + 1 + all_isf[pos + 1..].partition_point(|x| x.0.premz < up);
            mem_inc1.clear();
            mem_inc1.extend(
                all_isf[p0..p1]
                    .iter()
                    .map(|x| (*x, (x.1.premz - ann_p.premz).abs() < 0.01)),
            );
            while let Some((max_sum, _)) = mem_inc1
                .iter()
                .zip(&inc[p0..])
                .filter(|(((_, a_p), i1), i0)| {
                    **i0 && *i1 && a_p.rt <= ann_p.rt && ann_p.rt < a_p.rt + rt_shift
                })
                .map(|((aa_p, _), _)| {
                    let rt_up = aa_p.1.rt + rt_shift;
                    (
                        aa_p,
                        mem_inc1
                            .iter()
                            .zip(&inc[p0..])
                            .filter(|(((_, x), i1), i0)| {
                                **i0 && *i1 && aa_p.1.rt <= x.rt && x.rt < rt_up
                            })
                            .map(|(((_, x), _), _)| x.auc)
                            .sum::<f32>(),
                    )
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            {
                mem.push({
                    let rt_up = max_sum.1.rt + rt_shift;
                    let aa: Vec<_> = mem_inc1
                        .iter()
                        .zip(&mut inc[p0..])
                        .filter(|(((_, x), i1), i0)| {
                            **i0 && *i1 && max_sum.1.rt <= x.rt && x.rt < rt_up
                        })
                        .map(|((x, _), i0)| {
                            *i0 = false;
                            *x
                        })
                        .collect();
                    (median(aa.iter().map(|x| x.1.premz).collect()), aa)
                });
            }
        }
    }
    mem.sort_unstable_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    mem.into_iter().map(|x| x.1).collect()
}
fn get_ref_pt<'a>(all_dat: &'a [Ann], param_t: &crate::Param) -> Vec<Vec<&'a Ann<'a>>> {
    let &crate::Param {
        rt_shift, mz_shift, ..
    } = param_t;
    let mut mem: Vec<(f32, Vec<&Ann>)> = Vec::new();
    let mut inc: Vec<bool> = all_dat.iter().map(|x| !x.ent_l.is_empty()).collect();
    let mut all_dat_int: Vec<(usize, &Ann)> = all_dat.iter().enumerate().collect();
    all_dat_int.sort_unstable_by(|y, x| x.1.auc.partial_cmp(&y.1.auc).unwrap());
    for (pos, dat) in all_dat_int {
        if inc[pos] {
            let lo = dat.premz - mz_shift;
            let up = dat.premz + mz_shift;
            let p0 = all_dat[..pos].partition_point(|x| x.premz < lo);
            let p1 = pos + 1 + all_dat[pos + 1..].partition_point(|x| x.premz < up);
            let mem_: &[Ann] = &all_dat[p0..p1];
            while let Some((max_sum, _)) = mem_
                .iter()
                .zip(&inc[p0..])
                .filter(|(ann_, i0)| **i0 && ann_.rt <= dat.rt && dat.rt < ann_.rt + rt_shift)
                .map(|(ann_, _)| {
                    let rt_up = ann_.rt + rt_shift;
                    (
                        ann_,
                        mem_.iter()
                            .zip(&inc[p0..])
                            .filter(|(x, i0)| **i0 && ann_.rt <= x.rt && x.rt < rt_up)
                            .map(|x| x.0.auc)
                            .sum::<f32>(),
                    )
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            {
                mem.push({
                    let rt_up = max_sum.rt + rt_shift;
                    let aa: Vec<_> = mem_
                        .iter()
                        .zip(&mut inc[p0..])
                        .filter(|(x, i0)| **i0 && max_sum.rt <= x.rt && x.rt < rt_up)
                        .map(|(x, i0)| {
                            *i0 = false;
                            x
                        })
                        .collect();
                    (median(aa.iter().map(|x| x.premz).collect()), aa)
                });
            }
        }
    }
    mem.sort_unstable_by(|x, y| x.0.partial_cmp(&y.0).unwrap());
    mem.into_iter().map(|x| x.1).collect()
}
struct Iddat {
    inchikey: String,
    acc_hmdb: String,
    name_hmdb: String,
    super_hmdb: String,
    class_hmdb: String,
    sub_hmdb: String,
    lm_id: String,
    name_lm: String,
    abb_lm: String,
    core_lm: String,
    main_lm: String,
    sub_lm: String,
    abb_c_lm: String,
    formu: String,
}
fn read_inchik() -> Vec<Iddat> {
    let file_path = std::env::current_exe()
        .unwrap()
        .parent()
        .unwrap()
        .join("libs")
        .join("inchik.txt");
    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .trim(csv::Trim::All)
        .from_path(file_path)
        .unwrap();
    let mut iddat_d: Vec<Iddat> = rdr
        .into_records()
        .map(std::result::Result::unwrap)
        .map(|x| Iddat {
            inchikey: x[0].to_string(),
            acc_hmdb: x[1].to_string(),
            name_hmdb: x[2].to_string(),
            super_hmdb: x[3].to_string(),
            class_hmdb: x[4].to_string(),
            sub_hmdb: x[5].to_string(),
            lm_id: x[6].to_string(),
            name_lm: x[7].to_string(),
            abb_lm: x[8].to_string(),
            core_lm: x[9].to_string(),
            main_lm: x[10].to_string(),
            sub_lm: x[11].to_string(),
            abb_c_lm: x[12].to_string(),
            formu: x[13].to_string(),
        })
        .collect();
    iddat_d.sort_unstable_by(|a, b| a.inchikey.cmp(&b.inchikey));
    iddat_d
}

fn read_name_f() -> Vec<(String, String)> {
    let file_path = std::env::current_exe()
        .unwrap()
        .parent()
        .unwrap()
        .join("libs")
        .join("name_formu.txt");
    let rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path(file_path)
        .unwrap();
    let mut name_f: Vec<_> = rdr
        .into_records()
        .map(std::result::Result::unwrap)
        .map(|x| (x[0].to_string(), x[1].to_string()))
        .collect();
    name_f.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    name_f
}
fn median(mut num_l: Vec<f32>) -> f32 {
    let pos = num_l.len() / 2;
    num_l.select_nth_unstable_by(pos, |a, b| a.partial_cmp(b).unwrap());
    if num_l.len() % 2 == 1 {
        num_l[pos]
    } else {
        num_l[pos].midpoint(num_l[..pos].iter().copied().reduce(f32::max).unwrap())
    }
}
fn quantile(mut num_l: Vec<f32>) -> [f32; 3] {
    let pos = num_l.len() / 2;
    num_l.select_nth_unstable_by(pos, |a, b| a.partial_cmp(b).unwrap());
    [
        num_l[..pos.max(1)]
            .iter()
            .copied()
            .reduce(f32::min)
            .unwrap(),
        if num_l.len() % 2 == 1 {
            num_l[pos]
        } else {
            num_l[pos].midpoint(num_l[..pos].iter().copied().reduce(f32::max).unwrap())
        },
        num_l[pos..].iter().copied().reduce(f32::max).unwrap(),
    ]
}
fn write_spec_ent(
    [bufw_e, bufw_r, bufw_l]: [&mut BufWriter<File>; 3],
    ms2tol: f32,
    mem: &[&Ann],
    bns: &[&str],
    ispos: bool,
    all_spec: &[Vec<Spec>],
) -> io::Result<()> {
    let mut top: Vec<(&crate::Lib, &Ann)> = mem
        .iter()
        .flat_map(|ann_| {
            ann_.ent_l
                .iter()
                .filter(|d| !d.ent.inchik.trim().is_empty())
                .map(|d| (d, *ann_))
        })
        .collect();
    top.sort_unstable_by_key(|x| &x.0.ent.inchik);
    let mut top_one = Vec::<(&crate::Lib, &Ann)>::new();
    for d_ann in top {
        if top_one
            .last()
            .is_none_or(|x| x.0.ent.inchik != d_ann.0.ent.inchik)
        {
            top_one.push(d_ann);
        } else if top_one.last().unwrap().0.dotp < d_ann.0.dotp {
            *top_one.last_mut().unwrap() = d_ann;
        }
    }
    let write_ent =
        |bufw: &mut BufWriter<File>, d: &crate::Lib, ann_: &Ann, spec: &Spec| -> io::Result<()> {
            writeln!(bufw, "NAME: {}", d.ent.name)?;
            writeln!(bufw, "INCHIKEY: {}", d.ent.inchik)?;
            writeln!(
                bufw,
                "PRECURSORTYPE: [{}]{}",
                d.ent.adduct,
                if ispos { '+' } else { '-' }
            )?;
            writeln!(bufw, "PRECURSORMZ: {:.4}", d.ent.mmass)?;
            writeln!(bufw, "SAMPLE: {}", bns[ann_.nn])?;
            if spec.ce > 0.0 {
                writeln!(bufw, "COLLISIONENERGY: {}", spec.ce)?;
            }
            writeln!(bufw, "SCORE: {:.2}", d.dotp)?;
            writeln!(bufw, "SHAPE: {:.2}", ann_.shape)?;
            writeln!(bufw, "MASS_DIFF(LIB-EXP): {:.4}", d.ent.mmass - ann_.premz)?;
            writeln!(bufw, "RETENTIONTIME: {:.2}", ann_.rt)?;
            Ok(())
        };
    let write_exp = |bufw: &mut BufWriter<File>, d: &crate::Lib, spec: &Spec| -> io::Result<()> {
        writeln!(bufw, "Num Peaks: {}", spec.mz_i_l.len())?;
        let mut lib_mz_l: Vec<_> = d.ent.mz_i.iter().map(|x| x.0).collect();
        lib_mz_l.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        for (mz, i) in &spec.mz_i_l {
            let pos0 = mz - ms2tol;
            let pos0 = lib_mz_l.partition_point(|&x| x < pos0);
            writeln!(
                bufw,
                "{mz:.4} {i:.2}{}",
                if lib_mz_l[pos0..].first().is_some_and(|&x| x < mz + ms2tol) {
                    " *"
                } else {
                    ""
                }
            )?;
        }
        writeln!(bufw)?;
        Ok(())
    };
    let write_lib = |bufw: &mut BufWriter<File>, d: &crate::Lib, spec: &Spec| -> io::Result<()> {
        writeln!(bufw, "Num Peaks: {}", d.ent.mz_i.len())?;
        let mut lib_mz_l: Vec<_> = spec.mz_i_l.iter().map(|x| x.0).collect();
        lib_mz_l.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        for (mz, i) in &d.ent.mz_i {
            let pos0 = mz - ms2tol;
            let pos0 = lib_mz_l.partition_point(|&x| x < pos0);
            writeln!(
                bufw,
                "{mz:.4} {i:.2}{}",
                if lib_mz_l[pos0..].first().is_some_and(|&x| x < mz + ms2tol) {
                    " *"
                } else {
                    ""
                }
            )?;
        }
        writeln!(bufw)?;
        Ok(())
    };
    for ann_ in mem {
        for d in &ann_.ent_l {
            let pos0 = all_spec[ann_.nn].partition_point(|y| y.ms1mz < d.ms1mz);
            let spec = all_spec[ann_.nn][pos0..]
                .iter()
                .find(|y| y.ms1mz == d.ms1mz && y.rt == d.rt)
                .unwrap();
            write_ent(bufw_e, d, ann_, spec)?;
            write_exp(bufw_e, d, spec)?;
            write_ent(bufw_l, d, ann_, spec)?;
            write_lib(bufw_l, d, spec)?;
        }
    }
    for (d, ann_) in top_one {
        let pos0 = all_spec[ann_.nn].partition_point(|y| y.ms1mz < d.ms1mz);
        let spec = all_spec[ann_.nn][pos0..]
            .iter()
            .find(|y| y.ms1mz == d.ms1mz && y.rt == d.rt)
            .unwrap();
        write_ent(bufw_r, d, ann_, spec)?;
        write_exp(bufw_r, d, spec)?;
    }
    Ok(())
}
