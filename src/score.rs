use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufReader, BufWriter, Write};
fn cos_sim(l1l2: &[(f32, f32)]) -> f32 {
    let l1sum = l1l2.iter().map(|x| x.0).sum::<f32>();
    let l2sum = l1l2.iter().map(|x| x.1).sum::<f32>();
    if l1sum <= 0. || l2sum <= 0. {
        return 0.;
    }
    let numer = l1l2.iter().map(|x| (x.0 * x.1).sqrt()).sum::<f32>();
    let denom = (l1sum * l2sum).sqrt();
    numer / denom
}

struct Peak {
    mz: f32,
    rt: f32,
    sc: f32,
    coef: f32,
    shape: f32,
    smooth: f32,
}

fn ms1feat(bn: &str) -> Result<Vec<Peak>, Box<dyn Error>> {
    use crate::common::{unpack_f32, unpack_u32};
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("ms1f_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path)?);
    let len0 = unpack_u32(buf)?;
    (0..len0)
        .map(|_| {
            Ok(Peak {
                mz: unpack_f32(buf)?,
                rt: unpack_f32(buf)?,
                sc: unpack_f32(buf)?,
                coef: unpack_f32(buf)?,
                shape: unpack_f32(buf)?,
                smooth: unpack_f32(buf)?,
            })
        })
        .collect()
}

use super::Spec;
use crate::Param;

fn filter_ms2(ms2_scans: Vec<crate::Msms>, param_t: &Param) -> Vec<Spec> {
    ms2_scans
        .into_iter()
        .filter_map(|x| {
            let mut mz_i_l = x
                .mz_i_l
                .into_iter()
                .filter(|y| y.0 < x.ms1mz - 0.3)
                .collect::<Vec<_>>();
            mz_i_l.sort_unstable_by(|b, a| a.1.partial_cmp(&b.1).unwrap());
            mz_i_l.truncate(u8::MAX.into());
            let mz_i_all = mz_i_l.clone();
            if let Some(icut) = param_t.icut {
                if let Some(i) = mz_i_l.iter().position(|&x| x.1 < icut) {
                    mz_i_l.truncate(i);
                }
            } else {
                mz_i_l.truncate(param_t.match_n_fragments);
            }
            mz_i_l.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            (!mz_i_l.is_empty()).then_some(Spec {
                ms1mz: x.ms1mz,
                rt: x.rt,
                mz_i_l,
                mz_i_all,
                ce: x.ce,
            })
        })
        .collect()
}
fn ms2_feat_assign<'a>(
    ms2_sc: &'a [Spec],
    ms1peaks: &'a [Peak],
    param_t: &Param,
) -> (Vec<(Vec<&'a Spec>, &'a Peak)>, Vec<&'a Spec>) {
    let mut attached = Vec::<usize>::new();
    let mut cand = Vec::<(usize, &Spec, f32)>::new();
    let ms2_wp: Vec<(&Peak, Vec<(usize, &Spec)>)> = ms1peaks
        .iter()
        .filter_map(|ms1| {
            let pos0 = ms1.mz - param_t.ms1ms2;
            let pos0 = ms2_sc.partition_point(|x| x.ms1mz < pos0);
            let up = ms1.mz + param_t.ms1ms2;
            cand.clear();
            cand.extend(
                (pos0..)
                    .zip(ms2_sc[pos0..].iter().take_while(|x| x.ms1mz < up))
                    .map(|x| (x.0, x.1, (ms1.rt - x.1.rt).abs())),
            );
            attached.extend(cand.iter().filter(|x| x.2 < ms1.sc * 1.5).map(|x| x.0));
            let spec_i: Vec<(usize, &Spec)> = if cand.iter().any(|x| x.2 < ms1.sc) {
                cand.iter()
                    .filter(|x| x.2 < ms1.sc)
                    .map(|x| (x.0, x.1))
                    .collect()
            } else {
                cand.iter()
                    .filter(|x| x.2 < ms1.sc * 1.5)
                    .map(|x| (x.0, x.1))
                    .collect()
            };
            (!spec_i.is_empty()).then_some((ms1, spec_i))
        })
        .collect();
    let ms2_wp_ = ms2_wp
        .iter()
        .map(|(ms1, spec_i)| (spec_i.iter().map(|x| x.1).collect(), *ms1))
        .collect();
    let ms2_no_p = if param_t.features_only {
        Vec::new()
    } else {
        attached.sort_unstable();
        ms2_sc
            .iter()
            .enumerate()
            .filter(|(spec_i, _)| attached.binary_search(spec_i).is_err())
            .map(|x| x.1)
            .collect()
    };
    (ms2_wp_, ms2_no_p)
}
fn isf_match(
    ms2_wp: &[(Vec<&Spec>, &Peak)],
    &Param {
        ms2tol,
        isf_rt_diff,
        isf_p_diff,
        ..
    }: &Param,
) -> Vec<Vec<(usize, u8)>> {
    let mut isf_sc = vec![Vec::new(); ms2_wp.len()];
    for (ii, (ms2_vec, p)) in ms2_wp.iter().enumerate().rev() {
        for ms2 in ms2_vec {
            let bisect = |a| ms2.mz_i_l.partition_point(|x| x.0 < a);
            let pos0 = ms2.ms1mz - isf_p_diff;
            let pos0 = ms2_wp.partition_point(|x| x.1.mz < pos0);
            for ((i_ms2_vec, i_p), isf_sc_jj) in ms2_wp[..pos0]
                .iter()
                .zip(isf_sc.iter_mut())
                .filter(|((_, i_p), _)| (p.rt - i_p.rt).abs() < isf_rt_diff)
            {
                let pos0 = bisect(i_p.mz - ms2tol);
                let i_p_mz_tol = i_p.mz + ms2tol;
                if 0.1 * ms2.mz_i_l.iter().map(|x| x.1).reduce(f32::max).unwrap()
                    < ms2.mz_i_l[pos0..]
                        .iter()
                        .take_while(|x| x.0 < i_p_mz_tol)
                        .map(|x| x.1)
                        .reduce(f32::max)
                        .unwrap_or(0.)
                {
                    for i_ms2 in i_ms2_vec {
                        let m_peaks = u8::try_from(
                            i_ms2
                                .mz_i_all
                                .iter()
                                .filter(|x| x.0 < i_ms2.ms1mz + ms2tol)
                                .take(10)
                                .filter(|(f_mz, _)| {
                                    let pos0 = bisect(f_mz - ms2tol);
                                    ms2.mz_i_l[pos0..]
                                        .first()
                                        .is_some_and(|x| x.0 < f_mz + ms2tol)
                                })
                                .count(),
                        )
                        .unwrap();
                        if m_peaks > 1 {
                            if isf_sc_jj.last().is_none_or(|&(x, _)| x != ii) {
                                isf_sc_jj.push((ii, m_peaks));
                            } else if isf_sc_jj.last().is_some_and(|x| x.0 == ii && x.1 < m_peaks) {
                                isf_sc_jj.last_mut().unwrap().1 = m_peaks;
                            }
                        }
                    }
                }
            }
        }
    }
    isf_sc
}
struct ScMsms<'a, 'b> {
    cs: f32,
    ent: &'a super::readlib::Ent,
    m_peaks: u8,
    ma_p: Vec<(f32, f32, f32, f32)>,
    ms2: &'b Spec,
}
fn cpd_match<'a, 'b>(
    score_ent: &mut Vec<ScMsms<'a, 'b>>,
    (pos0, pos1): (usize, f32),
    ms2: &'b Spec,
    lib_ent: &'a [super::readlib::Ent],
    &Param {
        ms2tol,
        ms2_score,
        min_peaks,
        rt_shift,
        chimeric,
        top_only,
        ..
    }: &Param,
) {
    let mut xfrag = Vec::<usize>::new();
    let mut ent_ms2_i = Vec::<(f32, f32)>::new();
    for ent in lib_ent[pos0..]
        .iter()
        .take_while(|x| x.mmass < pos1)
        .filter(|x| x.rt.is_none_or(|rt| (rt - ms2.rt).abs() < rt_shift))
    {
        xfrag.clear();
        ent_ms2_i.clear();
        let mut ma_p = Vec::with_capacity(ent.mz_i.len());
        for (f_mz, f_i) in ent
            .mz_i
            .iter()
            .filter(|x| 0.1 < f32::from(ent.charge).mul_add(ms2.ms1mz, -x.0))
        {
            let pos0 = f_mz - ms2tol;
            let pos0 = ms2.mz_i_l.partition_point(|x| x.0 < pos0);
            let mz_up = f_mz + ms2tol;
            let mz_i = (pos0..)
                .zip(ms2.mz_i_l[pos0..].iter().take_while(|x| x.0 < mz_up))
                .map(|(x, y)| {
                    xfrag.push(x);
                    y
                })
                .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
                .unwrap_or(&(0., 0.));
            ent_ms2_i.push((*f_i, mz_i.1));
            if mz_i.1 > 0. {
                ma_p.push((*f_mz, *f_i, mz_i.0, mz_i.1));
            }
        }
        let m_peaks = u8::try_from(ent_ms2_i.iter().filter(|x| x.1 > 0.).count()).unwrap();
        if m_peaks >= min_peaks {
            if !chimeric {
                xfrag.sort_unstable();
                ent_ms2_i.extend(
                    ms2.mz_i_l
                        .iter()
                        .enumerate()
                        .filter(|(x, _)| xfrag.binary_search(x).is_err())
                        .take_while(|x| 0.1 < f32::from(ent.charge).mul_add(ms2.ms1mz, -x.1.0))
                        .map(|(_, (_, f_i))| (0., *f_i)),
                );
            }
            let cs = cos_sim(&ent_ms2_i);
            if ms2_score < cs {
                if score_ent.is_empty() || !top_only {
                    score_ent.push(ScMsms {
                        cs,
                        ent,
                        m_peaks,
                        ma_p,
                        ms2,
                    });
                } else if top_only && score_ent.first().is_some_and(|x| x.cs < cs) {
                    score_ent[0] = ScMsms {
                        cs,
                        ent,
                        m_peaks,
                        ma_p,
                        ms2,
                    };
                }
            }
        }
    }
}

fn print_plot(
    buf: &mut BufWriter<File>,
    feat: Option<&Peak>,
    (chrom, rt_l, rt_r): (&[(f32, f32)], usize, usize),
    ScMsms {
        cs,
        ent,
        m_peaks,
        ma_p,
        ms2,
    }: &ScMsms,
) -> Result<(), Box<dyn Error>> {
    buf.write_all(ent.name.as_bytes())?;
    buf.write_all(b"\0")?;
    let (premz, rt, shape, smooth) =
        feat.map_or((0., 0., 0., 0.), |x| (x.mz, x.rt, x.shape, x.smooth));
    buf.write_all(&premz.to_le_bytes())?;
    buf.write_all(&shape.to_le_bytes())?;
    buf.write_all(&smooth.to_le_bytes())?;
    buf.write_all(&rt.to_le_bytes())?;
    buf.write_all(&ms2.ms1mz.to_le_bytes())?;
    buf.write_all(&ms2.rt.to_le_bytes())?;
    buf.write_all(&u16::try_from(rt_l + 1)?.to_le_bytes())?;
    buf.write_all(&u16::try_from(rt_r)?.to_le_bytes())?;
    buf.write_all(&u16::try_from(chrom.len())?.to_le_bytes())?;
    for (x, y) in chrom {
        buf.write_all(&x.to_le_bytes())?;
        buf.write_all(&y.to_le_bytes())?;
    }
    buf.write_all(&ent.mmass.to_le_bytes())?;
    buf.write_all(&cs.to_le_bytes())?;
    buf.write_all(&u8::try_from(ms2.mz_i_all.len())?.to_le_bytes())?;
    for (x, y) in &ms2.mz_i_all {
        buf.write_all(&x.to_le_bytes())?;
        buf.write_all(&y.to_le_bytes())?;
    }
    buf.write_all(&ms2.ce.to_le_bytes())?;
    buf.write_all(&u8::try_from(ent.mz_i.len())?.to_le_bytes())?;
    for (x, y) in &ent.mz_i {
        buf.write_all(&x.to_le_bytes())?;
        buf.write_all(&y.to_le_bytes())?;
    }
    buf.write_all(&m_peaks.to_le_bytes())?;
    for x in ma_p {
        buf.write_all(&x.0.to_le_bytes())?;
        buf.write_all(&x.1.to_le_bytes())?;
    }
    for x in ma_p {
        buf.write_all(&x.2.to_le_bytes())?;
        buf.write_all(&x.3.to_le_bytes())?;
    }
    Ok(())
}
fn peak_grouping(ms1peaks: &[Peak], isf_rt_diff: f32) -> Vec<Vec<&Peak>> {
    let find_mplus = |iso: &Peak, pos: usize| -> Option<(usize, &Peak)> {
        let lo = iso.mz + 1.003_355 - 0.003;
        let up = iso.mz + 1.003_355 + 0.003;
        (pos + 1..)
            .zip(&ms1peaks[pos + 1..])
            .skip_while(|(_, x)| x.mz < lo)
            .take_while(|(_, x)| x.mz < up)
            .find(|(_, x)| {
                (x.rt - iso.rt).abs() < isf_rt_diff
                    && x.coef < iso.coef
                    && (x.sc - iso.sc).abs() < 1.01f32.max(iso.sc * 0.11)
            })
    };
    let mut peak_g: Vec<Vec<&Peak>> = Vec::new();
    for (pos, mono) in ms1peaks.iter().enumerate() {
        let Some((pos, iso)) = find_mplus(mono, pos) else {
            continue;
        };
        let mut peak_g1: Vec<&Peak> = vec![mono, iso];
        let Some((pos, iso)) = find_mplus(iso, pos) else {
            peak_g.push(peak_g1);
            continue;
        };
        peak_g1.push(iso);
        let Some((_, iso)) = find_mplus(iso, pos) else {
            peak_g.push(peak_g1);
            continue;
        };
        peak_g1.push(iso);
        peak_g.push(peak_g1);
    }
    peak_g
}
const DIS_WIDTH: f32 = 0.3;
pub fn score<'a>(
    nn: usize,
    param_t: &Param,
    bn: &str,
    lib_ent: &'a [super::readlib::Ent],
    ms1_scans: &[crate::Ms],
    ms2_scans: Vec<crate::Msms>,
) -> Result<(Vec<super::Ann<'a>>, Vec<Spec>), Box<dyn Error>> {
    let ms1peaks = ms1feat(bn)?;
    let ms2_sc = filter_ms2(ms2_scans, param_t);
    print_u(nn, bn, &ms2_sc, param_t.ispos)?;
    let (ms2_wp, ms2_no_p) = ms2_feat_assign(&ms2_sc, &ms1peaks, param_t);
    let isf_sc = isf_match(&ms2_wp, param_t);
    let peak_g = peak_grouping(&ms1peaks, param_t.isf_rt_diff);
    let mut iso_mono: Vec<(&Peak, &Peak)> = peak_g
        .iter()
        .flat_map(|peak_g0| peak_g0[1..].iter().map(|peak| (*peak, peak_g0[0])))
        .collect();
    iso_mono.sort_by(|x, y| x.0.mz.partial_cmp(&y.0.mz).unwrap());
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("plot_{bn}.bin"));
    let mut bufp = BufWriter::new(File::create(file_path)?);
    let mut all_dat = Vec::new();
    let mut score_ent: Vec<ScMsms> = Vec::new();
    for (spec_feat, isf_sc_) in ms2_wp.iter().zip(isf_sc) {
        let (ms2_l, feat): &(Vec<&Spec>, &Peak) = spec_feat;
        let rt_l = feat.sc.mul_add(-param_t.i_rt, feat.rt);
        let rt_r = feat.sc.mul_add(param_t.i_rt, feat.rt);
        let chrom = crate::common::get_chrom(
            (feat.mz, feat.rt, feat.sc.mul_add(param_t.i_rt, DIS_WIDTH)),
            ms1_scans,
            param_t.i_mz,
        );
        let pos0 = chrom.partition_point(|x| x.0 < rt_l);
        let pos1 = pos0 + chrom[pos0..].partition_point(|x| x.0 < rt_r);
        let s_n = {
            let pos0_ = feat.rt - feat.sc;
            let pos0_ = pos0 + chrom[pos0..pos1].partition_point(|x| x.0 < pos0_);
            let pos1_ = feat.rt + feat.sc;
            let pos1_ = pos0_ + chrom[pos0_..pos1].partition_point(|x| x.0 < pos1_);
            chrom[pos0_..pos1_]
                .iter()
                .map(|x| x.1)
                .reduce(f32::max)
                .unwrap()
                / (chrom[pos0..pos0_]
                    .iter()
                    .chain(&chrom[pos1_..pos1])
                    .map(|x| x.1)
                    .sum::<f32>()
                    / (pos0_ - pos0 + pos1 - pos1_).max(1) as f32)
                    .max(1.)
        };
        if s_n < param_t.s_n1 {
            continue;
        }
        let ms1_auc = chrom[pos0..]
            .iter()
            .zip(&chrom[pos0 + 1..pos1])
            .map(|((rt0, i0), (rt1, i1))| (i0 + i1) * (rt1 - rt0))
            .sum::<f32>()
            * 30.;
        let ms1tol = match param_t.ms1tol_u.as_str() {
            "ppm" => param_t.ms1tol * feat.mz / 1e6,
            "m/z" => param_t.ms1tol,
            _ => panic!(),
        };
        let pos = feat.mz - ms1tol;
        score_ent.clear();
        for ms2 in ms2_l {
            cpd_match(
                &mut score_ent,
                (lib_ent.partition_point(|x| x.mmass < pos), feat.mz + ms1tol),
                ms2,
                lib_ent,
                param_t,
            );
        }
        for tm in &score_ent {
            print_plot(&mut bufp, Some(feat), (&chrom, pos0, pos1), tm)?;
        }
        if !score_ent.is_empty() || !isf_sc_.is_empty() {
            let mono: Option<super::Par> = {
                iso_mono[iso_mono.partition_point(|x| x.0.mz < feat.mz)..]
                    .iter()
                    .take_while(|x| x.0.mz == feat.mz)
                    .find(|x| x.0.rt == feat.rt)
                    .map(|m| super::Par {
                        pmz: m.1.mz,
                        prt: m.1.rt,
                        dotp: 0.,
                        m_peaks: 0,
                    })
            };
            all_dat.push({
                super::Ann {
                    nn,
                    premz: feat.mz,
                    rt: feat.rt,
                    auc: ms1_auc,
                    feat: true,
                    s_n,
                    shape: feat.shape,
                    ent_l: score_ent
                        .iter()
                        .map(|x| super::Lib {
                            dotp: x.cs,
                            m_peaks: x.m_peaks,
                            ent: x.ent,
                            ms1mz: x.ms2.ms1mz,
                            rt: x.ms2.rt,
                        })
                        .collect(),
                    par_l: isf_sc_
                        .into_iter()
                        .map(|(ii, m_peaks)| {
                            let feat1: &Peak = ms2_wp[ii].1;
                            super::Par {
                                dotp: 0.,
                                pmz: feat1.mz,
                                prt: feat1.rt,
                                m_peaks,
                            }
                        })
                        .collect(),
                    mono,
                }
            });
        }
    }
    for ms2 in ms2_no_p {
        let rt_l = ms2.rt - param_t.impute_width;
        let rt_r = ms2.rt + param_t.impute_width;
        let chrom = crate::common::get_chrom(
            (ms2.ms1mz, ms2.rt, param_t.impute_width + DIS_WIDTH),
            ms1_scans,
            param_t.i_mz,
        );
        let pos0 = chrom.partition_point(|x| x.0 < rt_l);
        let pos1 = pos0 + chrom[pos0..].partition_point(|x| x.0 < rt_r);
        let ms1_auc = chrom[pos0..]
            .iter()
            .zip(&chrom[pos0 + 1..pos1])
            .map(|((rt0, i0), (rt1, i1))| (i0 + i1) * (rt1 - rt0))
            .sum::<f32>()
            * 30.;
        let ms1tol = match param_t.ms1tol_u.as_str() {
            "ppm" => param_t.ms1tol * ms2.ms1mz / 1e6,
            "m/z" => param_t.ms1tol,
            _ => panic!(),
        };
        let pos = ms2.ms1mz - ms1tol;
        score_ent.clear();
        cpd_match(
            &mut score_ent,
            (
                lib_ent.partition_point(|x| x.mmass < pos),
                ms2.ms1mz + ms1tol,
            ),
            ms2,
            lib_ent,
            param_t,
        );
        for tm in &score_ent {
            print_plot(&mut bufp, None, (&chrom, pos0, pos1), tm)?;
        }
        if !score_ent.is_empty() {
            all_dat.push({
                super::Ann {
                    nn,
                    premz: ms2.ms1mz,
                    rt: ms2.rt,
                    auc: ms1_auc,
                    feat: false,
                    s_n: 0.,
                    shape: 0.,
                    ent_l: score_ent
                        .iter()
                        .map(|x| super::Lib {
                            dotp: x.cs,
                            m_peaks: x.m_peaks,
                            ent: x.ent,
                            ms1mz: x.ms2.ms1mz,
                            rt: x.ms2.rt,
                        })
                        .collect(),
                    par_l: Vec::new(),
                    mono: None,
                }
            });
        }
    }
    Ok((all_dat, ms2_sc))
}

fn print_u(filei: usize, bn: &str, ms2_sc: &[Spec], ispos: bool) -> io::Result<()> {
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("u_{bn}.txt"));
    let mut buf = BufWriter::new(File::create(file_path)?);
    for (c, ms2) in ms2_sc.iter().enumerate() {
        writeln!(buf, "NAME: {filei}_{c}")?;
        writeln!(buf, "PRECURSORMZ: {:.4}", ms2.ms1mz)?;
        writeln!(
            buf,
            "PRECURSORTYPE: [unknown]{}",
            if ispos { "+" } else { "-" }
        )?;
        writeln!(buf, "RETENTIONTIME: {:.2}", ms2.rt)?;
        if ms2.ce > 0. {
            writeln!(buf, "COLLISIONENERGY: {:.1}", ms2.ce)?;
        }
        writeln!(buf, "Num Peaks: {}", ms2.mz_i_all.len().min(super::PRINT_N))?;
        for (mz, i) in ms2.mz_i_all.iter().take(super::PRINT_N) {
            writeln!(buf, "{mz:.4} {i:.2}")?;
        }
        writeln!(buf)?;
    }
    Ok(())
}
