use std::fs::File;
use std::io::{BufWriter, Write};

struct MzRtScC {
    mz: f32,
    rt: f32,
    sc: f32,
    coef: f32,
    shape: f32,
    smooth: f32,
}
fn findridge(
    peak_list: &mut Vec<MzRtScC>,
    (rt_all, ms1_scans): (&[f32], &[crate::Ms]),
    (rt_mz_i_l, mz_range): (&[Eic], (f32, f32)),
    (wave_scs, wave_sqrt): (&[f32], &[f32]),
    msms_rt: &mut [f32],
    p: &crate::Param,
) {
    struct RtScC {
        rt: f32,
        sc: f32,
        coef: f32,
    }
    let bisect = |a| rt_all.partition_point(|&x| x < a);
    let mut eic_p = vec![(0f32, 0f32); rt_all.len()];
    for &(rt, mz_i) in rt_mz_i_l {
        eic_p[rt] = mz_i;
    }
    let noise = {
        let mut i_l: Vec<_> = rt_mz_i_l.iter().map(|x| x.1.1).collect();
        let i = i_l.len() / 20;
        i_l.select_nth_unstable_by(i, |a, b| a.partial_cmp(b).unwrap());
        i_l[i]
    };

    msms_rt.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let mut prev = true;
    let (eic_rt, mut calc_p): (Vec<_>, Vec<_>) = rt_all
        .iter()
        .zip(&rt_all[1..])
        .filter_map(|(rt0, rt1)| {
            let midrt = rt0.midpoint(*rt1);
            let pos0 = midrt - p.rt_search;
            let pos0 = msms_rt.partition_point(|&x| x < pos0);
            let c = msms_rt[pos0..]
                .first()
                .is_some_and(|&x| x < midrt + p.rt_search);
            (c || prev).then(|| {
                prev = c;
                (midrt, c)
            })
        })
        .unzip();
    calc_p[0] = false;
    *calc_p.last_mut().unwrap() = false;

    for x in &mut eic_p {
        x.1 = 0f32.max(x.1 - noise);
    }

    let mut coefs = vec![0f32; eic_rt.len() * wave_scs.len()];
    let mut int_i = Vec::new();
    for ((wave_s, w_sqrt), coef_xx) in wave_scs
        .iter()
        .zip(wave_sqrt)
        .zip(coefs.chunks_exact_mut(eic_rt.len()))
    {
        let wave_s_ = wave_s.recip();
        for ((yy, wave_loc), _) in coef_xx
            .iter_mut()
            .zip(&eic_rt)
            .zip(&calc_p)
            .filter(|x| *x.1)
        {
            let pos0 = bisect(wave_loc - wave_s);
            let pos1 = bisect(wave_loc + wave_s);
            let rt_ = &rt_all[pos0..pos1];
            if rt_.len() < 2 {
                continue;
            }
            int_i.clear();
            int_i.extend(
                rt_.iter()
                    .map(|rt0| ((rt0 - wave_loc) * wave_s_).powi(2))
                    .zip(&eic_p[pos0..])
                    .map(|(tsig2, i0)| i0.1 * (-tsig2 / 2.).exp() * (1. - tsig2)),
            );
            *yy = rt_
                .iter()
                .zip(&int_i)
                .zip(&rt_[1..])
                .zip(&int_i[1..])
                .map(|(((rt0, i0), rt1), i1)| (i0 + i1) * (rt1 - rt0))
                .sum::<f32>()
                / w_sqrt;
        }
    }
    let mut local_max = Vec::new();
    for (wave_s, coef_xx) in wave_scs.iter().zip(coefs.chunks_exact_mut(eic_rt.len())) {
        let mut l_max = Vec::new();
        loop {
            let (max_i, &max_coef) = coef_xx
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .unwrap();
            if max_coef <= 0. {
                break;
            }
            if coef_xx[max_i - 1] <= 0. || coef_xx[max_i + 1] <= 0. {
                coef_xx[max_i] = 0.;
                continue;
            }
            l_max.push((max_i, max_coef));
            let lo = eic_rt[max_i] - wave_s;
            let up = eic_rt[max_i] + wave_s;
            let lo = eic_rt[..max_i].partition_point(|x| *x <= lo);
            let up = max_i + eic_rt[max_i..].partition_point(|x| *x < up);
            coef_xx[lo..up].fill(0.);
        }
        local_max.push(l_max);
    }
    let mut ridgels: Vec<Vec<RtScC>> = local_max[0]
        .iter()
        .map(|&(max_i, coef)| {
            vec![RtScC {
                rt: eic_rt[max_i],
                sc: wave_scs[0],
                coef,
            }]
        })
        .collect();
    for (xx, scale_coef) in local_max[1..].iter().enumerate() {
        for &(max_i, coef) in scale_coef {
            let rt = eic_rt[max_i];
            if let Some(rl) = ridgels.iter_mut().find(|rl| {
                let rl_last = rl.last().unwrap();
                rl_last.sc == wave_scs[xx]
                    && ((rl_last.rt - rt).abs() < 0.01
                        || bisect(rt).abs_diff(bisect(rl_last.rt)) < 2)
            }) {
                rl.push(RtScC {
                    rt,
                    sc: wave_scs[xx + 1],
                    coef,
                });
            } else {
                ridgels.push(vec![RtScC {
                    rt,
                    sc: wave_scs[xx + 1],
                    coef,
                }]);
            }
        }
    }
    let rlen = 3;
    let iter = ridgels
        .into_iter()
        .filter(|x| x.len() >= rlen)
        .filter_map(|rd| {
            let (p_loc_i, p_loc) = rd
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.coef.partial_cmp(&b.1.coef).unwrap())
                .unwrap();
            if p_loc.sc * 2.0 < p.peak_w.0 || p.peak_w.1 < p_loc.sc * 2.0 {
                return None;
            }
            if p_loc_i + rlen > rd.len() {
                return None;
            }
            let pos0 = bisect(p_loc.rt - p_loc.sc);
            let pos1 = bisect(p_loc.rt + p_loc.sc);
            if pos0 == 0 || pos1 == rt_all.len() {
                return None;
            }
            let peakmz: f32 = eic_p[pos0..pos1].iter().map(|x| x.0 * x.1).sum::<f32>()
                / eic_p[pos0..pos1].iter().map(|x| x.1).sum::<f32>();
            if peakmz < mz_range.0 + 0.005 || mz_range.1 - 0.005 < peakmz {
                return None;
            }
            let mut chrom =
                crate::common::get_chrom((peakmz, p_loc.rt, 3. * p_loc.sc), ms1_scans, p.i_mz);
            let leftbd = p_loc.rt - p_loc.sc;
            let rightbd = p_loc.rt + p_loc.sc;
            let pos0 = chrom.partition_point(|x| x.0 < leftbd);
            let pos1 = chrom.partition_point(|x| x.0 < rightbd);
            let smooth = acf1({
                let mut keep = vec![true; chrom.len()];
                for (i, c3) in keep[1..].iter_mut().zip(chrom.windows(3)) {
                    *i = c3.iter().any(|x| x.1 > 0.);
                }
                chrom
                    .iter()
                    .zip(keep)
                    .filter(|x| x.1)
                    .map(|x| x.0.1)
                    .collect()
            })
            .unwrap();
            if smooth < 0. {
                return None;
            }

            chrom.truncate(pos1);
            chrom.drain(..pos0);
            int_i.clear();
            int_i.extend(chrom.iter().map(|(rt0, _)| {
                let tsig2 = ((rt0 - p_loc.rt) / p_loc.sc).powi(2);
                (1. - tsig2) * (-tsig2 / 2.).exp()
            }));
            {
                let ave = int_i.iter().sum::<f32>() / int_i.len() as f32;
                for i in &mut int_i {
                    *i -= ave;
                }
                let ave = chrom.iter().map(|x| x.1).sum::<f32>() / chrom.len() as f32;
                for c in &mut chrom {
                    c.1 -= ave;
                }
            }
            let mut a_dot_b = 0.;
            let mut a_mag = 0.;
            let mut b_mag = 0.;
            for (((c0, c1), &i0), &i1) in chrom.iter().zip(&chrom[1..]).zip(&int_i).zip(&int_i[1..])
            {
                let d = c1.0 - c0.0;
                a_dot_b = c0.1.mul_add(i0, c1.1 * i1).mul_add(d, a_dot_b);
                a_mag = i0.mul_add(i0, i1 * i1).mul_add(d, a_mag);
                b_mag = c0.1.mul_add(c0.1, c1.1 * c1.1).mul_add(d, b_mag);
            }
            let shape = a_dot_b / (a_mag * b_mag).sqrt();
            (shape > p.peak_shape && smooth > p.sn_score).then_some(MzRtScC {
                mz: peakmz,
                rt: p_loc.rt,
                sc: p_loc.sc,
                coef: p_loc.coef,
                shape,
                smooth,
            })
        });
    peak_list.extend(iter);
}

type Eic = (usize, (f32, f32));

const MZ_S: f32 = 0.007;
pub fn cwt(
    bn: &str,
    ms1_scans: &[crate::Ms],
    ms2_scans: &[crate::Msms],
    lib_ent: &[crate::readlib::Ent],
    p: &crate::Param,
) -> std::io::Result<()> {
    let mut wave_scs = vec![0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09];
    wave_scs.extend((0..40).map(|i| 0.1 * 1.1f32.powi(i)));
    wave_scs.truncate(wave_scs.partition_point(|&x| x < p.peak_w.1 / 2.) + 3);
    let wave_sqrt: Vec<f32> = wave_scs.iter().map(|x| x.sqrt()).collect();
    let rt_all: Vec<_> = ms1_scans.iter().map(|x| x.rt).collect();
    let mut mz0 = ms2_scans[0].ms1mz - 0.02;
    let mut peak_list = Vec::<MzRtScC>::new();
    let mut rt_mz_i_l = Vec::new();
    let mut msms_rt = Vec::new();
    while mz0 < ms2_scans.last().unwrap().ms1mz {
        let mz1 = MZ_S.mul_add(3., mz0);
        let lo = mz0 - p.ms1ms2;
        let up = mz1 + p.ms1ms2;
        msms_rt.clear();
        msms_rt.extend(
            ms2_scans[ms2_scans.partition_point(|x| x.ms1mz < lo)..]
                .iter()
                .take_while(|x| x.ms1mz <= up)
                .map(|x| x.rt),
        );
        if !msms_rt.is_empty() && {
            let lo = mz0 - MZ_S;
            let pos0 = lib_ent.partition_point(|x| x.mmass < lo);
            lib_ent[pos0..]
                .first()
                .is_some_and(|x| x.mmass <= mz1 + MZ_S)
        } {
            rt_mz_i_l.clear();
            rt_mz_i_l.extend(ms1_scans.iter().enumerate().filter_map(|(rt_i, ms1_sc)| {
                let pos0 = ms1_sc.mz_i.partition_point(|x| x.0 < mz0);
                ms1_sc.mz_i[pos0..]
                    .iter()
                    .take_while(|x| x.0 <= mz1)
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .map(|x| (rt_i, *x))
            }));
            if rt_mz_i_l.len() > 2 {
                findridge(
                    &mut peak_list,
                    (&rt_all, ms1_scans),
                    (&rt_mz_i_l, (mz0, mz1)),
                    (&wave_scs, &wave_sqrt),
                    &mut msms_rt,
                    p,
                );
            }
        }
        mz0 += MZ_S;
    }
    peak_list.sort_unstable_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());

    let mut inc = vec![true; peak_list.len()];
    let mut pl_int: Vec<(usize, &MzRtScC)> = peak_list.iter().enumerate().collect();
    pl_int.sort_unstable_by(|y, x| x.1.coef.partial_cmp(&y.1.coef).unwrap());
    for (pos, peak) in pl_int {
        if inc[pos] {
            let lo = peak.mz - 0.009;
            let up = peak.mz + 0.009;
            peak_list[..pos]
                .iter()
                .enumerate()
                .rev()
                .take_while(|(_, x)| lo < x.mz)
                .chain((pos + 1..).zip(peak_list[pos + 1..].iter().take_while(|x| x.mz < up)))
                .filter(|(_, x)| (x.rt - peak.rt).abs() < x.sc + peak.sc)
                .for_each(|(j, _)| inc[j] = false);
        }
    }
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("ms1f_{bn}.bin"));
    let mut bufb = BufWriter::new(File::create(file_path)?);
    let len0 = inc.iter().filter(|x| **x).count();
    bufb.write_all(&u32::try_from(len0).unwrap().to_le_bytes())?;
    for (_, x) in inc.into_iter().zip(peak_list).filter(|x| x.0) {
        bufb.write_all(&x.mz.to_le_bytes())?;
        bufb.write_all(&x.rt.to_le_bytes())?;
        bufb.write_all(&x.sc.to_le_bytes())?;
        bufb.write_all(&x.coef.to_le_bytes())?;
        bufb.write_all(&x.shape.to_le_bytes())?;
        bufb.write_all(&x.smooth.to_le_bytes())?;
    }
    Ok(())
}
fn acf1(mut data: Vec<f32>) -> Option<f32> {
    if data.len() < 2 {
        return None; // Cannot compute autocorrelation with less than 2 points
    }
    let mean = data.iter().sum::<f32>() / data.len() as f32;
    for d in &mut data {
        *d -= mean;
    }
    let variance: f32 = data.iter().map(|&x| x * x).sum();
    if variance == 0.0 {
        return Some(0.0); // Or handle as NaN if variance is zero
    }
    let covariance: f32 = data.windows(2).map(|c| c[0] * c[1]).sum();
    Some(covariance / variance)
}
