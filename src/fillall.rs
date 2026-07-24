use std::error::Error;
pub fn fill(param_t: &crate::Param, tab_name: &str) -> Result<(), Box<dyn Error>> {
    let crate::Param { mzml_fs, .. } = param_t;
    let &crate::Param {
        impute_width, i_mz, ..
    } = param_t;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path(format!("Report_{tab_name}.csv"))?;
    let headers1 = rdr.records().next().unwrap()?;
    let headers: csv::StringRecord = rdr.records().next().unwrap()?;
    let starti = headers.iter().position(|x| x == "%detected").unwrap() + 1;
    let featmzi = headers.iter().position(|x| x == "feature_m/z").unwrap();
    let mrti = headers.iter().position(|x| x == "Median RT").unwrap();
    let pos: csv::Position = rdr.records().reader().position().clone();
    let mz_rt_l: Vec<(f32, f32)> = rdr
        .records()
        .map(std::result::Result::unwrap)
        .map(|x| (x[featmzi].parse().unwrap(), x[mrti].parse().unwrap()))
        .collect();

    let qmat: Vec<Vec<f32>> = mzml_fs
        .iter()
        .enumerate()
        .map(|(ii, mzml_f)| {
            rdr.records().reader_mut().seek(pos.clone())?;
            let bn = mzml_f.file_name().unwrap().to_str().unwrap();
            let ms1_scans = crate::common::get_ms1(bn)?;
            Ok(rdr
                .records()
                .zip(&mz_rt_l)
                .map(|(rec, (med_mz, med_rt))| {
                    rec.as_ref().unwrap()[starti + ii]
                        .parse::<f32>()
                        .unwrap_or_else(|_| {
                            let p_maxi = crate::common::get_chrom(
                                (*med_mz, *med_rt, impute_width),
                                &ms1_scans,
                                i_mz,
                            );
                            p_maxi
                                .iter()
                                .zip(p_maxi.iter().skip(1))
                                .map(|((rt0, i0), (rt1, i1))| (i0 + i1) * (rt1 - rt0))
                                .sum::<f32>()
                                * 30.
                        })
                })
                .collect())
        })
        .collect::<Result<_, Box<dyn Error>>>()?;
    let mut wtr = csv::WriterBuilder::new().from_path(format!("Report_{tab_name}_fill.csv"))?;
    wtr.write_record(&headers1)?;
    wtr.write_record(&headers)?;
    rdr.records().reader_mut().seek(pos)?;
    let mut qrow = Vec::with_capacity(qmat.len());
    for (jj, rec) in rdr.records().enumerate() {
        qrow.clear();
        qrow.extend(qmat.iter().map(|x| format!("{:.1}", x[jj])));
        let line: &csv::StringRecord = rec.as_ref().unwrap();
        wtr.write_record(
            line.iter()
                .take(starti)
                .chain(qrow.iter().map(std::string::String::as_str))
                .chain(line.iter().skip(starti + mzml_fs.len())),
        )?;
    }
    Ok(())
}
