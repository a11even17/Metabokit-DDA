use std::error::Error;
use std::fs::File;
use std::io;
use std::io::{BufRead, BufReader, Read};

fn pos_or_neg(mzml_fs: &[std::path::PathBuf]) -> Result<bool, Box<dyn Error>> {
    let line = &mut String::new();
    let mut reader = BufReader::new(File::open(&mzml_fs[0])?);
    while reader.read_line(line)? != 0 {
        let line0 = line.trim_start();
        if line0.starts_with(r#"<cvParam cvRef="MS" accession="MS:1000130"#) {
            return Ok(true);
        } else if line0.starts_with(r#"<cvParam cvRef="MS" accession="MS:1000129"#) {
            return Ok(false);
        }
        line.clear();
    }
    panic!("scan mode");
}
pub fn read_param() -> Result<crate::Param, Box<dyn Error>> {
    let param = std::fs::read_to_string("param.txt")?;
    let value = param.parse::<toml::Table>()?;

    let mzml_fs: Vec<_> = if let Ok(mut rdr) = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .delimiter(b'\t')
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path("file_order.txt")
    {
        let headers = rdr.records().next().unwrap()?;
        rdr.into_records()
            .map(std::result::Result::unwrap)
            .map(|x| std::path::Path::new(&headers[0]).join(&x[0]))
            .collect()
    } else {
        glob::glob(value["mzML_files"].as_str().unwrap())?
            .filter_map(Result::ok)
            .collect()
    };
    assert!(!mzml_fs.is_empty(), "mzML files not found");

    let mut lo_hi = value["length_of_ion_chromatogram"]
        .as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_float().unwrap() as f32);
    let ms1tol_unit = value["ms1tol"].as_array().unwrap();
    Ok(crate::Param {
        ispos: pos_or_neg(&mzml_fs).unwrap(),
        mzml_fs,
        lib: value["library"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().to_string())
            .collect(),
        ms1tol: ms1tol_unit[0].as_float().unwrap() as f32,
        ms1tol_u: ms1tol_unit[1].as_str().unwrap().to_string(),
        ms2tol: value["ms2tol"].as_float().unwrap() as f32,
        mz_shift: value["mz_shift"].as_float().unwrap() as f32,
        isf_rt_diff: 0.04,
        isf_p_diff: value["ISF_parent_mass_diff"].as_float().unwrap() as f32,
        ms2_score: value["MS2_score"].as_float().unwrap() as f32,
        min_peaks: u8::try_from(value["min_peaks"].as_integer().unwrap())?,
        rt_shift: value["RT_shift"].as_float().unwrap() as f32,
        features_only: value["features_only"].as_bool().unwrap(),
        peak_w: (lo_hi.next().unwrap(), lo_hi.next().unwrap()),
        num_t: usize::try_from(value["num_threads"].as_integer().unwrap())?,
        peak_shape: value["peak_shape"].as_float().unwrap() as f32,
        sn_score: value["sn_score"].as_float().unwrap() as f32,
        s_n1: value["S_N_1"].as_float().unwrap() as f32,
        top_only: value["top_scoring_only"].as_bool().unwrap(),
        icut: value
            .get("intensity_cutoff")
            .map(|x| x.as_float().unwrap() as f32),
        chimeric: value["chimeric_spectra"].as_bool().unwrap(),
        rt_search: 0.25,
        ms1ms2: value["MS1_MS2_pair"].as_float().unwrap() as f32,
        impute_width: value["impute_width"].as_float().unwrap() as f32,
        ex_add: value["exclude_adduct"]
            .as_array()
            .unwrap()
            .iter()
            .map(|x| x.as_str().unwrap().to_string())
            .collect(),
        match_n_fragments: usize::try_from(value["match_n_fragments"].as_integer().unwrap())?,
        i_rt: value["integration_RT"].as_float().unwrap() as f32,
        i_mz: value["integration_mz"].as_float().unwrap() as f32,
    })
}

macro_rules! unpack {
    ($sn:ident, $sn1:ident) => {
        pub fn $sn1(file: &mut BufReader<File>) -> io::Result<$sn> {
            let mut buffer = [0; std::mem::size_of::<$sn>()];
            file.read_exact(&mut buffer)?;
            Ok($sn::from_le_bytes(buffer))
        }
    };
}
unpack!(f32, unpack_f32);
unpack!(u32, unpack_u32);
unpack!(u8, unpack_u8);

pub fn get_ms1(bn: &str) -> Result<Vec<crate::Ms>, Box<dyn Error>> {
    let mut ms1_scans = Vec::new();
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("ms1_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path)?);
    buf.skip_until(b'\0')?;
    while let Ok(rt) = unpack_f32(buf) {
        let len0 = unpack_u32(buf)?;
        ms1_scans.push(crate::Ms {
            rt,
            mz_i: (0..len0)
                .map(|_| (unpack_f32(buf).unwrap(), unpack_f32(buf).unwrap()))
                .collect(),
        });
    }
    Ok(ms1_scans)
}
#[must_use]
pub fn get_chrom(cpd: (f32, f32, f32), scans: &[crate::Ms], tol: f32) -> Vec<(f32, f32)> {
    let (mz, rt, sc) = cpd;
    let rt_l = rt - sc;
    let rt_r = rt + sc;
    let mz0 = mz - tol;
    let mz1 = mz + tol;
    scans[scans.partition_point(|x| x.rt < rt_l)..]
        .iter()
        .take_while(|x| x.rt < rt_r)
        .map(|sc| {
            let pos0 = sc.mz_i.partition_point(|x| x.0 < mz0);
            (
                sc.rt,
                sc.mz_i[pos0..]
                    .iter()
                    .take_while(|x| x.0 < mz1)
                    .map(|x| x.1)
                    .reduce(f32::max)
                    .unwrap_or(0.0),
            )
        })
        .collect()
}
