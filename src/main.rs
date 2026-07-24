mod align;
mod common;
mod cwt;
mod fillall;
mod parse;
mod readlib;
mod score;
const MISCDIR: &str = "misc";
const PRINT_N: usize = 20;
use rayon::prelude::*;
use std::error::Error;
struct Param {
    mzml_fs: Vec<std::path::PathBuf>,
    ispos: bool,
    lib: Vec<String>,
    ms1tol: f32,
    ms1tol_u: String,
    ms2tol: f32,
    mz_shift: f32,
    isf_rt_diff: f32,
    isf_p_diff: f32,
    ms2_score: f32,
    min_peaks: u8,
    rt_shift: f32,
    features_only: bool,
    peak_w: (f32, f32),
    num_t: usize,
    peak_shape: f32,
    sn_score: f32,
    s_n1: f32,
    top_only: bool,
    icut: Option<f32>,
    chimeric: bool,
    rt_search: f32,
    ms1ms2: f32,
    impute_width: f32,
    ex_add: Vec<String>,
    match_n_fragments: usize,
    i_rt: f32,
    i_mz: f32,
}

struct Msms {
    ms1mz: f32,
    rt: f32,
    mz_i_l: Vec<(f32, f32)>,
    ce: f32,
}
struct Ms {
    rt: f32,
    mz_i: Vec<(f32, f32)>,
}
#[derive(Clone)]
struct Par {
    dotp: f32,
    m_peaks: u8,
    pmz: f32,
    prt: f32,
}
#[derive(Clone)]
struct Lib<'a> {
    dotp: f32,
    m_peaks: u8,
    ent: &'a readlib::Ent,
    ms1mz: f32,
    rt: f32,
}
struct Ann<'a> {
    nn: usize,
    premz: f32,
    rt: f32,
    auc: f32,
    feat: bool,
    s_n: f32,
    shape: f32,
    ent_l: Vec<Lib<'a>>,
    par_l: Vec<Par>,
    mono: Option<Par>,
}
struct Spec {
    ms1mz: f32,
    rt: f32,
    mz_i_l: Vec<(f32, f32)>,
    mz_i_all: Vec<(f32, f32)>,
    ce: f32,
}
fn main() -> Result<(), Box<dyn Error>> {
    let param_t = common::read_param()?;
    let lib_ent = readlib::get_cpds(&param_t);
    let _ = std::fs::remove_dir_all(MISCDIR);
    std::fs::create_dir(MISCDIR)?;
    print_param(&param_t)?;
    rayon::ThreadPoolBuilder::new()
        .num_threads(param_t.num_t)
        .build_global()?;
    let (time_stamp, (all_dat, all_spec)): (Vec<String>, (Vec<Vec<_>>, Vec<Vec<_>>)) = param_t
        .mzml_fs
        .par_iter()
        .enumerate()
        .map(|(i, mzml_f)| {
            let bn = mzml_f.file_name().unwrap().to_str().unwrap();
            println!("{bn}");
            let (ms1_scans, mut ms2_scans, ts) = parse::mzml(mzml_f).unwrap();
            ms2_scans.sort_unstable_by(|a, b| a.ms1mz.partial_cmp(&b.ms1mz).unwrap());
            if ms2_scans.is_empty() {
                println!("No MSMS in {bn}");
                ms2_scans.push(Msms {
                    ms1mz: 0.0,
                    rt: 0.0,
                    mz_i_l: Vec::new(),
                    ce: 0.0,
                });
            }
            cwt::cwt(bn, &ms1_scans, &ms2_scans, &lib_ent, &param_t).unwrap();
            (
                ts,
                score::score(i, &param_t, bn, &lib_ent, &ms1_scans, ms2_scans).unwrap(),
            )
        })
        .unzip();
    align::align(&param_t, all_dat, &all_spec, &time_stamp)?;
    fillall::fill(&param_t, "RTseparated")?;
    fillall::fill(&param_t, "by_ID")?;
    Ok(())
}
fn print_param(param_t: &Param) -> Result<(), Box<dyn Error>> {
    use std::io::Write;
    let file_path = std::path::Path::new(MISCDIR).join("param.bin");
    let mut bufb = std::io::BufWriter::new(std::fs::File::create(file_path)?);
    bufb.write_all(&param_t.ms1ms2.to_le_bytes())?;
    bufb.write_all(&param_t.i_rt.to_le_bytes())?;
    bufb.write_all(&param_t.i_mz.to_le_bytes())?;
    Ok(())
}
