
#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    std::env::set_current_dir(std::env::current_exe().unwrap().parent().unwrap()).unwrap();
    tauri::Builder::default()
        .plugin(tauri_plugin_opener::init())
        .invoke_handler(tauri::generate_handler![
            get_mzml_list,
            read_param,
            ms1feat,
            get_spec,
            get_mirror,
            get_ms1,
            read_report,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
use std::error::Error;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path, PathBuf};

#[derive(serde::Serialize)]
struct Peak {
    ms1mz: f32,
    rt: f32,
    sc: f32,
    shape: f32,
    smooth: f32,
    name_l: Vec<String>,
}

#[derive(Clone, serde::Serialize)]
struct Msms {
    ms1mz: f32,
    rt: f32,
    mz_i_l: Vec<(f32, f32)>,
    ce: f32,
}
const MISCDIR: &str = "misc";

#[tauri::command]
fn get_mzml_list() -> Vec<String> {
    let file_path: PathBuf = [MISCDIR, "ms1f_*.bin"].iter().collect();
    glob::glob(file_path.to_str().unwrap())
        .unwrap()
        .filter_map(Result::ok)
        .map(|x| {
            x.file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .strip_prefix("ms1f_")
                .unwrap()
                .to_string()
        })
        .collect()
}
#[derive(serde::Serialize)]
struct Param {
    ms1ms2: f32,
    i_rt: f32,
    i_mz: f32,
}
#[tauri::command]
fn read_param() -> Param {
    let file_path = Path::new(MISCDIR).join("param.bin");
    let bufw = &mut BufReader::new(File::open(file_path).unwrap());
    Param {
        ms1ms2: unpack_f32(bufw).unwrap(),
        i_rt: unpack_f32(bufw).unwrap(),
        i_mz: unpack_f32(bufw).unwrap(),
    }
}

macro_rules! unpack {
    ($sn:ident, $sn1:ident) => {
        fn $sn1(file: &mut BufReader<File>) -> io::Result<$sn> {
            let mut buffer = [0; std::mem::size_of::<$sn>()];
            file.read_exact(&mut buffer)?;
            Ok($sn::from_le_bytes(buffer))
        }
    };
}
unpack!(f32, unpack_f32);
unpack!(u32, unpack_u32);
unpack!(u16, unpack_u16);
unpack!(u8, unpack_u8);

fn unpack_f32_2(file: &mut BufReader<File>) -> io::Result<(f32, f32)> {
    let mut buffer = [0; std::mem::size_of::<f32>()];
    file.read_exact(&mut buffer)?;
    let a = f32::from_le_bytes(buffer);
    file.read_exact(&mut buffer)?;
    Ok((a, f32::from_le_bytes(buffer)))
}

#[tauri::command]
fn get_spec(bn: &str, ms1mz: f32, ms1rt: f32, rtwid: f32, ms1ms2: f32) -> Vec<Msms> {
    let mut ms2_scans = Vec::new();
    let file_path = Path::new(MISCDIR).join(format!("ms2_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path).unwrap());
    while let Ok(mz) = unpack_f32(buf) {
        let rt = unpack_f32(buf).unwrap();
        let ce = unpack_f32(buf).unwrap();
        let len0 = unpack_u32(buf).unwrap();
        if (mz - ms1mz).abs() < ms1ms2 && (rt - ms1rt).abs() < rtwid {
            ms2_scans.push(Msms {
                ms1mz: mz,
                rt,
                mz_i_l: (0..len0)
                    .map(|_| unpack_f32_2(buf))
                    .collect::<Result<_, _>>()
                    .unwrap(),
                ce,
            });
        } else {
            buf.seek_relative(i64::from(len0) * 8).unwrap();
        }
    }
    ms2_scans
}
#[tauri::command]
fn get_ms1(bn: &str, ms1mz: f32, ms1rt: f32, rtwid: f32, i_mz: f32) -> Vec<(f32, f32)> {
    let mut ms1_scans = Vec::new();
    let file_path = Path::new(MISCDIR).join(format!("ms1_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path).unwrap());
    buf.skip_until(b'\0').unwrap();
    while let Ok(rt) = unpack_f32(buf) {
        let len0 = unpack_u32(buf).unwrap();
        if rt < ms1rt - rtwid {
            buf.seek_relative(i64::from(len0) * 8).unwrap();
        } else if rt < ms1rt + rtwid {
            let pos0 = ms1mz - i_mz;
            let pos1 = ms1mz + i_mz;
            let mut ii = 0;
            ms1_scans.push((
                rt,
                (0..len0)
                    .map(|i| {
                        ii = i;
                        unpack_f32_2(buf).unwrap()
                    })
                    .skip_while(|x| x.0 <= pos0)
                    .take_while(|x| x.0 < pos1)
                    .map(|x| x.1)
                    .reduce(f32::max)
                    .unwrap_or(0.),
            ));
            if len0 > 1 {
                buf.seek_relative(i64::from(len0 - ii - 1) * 8).unwrap();
            }
        } else {
            break;
        }
    }
    ms1_scans
}

#[tauri::command]
fn ms1feat(bn: &str, ms1ms2: f32) -> Vec<Peak> {
    let spec_l = get_all_spec(bn).unwrap();
    let matched_l = get_all_mirror(bn).unwrap();
    let file_path = Path::new(MISCDIR).join(format!("ms1f_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path).unwrap());
    let len0 = unpack_u32(buf).unwrap();
    (0..len0)
        .filter_map(|_| {
            let ms1mz = unpack_f32(buf).unwrap();
            let rt = unpack_f32(buf).unwrap();
            let p0 = ms1mz - 0.0001;
            let p0 = matched_l.partition_point(|x| x.mz < p0);
            let p1 = ms1mz + 0.0001;
            let sc = unpack_f32(buf).unwrap();
            buf.seek_relative(4).unwrap();
            let shape = unpack_f32(buf).unwrap();
            let smooth = unpack_f32(buf).unwrap();
            let pos0 = ms1mz - ms1ms2;
            let pos0 = spec_l.partition_point(|x| x.ms1mz < pos0);
            let pos1 = ms1mz + ms1ms2;
            (spec_l[pos0..]
                .iter()
                .take_while(|x| x.ms1mz < pos1)
                .any(|x| (x.rt - rt).abs() < sc * 2.))
            .then(|| {
                let mut name_l: Vec<&Matched> = matched_l[p0..]
                    .iter()
                    .take_while(|x| x.mz < p1)
                    .filter(|x| (rt - x.rt).abs() < 0.0001)
                    .collect();
                name_l.sort_unstable_by(|y, x| x.dotp.partial_cmp(&y.dotp).unwrap());
                Peak {
                    ms1mz,
                    rt,
                    name_l: name_l.into_iter().map(|x| x.name.clone()).collect(),
                    shape,
                    smooth,
                    sc,
                }
            })
        })
        .collect()
}
fn unpack_string(file: &mut BufReader<File>) -> Result<String, Box<dyn Error>> {
    let mut str_buf = Vec::new();
    file.read_until(b'\0', &mut str_buf)?;
    str_buf.pop().ok_or("unpack_string")?;
    Ok(String::from_utf8(str_buf)?)
}

struct Matched {
    mz: f32,
    rt: f32,
    name: String,
    dotp: f32,
}
fn get_all_spec(bn: &str) -> io::Result<Vec<Msms>> {
    let mut ms2_scans = Vec::new();
    let file_path = Path::new(MISCDIR).join(format!("ms2_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path)?);
    while let Ok(ms1mz) = unpack_f32(buf) {
        let rt = unpack_f32(buf)?;
        let ce = unpack_f32(buf)?;
        let len0 = unpack_u32(buf)?;
        buf.seek_relative(i64::from(len0) * 8)?;
        ms2_scans.push(Msms {
            ms1mz,
            rt,
            ce,
            mz_i_l: Vec::new(),
        });
    }
    ms2_scans.sort_unstable_by(|x, y| x.ms1mz.partial_cmp(&y.ms1mz).unwrap());
    Ok(ms2_scans)
}
fn get_all_mirror(bn: &str) -> io::Result<Vec<Matched>> {
    let file_path = Path::new(MISCDIR).join(format!("plot_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path)?);
    let mut ma_vec = Vec::new();
    while let Ok(name) = unpack_string(buf) {
        let premz = unpack_f32(buf)?;
        buf.seek_relative(8)?;
        let peakrt = unpack_f32(buf)?;
        buf.seek_relative(12)?;
        let len0 = unpack_u16(buf)?;
        buf.seek_relative(i64::from(len0) * 8)?;
        buf.seek_relative(4)?;
        let dotp = unpack_f32(buf)?;
        let len0 = unpack_u8(buf)?;
        buf.seek_relative(i64::from(len0) * 8)?;
        buf.seek_relative(4)?;
        let len0 = unpack_u8(buf)?;
        buf.seek_relative(i64::from(len0) * 8)?;
        let len0 = unpack_u8(buf)?;
        buf.seek_relative(i64::from(len0) * 16)?;
        if premz <= 0. {
            break;
        }
        ma_vec.push(Matched {
            mz: premz,
            rt: peakrt,
            name,
            dotp,
        });
    }
    ma_vec.sort_unstable_by(|x, y| x.mz.partial_cmp(&y.mz).unwrap());
    Ok(ma_vec)
}

#[tauri::command]
fn read_report() -> Vec<Peak> {
    let Ok(mut rdr) = csv::ReaderBuilder::new()
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path("Report_RTseparated.csv")
    else {
        return Vec::new();
    };
    rdr.records().next().unwrap().unwrap();
    let headers: csv::StringRecord = rdr.records().next().unwrap().unwrap();
    let featmzi = headers.iter().position(|x| x == "feature_m/z").unwrap();
    let mrti = headers.iter().position(|x| x == "Median RT").unwrap();
    let namei = headers.iter().position(|x| x == "name").unwrap();
    rdr.into_records()
        .map(std::result::Result::unwrap)
        .filter(|x| !x[namei].starts_with("ISF of "))
        .map(|x| Peak {
            ms1mz: x[featmzi].parse().unwrap(),
            rt: x[mrti].parse().unwrap(),
            name_l: x[namei]
                .split(" --- ")
                .filter(|x| !x.is_empty())
                .map(std::string::ToString::to_string)
                .collect(),
            shape: 0.,
            smooth: 0.,
            sc: 0.,
        })
        .collect()
}
#[derive(serde::Serialize)]
struct Mirror {
    name: String,
    ms1mz: f32,
    lib_mass: f32,
    dotp: f32,
    ce: f32,
    specmz: f32,
    specrt: f32,
    exp_mz_i_l: Vec<(f32, f32)>,
    lib_mz_i_l: Vec<(f32, f32)>,
    m_exp_mz_i_l: Vec<(f32, f32)>,
    m_lib_mz_i_l: Vec<(f32, f32)>,
}
#[tauri::command]
fn get_mirror(bn: &str, ms1mz: f32, ms1rt: f32) -> Vec<Mirror> {
    let file_path = Path::new(MISCDIR).join(format!("plot_{bn}.bin"));
    let buf = &mut BufReader::new(File::open(file_path).unwrap());
    let mut ma_vec = Vec::new();
    while let Ok(name) = unpack_string(buf) {
        let premz = unpack_f32(buf).unwrap();
        buf.seek_relative(8).unwrap();
        let peakrt = unpack_f32(buf).unwrap();
        let specmz = unpack_f32(buf).unwrap();
        let specrt = unpack_f32(buf).unwrap();
        buf.seek_relative(4).unwrap();
        let len0 = unpack_u16(buf).unwrap();
        buf.seek_relative(i64::from(len0) * 8).unwrap();
        let lib_mass = unpack_f32(buf).unwrap();
        let dotp = unpack_f32(buf).unwrap();
        let len0 = unpack_u8(buf).unwrap();
        let exp_mz_i_l: Vec<_> = (0..len0)
            .map(|_| unpack_f32_2(buf))
            .collect::<Result<_, _>>()
            .unwrap();
        let ce = unpack_f32(buf).unwrap();
        let len0 = unpack_u8(buf).unwrap();
        let mut lib_mz_i_l: Vec<_> = (0..len0)
            .map(|_| unpack_f32_2(buf))
            .collect::<Result<_, _>>()
            .unwrap();
        let len0 = unpack_u8(buf).unwrap();
        let mut m_lib_mz_i_l: Vec<_> = (0..len0)
            .map(|_| unpack_f32_2(buf))
            .collect::<Result<_, _>>()
            .unwrap();
        let m_exp_mz_i_l: Vec<_> = (0..len0)
            .map(|_| unpack_f32_2(buf))
            .collect::<Result<_, _>>()
            .unwrap();
        if premz <= 0. {
            break;
        }
        if (premz - ms1mz).abs() < 0.0001 && (peakrt - ms1rt).abs() < 0.0001 {
            let maxi = lib_mz_i_l
                .iter()
                .find(|x| x.0 < premz - 0.1)
                .map_or(lib_mz_i_l[0].1, |x| x.1);
            for mz_i in &mut lib_mz_i_l {
                mz_i.1 /= maxi;
            }
            for mz_i in &mut m_lib_mz_i_l {
                mz_i.1 /= maxi;
            }
            let mut sorted: Vec<_> = m_exp_mz_i_l.into_iter().zip(m_lib_mz_i_l).collect();
            sorted.sort_unstable_by(|a, b| a.0 .0.partial_cmp(&b.0 .0).unwrap());
            ma_vec.push(Mirror {
                ms1mz: premz,
                lib_mass,
                name,
                dotp,
                ce,
                specmz,
                specrt,
                exp_mz_i_l,
                lib_mz_i_l,
                m_exp_mz_i_l: sorted.iter().map(|x| x.0).collect(),
                m_lib_mz_i_l: sorted.iter().map(|x| x.1).collect(),
            });
        }
    }
    ma_vec.sort_unstable_by(|y, x| x.dotp.partial_cmp(&y.dotp).unwrap());
    ma_vec.truncate(38);
    ma_vec
}
