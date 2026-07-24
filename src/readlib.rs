use regex::Regex;
use std::cmp::Ordering;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub struct Ent {
    pub mmass: f32,
    pub name: String,
    pub mz_i: Vec<(f32, f32)>,
    pub adduct: String,
    pub charge: i8,
    pub rt: Option<f32>,
    pub inchik: String,
    pub formu: String,
}

fn unpack_string(file: &mut BufReader<File>, len0: usize) -> Result<String, Box<dyn Error>> {
    let mut buffer = vec![0; len0];
    file.read_exact(&mut buffer)?;
    Ok(String::from_utf8(buffer)?)
}
fn read_bin(ispos: bool, mut infile: &str, ent_vec: &mut Vec<Ent>) -> Result<(), Box<dyn Error>> {
    let atlas_f = infile == "Atlas_filtered";
    if atlas_f {
        infile = "Atlas";
    }
    let infile = std::env::current_exe()?
        .parent()
        .unwrap()
        .join("libs")
        .join(format!("{infile}_{}", if ispos { "pos" } else { "neg" }));
    let re = Regex::new(r"\d\d?:\d\d?").unwrap();
    let atlas_f = infile.file_stem().unwrap().to_str().unwrap() == "Atlas_pos" && atlas_f;
    print!("\"{}\" ", infile.file_stem().unwrap().display());
    let buf = &mut BufReader::new(File::open(infile)?);
    let current_sz = ent_vec.len();
    while let Ok(mmass) = crate::common::unpack_f32(buf) {
        let mut buffer = [0; std::mem::size_of::<i8>()];
        buf.read_exact(&mut buffer)?;
        let charge = i8::from_le_bytes(buffer);
        let len0 = crate::common::unpack_u8(buf)?;
        let name = unpack_string(buf, len0.into())?;
        let len0 = crate::common::unpack_u8(buf)?;
        let adduct = unpack_string(buf, len0.into())?;
        let len0 = crate::common::unpack_u8(buf)?;
        let inchik = unpack_string(buf, len0.into())?;
        let len0 = crate::common::unpack_u8(buf)?;
        let mz = (0..len0)
            .map(|_| crate::common::unpack_f32(buf))
            .collect::<Result<Vec<f32>, _>>()?;
        let i = (0..len0)
            .map(|_| crate::common::unpack_f32(buf))
            .collect::<Result<Vec<f32>, _>>()?;
        if atlas_f
            && !["CAR ", "VAE ", "CE ", "CL "]
                .iter()
                .any(|x| name.starts_with(x))
        {
            let matches: Vec<(u8, u8)> = re
                .find_iter(&name)
                .map(|x| x.as_str().split_once(':').unwrap())
                .map(|x| (x.0.parse().unwrap(), x.1.parse().unwrap()))
                .collect();
            let fa_filter = || matches.iter().any(|x| !(12..=26).contains(&x.0) || x.1 > 6);
            match matches.len().cmp(&1) {
                Ordering::Greater => {
                    if fa_filter() {
                        continue;
                    }
                }
                Ordering::Equal => {
                    if name.starts_with("LP") {
                        if fa_filter() {
                            continue;
                        }
                    } else if !(24..=52).contains(&matches[0].0) {
                        continue;
                    }
                }
                Ordering::Less => {}
            }
        }
        ent_vec.push(Ent {
            mmass,
            mz_i: {
                if atlas_f {
                    mz.into_iter()
                        .zip(i)
                        .filter(|x| x.0 < mmass - 18.1)
                        .collect()
                } else {
                    mz.into_iter().zip(i).collect()
                }
            },
            name,
            formu: String::new(),
            adduct,
            charge,
            rt: None,
            inchik,
        });
    }
    println!("{}", ent_vec.len() - current_sz);
    Ok(())
}

fn read_lib(x: &str, ent_vec: &mut Vec<Ent>) {
    let libpath = x.strip_prefix("user ").unwrap().trim();
    print!("{libpath:?} ");
    let mut cursor = BufReader::new(File::open(libpath).unwrap());
    let mut name = String::new();
    let mut formu = String::new();
    let mut inchik = String::new();
    let mut adduct = String::new();
    let mut rt = None;
    let mut mmass = 0.0;
    let mut charge = 1i8;
    let re = Regex::new(r"\[(.*)\](.)").unwrap();
    let mut line = String::new();
    let mut line0 = String::new();
    let mut mz_i: Vec<(f32, f32, bool)> = Vec::new();
    let current_sz = ent_vec.len();
    while cursor.read_line(&mut line).unwrap() > 0 {
        let Some((lsp0, lsp1)) = line.split_once(": ") else {
            line.clear();
            continue;
        };
        let lsp1 = lsp1.trim();
        match lsp0.trim().to_uppercase().as_str() {
            "NAME" => name = lsp1.to_string(),
            "FORMULA" => formu = lsp1.to_string(),
            "PRECURSORMZ" => mmass = lsp1.parse().unwrap(),
            "PRECURSORTYPE" | "PRECURSOR_TYPE" => {
                (adduct, charge) = re.captures(lsp1).map_or_else(
                    || (lsp1.to_string(), 1),
                    |caps| (caps[1].to_string(), caps[2].parse().unwrap_or(1)),
                );
                if !name.contains(&adduct) {
                    name = format!("{name} {adduct}");
                }
            }
            "RETENTIONTIME" => rt = lsp1.parse().ok(),
            "INCHIKEY" => inchik = lsp1.to_string(),
            "NUM PEAKS" => {
                use std::fmt::Write;
                write!(&mut name, " ({libpath})").unwrap();
                line0.clear();
                while cursor.read_line(&mut line0).unwrap() > 0 {
                    if line0.trim().is_empty() {
                        break;
                    }
                    let mut iter = line0.split_whitespace();
                    mz_i.push((
                        iter.next().unwrap().parse().unwrap(),
                        iter.next().unwrap().parse().unwrap(),
                        iter.next().is_some_and(|x| x == "*"),
                    ));
                    line0.clear();
                }
                let mut mz_i: Vec<_> = if mz_i.iter().any(|x| x.2) {
                    mz_i.drain(..).filter(|x| x.2).map(|x| (x.0, x.1)).collect()
                } else {
                    mz_i.drain(..).map(|x| (x.0, x.1)).collect()
                };
                mz_i.sort_unstable_by(|b, a| a.1.partial_cmp(&b.1).unwrap());
                ent_vec.push(Ent {
                    mmass,
                    name: std::mem::take(&mut name),
                    formu: std::mem::take(&mut formu),
                    mz_i,
                    adduct: std::mem::take(&mut adduct),
                    charge,
                    rt,
                    inchik: std::mem::take(&mut inchik),
                });
            }
            _ => (),
        }
        line.clear();
    }
    println!("{}", ent_vec.len() - current_sz);
}

fn read_mz_rt(x: &str, ent_vec: &mut Vec<Ent>) {
    let libpath = x.strip_prefix("csv ").unwrap().trim();
    print!("{libpath:?} ");
    let rdr = csv::ReaderBuilder::new()
        .comment(Some(b'#'))
        .has_headers(true)
        .trim(csv::Trim::All)
        .from_path(libpath)
        .unwrap();
    let current_sz = ent_vec.len();
    ent_vec.extend(
        rdr.into_records()
            .map(std::result::Result::unwrap)
            .map(|x| {
                let mut name = x[0].to_string();
                let adduct = x[1].to_string();
                if !name.contains(&adduct) {
                    name = format!("{name} {adduct}");
                }
                Ent {
                    mmass: x[2].parse().unwrap(),
                    name,
                    formu: String::new(),
                    mz_i: Vec::new(),
                    adduct,
                    charge: 1,
                    rt: x[3].parse().ok(),
                    inchik: String::new(),
                }
            }),
    );
    println!("{}", ent_vec.len() - current_sz);
}

#[must_use]
pub fn get_cpds(param_t: &crate::Param) -> Vec<Ent> {
    let &crate::Param { ispos, .. } = param_t;
    let crate::Param { lib: lib_types, .. } = param_t;
    let mut lib_ent = Vec::new();
    for g_lib in [
        "nist",
        "Atlas",
        "Atlas_filtered",
        "MSDIAL",
        "hmdb",
        "sling",
        "MassBank",
        "FiehnHILIC",
    ] {
        if lib_types
            .iter()
            .any(|x| x.to_lowercase() == g_lib.to_lowercase())
        {
            read_bin(ispos, g_lib, &mut lib_ent).unwrap();
        }
    }
    for x in lib_types.iter().filter(|x| x.starts_with("user ")) {
        read_lib(x, &mut lib_ent);
    }
    for x in lib_types.iter().filter(|x| x.starts_with("csv ")) {
        read_mz_rt(x, &mut lib_ent);
    }
    let rm_i: Vec<_> = lib_ent
        .iter()
        .enumerate()
        .filter(|(_, ent)| param_t.ex_add.contains(&ent.adduct))
        .map(|(i, _)| i)
        .collect();
    for i in rm_i.into_iter().rev() {
        lib_ent.swap_remove(i);
    }
    println!("total: {}", lib_ent.len());
    for x in &mut lib_ent {
        x.mz_i = x
            .mz_i
            .drain(..)
            .filter(|y| y.0 < x.mmass - 0.3)
            .take(param_t.match_n_fragments)
            .collect();
    }
    lib_ent.sort_unstable_by(|a, b| a.mmass.partial_cmp(&b.mmass).unwrap());
    lib_ent
}
