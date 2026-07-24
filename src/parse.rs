use quick_xml::events::Event;
use quick_xml::reader::Reader;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::str;
type MsTs = (Vec<crate::Ms>, Vec<crate::Msms>, String);
pub fn mzml(mzml_f: &std::path::Path) -> Result<MsTs, Box<dyn Error>> {
    let bn = mzml_f.file_name().unwrap().to_str().unwrap();
    let mut reader = Reader::from_file(mzml_f)?;
    reader.config_mut().trim_text(true);
    let mut stack = Vec::<Vec<u8>>::new();
    let mut buf = Vec::new();
    let mut mslevel = Vec::<u8>::new();
    let mut rt = f32::NAN;
    let mut ms1mz = f32::NAN;
    let mut ce = 0f32;
    let mut mz_l = Vec::<f32>::new();
    let mut i_l = Vec::<f32>::new();
    let mut buf0 = Vec::<u8>::new();
    let mut buf1 = Vec::<u8>::new();
    let mut zlibc = None;
    let mut pre64 = None;
    let mut mz_arr = None;
    let mut profile_m = true;
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("ms1_{bn}.bin"));
    let mut bw1 = BufWriter::new(File::create(file_path)?);
    let file_path = std::path::Path::new(crate::MISCDIR).join(format!("ms2_{bn}.bin"));
    let mut bw2 = BufWriter::new(File::create(file_path)?);
    let mut ms1_sc = Vec::new();
    let mut ms2_sc = Vec::new();
    let mut ts = String::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(e)) => {
                stack.push(e.local_name().as_ref().to_vec());
                if e.local_name().as_ref() == b"run" {
                    if let Ok(Some(value)) = e.try_get_attribute("startTimeStamp") {
                        bw1.write_all(&value.value)?;
                        bw1.write_all(b"\0")?;
                        ts = String::from_utf8(value.value.to_vec())?;
                    }
                    break;
                }
            }
            Ok(Event::End(_)) => {
                stack.pop();
            }
            _ => (),
        }
        buf.clear();
    }
    loop {
        match reader.read_event_into(&mut buf) {
            Err(e) => panic!("Error at position {}: {e:?}", reader.buffer_position()),
            Ok(Event::Eof) => break,
            Ok(Event::Start(e)) => {
                match e.local_name().as_ref() {
                    b"spectrum" => {
                        profile_m = true;
                        mslevel.clear();
                    }
                    b"binaryDataArray" => {
                        zlibc = None;
                        pre64 = None;
                        mz_arr = None;
                    }
                    _ => (),
                }
                stack.push(e.local_name().as_ref().to_vec());
            }
            Ok(Event::End(e)) => {
                if e.local_name().as_ref() == b"spectrum" {
                    assert!(mslevel.is_empty() || !profile_m, "profile mode detected");
                    match mslevel.as_slice() {
                        b"1" if rt.is_finite() => {
                            let ms1 = crate::Ms {
                                rt,
                                mz_i: mz_l
                                    .drain(..)
                                    .zip(i_l.drain(..))
                                    .filter(|x| x.1 > 0.)
                                    .collect(),
                            };
                            bw1.write_all(&rt.to_le_bytes())?;
                            bw1.write_all(&u32::try_from(ms1.mz_i.len())?.to_le_bytes())?;
                            for (x, y) in &ms1.mz_i {
                                bw1.write_all(&x.to_le_bytes())?;
                                bw1.write_all(&y.to_le_bytes())?;
                            }
                            ms1_sc.push(ms1);
                        }
                        b"2" if ms1mz.is_finite() && rt.is_finite() => {
                            let ms2 = crate::Msms {
                                ms1mz,
                                rt,
                                mz_i_l: mz_l
                                    .drain(..)
                                    .zip(i_l.drain(..))
                                    .filter(|x| x.1 > 0.)
                                    .collect(),
                                ce,
                            };
                            bw2.write_all(&ms1mz.to_le_bytes())?;
                            bw2.write_all(&rt.to_le_bytes())?;
                            bw2.write_all(&ce.to_le_bytes())?;
                            bw2.write_all(&u32::try_from(ms2.mz_i_l.len())?.to_le_bytes())?;
                            for (x, y) in &ms2.mz_i_l {
                                bw2.write_all(&x.to_le_bytes())?;
                                bw2.write_all(&y.to_le_bytes())?;
                            }
                            ms2_sc.push(ms2);
                        }
                        _ => (),
                    }
                    ce = 0.0;
                    ms1mz = f32::NAN;
                    rt = f32::NAN;
                }
                stack.pop();
            }
            Ok(Event::Empty(e)) if e.local_name().as_ref() == b"cvParam" => {
                let Ok(Some(accession)) = e.try_get_attribute("accession") else {
                    continue;
                };
                match accession.value.as_ref() {
                    b"MS:1000511" => {
                        if let Ok(Some(value)) = e.try_get_attribute("value") {
                            mslevel = value.value.to_vec();
                        }
                    }
                    b"MS:1000016" => {
                        if let Ok(Some(value)) = e.try_get_attribute("value") {
                            rt = str::from_utf8(&value.value)?.parse()?;
                        }
                    }
                    b"MS:1000744" => {
                        if let Ok(Some(value)) = e.try_get_attribute("value") {
                            ms1mz = str::from_utf8(&value.value)?.parse()?;
                        }
                    }
                    b"MS:1000045" => {
                        if let Ok(Some(value)) = e.try_get_attribute("value") {
                            ce = str::from_utf8(&value.value)?.parse()?;
                        }
                    }
                    b"MS:1000127" => profile_m = false,
                    b"MS:1000523" => pre64 = Some(true),
                    b"MS:1000521" => pre64 = Some(false),
                    b"MS:1000574" => zlibc = Some(true),
                    b"MS:1000576" => zlibc = Some(false),
                    b"MS:1000514" => mz_arr = Some(true),
                    b"MS:1000515" => mz_arr = Some(false),
                    _ => {}
                }
            }
            Ok(Event::Text(e))
                if stack[stack.len() - 1] == b"binary"
                    && stack[3] == b"spectrumList"
                    && !mslevel.is_empty() =>
            {
                let arr_l = if mz_arr.expect("array type not set") {
                    &mut mz_l
                } else {
                    &mut i_l
                };
                decode_bin(
                    &e.decode()?,
                    zlibc.expect("zlib not set"),
                    pre64.expect("pre64 not set"),
                    arr_l,
                    &mut buf0,
                    &mut buf1,
                )?;
            }
            _ => (),
        }
        buf.clear();
    }
    Ok((ms1_sc, ms2_sc, ts))
}
fn decode_bin(
    bin: &str,
    zlibc: bool,
    pre64: bool,
    arr_l: &mut Vec<f32>,
    buf0: &mut Vec<u8>,
    buf1: &mut Vec<u8>,
) -> std::io::Result<()> {
    let mut wrapped_reader = bin.as_bytes();
    let mut decoder = base64::read::DecoderReader::new(
        &mut wrapped_reader,
        &base64::engine::general_purpose::STANDARD,
    );
    buf0.clear();
    decoder.read_to_end(buf0)?;
    let buf2 = if zlibc {
        buf1.clear();
        flate2::read::ZlibDecoder::new(buf0.as_slice()).read_to_end(buf1)?;
        buf1
    } else {
        buf0
    };
    arr_l.clear();
    if pre64 {
        arr_l.extend(
            buf2.chunks_exact(std::mem::size_of::<f64>())
                .map(|s| f64::from_le_bytes(s.try_into().unwrap()) as f32),
        );
    } else {
        arr_l.extend(
            buf2.chunks_exact(std::mem::size_of::<f32>())
                .map(|s| f32::from_le_bytes(s.try_into().unwrap())),
        );
    }
    Ok(())
}
