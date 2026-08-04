//! Columnar scan storage and its on-disk cache.
//!
//! # Why this shape
//!
//! 0.1 modelled a run as `Vec<Ms>` where `Ms { rt: f32, mz_i: Vec<(f32, f32)> }`.
//! For a typical 40-minute DDA run that is one heap allocation *per scan*
//! (tens of thousands), each holding interleaved `(mz, intensity)` pairs. Two
//! costs follow:
//!
//! * **Allocator pressure and fragmentation.** Every scan is a separate malloc
//!   whose size is only known after decoding, so the vectors grow by doubling
//!   and leave holes behind.
//! * **Cache misses on the hot path.** Every m/z lookup is a binary search, and
//!   with interleaved pairs each probe pulls an intensity it does not need into
//!   L1, halving the useful density of every cache line. Feature detection does
//!   millions of these.
//!
//! Here the whole run is four allocations. m/z and intensity live in separate
//! columns, so a binary search touches only m/z, and the compiler can
//! auto-vectorise the intensity reductions.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use memmap2::Mmap;

use crate::error::{Error, IoContext, Result};

// The cache stores raw little-endian f32 columns and maps them without a copy.
#[cfg(target_endian = "big")]
compile_error!("metabokit-core assumes a little-endian target (x86_64 / aarch64)");

/// MS1 scans in structure-of-arrays form. Scans are in acquisition order, so
/// `rt` is ascending; within a scan, m/z is ascending (guaranteed by mzML).
#[derive(Debug, Clone)]
pub struct Ms1Set {
    rt: Vec<f32>,
    /// `len() + 1` entries; scan `i` spans `off[i]..off[i + 1]`.
    off: Vec<u32>,
    mz: Vec<f32>,
    inten: Vec<f32>,
}

impl Default for Ms1Set {
    fn default() -> Self {
        Ms1Set {
            rt: Vec::new(),
            off: vec![0],
            mz: Vec::new(),
            inten: Vec::new(),
        }
    }
}

impl Ms1Set {
    pub fn new() -> Self {
        Self::default()
    }

    /// Pre-size the columns. Callers estimate from the file length; being
    /// wrong only costs a realloc, being right removes all of them.
    pub fn with_capacity(scans: usize, points: usize) -> Self {
        let mut off = Vec::with_capacity(scans + 1);
        off.push(0);
        Ms1Set {
            rt: Vec::with_capacity(scans),
            off,
            mz: Vec::with_capacity(points),
            inten: Vec::with_capacity(points),
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.rt.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.rt.is_empty()
    }

    #[inline]
    pub fn point_count(&self) -> usize {
        self.mz.len()
    }

    /// Retention times, ascending.
    #[inline]
    pub fn rts(&self) -> &[f32] {
        &self.rt
    }

    #[inline]
    pub fn rt_at(&self, i: usize) -> f32 {
        self.rt[i]
    }

    /// `(mz, intensity)` columns for one scan.
    #[inline]
    pub fn scan(&self, i: usize) -> (&[f32], &[f32]) {
        let a = self.off[i] as usize;
        let b = self.off[i + 1] as usize;
        (&self.mz[a..b], &self.inten[a..b])
    }

    /// Append a scan. `mz` must be ascending.
    pub fn push_scan(&mut self, rt: f32, mz: &[f32], inten: &[f32]) {
        debug_assert_eq!(mz.len(), inten.len());
        self.rt.push(rt);
        self.mz.extend_from_slice(mz);
        self.inten.extend_from_slice(inten);
        self.off.push(self.mz.len() as u32);
    }

    pub fn shrink_to_fit(&mut self) {
        self.rt.shrink_to_fit();
        self.off.shrink_to_fit();
        self.mz.shrink_to_fit();
        self.inten.shrink_to_fit();
    }

    /// Bytes of heap held by this set. Surfaced in the UI so a user can see
    /// what a given run costs.
    pub fn heap_bytes(&self) -> usize {
        self.rt.capacity() * 4
            + self.off.capacity() * 4
            + self.mz.capacity() * 4
            + self.inten.capacity() * 4
    }

    /// Index of the first scan at or after `rt`.
    #[inline]
    pub fn scan_at_or_after(&self, rt: f32) -> usize {
        self.rt.partition_point(|&x| x < rt)
    }

    /// Extracted ion chromatogram.
    ///
    /// One point per MS1 scan whose retention time falls in
    /// `[rt - half_rt, rt + half_rt)`, carrying the maximum intensity found in
    /// `[mz - tol, mz + tol)`; scans with nothing in the window contribute a
    /// zero so the trace stays evenly sampled.
    ///
    /// `out` is cleared and reused — callers hold one buffer for a whole
    /// sample rather than allocating per peak.
    pub fn xic(&self, mz: f32, rt: f32, half_rt: f32, tol: f32, out: &mut Vec<(f32, f32)>) {
        out.clear();
        let lo_mz = mz - tol;
        let hi_mz = mz + tol;
        let hi_rt = rt + half_rt;
        let start = self.scan_at_or_after(rt - half_rt);
        for i in start..self.rt.len() {
            let scan_rt = self.rt[i];
            if scan_rt >= hi_rt {
                break;
            }
            let (mzs, ints) = self.scan(i);
            let p = mzs.partition_point(|&x| x < lo_mz);
            let mut best = 0.0f32;
            for k in p..mzs.len() {
                if mzs[k] >= hi_mz {
                    break;
                }
                let v = ints[k];
                if v > best {
                    best = v;
                }
            }
            out.push((scan_rt, best));
        }
    }

    /// Per-scan maximum in `[lo, hi]`, used to build the coarse m/z-slice
    /// traces that feature detection runs on. Only scans with a hit are
    /// emitted, as `(scan index, mz, intensity)`.
    pub fn slice_trace(&self, lo: f32, hi: f32, out: &mut Vec<(u32, f32, f32)>) {
        out.clear();
        for i in 0..self.rt.len() {
            let (mzs, ints) = self.scan(i);
            let p = mzs.partition_point(|&x| x < lo);
            let mut found = false;
            let mut best_i = 0.0f32;
            let mut best_mz = 0.0f32;
            for k in p..mzs.len() {
                if mzs[k] > hi {
                    break;
                }
                // `>=` so ties resolve to the last point, matching the
                // `max_by` the 0.1 pipeline used here.
                if !found || ints[k] >= best_i {
                    found = true;
                    best_i = ints[k];
                    best_mz = mzs[k];
                }
            }
            if found {
                out.push((i as u32, best_mz, best_i));
            }
        }
    }
}

/// MS2 scans, columnar. Sorted by precursor m/z after [`Ms2Set::sort_by_precursor`].
#[derive(Debug, Clone)]
pub struct Ms2Set {
    prec_mz: Vec<f32>,
    rt: Vec<f32>,
    ce: Vec<f32>,
    off: Vec<u32>,
    mz: Vec<f32>,
    inten: Vec<f32>,
}

impl Default for Ms2Set {
    fn default() -> Self {
        Ms2Set {
            prec_mz: Vec::new(),
            rt: Vec::new(),
            ce: Vec::new(),
            // Must never be empty: `scan()` reads `off[i + 1]`.
            off: vec![0],
            mz: Vec::new(),
            inten: Vec::new(),
        }
    }
}

impl Ms2Set {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.rt.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.rt.is_empty()
    }

    #[inline]
    pub fn precursors(&self) -> &[f32] {
        &self.prec_mz
    }

    #[inline]
    pub fn prec_mz(&self, i: usize) -> f32 {
        self.prec_mz[i]
    }

    #[inline]
    pub fn rt_at(&self, i: usize) -> f32 {
        self.rt[i]
    }

    #[inline]
    pub fn ce_at(&self, i: usize) -> f32 {
        self.ce[i]
    }

    #[inline]
    pub fn scan(&self, i: usize) -> (&[f32], &[f32]) {
        let a = self.off[i] as usize;
        let b = self.off[i + 1] as usize;
        (&self.mz[a..b], &self.inten[a..b])
    }

    pub fn push_scan(&mut self, prec_mz: f32, rt: f32, ce: f32, mz: &[f32], inten: &[f32]) {
        debug_assert_eq!(mz.len(), inten.len());
        self.prec_mz.push(prec_mz);
        self.rt.push(rt);
        self.ce.push(ce);
        self.mz.extend_from_slice(mz);
        self.inten.extend_from_slice(inten);
        self.off.push(self.mz.len() as u32);
    }

    pub fn heap_bytes(&self) -> usize {
        (self.prec_mz.capacity()
            + self.rt.capacity()
            + self.ce.capacity()
            + self.off.capacity()
            + self.mz.capacity()
            + self.inten.capacity())
            * 4
    }

    /// Reorder scans by ascending precursor m/z.
    ///
    /// Everything downstream binary-searches this axis. Sorting columns means
    /// permuting them all, which is a single pass per column plus one gather
    /// over the fragment arena — still far cheaper than the pointer-chasing
    /// sort of a `Vec<Msms>` that owns a `Vec` per scan.
    pub fn sort_by_precursor(&mut self) {
        let n = self.len();
        if n < 2 {
            return;
        }
        let mut order: Vec<u32> = (0..n as u32).collect();
        order.sort_unstable_by(|&a, &b| {
            let (x, y) = (self.prec_mz[a as usize], self.prec_mz[b as usize]);
            x.partial_cmp(&y)
                .unwrap_or(std::cmp::Ordering::Equal)
                // Keep the sort total so equal precursors stay in acquisition
                // order; results must not depend on sort implementation.
                .then_with(|| a.cmp(&b))
        });

        self.prec_mz = gather(&self.prec_mz, &order);
        self.rt = gather(&self.rt, &order);
        self.ce = gather(&self.ce, &order);

        let mut mz = Vec::with_capacity(self.mz.len());
        let mut inten = Vec::with_capacity(self.inten.len());
        let mut off = Vec::with_capacity(n + 1);
        off.push(0u32);
        for &o in &order {
            let a = self.off[o as usize] as usize;
            let b = self.off[o as usize + 1] as usize;
            mz.extend_from_slice(&self.mz[a..b]);
            inten.extend_from_slice(&self.inten[a..b]);
            off.push(mz.len() as u32);
        }
        self.mz = mz;
        self.inten = inten;
        self.off = off;
    }

    /// Index of the first scan with precursor m/z at or above `mz`.
    #[inline]
    pub fn prec_at_or_after(&self, mz: f32) -> usize {
        self.prec_mz.partition_point(|&x| x < mz)
    }
}

fn gather<T: Copy>(src: &[T], order: &[u32]) -> Vec<T> {
    let mut out = Vec::with_capacity(order.len());
    for &i in order {
        out.push(src[i as usize]);
    }
    out
}

// ---------------------------------------------------------------------------
// On-disk cache
// ---------------------------------------------------------------------------

const HEADER_LEN: usize = 32;
const MS1_MAGIC: &[u8; 8] = b"MKMS1\x02\0\0";
const MS2_MAGIC: &[u8; 8] = b"MKMS2\x02\0\0";
const FEATURE_MAGIC: &[u8; 8] = b"MKFEAT\x02\0";

fn pad4(n: usize) -> usize {
    (n + 3) & !3
}

fn write_header(
    w: &mut BufWriter<File>,
    magic: &[u8; 8],
    n_scans: u32,
    n_points: u32,
    meta: &[u8],
) -> std::io::Result<()> {
    w.write_all(magic)?;
    w.write_all(&n_scans.to_le_bytes())?;
    w.write_all(&n_points.to_le_bytes())?;
    w.write_all(&(meta.len() as u32).to_le_bytes())?;
    w.write_all(&[0u8; 12])?; // reserved, brings the header to 32 bytes
    w.write_all(meta)?;
    let pad = pad4(meta.len()) - meta.len();
    if pad > 0 {
        w.write_all(&[0u8; 3][..pad])?;
    }
    Ok(())
}

fn write_f32s(w: &mut BufWriter<File>, xs: &[f32]) -> std::io::Result<()> {
    // `f32` and `u8` are both `Pod`, so this is a borrow, not a conversion.
    w.write_all(bytemuck::cast_slice(xs))
}

fn write_u32s(w: &mut BufWriter<File>, xs: &[u32]) -> std::io::Result<()> {
    w.write_all(bytemuck::cast_slice(xs))
}

/// Persist MS1 scans plus the run's acquisition timestamp.
pub fn write_ms1_cache(path: &Path, set: &Ms1Set, timestamp: &str) -> Result<()> {
    let file = File::create(path).at(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    let meta = timestamp.as_bytes();
    (|| -> std::io::Result<()> {
        write_header(
            &mut w,
            MS1_MAGIC,
            set.len() as u32,
            set.point_count() as u32,
            meta,
        )?;
        write_f32s(&mut w, &set.rt)?;
        write_u32s(&mut w, &set.off)?;
        write_f32s(&mut w, &set.mz)?;
        write_f32s(&mut w, &set.inten)?;
        w.flush()
    })()
    .at(path)
}

pub fn write_ms2_cache(path: &Path, set: &Ms2Set) -> Result<()> {
    let file = File::create(path).at(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    (|| -> std::io::Result<()> {
        write_header(
            &mut w,
            MS2_MAGIC,
            set.len() as u32,
            set.mz.len() as u32,
            &[],
        )?;
        write_f32s(&mut w, &set.prec_mz)?;
        write_f32s(&mut w, &set.rt)?;
        write_f32s(&mut w, &set.ce)?;
        write_u32s(&mut w, &set.off)?;
        write_f32s(&mut w, &set.mz)?;
        write_f32s(&mut w, &set.inten)?;
        w.flush()
    })()
    .at(path)
}

/// Persist the compact per-sample feature table used by the visualizer. Six
/// column writes avoid per-feature records and make the result directly
/// memory-mappable.
pub fn write_feature_cache(path: &Path, peaks: &[crate::features::Peak]) -> Result<()> {
    let file = File::create(path).at(path)?;
    let mut w = BufWriter::with_capacity(1 << 20, file);
    (|| -> std::io::Result<()> {
        write_header(&mut w, FEATURE_MAGIC, peaks.len() as u32, 0, &[])?;
        let fields: [fn(&crate::features::Peak) -> f32; 6] = [
            |p: &crate::features::Peak| p.mz,
            |p: &crate::features::Peak| p.rt,
            |p: &crate::features::Peak| p.half_width,
            |p: &crate::features::Peak| p.coef,
            |p: &crate::features::Peak| p.shape,
            |p: &crate::features::Peak| p.smooth,
        ];
        let mut column = Vec::with_capacity(peaks.len());
        for field in fields {
            column.clear();
            column.extend(peaks.iter().map(field));
            write_f32s(&mut w, &column)?;
        }
        w.flush()
    })()
    .at(path)
}

/// A memory-mapped MS1 cache.
///
/// Gap filling and the visualizer both need random access to every sample's
/// MS1 data. Reading them all into RAM would undo the memory work done
/// elsewhere; mapping them lets the OS page cache decide what stays resident,
/// so the process holds addresses rather than data.
pub struct Ms1View {
    map: Mmap,
    n_scans: usize,
    n_points: usize,
    rt_at: usize,
    off_at: usize,
    mz_at: usize,
    inten_at: usize,
    timestamp: String,
}

#[inline]
fn f32s(bytes: &[u8]) -> &[f32] {
    // Validated at open(): the mapping is page-aligned and every offset is a
    // multiple of 4, so this cast always succeeds.
    bytemuck::try_cast_slice(bytes).unwrap_or(&[])
}

#[inline]
fn u32s(bytes: &[u8]) -> &[u32] {
    bytemuck::try_cast_slice(bytes).unwrap_or(&[])
}

fn map_file(path: &Path) -> Result<Mmap> {
    let file = File::open(path).at(path)?;
    // SAFETY: the cache is written by this process into its own output
    // directory. If something external truncates it mid-run the mapping can
    // fault; that is the standard, accepted mmap trade-off and the reason the
    // cache lives under `misc/` rather than anywhere user-managed.
    unsafe { Mmap::map(&file) }.at(path)
}

fn read_header(map: &[u8], magic: &[u8; 8], path: &Path) -> Result<(usize, usize, usize)> {
    if map.len() < HEADER_LEN || map[..8] != magic[..] {
        return Err(Error::io(
            std::io::Error::new(std::io::ErrorKind::InvalidData, "bad cache header"),
            path,
        ));
    }
    let n_scans = u32::from_le_bytes([map[8], map[9], map[10], map[11]]) as usize;
    let n_points = u32::from_le_bytes([map[12], map[13], map[14], map[15]]) as usize;
    let meta_len = u32::from_le_bytes([map[16], map[17], map[18], map[19]]) as usize;
    Ok((n_scans, n_points, HEADER_LEN + pad4(meta_len)))
}

impl Ms1View {
    pub fn open(path: &Path) -> Result<Ms1View> {
        let map = map_file(path)?;
        let (n_scans, n_points, body) = read_header(&map, MS1_MAGIC, path)?;
        let meta_len = u32::from_le_bytes([map[16], map[17], map[18], map[19]]) as usize;
        if map.len() < HEADER_LEN + meta_len {
            return Err(Error::io(
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "cache truncated"),
                path,
            ));
        }
        let timestamp =
            String::from_utf8_lossy(&map[HEADER_LEN..HEADER_LEN + meta_len]).into_owned();

        let rt_at = body;
        let off_at = rt_at + n_scans * 4;
        let mz_at = off_at + (n_scans + 1) * 4;
        let inten_at = mz_at + n_points * 4;
        let end = inten_at + n_points * 4;
        if map.len() < end {
            return Err(Error::io(
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "cache truncated"),
                path,
            ));
        }
        Ok(Ms1View {
            map,
            n_scans,
            n_points,
            rt_at,
            off_at,
            mz_at,
            inten_at,
            timestamp,
        })
    }

    pub fn timestamp(&self) -> &str {
        &self.timestamp
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n_scans
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n_scans == 0
    }

    #[inline]
    pub fn rts(&self) -> &[f32] {
        f32s(&self.map[self.rt_at..self.rt_at + self.n_scans * 4])
    }

    #[inline]
    fn offsets(&self) -> &[u32] {
        u32s(&self.map[self.off_at..self.off_at + (self.n_scans + 1) * 4])
    }

    #[inline]
    fn mzs(&self) -> &[f32] {
        f32s(&self.map[self.mz_at..self.mz_at + self.n_points * 4])
    }

    #[inline]
    fn intens(&self) -> &[f32] {
        f32s(&self.map[self.inten_at..self.inten_at + self.n_points * 4])
    }

    /// See [`Ms1Set::xic`].
    pub fn xic(&self, mz: f32, rt: f32, half_rt: f32, tol: f32, out: &mut Vec<(f32, f32)>) {
        out.clear();
        let rts = self.rts();
        let offs = self.offsets();
        let mzs = self.mzs();
        let ints = self.intens();
        if rts.is_empty() || offs.len() != rts.len() + 1 {
            return;
        }
        let lo_mz = mz - tol;
        let hi_mz = mz + tol;
        let hi_rt = rt + half_rt;
        let start = rts.partition_point(|&x| x < rt - half_rt);
        for i in start..rts.len() {
            let scan_rt = rts[i];
            if scan_rt >= hi_rt {
                break;
            }
            let a = offs[i] as usize;
            let b = offs[i + 1] as usize;
            let m = &mzs[a..b];
            let v = &ints[a..b];
            let p = m.partition_point(|&x| x < lo_mz);
            let mut best = 0.0f32;
            for k in p..m.len() {
                if m[k] >= hi_mz {
                    break;
                }
                if v[k] > best {
                    best = v[k];
                }
            }
            out.push((scan_rt, best));
        }
    }
}

/// A memory-mapped MS2 cache. Used by the visualizer to pull the spectra
/// behind a single feature without re-reading the mzML.
pub struct Ms2View {
    map: Mmap,
    n_scans: usize,
    n_points: usize,
    prec_at: usize,
    rt_at: usize,
    ce_at: usize,
    off_at: usize,
    mz_at: usize,
    inten_at: usize,
}

/// Memory-mapped detected features for one sample. The six slices borrow the
/// cache directly; opening a hundred-thousand-feature map allocates nothing.
pub struct FeatureView {
    map: Mmap,
    len: usize,
}

impl FeatureView {
    pub fn open(path: &Path) -> Result<FeatureView> {
        let map = map_file(path)?;
        let (len, _, body) = read_header(&map, FEATURE_MAGIC, path)?;
        let end = body + len * 6 * 4;
        if map.len() < end {
            return Err(Error::io(
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "cache truncated"),
                path,
            ));
        }
        Ok(FeatureView { map, len })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    fn column(&self, index: usize) -> &[f32] {
        let start = HEADER_LEN + index * self.len * 4;
        f32s(&self.map[start..start + self.len * 4])
    }

    pub fn mzs(&self) -> &[f32] {
        self.column(0)
    }

    pub fn rts(&self) -> &[f32] {
        self.column(1)
    }

    pub fn half_widths(&self) -> &[f32] {
        self.column(2)
    }

    pub fn coefficients(&self) -> &[f32] {
        self.column(3)
    }

    pub fn shapes(&self) -> &[f32] {
        self.column(4)
    }

    pub fn smoothness(&self) -> &[f32] {
        self.column(5)
    }
}

impl Ms2View {
    pub fn open(path: &Path) -> Result<Ms2View> {
        let map = map_file(path)?;
        let (n_scans, n_points, body) = read_header(&map, MS2_MAGIC, path)?;
        let prec_at = body;
        let rt_at = prec_at + n_scans * 4;
        let ce_at = rt_at + n_scans * 4;
        let off_at = ce_at + n_scans * 4;
        let mz_at = off_at + (n_scans + 1) * 4;
        let inten_at = mz_at + n_points * 4;
        let end = inten_at + n_points * 4;
        if map.len() < end {
            return Err(Error::io(
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "cache truncated"),
                path,
            ));
        }
        Ok(Ms2View {
            map,
            n_scans,
            n_points,
            prec_at,
            rt_at,
            ce_at,
            off_at,
            mz_at,
            inten_at,
        })
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n_scans
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n_scans == 0
    }

    #[inline]
    pub fn precursors(&self) -> &[f32] {
        f32s(&self.map[self.prec_at..self.prec_at + self.n_scans * 4])
    }

    #[inline]
    pub fn rts(&self) -> &[f32] {
        f32s(&self.map[self.rt_at..self.rt_at + self.n_scans * 4])
    }

    #[inline]
    pub fn ces(&self) -> &[f32] {
        f32s(&self.map[self.ce_at..self.ce_at + self.n_scans * 4])
    }

    pub fn scan(&self, i: usize) -> (&[f32], &[f32]) {
        let offs = u32s(&self.map[self.off_at..self.off_at + (self.n_scans + 1) * 4]);
        let a = offs[i] as usize;
        let b = offs[i + 1] as usize;
        let mzs = f32s(&self.map[self.mz_at..self.mz_at + self.n_points * 4]);
        let ints = f32s(&self.map[self.inten_at..self.inten_at + self.n_points * 4]);
        (&mzs[a..b], &ints[a..b])
    }
}

/// Cache file name for a sample. `stem` is the mzML file name without its
/// extension.
pub fn ms1_cache_name(stem: &str) -> String {
    format!("ms1_{stem}.mkc")
}

pub fn ms2_cache_name(stem: &str) -> String {
    format!("ms2_{stem}.mkc")
}

pub fn feature_cache_name(stem: &str) -> String {
    format!("features_{stem}.mkc")
}

#[cfg(test)]
mod feature_cache_tests {
    use super::*;

    #[test]
    fn feature_cache_round_trips_without_heap_decoding() {
        let path = std::env::temp_dir().join(format!(
            "metabokit-feature-cache-{}-{}.mkc",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let peaks = [
            crate::features::Peak {
                mz: 101.25,
                rt: 2.5,
                half_width: 0.12,
                coef: 42.0,
                shape: 0.91,
                smooth: 0.82,
            },
            crate::features::Peak {
                mz: 250.75,
                rt: 7.25,
                half_width: 0.2,
                coef: 18.0,
                shape: 0.88,
                smooth: 0.77,
            },
        ];
        write_feature_cache(&path, &peaks).unwrap();
        let view = FeatureView::open(&path).unwrap();
        assert_eq!(view.len(), 2);
        assert_eq!(view.mzs(), &[101.25, 250.75]);
        assert_eq!(view.rts(), &[2.5, 7.25]);
        assert_eq!(view.half_widths(), &[0.12, 0.2]);
        assert_eq!(view.coefficients(), &[42.0, 18.0]);
        assert_eq!(view.shapes(), &[0.91, 0.88]);
        assert_eq!(view.smoothness(), &[0.82, 0.77]);
        drop(view);
        let _ = std::fs::remove_file(path);
    }
}
