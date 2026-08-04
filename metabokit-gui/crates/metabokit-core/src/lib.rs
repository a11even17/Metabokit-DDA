//! MetaboKit-DDA processing engine.
//!
//! A rewrite of the 0.1 command-line pipeline as a library the GUI drives.
//! The science is unchanged — mzML parsing, CWT feature detection, spectral
//! library scoring, in-source-fragment and isotope detection, cross-sample
//! alignment, reporting and gap filling — but the execution model is not.
//!
//! # What changed and why
//!
//! **Nothing panics.** 0.1 reached for `unwrap()`, `expect()` and `panic!` in
//! the parser, the library reader and the aligner. Inside a rayon worker that
//! aborts the process, which is survivable for a CLI and unacceptable for a
//! GUI. Every fallible path returns [`error::Error`] with the offending file
//! named.
//!
//! **Data is columnar.** Scans, spectra and library entries were vectors of
//! structs that each owned their own heap allocations. They are now parallel
//! arrays: a run costs a handful of allocations instead of hundreds of
//! thousands, binary searches touch only the column they search, and the
//! numeric loops vectorise.
//!
//! **Parallelism moved down a level.** 0.1 parallelised across input files, so
//! peak memory scaled with core count. Feature detection now parallelises
//! across m/z slices, and samples are processed in bounded batches.
//!
//! **Large data is memory-mapped.** Per-sample scan caches and the built-in
//! binary libraries are mapped rather than read, so the OS page cache — not
//! the process heap — decides what stays resident.
//!
//! # Entry point
//!
//! [`pipeline::run`] executes a complete analysis, reporting through a
//! [`progress::Reporter`] and stopping at the next checkpoint when a
//! [`progress::Cancel`] token is tripped.

pub mod align;
pub mod discover;
pub mod error;
pub mod features;
pub mod fill;
pub mod library;
pub mod mzml;
pub mod params;
pub mod pipeline;
pub mod progress;
pub mod report;
pub mod scans;
pub mod score;

pub use discover::{DatasetScan, Note, NoteLevel};
pub use error::{Error, Result};
pub use params::{LibrarySource, Params, Polarity, Problem, TolUnit, BUILTIN_LIBRARIES};
pub use pipeline::{run, RunOutcome, SampleStat};
pub use progress::{Cancel, Event, Level, Reporter, Silent, Stage};
pub use report::ReportSummary;

/// Crate version, surfaced in the UI's about panel.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defaults_are_invalid_until_configured() {
        let p = Params::default();
        let problems = p.validate();
        assert!(
            problems.iter().any(|x| x.fatal),
            "empty params must not run"
        );
    }

    #[test]
    fn ppm_tolerance_scales_with_mass() {
        let narrow = TolUnit::Ppm.absolute(10.0, 100.0);
        let wide = TolUnit::Ppm.absolute(10.0, 1000.0);
        assert!(wide > narrow);
        assert_eq!(TolUnit::Mz.absolute(0.01, 500.0), 0.01);
    }

    #[test]
    fn median_handles_both_parities() {
        let mut odd = [3.0f32, 1.0, 2.0];
        assert_eq!(align::median(&mut odd), 2.0);
        let mut even = [4.0f32, 1.0, 3.0, 2.0];
        assert_eq!(align::median(&mut even), 2.5);
        let mut empty: [f32; 0] = [];
        assert_eq!(align::median(&mut empty), 0.0);
    }

    #[test]
    fn min_median_max_spans_the_data() {
        let mut v = [5.0f32, 1.0, 3.0];
        let [lo, mid, hi] = align::min_median_max(&mut v);
        assert_eq!(lo, 1.0);
        assert_eq!(mid, 3.0);
        assert_eq!(hi, 5.0);
    }

    #[test]
    fn scan_set_extracts_a_chromatogram() {
        let mut ms1 = scans::Ms1Set::new();
        ms1.push_scan(1.0, &[100.0, 200.0], &[10.0, 20.0]);
        ms1.push_scan(2.0, &[100.0, 200.0], &[30.0, 40.0]);
        ms1.push_scan(9.0, &[100.0], &[99.0]);

        let mut out = Vec::new();
        ms1.xic(100.0, 1.5, 1.0, 0.01, &mut out);
        assert_eq!(out, vec![(1.0, 10.0), (2.0, 30.0)]);

        // A window with no matching m/z still yields an evenly sampled trace.
        ms1.xic(555.0, 1.5, 1.0, 0.01, &mut out);
        assert_eq!(out, vec![(1.0, 0.0), (2.0, 0.0)]);
    }

    #[test]
    fn ms2_sorting_permutes_fragments_with_their_scan() {
        let mut ms2 = scans::Ms2Set::new();
        ms2.push_scan(300.0, 1.0, 20.0, &[10.0, 11.0], &[1.0, 2.0]);
        ms2.push_scan(100.0, 2.0, 25.0, &[30.0], &[3.0]);
        ms2.sort_by_precursor();

        assert_eq!(ms2.prec_mz(0), 100.0);
        assert_eq!(ms2.scan(0), (&[30.0f32][..], &[3.0f32][..]));
        assert_eq!(ms2.prec_mz(1), 300.0);
        assert_eq!(ms2.scan(1), (&[10.0f32, 11.0][..], &[1.0f32, 2.0][..]));
    }

    #[test]
    fn string_arena_round_trips() {
        let mut arena = library::StrArena::default();
        let a = arena.push("alpha");
        let b = arena.push("beta");
        let empty = arena.push("");
        assert_eq!(arena.get(a), "alpha");
        assert_eq!(arena.get(b), "beta");
        assert_eq!(arena.get(empty), "");
    }

    #[test]
    fn cancellation_is_observed() {
        let cancel = Cancel::new();
        assert!(cancel.check().is_ok());
        cancel.cancel();
        assert!(matches!(cancel.check(), Err(Error::Cancelled)));
        cancel.reset();
        assert!(cancel.check().is_ok());
    }

    #[test]
    fn params_round_trip_through_json() {
        let mut p = Params::default();
        p.ms1_tol = 7.5;
        p.libraries = vec![LibrarySource::Builtin("hmdb".into())];
        let text = p.to_json().expect("serialise");
        let back = Params::from_json(&text).expect("deserialise");
        assert_eq!(back.ms1_tol, 7.5);
        assert_eq!(back.libraries.len(), 1);
    }
}
