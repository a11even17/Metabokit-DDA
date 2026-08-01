//! Progress reporting and cooperative cancellation.
//!
//! The engine knows nothing about Tauri. It emits `Event`s into a `Reporter`
//! the host installs, and polls a `Cancel` token at coarse checkpoints (never
//! inside an inner numeric loop, where the atomic load would cost more than the
//! work it guards).

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::error::{Error, Result};

#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    Info,
    Warn,
    Error,
}

/// A pipeline phase. Ordering matches execution order and drives the UI's
/// stage rail.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum Stage {
    Preparing,
    Library,
    /// Parsing, feature detection and scoring. These overlap across samples,
    /// so they are one stage rather than three that would appear to run
    /// backwards in the UI.
    Processing,
    Aligning,
    Reporting,
    GapFilling,
}

impl Stage {
    pub const ALL: [Stage; 6] = [
        Stage::Preparing,
        Stage::Library,
        Stage::Processing,
        Stage::Aligning,
        Stage::Reporting,
        Stage::GapFilling,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Stage::Preparing => "Preparing",
            Stage::Library => "Reading libraries",
            Stage::Processing => "Processing samples",
            Stage::Aligning => "Aligning samples",
            Stage::Reporting => "Writing reports",
            Stage::GapFilling => "Gap filling",
        }
    }
}

#[derive(Clone, Debug, serde::Serialize)]
#[serde(tag = "type", rename_all = "kebab-case")]
pub enum Event {
    /// A new phase began.
    Stage { stage: Stage, label: String },
    /// Fractional progress within the current phase, 0.0..=1.0.
    Progress { done: u64, total: u64 },
    /// A named per-sample sub-task started.
    Sample { name: String, index: usize, total: usize },
    Log { level: Level, message: String },
    /// A key/value fact worth surfacing (counts, timings, peak RSS).
    Metric { key: String, value: String },
}

/// Where the engine sends events. Implemented by the Tauri layer; `Silent` is
/// used by tests and headless runs.
pub trait Reporter: Send + Sync {
    fn emit(&self, event: Event);
}

pub struct Silent;

impl Reporter for Silent {
    fn emit(&self, _event: Event) {}
}

/// Convenience helpers so call sites stay short.
impl<'a> dyn Reporter + 'a {
    pub fn stage(&self, stage: Stage) {
        self.emit(Event::Stage {
            stage,
            label: stage.label().to_string(),
        });
    }

    pub fn info(&self, message: impl Into<String>) {
        self.emit(Event::Log {
            level: Level::Info,
            message: message.into(),
        });
    }

    pub fn warn(&self, message: impl Into<String>) {
        self.emit(Event::Log {
            level: Level::Warn,
            message: message.into(),
        });
    }

    pub fn error(&self, message: impl Into<String>) {
        self.emit(Event::Log {
            level: Level::Error,
            message: message.into(),
        });
    }

    pub fn metric(&self, key: impl Into<String>, value: impl Into<String>) {
        self.emit(Event::Metric {
            key: key.into(),
            value: value.into(),
        });
    }

    pub fn progress(&self, done: u64, total: u64) {
        self.emit(Event::Progress { done, total });
    }
}

/// Cooperative cancellation token, cheap to clone across rayon workers.
#[derive(Clone, Default)]
pub struct Cancel(Arc<AtomicBool>);

impl Cancel {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    pub fn reset(&self) {
        self.0.store(false, Ordering::Relaxed);
    }

    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }

    /// Checkpoint. Call between units of work, not inside them.
    #[inline]
    pub fn check(&self) -> Result<()> {
        if self.is_cancelled() {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
}
