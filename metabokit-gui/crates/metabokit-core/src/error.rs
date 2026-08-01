//! A concrete error type.
//!
//! The 0.1 CLI used `Box<dyn Error>` everywhere and `unwrap()`/`panic!` in the
//! hot paths. A GUI cannot afford either: a panic inside a worker thread takes
//! the whole run down with no message the user can act on, and boxed errors
//! allocate on every `?` in code that runs millions of times.

use std::fmt;
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub enum Error {
    /// I/O failure, carrying the path so the UI can name the offending file.
    Io(std::io::Error, Option<PathBuf>),
    Csv(csv::Error),
    /// Malformed or unsupported mzML.
    Mzml { file: PathBuf, detail: String },
    /// A binary data array could not be base64/zlib decoded.
    Decode(String),
    /// Invalid or contradictory run parameters.
    Param(String),
    /// A spectral library could not be read.
    Library { source: String, detail: String },
    /// The user asked to stop.
    Cancelled,
}

impl Error {
    pub fn io(err: std::io::Error, path: impl AsRef<Path>) -> Self {
        Error::Io(err, Some(path.as_ref().to_path_buf()))
    }

    pub fn param(msg: impl Into<String>) -> Self {
        Error::Param(msg.into())
    }

    pub fn mzml(file: impl AsRef<Path>, detail: impl Into<String>) -> Self {
        Error::Mzml {
            file: file.as_ref().to_path_buf(),
            detail: detail.into(),
        }
    }

    pub fn library(source: impl Into<String>, detail: impl Into<String>) -> Self {
        Error::Library {
            source: source.into(),
            detail: detail.into(),
        }
    }

    pub fn is_cancelled(&self) -> bool {
        matches!(self, Error::Cancelled)
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e, Some(p)) => write!(f, "{}: {e}", p.display()),
            Error::Io(e, None) => write!(f, "{e}"),
            Error::Csv(e) => write!(f, "csv: {e}"),
            Error::Mzml { file, detail } => {
                write!(f, "{}: {detail}", file.display())
            }
            Error::Decode(d) => write!(f, "binary array: {d}"),
            Error::Param(d) => write!(f, "parameter: {d}"),
            Error::Library { source, detail } => write!(f, "library {source}: {detail}"),
            Error::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e, _) => Some(e),
            Error::Csv(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e, None)
    }
}

impl From<csv::Error> for Error {
    fn from(e: csv::Error) -> Self {
        Error::Csv(e)
    }
}

pub type Result<T> = std::result::Result<T, Error>;

/// Attach a path to an `io::Result`.
pub trait IoContext<T> {
    fn at(self, path: impl AsRef<Path>) -> Result<T>;
}

impl<T> IoContext<T> for std::io::Result<T> {
    fn at(self, path: impl AsRef<Path>) -> Result<T> {
        self.map_err(|e| Error::io(e, path))
    }
}
