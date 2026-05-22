use rayon::ThreadPoolBuildError;
use std::error::Error;
use std::fmt;
use std::path::PathBuf;

#[derive(Debug)]
pub enum CompareError {
    InvalidThreadCount {
        threads: usize,
    },
    ThreadPoolBuild {
        threads: usize,
        source: ThreadPoolBuildError,
    },
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    FastaParse {
        path: PathBuf,
        source: seq_io::fasta::Error,
    },
    HeaderUtf8 {
        path: PathBuf,
        source: std::str::Utf8Error,
    },
    InvalidPath {
        path: PathBuf,
    },
    EmptyAlignment {
        path: PathBuf,
    },
    EmptySequence {
        path: PathBuf,
        sample: String,
    },
    VariableLength {
        path: PathBuf,
        sample: String,
        expected: usize,
        observed: usize,
    },
    AlignmentTooLong {
        path: PathBuf,
        length: usize,
    },
    DuplicateSample {
        path: PathBuf,
        sample: String,
    },
    SampleSetMismatch {
        path: PathBuf,
        missing: Vec<String>,
        extra: Vec<String>,
    },
}

impl fmt::Display for CompareError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CompareError::InvalidThreadCount { threads } => {
                write!(f, "thread count must be at least 1, got {threads}")
            }
            CompareError::ThreadPoolBuild { threads, source } => {
                write!(
                    f,
                    "failed to build Rayon thread pool with {threads} threads: {source}"
                )
            }
            CompareError::Io { path, source } => {
                write!(f, "failed to read alignment '{}': {source}", path.display())
            }
            CompareError::FastaParse { path, source } => {
                write!(
                    f,
                    "failed to parse FASTA alignment '{}': {source}",
                    path.display()
                )
            }
            CompareError::HeaderUtf8 { path, source } => {
                write!(
                    f,
                    "sample header in alignment '{}' is not valid UTF-8: {source}",
                    path.display()
                )
            }
            CompareError::InvalidPath { path } => {
                write!(f, "alignment path has no filename: '{}'", path.display())
            }
            CompareError::EmptyAlignment { path } => {
                write!(
                    f,
                    "alignment '{}' contains no FASTA records",
                    path.display()
                )
            }
            CompareError::EmptySequence { path, sample } => {
                write!(
                    f,
                    "sample '{sample}' in alignment '{}' has a zero-length sequence",
                    path.display()
                )
            }
            CompareError::VariableLength {
                path,
                sample,
                expected,
                observed,
            } => {
                write!(
                    f,
                    "alignment '{}' has variable sequence lengths: sample '{sample}' has length {observed}, expected {expected}",
                    path.display()
                )
            }
            CompareError::AlignmentTooLong { path, length } => {
                write!(
                    f,
                    "alignment '{}' has {length} columns, exceeding the {} column limit",
                    path.display(),
                    u32::MAX
                )
            }
            CompareError::DuplicateSample { path, sample } => {
                write!(
                    f,
                    "alignment '{}' contains duplicate sample header '{sample}'",
                    path.display()
                )
            }
            CompareError::SampleSetMismatch {
                path,
                missing,
                extra,
            } => {
                write!(
                    f,
                    "alignment '{}' sample set does not match the first alignment",
                    path.display()
                )?;
                if !missing.is_empty() {
                    write!(f, "; missing samples: {}", missing.join(", "))?;
                }
                if !extra.is_empty() {
                    write!(f, "; extra samples: {}", extra.join(", "))?;
                }
                Ok(())
            }
        }
    }
}

impl Error for CompareError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CompareError::ThreadPoolBuild { source, .. } => Some(source),
            CompareError::Io { source, .. } => Some(source),
            CompareError::FastaParse { source, .. } => Some(source),
            CompareError::HeaderUtf8 { source, .. } => Some(source),
            _ => None,
        }
    }
}
