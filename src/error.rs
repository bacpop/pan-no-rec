use std::error::Error;
use std::fmt;
use std::path::PathBuf;

#[derive(Debug)]
pub enum CliError {
    ReadMsaList(MsaListError),
    WriteStdout {
        source: std::io::Error,
    },
    ThreadPoolBuild {
        threads: usize,
        source: ThreadPoolBuildError,
    },
    Compare(CompareError),
}

impl fmt::Display for CliError {
    // Formats CLI errors as user-facing messages.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CliError::ReadMsaList(error) => write!(f, "{error}"),
            CliError::WriteStdout { source } => {
                write!(f, "failed to write recombination table to stdout: {source}")
            }
            CliError::ThreadPoolBuild { threads, source } => {
                write!(
                    f,
                    "failed to initialize Rayon global thread pool with {threads} threads: {source}"
                )
            }
            CliError::Compare(error) => write!(f, "{error}"),
        }
    }
}

impl Error for CliError {
    // Exposes wrapped errors for standard error chaining.
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CliError::ReadMsaList(error) => Some(error),
            CliError::WriteStdout { source } => Some(source),
            CliError::ThreadPoolBuild { source, .. } => Some(source),
            CliError::Compare(error) => Some(error),
        }
    }
}

impl From<MsaListError> for CliError {
    // Converts MSA-list failures into the CLI error type.
    fn from(error: MsaListError) -> Self {
        CliError::ReadMsaList(error)
    }
}

impl From<CompareError> for CliError {
    // Converts recombination failures into the CLI error type.
    fn from(error: CompareError) -> Self {
        CliError::Compare(error)
    }
}

#[derive(Debug)]
pub enum CompareError {
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
    // Formats comparison errors as user-facing messages.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
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
    // Exposes wrapped parsing and I/O errors.
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CompareError::Io { source, .. } => Some(source),
            CompareError::FastaParse { source, .. } => Some(source),
            CompareError::HeaderUtf8 { source, .. } => Some(source),
            _ => None,
        }
    }
}
