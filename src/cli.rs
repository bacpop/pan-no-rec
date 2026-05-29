use clap::Parser;
use pangenome_recombination::io::{MsaListError, read_msa_list, write_recombination_table};
use pangenome_recombination::{CompareError, infer_recombination_presence};
use rayon::{ThreadPoolBuildError, ThreadPoolBuilder};
use std::error::Error;
use std::fmt;
use std::io;
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Debug, Parser)]
#[command(version, about = "Compare pangenome multiple-sequence alignments")]
struct Args {
    #[arg(long, value_name = "PATH")]
    msa_list: PathBuf,

    #[arg(long, default_value_t = 1, value_parser = parse_threads)]
    threads: usize,
}

// Runs the CLI and converts its result into a process exit code.
pub fn launch() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

// Parses arguments, runs recombination inference, and writes the TSV table.
fn run() -> Result<(), CliError> {
    let args = Args::parse();
    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .map_err(|source| CliError::ThreadPoolBuild {
            threads: args.threads,
            source,
        })?;
    let aln_paths = read_msa_list(&args.msa_list)?;
    let table = infer_recombination_presence(&aln_paths)?;
    write_recombination_table(&table, io::stdout().lock())
        .map_err(|source| CliError::WriteStdout { source })?;
    Ok(())
}

// Parses and validates a positive thread count for clap.
fn parse_threads(value: &str) -> Result<usize, String> {
    let threads = value
        .parse::<usize>()
        .map_err(|source| format!("invalid thread count: {source}"))?;

    if threads == 0 {
        Err("thread count must be at least 1".to_string())
    } else {
        Ok(threads)
    }
}

#[derive(Debug)]
enum CliError {
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
