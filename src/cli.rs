use clap::Parser;
use pangenome_recombination::io::{MsaListError, read_msa_list, write_recombination_table};
use pangenome_recombination::{CompareError, infer_recombination_presence};
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

pub fn launch() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("{error}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), CliError> {
    let args = Args::parse();
    let aln_paths = read_msa_list(&args.msa_list)?;
    let table = infer_recombination_presence(&aln_paths, args.threads)?;
    write_recombination_table(&table, io::stdout().lock())
        .map_err(|source| CliError::WriteStdout { source })?;
    Ok(())
}

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
    WriteStdout { source: std::io::Error },
    Compare(CompareError),
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CliError::ReadMsaList(error) => write!(f, "{error}"),
            CliError::WriteStdout { source } => {
                write!(f, "failed to write recombination table to stdout: {source}")
            }
            CliError::Compare(error) => write!(f, "{error}"),
        }
    }
}

impl Error for CliError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            CliError::ReadMsaList(error) => Some(error),
            CliError::WriteStdout { source } => Some(source),
            CliError::Compare(error) => Some(error),
        }
    }
}

impl From<MsaListError> for CliError {
    fn from(error: MsaListError) -> Self {
        CliError::ReadMsaList(error)
    }
}

impl From<CompareError> for CliError {
    fn from(error: CompareError) -> Self {
        CliError::Compare(error)
    }
}
