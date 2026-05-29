use clap::{ArgGroup, Parser, ValueEnum};

use std::fmt;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[value(rename_all = "kebab-case")]
pub(crate) enum ParalogMode {
    First,
    Skip,
    Longest,
}

impl ParalogMode {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            ParalogMode::First => "first",
            ParalogMode::Skip => "skip",
            ParalogMode::Longest => "longest",
        }
    }
}

impl fmt::Display for ParalogMode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Debug, Parser)]
#[command(version, about = "Find recombination from pangenome alignments")]
#[command(group(
    ArgGroup::new("input")
        .required(true)
        .multiple(false)
        .args(["msa_list", "panaroo_dir"])
))]
pub(crate) struct Args {
    #[arg(long, value_name = "PATH")]
    pub msa_list: Option<PathBuf>,

    #[arg(long, value_name = "DIR")]
    pub panaroo_dir: Option<PathBuf>,

    #[arg(long, default_value_t = 1, value_parser = parse_threads)]
    pub threads: usize,

    #[arg(long, value_enum, default_value_t = ParalogMode::First)]
    pub paralog_mode: ParalogMode,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Don't show any messages
    #[arg(long, global = true)]
    pub quiet: bool,
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

/// Function to parse command line args into [`MainArgs`] struct
pub fn cli_args() -> Args {
    Args::parse()
}
