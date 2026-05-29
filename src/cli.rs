use clap::Parser;

use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(version, about = "Find recombination from pangenome alignments")]
struct Args {
    #[arg(long, value_name = "PATH")]
    pub msa_list: PathBuf,

    #[arg(long, default_value_t = 1, value_parser = parse_threads)]
    pub threads: usize,

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
