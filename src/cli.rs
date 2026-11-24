use clap::{Parser, ValueEnum};

use std::fmt;
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[value(rename_all = "kebab-case")]
pub(crate) enum ParalogMode {
    /// Take the first paralog sequence
    First,
    /// Skip any paralogs entirely
    Skip,
    /// Take the sequence with the fewest gaps
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
pub(crate) struct Args {
    /// Directory with panaroo output
    #[arg(long, value_name = "DIR")]
    pub panaroo_dir: PathBuf,

    /// Optional debug segment TSV matching stdout, or the joint TSV with --igrs
    #[arg(long, value_name = "FILE")]
    pub layout_segments_output: Option<PathBuf>,

    /// Directory with panIGRoo output and aligned_IGR_sequences
    #[arg(long, value_name = "DIR")]
    pub igrs: Option<PathBuf>,

    /// Output filename prefix for combined gene and IGR TSVs (required with --igrs)
    #[arg(long, value_name = "PREFIX")]
    pub output: Option<PathBuf>,

    /// Filter alignments above a certain entropy
    #[arg(long, value_name = "FLOAT", value_parser = parse_max_entropy)]
    pub max_entropy: Option<f64>,

    /// How to handle duplicated genes
    #[arg(long, value_enum, default_value_t = ParalogMode::First)]
    pub paralog_mode: ParalogMode,

    /// Location to write report on paralogs
    #[arg(long, value_name = "PATH", default_value = "paralogs.txt")]
    pub paralog_report: PathBuf,

    /// Include gaps as SNPs
    #[arg(long, default_value_t = false)]
    pub gaps: bool,

    /// CPU threads to use
    #[arg(long, default_value_t = 1, value_parser = parse_threads)]
    pub threads: usize,

    /// Show progress messages
    #[arg(short, long, global = true)]
    pub verbose: bool,

    /// Don't show any messages
    #[arg(long, global = true)]
    pub quiet: bool,

    /// Directory to write Cytoscape-compatible GraphML pruning snapshots
    #[arg(long, value_name = "DIR")]
    pub debug_graph_dir: Option<PathBuf>,

    /// Gene names to write GraphML pruning snapshots for
    #[arg(long, value_name = "GENE", value_delimiter = ',', value_parser = parse_debug_graph_gene)]
    pub debug_graph_genes: Vec<String>,

    /// Directory for recombination-free gene and core alignments (enables alignment output)
    #[arg(long, value_name = "DIR")]
    pub alignment_output: Option<PathBuf>,

    /// Minimum sample frequency for the concatenated core alignment
    #[arg(long, value_name = "FLOAT", default_value = "0.95", value_parser = parse_core_threshold)]
    pub core_threshold: f64,
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

// Parses an optional finite entropy threshold for filtering Panaroo alignments.
fn parse_max_entropy(value: &str) -> Result<f64, String> {
    let threshold = value
        .parse::<f64>()
        .map_err(|source| format!("invalid entropy threshold: {source}"))?;

    if threshold.is_finite() {
        Ok(threshold)
    } else {
        Err("entropy threshold must be finite".to_string())
    }
}

// Parses one non-empty gene name for graph debugging.
fn parse_debug_graph_gene(value: &str) -> Result<String, String> {
    let gene = value.trim();
    if gene.is_empty() {
        Err("debug graph gene name must not be empty".to_string())
    } else {
        Ok(gene.to_string())
    }
}

// Parses a core frequency threshold in the open-closed interval (0, 1].
fn parse_core_threshold(value: &str) -> Result<f64, String> {
    let threshold = value
        .parse::<f64>()
        .map_err(|source| format!("invalid core threshold: {source}"))?;

    if threshold.is_finite() && threshold > 0.0 && threshold <= 1.0 {
        Ok(threshold)
    } else {
        Err("core threshold must be greater than 0 and at most 1".to_string())
    }
}

/// Function to parse command line args into [`MainArgs`] struct
pub fn cli_args() -> Args {
    Args::parse()
}
