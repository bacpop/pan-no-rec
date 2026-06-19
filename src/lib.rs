use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::ThreadPoolBuilder;

use std::io::stdout;
use std::time::Instant;

mod cli;
use crate::cli::cli_args;

mod panaroo_io;

mod output;
use output::{write_paralog_report, write_recombination_table};

mod dists;
use dists::compare_loaded_alignments;

mod gene;
mod genome;
use crate::genome::load_genes;

mod model;

mod graph;
use graph::presence_table_from_pair_hits;

/// Create a progress bar for use in iterators
pub fn get_progress_bar(length: usize, percent: bool, quiet: bool) -> ProgressBar {
    if quiet {
        ProgressBar::hidden()
    } else {
        let style = if percent {
            ProgressStyle::with_template("{percent}% {bar:80.cyan/blue} eta:{eta}").unwrap()
        } else {
            ProgressStyle::with_template("{human_pos}/{human_len} {bar:80.cyan/blue} eta:{eta}")
                .unwrap()
        };
        ProgressBar::new(length as u64).with_style(style)
    }
}

#[doc(hidden)]
pub fn main() -> Result<()> {
    let args = cli_args();
    if args.quiet {
        simple_logger::init_with_level(log::Level::Error).unwrap();
    } else if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
        // simple_logger::init_with_level(log::Level::Trace).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    let print_success = true;
    ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .with_context(|| {
            format!(
                "failed to initialize Rayon global thread pool with {} threads",
                args.threads
            )
        })?;

    let start = Instant::now();

    log::info!("Reading Panaroo input files");
    let (sample_names, genes) = load_genes(
        &args.panaroo_dir,
        args.paralog_mode,
        args.max_entropy,
        args.quiet,
    )?;

    let (n_genes, n_samples) = genes.get_summary()?;
    log::info!("Read {n_genes} alignments and {n_samples} samples");

    let n_paralogs = write_paralog_report(&args.paralog_report, genes.gene_metadata())?;
    if n_paralogs > 1 {
        log::warn!(
            "{} alignments contained paralogs; wrote paralog report to '{}'; using paralog mode '{}'",
            n_paralogs,
            args.paralog_report.display(),
            args.paralog_mode
        );
    }

    log::info!("Running recombination detection: fitting pairwise distance models");
    let gene_hits = compare_loaded_alignments(n_samples, &genes, args.gaps, args.quiet);
    let gene_metadata = genes.into_gene_metadata();

    log::info!("Running recombination detection: using graphs to find genes");
    let rows = presence_table_from_pair_hits(n_samples, n_genes, &gene_hits, args.quiet);

    log::info!("Writing output");
    write_recombination_table(&sample_names, &gene_metadata, &rows, stdout().lock())
        .with_context(|| "failed to write recombination table to stdout")?;
    let end = Instant::now();

    log::info!("Complete");
    if print_success && !args.quiet {
        eprintln!(
            "🦘 pan-no-rec done in {}s",
            end.duration_since(start).as_secs()
        );
    }
    Ok(())
}
