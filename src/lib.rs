use anyhow::{Context, Result, bail};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::ThreadPoolBuilder;

use std::io::stdout;
use std::time::Instant;

mod cli;
use crate::cli::cli_args;

mod panaroo_io;
use panaroo_io::load_genes;

mod output;
use output::{write_paralog_report, write_recombination_table};

mod dists;
use dists::compare_loaded_alignments;

mod gene;
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
    let genes = load_genes(
        &args.panaroo_dir,
        args.paralog_mode,
        args.max_entropy,
        args.quiet,
    )?;

    let sample_names = genes.sample_names;
    let sequences = genes.gene_sequences;
    if sequences.is_empty() {
        bail!("No valid genes loaded");
    } else if sample_names.is_empty() {
        bail!("Alignments are empty");
    } else {
        log::info!(
            "Read {} alignments and {} samples",
            sequences.len(),
            sample_names.len()
        );
    }

    let n_paralogs = write_paralog_report(&args.paralog_report, &sequences)?;
    if n_paralogs > 1 {
        log::warn!(
            "{} alignments contained paralogs; wrote paralog report to '{}'; using paralog mode '{}'",
            n_paralogs,
            args.paralog_report.display(),
            args.paralog_mode
        );
    }

    log::info!("Running recombination detection: fitting pairwise distance models");
    let gene_hits =
        compare_loaded_alignments(sample_names.len(), &sequences, args.gaps, args.quiet);

    log::info!("Running recombination detection: using graphs to find genes");
    let rows =
        presence_table_from_pair_hits(sample_names.len(), &sequences, &gene_hits, args.quiet);

    log::info!("Writing output");
    write_recombination_table(&sample_names, &sequences, &rows, stdout().lock())
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
