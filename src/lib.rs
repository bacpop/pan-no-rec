use anyhow::{Context, Result, bail};
use hashbrown::HashMap;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::ThreadPoolBuilder;

use std::fs::File;
use std::io::{BufWriter, stdout};
use std::path::{Path, PathBuf};
use std::time::Instant;

mod cli;
use crate::cli::cli_args;

mod panaroo_io;
mod recombination_free;

mod layout_segments;
mod output;
use output::{order_rows_by_gene_names, write_paralog_report, write_recombination_table};

mod dists;
use dists::compare_loaded_alignments;

mod gene;
mod genome;
use crate::genome::load_genes;
use crate::genome::load_igrs;

mod model;

mod graph;
use graph::{GraphDebugConfig, presence_table_from_pair_hits};

pub mod blocks;
mod candidate;
pub mod plot;

fn build_graph_debug_config(args: &cli::Args) -> Result<Option<GraphDebugConfig>> {
    match (&args.debug_graph_dir, args.debug_graph_genes.is_empty()) {
        (Some(output_dir), false) => Ok(Some(GraphDebugConfig::new(
            output_dir.clone(),
            args.debug_graph_genes.clone(),
        )?)),
        (Some(_), true) => {
            bail!("--debug-graph-dir requires at least one --debug-graph-genes value")
        }
        (None, false) => bail!("--debug-graph-genes requires --debug-graph-dir"),
        (None, true) => Ok(None),
    }
}

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

fn graph_order_and_sizes(
    graph: &mut blocks::BlockGraph,
    sample_count: usize,
) -> Result<(Vec<String>, HashMap<String, usize>, blocks::LayoutResult)> {
    let layout = graph.layout_with_metadata(sample_count)?;
    let order = &layout.node_order;
    let names = order
        .iter()
        .map(|&node_id| {
            graph
                .node_name(node_id)
                .map(str::to_owned)
                .with_context(|| format!("Panaroo graph node {node_id} has no feature name"))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut sizes = HashMap::with_capacity(order.len());
    for &node_id in order {
        let name = graph
            .node_name(node_id)
            .with_context(|| format!("Panaroo graph node {node_id} has no feature name"))?;
        if let Some(size) = graph.node_size(node_id)? {
            // Keep emitted frequency rows valid even if paralogs inflate the GML size.
            sizes.insert(name.to_owned(), size.min(sample_count));
        }
    }
    Ok((names, sizes, layout))
}

fn export_layout_segments(
    path: &Path,
    graph: &blocks::BlockGraph,
    layout: &blocks::LayoutResult,
    genes: &[gene::GeneMetadata],
    rows: &[output::OutputRow],
) -> Result<()> {
    let names = rows
        .iter()
        .map(|row| {
            genes
                .get(row.gene_index)
                .map(|gene| gene.name().to_owned())
                .context("segment export references missing feature")
        })
        .collect::<Result<Vec<_>>>()?;
    layout_segments::write(path, graph, layout, &names)
}

fn prefixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut value = prefix.as_os_str().to_os_string();
    value.push(suffix);
    value.into()
}

fn write_table_file(
    path: &Path,
    sample_names: &[String],
    genes: &[gene::GeneMetadata],
    rows: &[output::OutputRow],
    sizes: &HashMap<String, usize>,
) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("failed to create TSV '{}'", path.display()))?;
    write_recombination_table(sample_names, genes, rows, Some(sizes), BufWriter::new(file))
        .with_context(|| format!("failed to write TSV '{}'", path.display()))
}

#[doc(hidden)]
pub fn main() -> Result<()> {
    let args = cli_args();
    let panaroo_graph_path = args.panaroo_dir.join("final_graph.gml");
    if args.layout_segments_output.is_some() {
        let graph = args
            .igrs
            .as_ref()
            .unwrap_or(&args.panaroo_dir)
            .join("final_graph.gml");
        if !graph.is_file() {
            bail!(
                "--layout-segments-output requires graph '{}'",
                graph.display()
            );
        }
    }
    if args.output.is_some() && args.igrs.is_none() {
        bail!("--output is only used with --igrs");
    }
    if args.igrs.is_some() && args.output.is_none() {
        bail!("--igrs requires --output PREFIX");
    }
    if args.alignment_output.is_some() && !panaroo_graph_path.is_file() {
        bail!(
            "--alignment-output requires Panaroo graph '{}'",
            panaroo_graph_path.display()
        );
    }
    if args.quiet {
        simple_logger::init_with_level(log::Level::Error).unwrap();
    } else if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
        // simple_logger::init_with_level(log::Level::Trace).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    let print_success = true;
    let graph_debug_config = build_graph_debug_config(&args)?;

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

    log::info!("Running gene recombination detection");
    let gene_hits = compare_loaded_alignments(n_samples, &genes, args.gaps, args.quiet);
    let gene_metadata = genes.into_gene_metadata();
    let gene_rows = presence_table_from_pair_hits(
        &sample_names,
        &gene_metadata,
        &gene_hits,
        graph_debug_config.as_ref(),
        args.quiet,
    )?;
    drop(gene_hits);

    if let (Some(igr_dir), Some(prefix)) = (args.igrs.as_deref(), args.output.as_deref()) {
        let igr_rtab = igr_dir.join("IGR_presence_absence.Rtab");
        let igr_samples_from_matrix = panaroo_io::read_rtab_sample_names(&igr_rtab)?;
        let gene_set: std::collections::HashSet<_> = sample_names.iter().collect();
        let igr_set: std::collections::HashSet<_> = igr_samples_from_matrix.iter().collect();
        if gene_set != igr_set {
            bail!("Panaroo and panIGRoo isolate names differ in their presence/absence matrices");
        }
        if !panaroo_graph_path.is_file() {
            bail!(
                "--igrs requires Panaroo graph '{}'",
                panaroo_graph_path.display()
            );
        }
        let igr_graph_path = igr_dir.join("final_graph.gml");
        if !igr_graph_path.is_file() {
            bail!(
                "--igrs requires panIGRoo graph '{}'",
                igr_graph_path.display()
            );
        }

        log::info!("Loading IGR alignments longer than 30 bp");
        let (igr_sample_names, igr_genome) = load_igrs(igr_dir, args.max_entropy, args.quiet)?;
        let igr_hits = compare_loaded_alignments(n_samples, &igr_genome, args.gaps, args.quiet);
        let igr_metadata = igr_genome.into_gene_metadata();
        let mut igr_rows = presence_table_from_pair_hits(
            &igr_sample_names,
            &igr_metadata,
            &igr_hits,
            None,
            args.quiet,
        )?;
        drop(igr_hits);
        if igr_sample_names != sample_names {
            let source_indices: HashMap<_, _> = igr_sample_names
                .iter()
                .enumerate()
                .map(|(index, name)| (name.as_str(), index))
                .collect();
            for row in &mut igr_rows {
                row.presence = sample_names
                    .iter()
                    .map(|sample| row.presence[source_indices[sample.as_str()]])
                    .collect();
            }
        }

        let mut panaroo_graph = blocks::load_gml(&panaroo_graph_path)?;
        let mut igr_graph = blocks::load_gml(&igr_graph_path)?;
        let (gene_order, gene_sizes, _) = graph_order_and_sizes(&mut panaroo_graph, n_samples)?;
        let (igr_order, igr_sizes, igr_layout) = graph_order_and_sizes(&mut igr_graph, n_samples)?;
        let ordered_gene_rows = order_rows_by_gene_names(&gene_rows, &gene_metadata, &gene_order)?;
        let ordered_igr_rows = order_rows_by_gene_names(&igr_rows, &igr_metadata, &igr_order)?;

        let mut combined_metadata = gene_metadata.clone();
        let gene_name_set: std::collections::HashSet<_> =
            gene_metadata.iter().map(|gene| gene.name()).collect();
        for gene in &igr_metadata {
            if gene_name_set.contains(gene.name()) {
                bail!("gene and IGR feature names collide at '{}'", gene.name());
            }
            combined_metadata.push(gene.clone());
        }
        let mut combined_rows = gene_rows.clone();
        combined_rows.extend(igr_rows.iter().map(|row| output::OutputRow {
            gene_index: row.gene_index + gene_metadata.len(),
            presence: row.presence.clone(),
        }));
        let combined_order = igr_order;
        let combined_rows =
            order_rows_by_gene_names(&combined_rows, &combined_metadata, &combined_order)?;

        if let Some(path) = args.layout_segments_output.as_deref() {
            export_layout_segments(
                path,
                &igr_graph,
                &igr_layout,
                &combined_metadata,
                &combined_rows,
            )?;
        }
        let genes_path = prefixed_path(prefix, ".genes.tsv");
        let igrs_path = prefixed_path(prefix, ".igrs.tsv");
        let joint_path = prefixed_path(prefix, ".tsv");
        write_table_file(
            &genes_path,
            &sample_names,
            &gene_metadata,
            &ordered_gene_rows,
            &gene_sizes,
        )?;
        write_table_file(
            &igrs_path,
            &sample_names,
            &igr_metadata,
            &ordered_igr_rows,
            &igr_sizes,
        )?;
        write_table_file(
            &joint_path,
            &sample_names,
            &combined_metadata,
            &combined_rows,
            &igr_sizes,
        )?;

        if let Some(output_dir) = args.alignment_output.as_deref() {
            let gene_names: std::collections::HashSet<_> =
                gene_metadata.iter().map(|g| g.name()).collect();
            let igr_names: std::collections::HashSet<_> =
                igr_metadata.iter().map(|g| g.name()).collect();
            let gene_paths: Vec<_> = panaroo_io::read_panaroo_dir(&args.panaroo_dir)?
                .into_iter()
                .filter_map(|path| match panaroo_io::gene_name_from_path(&path) {
                    Ok(name) if gene_names.contains(name.as_str()) => Some(Ok(path)),
                    Ok(_) => None,
                    Err(error) => Some(Err(error)),
                })
                .collect::<Result<_>>()?;
            let igr_paths: Vec<_> = panaroo_io::read_igr_dir(igr_dir)?
                .into_iter()
                .filter_map(|path| match panaroo_io::gene_name_from_path(&path) {
                    Ok(name) if igr_names.contains(name.as_str()) => Some(Ok(path)),
                    Ok(_) => None,
                    Err(error) => Some(Err(error)),
                })
                .collect::<Result<_>>()?;
            recombination_free::write_combined_recombination_free_alignments(
                &args.panaroo_dir,
                igr_dir,
                output_dir,
                &gene_paths,
                &igr_paths,
                &sample_names,
                &gene_metadata,
                &gene_rows,
                &igr_metadata,
                &igr_rows,
                &combined_metadata,
                &combined_rows,
                &mut panaroo_graph,
                &mut igr_graph,
                args.core_threshold,
            )?;
        }
        log::info!(
            "Wrote gene, IGR, and joint TSV outputs using prefix '{}'…",
            prefix.display()
        );
    } else {
        let mut gene_rows = gene_rows;
        let mut gene_sizes = None;
        if panaroo_graph_path.is_file() {
            let mut graph = blocks::load_gml(&panaroo_graph_path)?;
            let (order, sizes, layout) = graph_order_and_sizes(&mut graph, n_samples)?;
            gene_rows = order_rows_by_gene_names(&gene_rows, &gene_metadata, &order)?;
            gene_sizes = Some(sizes);
            if let Some(path) = args.layout_segments_output.as_deref() {
                export_layout_segments(path, &graph, &layout, &gene_metadata, &gene_rows)?;
            }
            if let Some(output_dir) = args.alignment_output.as_deref() {
                let alignment_paths = panaroo_io::read_panaroo_dir(&args.panaroo_dir)?;
                recombination_free::write_recombination_free_alignments(
                    &args.panaroo_dir,
                    output_dir,
                    &alignment_paths,
                    &sample_names,
                    &gene_metadata,
                    &gene_rows,
                    &graph,
                    args.core_threshold,
                )?;
            }
        } else {
            log::warn!(
                "Panaroo graph '{}' not found; preserving alignment order",
                panaroo_graph_path.display()
            );
        }
        write_recombination_table(
            &sample_names,
            &gene_metadata,
            &gene_rows,
            gene_sizes.as_ref(),
            stdout().lock(),
        )
        .with_context(|| "failed to write recombination table to stdout")?;
    }
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
