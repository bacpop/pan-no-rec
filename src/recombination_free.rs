use crate::blocks::BlockGraph;
use crate::gene::GeneMetadata;
use crate::output::OutputRow;
use crate::panaroo_io::{gene_name_from_path, normalize_sample_id, open_alignment_reader};
use anyhow::{Context, Result, bail};
use hashbrown::{HashMap, HashSet};
use seq_io::fasta::{Reader, Record};
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

const PER_GENE_DIR: &str = "recombination_free_aligned_genes";
const CORE_ALIGNMENT: &str = "recombination_free_core_gene_alignment.aln";
const CORE_HEADER: &str = "recombination_free_core_alignment_header.embl";

pub(crate) fn write_recombination_free_alignments(
    panaroo_dir: &Path,
    output_dir: &Path,
    alignment_paths: &[PathBuf],
    sample_names: &[String],
    genes: &[GeneMetadata],
    rows: &[OutputRow],
    graph: &BlockGraph,
    core_threshold: f64,
) -> Result<()> {
    ensure_output_outside_panaroo(panaroo_dir, output_dir)?;
    if output_dir.exists() && !output_dir.is_dir() {
        bail!(
            "alignment output path '{}' exists and is not a directory",
            output_dir.display()
        );
    }

    let mut row_by_gene = HashMap::with_capacity(rows.len());
    for row in rows {
        let gene = genes.get(row.gene_index).with_context(|| {
            format!(
                "recombination row references missing gene index {}",
                row.gene_index
            )
        })?;
        if row.presence.len() != sample_names.len() {
            bail!(
                "recombination row for '{}' has {} sample values; expected {}",
                gene.name(),
                row.presence.len(),
                sample_names.len()
            );
        }
        if row_by_gene.insert(gene.name(), row).is_some() {
            bail!(
                "multiple recombination rows found for gene '{}'",
                gene.name()
            );
        }
    }

    let mut recombinant_samples = HashMap::with_capacity(rows.len());
    for (&gene_name, row) in &row_by_gene {
        let samples: HashSet<_> = sample_names
            .iter()
            .zip(&row.presence)
            .filter_map(|(sample, &is_recombinant)| (is_recombinant != 0).then_some(sample.clone()))
            .collect();
        recombinant_samples.insert(gene_name.to_owned(), samples);
    }

    let per_gene_dir = output_dir.join(PER_GENE_DIR);
    fs::create_dir_all(&per_gene_dir).with_context(|| {
        format!(
            "failed to create recombination-free alignment directory '{}'",
            per_gene_dir.display()
        )
    })?;

    let mut gene_names_seen = HashSet::with_capacity(alignment_paths.len());
    for alignment_path in alignment_paths {
        let gene_name = gene_name_from_path(alignment_path)?;
        if !gene_names_seen.insert(gene_name.clone()) {
            bail!("multiple input alignments resolve to gene name '{gene_name}'");
        }
        if !row_by_gene.contains_key(gene_name.as_str()) {
            bail!("no recombination row found for input gene '{gene_name}'");
        }
        let filename = alignment_path.file_name().with_context(|| {
            format!(
                "alignment path has no filename: '{}'",
                alignment_path.display()
            )
        })?;
        let output_path = per_gene_dir.join(filename);
        let removed = recombinant_samples
            .get(gene_name.as_str())
            .is_some_and(|samples| !samples.is_empty());
        if removed {
            write_filtered_alignment(
                alignment_path,
                &output_path,
                &recombinant_samples[gene_name.as_str()],
            )?;
        } else {
            fs::copy(alignment_path, &output_path).with_context(|| {
                format!(
                    "failed to copy alignment '{}' to '{}'",
                    alignment_path.display(),
                    output_path.display()
                )
            })?;
        }
    }

    let mut core_gene_names = HashSet::new();
    for node in &graph.nodes {
        let size = graph
            .node_size(node.id)?
            .with_context(|| format!("Panaroo graph node {} has no size", node.id))?;
        if (size as f64) / (sample_names.len() as f64) >= core_threshold {
            let name = graph
                .node_name(node.id)
                .with_context(|| format!("Panaroo graph node {} has no gene name", node.id))?;
            core_gene_names.insert(name.to_owned());
        }
    }

    let mut path_by_name = HashMap::new();
    for path in alignment_paths {
        let name = gene_name_from_path(path)?;
        let filename = path
            .file_name()
            .with_context(|| format!("alignment path has no filename: '{}'", path.display()))?;
        let written = per_gene_dir.join(filename);
        if core_gene_names.contains(&name) && written.is_file() {
            path_by_name.insert(name, written);
        }
    }
    let mut core_alignments = Vec::new();
    for row in rows {
        let gene = genes
            .get(row.gene_index)
            .context("core alignment row references missing gene")?;
        if let Some(path) = path_by_name.remove(gene.name()) {
            core_alignments.push((gene.name().to_owned(), path));
        }
    }

    write_core_alignment(output_dir, &core_alignments)?;
    log::info!(
        "Wrote recombination-free gene alignments and core alignment to '{}'",
        output_dir.display()
    );
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_combined_recombination_free_alignments(
    panaroo_dir: &Path,
    igr_dir: &Path,
    output_dir: &Path,
    gene_paths: &[PathBuf],
    igr_paths: &[PathBuf],
    sample_names: &[String],
    genes: &[GeneMetadata],
    gene_rows: &[OutputRow],
    igrs: &[GeneMetadata],
    igr_rows: &[OutputRow],
    combined: &[GeneMetadata],
    combined_rows: &[OutputRow],
    gene_graph: &mut BlockGraph,
    igr_graph: &mut BlockGraph,
    core_threshold: f64,
) -> Result<()> {
    ensure_output_outside_panaroo(panaroo_dir, output_dir)?;
    ensure_output_outside_panaroo(igr_dir, output_dir)?;
    if output_dir.exists() && !output_dir.is_dir() {
        bail!(
            "alignment output path '{}' exists and is not a directory",
            output_dir.display()
        );
    }
    let gene_output = output_dir.join("recombination_free_aligned_genes");
    let igr_output = output_dir.join("recombination_free_aligned_IGRs");
    write_filtered_feature_set(&gene_output, gene_paths, sample_names, genes, gene_rows)?;
    write_filtered_feature_set(&igr_output, igr_paths, sample_names, igrs, igr_rows)?;

    let gene_core = core_feature_paths(
        gene_graph,
        gene_paths,
        &gene_output,
        sample_names.len(),
        core_threshold,
    )?;
    let igr_core = core_feature_paths(
        igr_graph,
        igr_paths,
        &igr_output,
        sample_names.len(),
        core_threshold,
    )?;
    write_core_alignment_to(output_dir, &gene_core, CORE_ALIGNMENT, CORE_HEADER)?;
    write_core_alignment_to(
        output_dir,
        &igr_core,
        "recombination_free_core_IGR_alignment.aln",
        "recombination_free_core_IGR_alignment_header.embl",
    )?;
    let mut available_paths = HashMap::new();
    for (paths, out_dir) in [(gene_paths, &gene_output), (igr_paths, &igr_output)] {
        for path in paths {
            let name = gene_name_from_path(path)?;
            let filename = path.file_name().context("alignment path has no filename")?;
            let written = out_dir.join(filename);
            if written.is_file() {
                available_paths.insert(name, written);
            }
        }
    }
    let mut joint_sizes = HashMap::new();
    for node in &igr_graph.nodes {
        if let (Some(name), Some(size)) =
            (igr_graph.node_name(node.id), igr_graph.node_size(node.id)?)
        {
            joint_sizes.insert(name.to_owned(), size);
        }
    }
    let mut joint_core = Vec::new();
    for row in combined_rows {
        let feature = combined
            .get(row.gene_index)
            .context("joint recombination row references missing feature")?;
        let Some(&size) = joint_sizes.get(feature.name()) else {
            continue;
        };
        if (size as f64) / (sample_names.len() as f64) < core_threshold {
            continue;
        }
        if let Some(path) = available_paths.remove(feature.name()) {
            joint_core.push((feature.name().to_owned(), path));
        }
    }
    write_core_alignment_to(
        output_dir,
        &joint_core,
        "recombination_free_core_genome.aln.fas",
        "recombination_free_core_genome_header.embl",
    )?;
    log::info!(
        "Wrote recombination-free gene, IGR, and core genome alignments to '{}'",
        output_dir.display()
    );
    Ok(())
}

fn core_feature_paths(
    graph: &mut BlockGraph,
    paths: &[PathBuf],
    output_dir: &Path,
    sample_count: usize,
    threshold: f64,
) -> Result<Vec<(String, PathBuf)>> {
    let order = graph.layout_core_genes(sample_count)?;
    let mut path_by_name = HashMap::new();
    for path in paths {
        let name = gene_name_from_path(path)?;
        let filename = path.file_name().context("alignment path has no filename")?;
        let written = output_dir.join(filename);
        if written.is_file() {
            path_by_name.insert(name, written);
        }
    }
    let mut core = Vec::new();
    for node in order {
        let Some(size) = graph.node_size(node)? else {
            continue;
        };
        if (size as f64) / (sample_count as f64) < threshold {
            continue;
        }
        let name = graph
            .node_name(node)
            .context("graph node has no feature name")?;
        if let Some(path) = path_by_name.remove(name) {
            core.push((name.to_owned(), path));
        }
    }
    Ok(core)
}

fn write_filtered_feature_set(
    output_dir: &Path,
    alignment_paths: &[PathBuf],
    sample_names: &[String],
    features: &[GeneMetadata],
    rows: &[OutputRow],
) -> Result<()> {
    let mut rows_by_name = HashMap::with_capacity(rows.len());
    for row in rows {
        let feature = features
            .get(row.gene_index)
            .context("output row references missing feature")?;
        if row.presence.len() != sample_names.len() {
            bail!(
                "recombination row for '{}' has an unexpected number of samples",
                feature.name()
            );
        }
        if rows_by_name.insert(feature.name(), row).is_some() {
            bail!(
                "multiple recombination rows found for feature '{}'",
                feature.name()
            );
        }
    }
    fs::create_dir_all(output_dir).with_context(|| {
        format!(
            "failed to create alignment directory '{}'",
            output_dir.display()
        )
    })?;
    let mut seen = HashSet::new();
    for path in alignment_paths {
        let name = gene_name_from_path(path)?;
        if !seen.insert(name.clone()) {
            bail!("multiple input alignments resolve to feature '{name}'");
        }
        let row = rows_by_name
            .get(name.as_str())
            .with_context(|| format!("no recombination row found for feature '{name}'"))?;
        let recombinant_samples: HashSet<_> = sample_names
            .iter()
            .zip(&row.presence)
            .filter_map(|(sample, &is_recombinant)| (is_recombinant != 0).then_some(sample.clone()))
            .collect();
        let filename = path.file_name().context("alignment path has no filename")?;
        let output_path = output_dir.join(filename);
        if recombinant_samples.is_empty() {
            fs::copy(path, &output_path)
                .with_context(|| format!("failed to copy '{}'", path.display()))?;
        } else {
            write_filtered_alignment(path, &output_path, &recombinant_samples)?;
        }
    }
    Ok(())
}

fn write_filtered_alignment(
    input_path: &Path,
    output_path: &Path,
    recombinant_samples: &HashSet<String>,
) -> Result<()> {
    let temporary_path = temporary_sibling(output_path);
    let input = open_alignment_reader(input_path)?;
    let mut reader = Reader::new(input);
    let mut kept_count = 0usize;
    {
        let file = File::create(&temporary_path).with_context(|| {
            format!(
                "failed to create temporary alignment '{}'",
                temporary_path.display()
            )
        })?;
        let mut writer = BufWriter::new(file);
        while let Some(record) = reader.next() {
            let record = record
                .with_context(|| format!("failed to parse alignment '{}'", input_path.display()))?;
            let id = record.id().with_context(|| {
                format!(
                    "alignment '{}' contains a non-UTF-8 FASTA identifier",
                    input_path.display()
                )
            })?;
            if recombinant_samples.contains(normalize_sample_id(id)) {
                continue;
            }
            record.write(&mut writer).with_context(|| {
                format!(
                    "failed to write filtered alignment '{}'",
                    temporary_path.display()
                )
            })?;
            kept_count += 1;
        }
        writer.flush().with_context(|| {
            format!(
                "failed to flush filtered alignment '{}'",
                temporary_path.display()
            )
        })?;
    }

    if kept_count == 0 {
        if temporary_path.exists() {
            fs::remove_file(&temporary_path)?;
        }
        if output_path.exists() {
            fs::remove_file(output_path).with_context(|| {
                format!(
                    "failed to remove empty alignment '{}'",
                    output_path.display()
                )
            })?;
        }
    } else {
        fs::rename(&temporary_path, output_path).with_context(|| {
            format!(
                "failed to move filtered alignment '{}' to '{}'",
                temporary_path.display(),
                output_path.display()
            )
        })?;
    }
    Ok(())
}

#[derive(Debug)]
struct GeneSequences {
    length: usize,
    by_sample: HashMap<String, Vec<u8>>,
}

fn read_best_gene_sequences(path: &Path) -> Result<GeneSequences> {
    let input = open_alignment_reader(path)?;
    let mut reader = Reader::new(input);
    let mut alignment_length = None;
    let mut by_sample: HashMap<String, (Vec<u8>, usize)> = HashMap::new();

    while let Some(record) = reader.next() {
        let record =
            record.with_context(|| format!("failed to parse alignment '{}'", path.display()))?;
        let id = record.id().with_context(|| {
            format!(
                "alignment '{}' contains a non-UTF-8 FASTA identifier",
                path.display()
            )
        })?;
        let id = id.strip_prefix("_R_").unwrap_or(id);
        let sample = normalize_sample_id(id).to_owned();
        let sequence = record.full_seq().into_owned();
        if let Some(expected) = alignment_length {
            if sequence.len() != expected {
                bail!(
                    "alignment '{}' has variable sequence lengths: found {}, expected {}",
                    path.display(),
                    sequence.len(),
                    expected
                );
            }
        } else {
            alignment_length = Some(sequence.len());
        }
        let gaps = sequence.iter().filter(|&&base| base == b'-').count();
        match by_sample.entry(sample) {
            hashbrown::hash_map::Entry::Vacant(entry) => {
                entry.insert((sequence, gaps));
            }
            hashbrown::hash_map::Entry::Occupied(mut entry) => {
                if gaps < entry.get().1 {
                    entry.insert((sequence, gaps));
                }
            }
        }
    }

    let length = alignment_length
        .with_context(|| format!("alignment '{}' contains no FASTA records", path.display()))?;
    Ok(GeneSequences {
        length,
        by_sample: by_sample
            .into_iter()
            .map(|(sample, (sequence, _))| (sample, sequence))
            .collect(),
    })
}

fn write_core_alignment(output_dir: &Path, core_alignments: &[(String, PathBuf)]) -> Result<()> {
    write_core_alignment_to(output_dir, core_alignments, CORE_ALIGNMENT, CORE_HEADER)
}

fn write_core_alignment_to(
    output_dir: &Path,
    core_alignments: &[(String, PathBuf)],
    alignment_filename: &str,
    header_filename: &str,
) -> Result<()> {
    let mut samples = BTreeSet::new();
    let mut gene_lengths = Vec::with_capacity(core_alignments.len());
    let mut total_length = 0usize;
    let mut retained_genes = Vec::new();

    // First pass determines the sorted sample set and concatenated dimensions.
    for (gene_name, path) in core_alignments {
        let records = read_best_gene_sequences(path)?;
        total_length = total_length
            .checked_add(records.length)
            .context("core alignment length overflows usize")?;
        samples.extend(records.by_sample.keys().cloned());
        gene_lengths.push(records.length);
        retained_genes.push((gene_name.as_str(), path.as_path()));
    }
    let samples: Vec<_> = samples.into_iter().collect();
    let sample_indices: HashMap<_, _> = samples
        .iter()
        .enumerate()
        .map(|(index, sample)| (sample.as_str(), index))
        .collect();
    let total_length_u64 = u64::try_from(total_length).context("core alignment is too large")?;
    let spool_bytes = total_length_u64
        .checked_mul(samples.len() as u64)
        .context("temporary core alignment size overflows u64")?;
    let spool_path = output_dir.join(format!(".{alignment_filename}.tmp"));
    let mut spool = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open(&spool_path)
        .with_context(|| {
            format!(
                "failed to create temporary core alignment '{}'",
                spool_path.display()
            )
        })?;
    write_gap_filled_spool(&mut spool, spool_bytes)?;

    // Second pass writes one gene at a time into disk-backed sample rows.
    let mut gene_offset = 0u64;
    for ((_, path), &gene_length) in retained_genes.iter().zip(&gene_lengths) {
        let records = read_best_gene_sequences(path)?;
        if records.length != gene_length {
            bail!(
                "alignment '{}' changed length while writing core alignment",
                path.display()
            );
        }
        for (sample, sequence) in records.by_sample {
            let sample_index = *sample_indices.get(sample.as_str()).with_context(|| {
                format!("sample '{sample}' disappeared from core alignment sample set")
            })? as u64;
            let sample_offset = sample_index
                .checked_mul(total_length_u64)
                .and_then(|offset| offset.checked_add(gene_offset))
                .context("core alignment file offset overflows u64")?;
            spool.seek(SeekFrom::Start(sample_offset))?;
            spool.write_all(&sequence)?;
        }
        gene_offset = gene_offset
            .checked_add(gene_length as u64)
            .context("core alignment file offset overflows u64")?;
    }
    spool.flush()?;

    let core_path = output_dir.join(alignment_filename);
    write_sample_major_fasta(&mut spool, &samples, total_length, &core_path)?;
    drop(spool);
    fs::remove_file(&spool_path).with_context(|| {
        format!(
            "failed to remove temporary core alignment '{}'",
            spool_path.display()
        )
    })?;
    write_core_header(
        output_dir.join(header_filename),
        &retained_genes,
        &gene_lengths,
    )?;
    Ok(())
}

fn write_gap_filled_spool(spool: &mut File, bytes: u64) -> Result<()> {
    spool.set_len(bytes)?;
    spool.seek(SeekFrom::Start(0))?;
    let gaps = [b'-'; 64 * 1024];
    let mut remaining = bytes;
    while remaining > 0 {
        let count = remaining.min(gaps.len() as u64) as usize;
        spool.write_all(&gaps[..count])?;
        remaining -= count as u64;
    }
    Ok(())
}

fn write_sample_major_fasta(
    spool: &mut File,
    samples: &[String],
    sequence_length: usize,
    output_path: &Path,
) -> Result<()> {
    let mut output = BufWriter::new(File::create(output_path).with_context(|| {
        format!(
            "failed to create core alignment '{}'",
            output_path.display()
        )
    })?);
    let mut buffer = [0u8; 60];
    let sequence_length_u64 =
        u64::try_from(sequence_length).context("core alignment is too large")?;
    for (sample_index, sample) in samples.iter().enumerate() {
        writeln!(output, ">{sample}")?;
        let row_offset = (sample_index as u64)
            .checked_mul(sequence_length_u64)
            .context("core alignment file offset overflows u64")?;
        spool.seek(SeekFrom::Start(row_offset))?;
        let mut remaining = sequence_length;
        while remaining > 0 {
            let count = remaining.min(buffer.len());
            spool.read_exact(&mut buffer[..count])?;
            output.write_all(&buffer[..count])?;
            output.write_all(b"\n")?;
            remaining -= count;
        }
    }
    output
        .flush()
        .with_context(|| format!("failed to flush core alignment '{}'", output_path.display()))?;
    Ok(())
}

fn write_core_header(path: PathBuf, genes: &[(&str, &Path)], lengths: &[usize]) -> Result<()> {
    let mut output = BufWriter::new(File::create(&path).with_context(|| {
        format!(
            "failed to create core alignment header '{}'",
            path.display()
        )
    })?);
    writeln!(output, "ID   Genome standard; DNA; PRO; 1234 BP.")?;
    writeln!(output, "XX")?;
    writeln!(output, "FH   Key             Location/Qualifiers")?;
    writeln!(output, "FH")?;

    let mut gene_start = 1usize;
    let mut gene_end = 0usize;
    for ((_, path), &length) in genes.iter().zip(lengths) {
        gene_end = gene_end
            .checked_add(length)
            .context("core feature coordinate overflows usize")?;
        let filename = path
            .file_name()
            .and_then(|name| name.to_str())
            .with_context(|| {
                format!(
                    "core alignment path has invalid filename '{}'",
                    path.display()
                )
            })?;
        // Match the legacy script's splitext(filename)[0] feature label.
        let label = filename.strip_suffix(".fas").unwrap_or(filename);
        writeln!(output, "FT   feature         {gene_start}..{gene_end}")?;
        writeln!(output, "FT                   /label={label}")?;
        writeln!(output, "FT                   /locus_tag={label}")?;
        gene_start = gene_start
            .checked_add(length)
            .context("core feature coordinate overflows usize")?;
    }

    writeln!(output, "XX")?;
    writeln!(
        output,
        "SQ   Sequence 1234 BP; 789 A; 1717 C; 1693 G; 691 T; 0 other;"
    )?;
    writeln!(output, "//")?;
    output
        .flush()
        .with_context(|| format!("failed to flush core alignment header '{}'", path.display()))?;
    Ok(())
}

fn temporary_sibling(path: &Path) -> PathBuf {
    let mut filename = OsString::from(".");
    filename.push(
        path.file_name()
            .unwrap_or_else(|| std::ffi::OsStr::new("alignment")),
    );
    filename.push(".tmp");
    path.with_file_name(filename)
}

fn ensure_output_outside_panaroo(panaroo_dir: &Path, output_dir: &Path) -> Result<()> {
    let input = fs::canonicalize(panaroo_dir).with_context(|| {
        format!(
            "failed to resolve Panaroo directory '{}'",
            panaroo_dir.display()
        )
    })?;
    let output = resolve_path(output_dir)?;
    if output.starts_with(&input) {
        bail!(
            "alignment output directory '{}' must not be inside the Panaroo directory '{}'",
            output_dir.display(),
            panaroo_dir.display()
        );
    }
    Ok(())
}

fn resolve_path(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir => {
                normalized.pop();
            }
            other => normalized.push(other.as_os_str()),
        }
    }

    let mut ancestor = normalized.as_path();
    let mut missing = Vec::new();
    while !ancestor.exists() {
        let name = ancestor
            .file_name()
            .with_context(|| format!("cannot resolve output path '{}'", path.display()))?;
        missing.push(name.to_os_string());
        ancestor = ancestor
            .parent()
            .with_context(|| format!("cannot resolve output path '{}'", path.display()))?;
    }
    let mut resolved = fs::canonicalize(ancestor)?;
    for component in missing.iter().rev() {
        resolved.push(component);
    }
    Ok(resolved)
}
