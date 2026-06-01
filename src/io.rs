use crate::cli::ParalogMode;
use crate::gene::Gene;
use crate::get_progress_bar;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use std::collections::{HashMap, HashSet, hash_map::Entry};
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub(crate) struct LoadedAlignments {
    pub(crate) sample_names: Vec<String>,
    pub(crate) genes: Vec<Gene>,
    pub(crate) paralogs: Vec<(String, usize)>,
}

#[derive(Debug)]
struct RawAlignment {
    gene_name: String,
    sample_names: Vec<String>,
    sequences: Vec<Vec<u8>>,
    alignment_len: usize,
    paralog: Option<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RecombinationRow {
    pub(crate) gene_index: usize,
    pub(crate) presence: Vec<u8>,
}

#[derive(Debug)]
struct RawSampleGroup {
    sample_name: String,
    records: Vec<RawSampleRecord>,
}

#[derive(Debug)]
struct RawSampleRecord {
    sequence: Vec<u8>,
    non_gap_count: usize,
}

// Reads an MSA list file and resolves its entries to paths.
pub fn read_msa_list(path: &Path) -> Result<Vec<PathBuf>> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read MSA list '{}'", path.display()))?;
    Ok(parse_msa_list(path, &contents))
}

// Finds Panaroo gene alignment files in a directory or its standard fallback.
pub fn read_panaroo_dir(path: &Path) -> Result<Vec<PathBuf>> {
    let mut paths = find_panaroo_alignments(path)?;
    if paths.is_empty() {
        let fallback = path.join("aligned_gene_sequences");
        paths = find_optional_panaroo_alignments(&fallback)?;
    }

    if paths.is_empty() {
        bail!(
            "found no Panaroo alignment files ending in .aln.fas in '{}' or '{}'",
            path.display(),
            path.join("aligned_gene_sequences").display()
        );
    }

    Ok(paths)
}

// Collects top-level Panaroo alignment files from a readable directory.
fn find_panaroo_alignments(dir: &Path) -> Result<Vec<PathBuf>> {
    collect_panaroo_alignments(dir, true)
}

// Collects Panaroo alignment files from an optional fallback directory.
fn find_optional_panaroo_alignments(dir: &Path) -> Result<Vec<PathBuf>> {
    collect_panaroo_alignments(dir, false)
}

fn collect_panaroo_alignments(dir: &Path, required: bool) -> Result<Vec<PathBuf>> {
    let entries = match fs::read_dir(dir) {
        Ok(entries) => entries,
        Err(error) if !required && error.kind() == std::io::ErrorKind::NotFound => {
            return Ok(Vec::new());
        }
        Err(error) => {
            return Err(error)
                .with_context(|| format!("failed to read Panaroo directory '{}'", dir.display()));
        }
    };

    let mut paths = Vec::new();
    for entry in entries {
        let entry = entry.with_context(|| {
            format!(
                "failed to read an entry in Panaroo directory '{}'",
                dir.display()
            )
        })?;
        let file_type = entry.file_type().with_context(|| {
            format!(
                "failed to read file type for Panaroo path '{}'",
                entry.path().display()
            )
        })?;
        if file_type.is_file()
            && entry
                .file_name()
                .to_str()
                .is_some_and(|name| name.ends_with(".aln.fas"))
        {
            paths.push(entry.path());
        }
    }

    paths.sort_unstable();
    Ok(paths)
}

// Writes a recombination presence table as tab-separated text.
pub fn write_recombination_table<W: Write>(
    sample_names: &[String],
    genes: &[Gene],
    rows: &[RecombinationRow],
    mut writer: W,
) -> Result<()> {
    write!(writer, "gene")?;
    for sample in sample_names {
        write!(writer, "\t{sample}")?;
    }
    writeln!(writer)?;

    for row in rows {
        let gene = genes.get(row.gene_index).with_context(|| {
            format!(
                "recombination table row references missing gene index {}",
                row.gene_index
            )
        })?;
        write!(writer, "{}", gene.name())?;
        for value in &row.presence {
            write!(writer, "\t{value}")?;
        }
        writeln!(writer)?;
    }

    Ok(())
}

// Writes per-gene paralog metadata as tab-separated text.
pub(crate) fn write_paralog_report(path: &Path, rows: &[(String, usize)]) -> Result<()> {
    let mut writer = File::create(path)
        .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;

    writeln!(writer, "gene\tparalog_samples")
        .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    for (gene, paralog_samples) in rows {
        writeln!(writer, "{gene}\t{paralog_samples}")
            .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    }

    Ok(())
}

// Parses list-file contents, ignoring blanks and comments.
fn parse_msa_list(list_path: &Path, contents: &str) -> Vec<PathBuf> {
    let base_dir = list_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    contents
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let path = PathBuf::from(line);
            if path.is_absolute() {
                path
            } else {
                base_dir.join(path)
            }
        })
        .collect()
}

// Loads all input alignments into genes and sorted sample names.
pub(crate) fn load_genes<P>(
    aln_paths: &[P],
    paralog_mode: ParalogMode,
    quiet: bool,
) -> Result<LoadedAlignments>
where
    P: AsRef<Path> + Sync,
{
    let mut all_samples = HashSet::new();

    let pbar = get_progress_bar(aln_paths.len(), false, quiet);
    let raw_alignments: Vec<RawAlignment> = aln_paths
        .par_iter()
        .progress_with(pbar)
        .map(|aln| parse_raw_alignment(aln.as_ref(), paralog_mode))
        .collect::<Result<Vec<_>>>()?;

    for raw in &raw_alignments {
        all_samples.extend(raw.sample_names.iter().cloned());
    }

    let mut sample_names: Vec<_> = all_samples.into_iter().collect();
    sample_names.sort_unstable();
    let sample_indices: HashMap<_, _> = sample_names
        .iter()
        .enumerate()
        .map(|(index, sample)| (sample.as_str(), index))
        .collect();

    let mut genes = Vec::with_capacity(raw_alignments.len());
    let mut paralogs = Vec::new();
    for raw in raw_alignments {
        if let Some(paralog_samples) = raw.paralog {
            paralogs.push((raw.gene_name.clone(), paralog_samples));
        }
        genes.push(build_gene(raw, &sample_indices));
    }

    Ok(LoadedAlignments {
        sample_names,
        genes,
        paralogs,
    })
}

// Converts a parsed raw alignment into the compact gene representation.
fn build_gene(raw: RawAlignment, sample_indices: &HashMap<&str, usize>) -> Gene {
    let global_sample_indices = raw
        .sample_names
        .iter()
        .map(|sample| {
            sample_indices
                .get(sample.as_str())
                .copied()
                .expect("raw alignment sample should exist in global sample index")
        })
        .collect();

    Gene::new(
        raw.gene_name,
        raw.alignment_len,
        global_sample_indices,
        raw.sequences,
    )
}

// Parses one FASTA alignment and validates its records.
fn parse_raw_alignment(path: &Path, paralog_mode: ParalogMode) -> Result<RawAlignment> {
    let path = path.to_path_buf();
    let gene_name = gene_name_from_path(&path)?;
    let reader = open_alignment_reader(&path)?;
    let mut reader = Reader::new(reader);
    let mut sample_groups: Vec<RawSampleGroup> = Vec::new();
    let mut sample_group_indices: HashMap<String, usize> = HashMap::new();
    let mut alignment_len = None;

    while let Some(record) = reader.next() {
        let record = record
            .with_context(|| format!("failed to parse FASTA alignment '{}'", path.display()))?;

        let sample = normalize_sample_id(record.id().with_context(|| {
            format!(
                "sample header in alignment '{}' is not valid UTF-8",
                path.display()
            )
        })?);

        let sequence = record.full_seq();
        let observed_len = sequence.len();

        if observed_len == 0 {
            bail!(
                "sample '{sample}' in alignment '{}' has a zero-length sequence",
                path.display()
            );
        }

        if observed_len > u32::MAX as usize {
            bail!(
                "alignment '{}' has {observed_len} columns, exceeding the {} column limit",
                path.display(),
                u32::MAX
            );
        }

        match alignment_len {
            Some(expected) if expected != observed_len => {
                bail!(
                    "alignment '{}' has variable sequence lengths: sample '{sample}' has length {observed_len}, expected {expected}",
                    path.display()
                );
            }
            Some(_) => {}
            None => alignment_len = Some(observed_len),
        }

        let sequence = sequence.into_owned();
        let record = RawSampleRecord {
            non_gap_count: non_gap_count(&sequence),
            sequence,
        };

        match sample_group_indices.entry(sample.clone()) {
            Entry::Occupied(entry) => {
                sample_groups[*entry.get()].records.push(record);
            }
            Entry::Vacant(entry) => {
                let group_index = sample_groups.len();
                entry.insert(group_index);
                sample_groups.push(RawSampleGroup {
                    sample_name: sample,
                    records: vec![record],
                });
            }
        }
    }

    let Some(alignment_len) = alignment_len else {
        bail!("alignment '{}' contains no FASTA records", path.display());
    };

    let affected_sample_count = sample_groups
        .iter()
        .filter(|group| group.records.len() > 1)
        .count();
    let (sample_names, sequences) = resolve_paralogs(sample_groups, paralog_mode);

    Ok(RawAlignment {
        gene_name,
        sample_names,
        sequences,
        alignment_len,
        paralog: (affected_sample_count > 0).then_some(affected_sample_count),
    })
}

// Resolves duplicate normalized sample IDs according to the selected mode.
fn resolve_paralogs(
    sample_groups: Vec<RawSampleGroup>,
    paralog_mode: ParalogMode,
) -> (Vec<String>, Vec<Vec<u8>>) {
    let mut sample_names = Vec::with_capacity(sample_groups.len());
    let mut sequences = Vec::with_capacity(sample_groups.len());

    for group in sample_groups {
        let Some(sequence) = resolve_sample_group(group.records, paralog_mode) else {
            continue;
        };

        sample_names.push(group.sample_name);
        sequences.push(sequence);
    }

    (sample_names, sequences)
}

// Selects the sequence to keep for a single normalized sample group.
fn resolve_sample_group(
    mut records: Vec<RawSampleRecord>,
    paralog_mode: ParalogMode,
) -> Option<Vec<u8>> {
    debug_assert!(!records.is_empty());

    match paralog_mode {
        ParalogMode::First => Some(records.remove(0).sequence),
        ParalogMode::Skip if records.len() > 1 => None,
        ParalogMode::Skip => Some(records.remove(0).sequence),
        ParalogMode::Longest => {
            let mut best_index = 0;
            let mut best_non_gap_count = records[0].non_gap_count;
            for (index, record) in records.iter().enumerate().skip(1) {
                if record.non_gap_count > best_non_gap_count {
                    best_index = index;
                    best_non_gap_count = record.non_gap_count;
                }
            }

            Some(records.swap_remove(best_index).sequence)
        }
    }
}

// Counts non-gap characters for longest-paralog resolution.
fn non_gap_count(sequence: &[u8]) -> usize {
    sequence.iter().filter(|&&base| base != b'-').count()
}

// Normalizes a FASTA record identifier to its sample name.
fn normalize_sample_id(record_id: &str) -> String {
    record_id
        .split_once(';')
        .map_or(record_id, |(sample, _)| sample)
        .to_owned()
}

// Opens plain or gzip-compressed alignment input for reading.
fn open_alignment_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path)
        .with_context(|| format!("failed to read alignment '{}'", path.display()))?;

    if has_extension(path, "gz") || has_extension(path, "bgz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

// Derives a gene name by stripping known alignment/compression suffixes.
fn gene_name_from_path(path: &Path) -> Result<String> {
    let Some(file_name) = path.file_name().and_then(|name| name.to_str()) else {
        bail!("alignment path has no filename: '{}'", path.display());
    };

    let mut name = file_name.to_owned();

    for suffix in ["gz", "bgz", "bz2", "xz", "zst"] {
        strip_extension(&mut name, suffix);
    }

    for suffix in ["aln.fas", "aln", "fa", "fasta", "fna", "fas"] {
        strip_extension(&mut name, suffix);
    }

    if name.is_empty() {
        bail!("alignment path has no filename: '{}'", path.display());
    }

    Ok(name)
}

// Checks a path extension without case sensitivity.
fn has_extension(path: &Path, extension: &str) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case(extension))
}

// Strips a matching suffix from a mutable file-name string.
fn strip_extension(name: &mut String, extension: &str) {
    let suffix = format!(".{extension}");
    if name
        .get(name.len().saturating_sub(suffix.len())..)
        .is_some_and(|candidate| candidate.eq_ignore_ascii_case(&suffix))
    {
        name.truncate(name.len() - suffix.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    // Writes a temporary alignment file and returns its path.
    fn write_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let path = dir.path().join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    // Parses a temporary alignment with the default paralog behavior.
    fn parse_default(path: &Path) -> Result<RawAlignment> {
        parse_raw_alignment(path, ParalogMode::First)
    }

    #[test]
    // Verifies valid equal-length FASTA records are accepted.
    fn parser_accepts_valid_equal_length_fasta() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1\nACGT\n>sample2 description\nTGCA\n",
        );

        let raw = parse_default(&path).unwrap();
        assert_eq!(raw.gene_name, "gene");
        assert_eq!(raw.alignment_len, 4);
        assert_eq!(raw.paralog, None);
        assert_eq!(raw.sample_names, vec!["sample1", "sample2"]);
        assert_eq!(raw.sequences.len(), 2);
    }

    #[test]
    // Verifies Panaroo compound suffixes are stripped from gene names.
    fn parser_strips_panaroo_alignment_suffix_from_gene_name() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln.fas", ">sample1\nACGT\n>sample2\nTGCA\n");

        let raw = parse_default(&path).unwrap();

        assert_eq!(raw.gene_name, "gene");
    }

    #[test]
    // Verifies sample IDs are truncated at the first semicolon.
    fn parser_normalizes_sample_ids_before_first_semicolon() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">cl7_bakta;4_0_1479\nACGT\n>cl3_bakta;2_1_12 description\nTGCA\n",
        );

        let raw = parse_default(&path).unwrap();

        assert_eq!(raw.sample_names, vec!["cl7_bakta", "cl3_bakta"]);
    }

    #[test]
    // Verifies duplicate normalized sample IDs keep the first sequence by default.
    fn parser_keeps_first_duplicate_normalized_sample_id_by_default() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample2\nCCCC\n>sample1;second\nTGCA\n",
        );

        let raw = parse_default(&path).unwrap();

        assert_eq!(raw.paralog, Some(1));
        assert_eq!(raw.sample_names, vec!["sample1", "sample2"]);
        assert_eq!(raw.sequences, vec![b"ACGT".to_vec(), b"CCCC".to_vec()]);
    }

    #[test]
    // Verifies skip mode removes duplicated samples from that gene.
    fn parser_skip_mode_removes_duplicated_samples() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample2\nCCCC\n>sample1;second\nTGCA\n>sample3\nGGGG\n",
        );

        let raw = parse_raw_alignment(&path, ParalogMode::Skip).unwrap();

        assert_eq!(raw.paralog, Some(1));
        assert_eq!(raw.sample_names, vec!["sample2", "sample3"]);
        assert_eq!(raw.sequences, vec![b"CCCC".to_vec(), b"GGGG".to_vec()]);
    }

    #[test]
    // Verifies longest mode keeps the duplicate with the most non-gap bytes.
    fn parser_longest_mode_keeps_duplicate_with_most_non_gap_bases() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;short\nAA--\n>sample2\nCCCC\n>sample1;long\nA-CG\n",
        );

        let raw = parse_raw_alignment(&path, ParalogMode::Longest).unwrap();

        assert_eq!(raw.paralog, Some(1));
        assert_eq!(raw.sample_names, vec!["sample1", "sample2"]);
        assert_eq!(raw.sequences, vec![b"A-CG".to_vec(), b"CCCC".to_vec()]);
    }

    #[test]
    // Verifies longest mode keeps the first duplicate when non-gap counts tie.
    fn parser_longest_mode_keeps_first_duplicate_on_equal_non_gap_count() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nA--C\n>sample1;second\nG--T\n",
        );

        let raw = parse_raw_alignment(&path, ParalogMode::Longest).unwrap();

        assert_eq!(raw.sample_names, vec!["sample1"]);
        assert_eq!(raw.sequences, vec![b"A--C".to_vec()]);
    }

    #[test]
    // Verifies malformed duplicate records are still validated before resolution.
    fn parser_rejects_variable_sequence_lengths_in_duplicate_records() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample1;second\nACG\n",
        );

        let error = parse_raw_alignment(&path, ParalogMode::First).unwrap_err();
        assert!(
            error.to_string().contains("variable sequence lengths"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies records with inconsistent sequence lengths are rejected.
    fn parser_rejects_variable_sequence_lengths() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\nACGT\n>sample2\nACG\n");

        let error = parse_default(&path).unwrap_err();
        let message = error.to_string();
        assert!(
            message.contains("variable sequence lengths"),
            "error: {message}"
        );
        assert!(message.contains("length 3"), "error: {message}");
        assert!(message.contains("expected 4"), "error: {message}");
    }

    #[test]
    // Verifies per-gene sample order is preserved independently.
    fn sample_order_can_differ_between_genes() {
        let dir = TempDir::new().unwrap();
        let gene_a = write_alignment(&dir, "gene_a.aln", ">s1\nAAAA\n>s2\nCCCC\n");
        let gene_b = write_alignment(&dir, "gene_b.aln", ">s2\nCCCC\n>s1\nAAAA\n");

        let loaded = load_genes(&[gene_a, gene_b], ParalogMode::First, true).unwrap();

        assert_eq!(loaded.sample_names, vec!["s1", "s2"]);
        assert_eq!(loaded.genes.len(), 2);
        assert_eq!(loaded.genes[0].sample_index(0), Some(0));
        assert_eq!(loaded.genes[0].sample_index(1), Some(1));
        assert_eq!(loaded.genes[1].sample_index(1), Some(0));
        assert_eq!(loaded.genes[1].sample_index(0), Some(1));
        assert_eq!(loaded.genes[0].snp_count(0, 1, false), (4, 4));
        assert_eq!(loaded.genes[1].snp_count(0, 1, false), (4, 4));
    }

    #[test]
    // Verifies genes may contain different subsets of samples.
    fn genes_may_have_different_sample_sets() {
        let dir = TempDir::new().unwrap();
        let gene_a = write_alignment(
            &dir,
            "gene_a.aln",
            ">s1;contig_a\nAAAA\n>s2;contig_a\nCCCC\n",
        );
        let gene_b = write_alignment(
            &dir,
            "gene_b.aln",
            ">s3;contig_b\nCCCC\n>s1;contig_b\nAAAA\n",
        );

        let loaded = load_genes(&[gene_a, gene_b], ParalogMode::First, true).unwrap();

        assert_eq!(loaded.sample_names, vec!["s1", "s2", "s3"]);
        assert_eq!(loaded.genes.len(), 2);
        assert_eq!(loaded.genes[0].sample_index(0), Some(0));
        assert_eq!(loaded.genes[0].sample_index(1), Some(1));
        assert_eq!(loaded.genes[1].sample_index(2), Some(0));
        assert_eq!(loaded.genes[1].sample_index(0), Some(1));
    }

    #[test]
    // Verifies load metadata records affected genes in input order.
    fn load_genes_collects_paralog_report_rows_in_input_order() {
        let dir = TempDir::new().unwrap();
        let gene_a = write_alignment(
            &dir,
            "gene_a.aln",
            ">s1;first\nAAAA\n>s2\nCCCC\n>s1;second\nTTTT\n",
        );
        let gene_clean = write_alignment(&dir, "gene_clean.aln", ">s1\nAAAA\n>s2\nCCCC\n");
        let gene_b = write_alignment(
            &dir,
            "gene_b.aln",
            concat!(
                ">s1;first\nAAAA\n",
                ">s2;first\nCCCC\n",
                ">s1;second\nTTTT\n",
                ">s2;second\nGGGG\n",
            ),
        );

        let loaded = load_genes(&[gene_a, gene_clean, gene_b], ParalogMode::First, true).unwrap();

        assert_eq!(
            loaded.paralogs,
            vec![("gene_a".to_string(), 1), ("gene_b".to_string(), 2)]
        );
    }

    #[test]
    // Verifies MSA lists ignore comments and blank lines.
    fn parse_msa_list_ignores_blank_lines_and_comments() {
        let observed = parse_msa_list(
            Path::new("fixtures/list.txt"),
            "\n  \n# comment\n  # indented comment\nalpha.aln\n\nbeta.aln\n",
        );

        assert_eq!(
            observed,
            vec![
                PathBuf::from("fixtures/alpha.aln"),
                PathBuf::from("fixtures/beta.aln")
            ]
        );
    }

    #[test]
    // Verifies relative MSA entries are resolved against the list directory.
    fn parse_msa_list_resolves_relative_paths_against_list_directory() {
        let observed = parse_msa_list(
            Path::new("data/lists/msa.txt"),
            "../gene.aln\nnested/gene.aln",
        );

        assert_eq!(
            observed,
            vec![
                PathBuf::from("data/lists/../gene.aln"),
                PathBuf::from("data/lists/nested/gene.aln")
            ]
        );
    }

    #[test]
    // Verifies absolute MSA entries are preserved unchanged.
    fn parse_msa_list_preserves_absolute_paths() {
        let absolute = std::env::current_dir().unwrap().join("gene.aln");
        let observed = parse_msa_list(Path::new("data/list.txt"), &absolute.to_string_lossy());

        assert_eq!(observed, vec![absolute]);
    }

    #[test]
    // Verifies Panaroo discovery uses sorted top-level .aln.fas files first.
    fn read_panaroo_dir_uses_sorted_top_level_aln_fas_files() {
        let dir = TempDir::new().unwrap();
        let first = dir.path().join("alpha.aln.fas");
        let second = dir.path().join("zeta.aln.fas");
        fs::write(&second, "").unwrap();
        fs::write(dir.path().join("ignored.fas"), "").unwrap();
        fs::write(&first, "").unwrap();

        let fallback = dir.path().join("aligned_gene_sequences");
        fs::create_dir(&fallback).unwrap();
        fs::write(fallback.join("fallback.aln.fas"), "").unwrap();

        let observed = read_panaroo_dir(dir.path()).unwrap();

        assert_eq!(observed, vec![first, second]);
    }

    #[test]
    // Verifies Panaroo discovery falls back to aligned_gene_sequences.
    fn read_panaroo_dir_falls_back_to_aligned_gene_sequences() {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("ignored.fas"), "").unwrap();
        let fallback = dir.path().join("aligned_gene_sequences");
        fs::create_dir(&fallback).unwrap();
        let gene = fallback.join("gene.aln.fas");
        fs::write(&gene, "").unwrap();

        let observed = read_panaroo_dir(dir.path()).unwrap();

        assert_eq!(observed, vec![gene]);
    }

    #[test]
    // Verifies Panaroo discovery errors when neither location contains alignments.
    fn read_panaroo_dir_rejects_empty_locations() {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("ignored.fas"), "").unwrap();

        let error = read_panaroo_dir(dir.path()).unwrap_err();

        assert!(
            error
                .to_string()
                .contains("found no Panaroo alignment files ending in .aln.fas"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies recombination tables are emitted as expected TSV.
    fn write_recombination_table_emits_tsv() {
        let sample_names = vec!["alpha".to_string(), "beta".to_string()];
        let genes = vec![
            Gene::new(
                "gene1".to_string(),
                1,
                vec![0, 1],
                vec![b"A".to_vec(), b"A".to_vec()],
            ),
            Gene::new(
                "gene2".to_string(),
                1,
                vec![0, 1],
                vec![b"A".to_vec(), b"A".to_vec()],
            ),
        ];
        let rows = vec![
            RecombinationRow {
                gene_index: 0,
                presence: vec![1, 0],
            },
            RecombinationRow {
                gene_index: 1,
                presence: vec![0, 1],
            },
        ];
        let mut output = Vec::new();

        write_recombination_table(&sample_names, &genes, &rows, &mut output).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "gene\talpha\tbeta\ngene1\t1\t0\ngene2\t0\t1\n"
        );
    }
}
