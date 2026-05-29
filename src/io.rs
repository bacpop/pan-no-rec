use crate::gene::Gene;
use crate::graph::RecombinationTable;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use seq_io::fasta::{Reader, Record};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
struct RawAlignment {
    gene_name: String,
    sample_names: Vec<String>,
    sequences: Vec<Vec<u8>>,
    alignment_len: usize,
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

    paths.sort();
    Ok(paths)
}

// Writes a recombination presence table as tab-separated text.
pub fn write_recombination_table<W: Write>(
    table: &RecombinationTable,
    mut writer: W,
) -> Result<()> {
    write!(writer, "gene")?;
    for sample in &table.sample_names {
        write!(writer, "\t{sample}")?;
    }
    writeln!(writer)?;

    for row in &table.rows {
        write!(writer, "{}", row.gene)?;
        for value in &row.presence {
            write!(writer, "\t{value}")?;
        }
        writeln!(writer)?;
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
pub(crate) fn load_genes<P>(aln_paths: &[P]) -> Result<(Vec<String>, Vec<Gene>)>
where
    P: AsRef<Path>,
{
    let mut genes = Vec::with_capacity(aln_paths.len());
    let mut all_samples = HashSet::new();

    for aln_path in aln_paths {
        let raw = parse_raw_alignment(aln_path.as_ref())?;
        all_samples.extend(raw.sample_names.iter().cloned());
        genes.push(build_gene(raw));
    }

    let mut sample_names: Vec<_> = all_samples.into_iter().collect();
    sample_names.sort();

    Ok((sample_names, genes))
}

// Converts a parsed raw alignment into the compact gene representation.
fn build_gene(raw: RawAlignment) -> Gene {
    Gene::new(
        raw.gene_name,
        raw.alignment_len,
        raw.sample_names,
        raw.sequences,
    )
}

// Parses one FASTA alignment and validates its records.
fn parse_raw_alignment(path: &Path) -> Result<RawAlignment> {
    let path = path.to_path_buf();
    let gene_name = gene_name_from_path(&path)?;
    let reader = open_alignment_reader(&path)?;
    let mut reader = Reader::new(reader);
    let mut sample_names = Vec::new();
    let mut sequences = Vec::new();
    let mut seen_samples = HashSet::new();
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

        if !seen_samples.insert(sample.clone()) {
            bail!(
                "alignment '{}' contains duplicate sample header '{sample}'",
                path.display()
            );
        }

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

        sample_names.push(sample);
        sequences.push(sequence.into_owned());
    }

    let Some(alignment_len) = alignment_len else {
        bail!("alignment '{}' contains no FASTA records", path.display());
    };

    Ok(RawAlignment {
        gene_name,
        sample_names,
        sequences,
        alignment_len,
    })
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
    use std::path::PathBuf;
    use tempfile::TempDir;

    // Writes a temporary alignment file and returns its path.
    fn write_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let path = dir.path().join(name);
        fs::write(&path, contents).unwrap();
        path
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

        let raw = parse_raw_alignment(&path).unwrap();
        assert_eq!(raw.gene_name, "gene");
        assert_eq!(raw.alignment_len, 4);
        assert_eq!(raw.sample_names, vec!["sample1", "sample2"]);
        assert_eq!(raw.sequences.len(), 2);
    }

    #[test]
    // Verifies Panaroo compound suffixes are stripped from gene names.
    fn parser_strips_panaroo_alignment_suffix_from_gene_name() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln.fas", ">sample1\nACGT\n>sample2\nTGCA\n");

        let raw = parse_raw_alignment(&path).unwrap();

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

        let raw = parse_raw_alignment(&path).unwrap();

        assert_eq!(raw.sample_names, vec!["cl7_bakta", "cl3_bakta"]);
    }

    #[test]
    // Verifies duplicate normalized sample IDs are rejected.
    fn parser_rejects_duplicate_normalized_sample_ids() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample1;second\nTGCA\n",
        );

        let error = parse_raw_alignment(&path).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("duplicate sample header 'sample1'"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies records with inconsistent sequence lengths are rejected.
    fn parser_rejects_variable_sequence_lengths() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\nACGT\n>sample2\nACG\n");

        let error = parse_raw_alignment(&path).unwrap_err();
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

        let (sample_names, genes) = load_genes(&[gene_a, gene_b]).unwrap();

        assert_eq!(sample_names, vec!["s1", "s2"]);
        assert_eq!(genes.len(), 2);
        assert_eq!(genes[0].sample_index("s1"), Some(0));
        assert_eq!(genes[0].sample_index("s2"), Some(1));
        assert_eq!(genes[1].sample_index("s2"), Some(0));
        assert_eq!(genes[1].sample_index("s1"), Some(1));
        assert_eq!(genes[0].snp_count(0, 1), 4);
        assert_eq!(genes[1].snp_count(0, 1), 4);
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

        let (sample_names, genes) = load_genes(&[gene_a, gene_b]).unwrap();

        assert_eq!(sample_names, vec!["s1", "s2", "s3"]);
        assert_eq!(genes.len(), 2);
        assert_eq!(genes[0].sample_index("s1"), Some(0));
        assert_eq!(genes[0].sample_index("s2"), Some(1));
        assert_eq!(genes[1].sample_index("s3"), Some(0));
        assert_eq!(genes[1].sample_index("s1"), Some(1));
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
        let table = RecombinationTable {
            sample_names: vec!["alpha".to_string(), "beta".to_string()],
            rows: vec![
                crate::graph::RecombinationRow {
                    gene: "gene1".to_string(),
                    presence: vec![1, 0],
                },
                crate::graph::RecombinationRow {
                    gene: "gene2".to_string(),
                    presence: vec![0, 1],
                },
            ],
        };
        let mut output = Vec::new();

        write_recombination_table(&table, &mut output).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "gene\talpha\tbeta\ngene1\t1\t0\ngene2\t0\t1\n"
        );
    }
}
