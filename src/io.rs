use crate::RecombinationTable;
use crate::error::CompareError;
use crate::gene::Gene;
use flate2::read::MultiGzDecoder;
use seq_io::fasta::{Reader, Record};
use std::collections::HashSet;
use std::error::Error;
use std::fmt;
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

#[derive(Debug)]
pub enum MsaListError {
    Read {
        path: PathBuf,
        source: std::io::Error,
    },
}

impl fmt::Display for MsaListError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MsaListError::Read { path, source } => {
                write!(f, "failed to read MSA list '{}': {source}", path.display())
            }
        }
    }
}

impl Error for MsaListError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            MsaListError::Read { source, .. } => Some(source),
        }
    }
}

pub fn read_msa_list(path: &Path) -> Result<Vec<PathBuf>, MsaListError> {
    let contents = fs::read_to_string(path).map_err(|source| MsaListError::Read {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(parse_msa_list(path, &contents))
}

pub fn write_recombination_table<W: Write>(
    table: &RecombinationTable,
    mut writer: W,
) -> Result<(), std::io::Error> {
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

pub(crate) fn load_genes<P>(aln_paths: &[P]) -> Result<(Vec<String>, Vec<Gene>), CompareError>
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

fn build_gene(raw: RawAlignment) -> Gene {
    Gene::new(
        raw.gene_name,
        raw.alignment_len,
        raw.sample_names,
        raw.sequences,
    )
}

fn parse_raw_alignment(path: &Path) -> Result<RawAlignment, CompareError> {
    let path = path.to_path_buf();
    let gene_name = gene_name_from_path(&path)?;
    let reader = open_alignment_reader(&path)?;
    let mut reader = Reader::new(reader);
    let mut sample_names = Vec::new();
    let mut sequences = Vec::new();
    let mut seen_samples = HashSet::new();
    let mut alignment_len = None;

    while let Some(record) = reader.next() {
        let record = record.map_err(|source| CompareError::FastaParse {
            path: path.clone(),
            source,
        })?;

        let sample =
            normalize_sample_id(record.id().map_err(|source| CompareError::HeaderUtf8 {
                path: path.clone(),
                source,
            })?);

        if !seen_samples.insert(sample.clone()) {
            return Err(CompareError::DuplicateSample { path, sample });
        }

        let sequence = record.full_seq();
        let observed_len = sequence.len();

        if observed_len == 0 {
            return Err(CompareError::EmptySequence { path, sample });
        }

        if observed_len > u32::MAX as usize {
            return Err(CompareError::AlignmentTooLong {
                path,
                length: observed_len,
            });
        }

        match alignment_len {
            Some(expected) if expected != observed_len => {
                return Err(CompareError::VariableLength {
                    path,
                    sample,
                    expected,
                    observed: observed_len,
                });
            }
            Some(_) => {}
            None => alignment_len = Some(observed_len),
        }

        sample_names.push(sample);
        sequences.push(sequence.into_owned());
    }

    let alignment_len =
        alignment_len.ok_or_else(|| CompareError::EmptyAlignment { path: path.clone() })?;

    Ok(RawAlignment {
        gene_name,
        sample_names,
        sequences,
        alignment_len,
    })
}

fn normalize_sample_id(record_id: &str) -> String {
    record_id
        .split_once(';')
        .map_or(record_id, |(sample, _)| sample)
        .to_owned()
}

fn open_alignment_reader(path: &Path) -> Result<Box<dyn Read>, CompareError> {
    let file = File::open(path).map_err(|source| CompareError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    if has_extension(path, "gz") || has_extension(path, "bgz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

fn gene_name_from_path(path: &Path) -> Result<String, CompareError> {
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| CompareError::InvalidPath {
            path: path.to_path_buf(),
        })?;

    let mut name = file_name.to_owned();

    for suffix in ["gz", "bgz", "bz2", "xz", "zst"] {
        strip_extension(&mut name, suffix);
    }

    for suffix in ["aln", "fa", "fasta", "fna", "fas"] {
        strip_extension(&mut name, suffix);
    }

    if name.is_empty() {
        return Err(CompareError::InvalidPath {
            path: path.to_path_buf(),
        });
    }

    Ok(name)
}

fn has_extension(path: &Path, extension: &str) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case(extension))
}

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

    fn write_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let path = dir.path().join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
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
    fn parser_rejects_duplicate_normalized_sample_ids() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample1;second\nTGCA\n",
        );

        let error = parse_raw_alignment(&path).unwrap_err();

        match error {
            CompareError::DuplicateSample { sample, .. } => {
                assert_eq!(sample, "sample1");
            }
            other => panic!("expected duplicate sample, got {other:?}"),
        }
    }

    #[test]
    fn parser_rejects_variable_sequence_lengths() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\nACGT\n>sample2\nACG\n");

        let error = parse_raw_alignment(&path).unwrap_err();
        assert!(matches!(error, CompareError::VariableLength { .. }));
    }

    #[test]
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
    fn parse_msa_list_preserves_absolute_paths() {
        let absolute = std::env::current_dir().unwrap().join("gene.aln");
        let observed = parse_msa_list(Path::new("data/list.txt"), &absolute.to_string_lossy());

        assert_eq!(observed, vec![absolute]);
    }

    #[test]
    fn write_recombination_table_emits_tsv() {
        let table = RecombinationTable {
            sample_names: vec!["alpha".to_string(), "beta".to_string()],
            rows: vec![
                crate::RecombinationRow {
                    gene: "gene1".to_string(),
                    presence: vec![1, 0],
                },
                crate::RecombinationRow {
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
