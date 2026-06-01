use crate::cli::ParalogMode;
use crate::gene::{Gene, SampleBases};
use crate::get_progress_bar;
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use hashbrown::{HashMap, hash_map::Entry};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use seq_io::fasta::{Reader, Record};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

#[derive(Debug)]
pub(crate) struct LoadedAlignments {
    pub(crate) sample_names: Vec<String>,
    pub(crate) gene_sequences: Vec<Gene>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct OutputRow {
    pub(crate) gene_index: usize,
    pub(crate) presence: Vec<u8>,
}

// Finds Panaroo gene alignment files in the standard alignment directory.
pub fn read_panaroo_dir(path: &Path) -> Result<Vec<PathBuf>> {
    let alignment_dir = path.join("aligned_gene_sequences");
    let paths = collect_panaroo_alignments(&alignment_dir)?;

    if paths.is_empty() {
        bail!(
            "found no Panaroo alignment files ending in .aln.fas in '{}'",
            alignment_dir.display()
        );
    }

    Ok(paths)
}

fn collect_panaroo_alignments(dir: &Path) -> Result<Vec<PathBuf>> {
    let entries = fs::read_dir(dir).with_context(|| {
        format!(
            "failed to read Panaroo alignment directory '{}'",
            dir.display()
        )
    })?;

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
    rows: &[OutputRow],
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
pub(crate) fn write_paralog_report(path: &Path, genes: &[Gene]) -> Result<usize> {
    let paralog_rows: Vec<_> = genes
        .iter()
        .filter_map(|gene| {
            gene.paralog_count()
                .map(|paralog_count| (gene.name(), paralog_count))
        })
        .collect();

    if paralog_rows.is_empty() {
        return Ok(0);
    }

    let mut writer = File::create(path)
        .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;

    writeln!(writer, "gene\tparalog_samples")
        .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    for (gene_name, paralog_count) in &paralog_rows {
        writeln!(writer, "{gene_name}\t{paralog_count}")
            .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    }

    Ok(paralog_rows.len())
}

// Loads all Panaroo alignments into genes using the Rtab header sample order.
pub(crate) fn load_genes(
    panaroo_dir: &Path,
    paralog_mode: ParalogMode,
    max_entropy: Option<f64>,
    quiet: bool,
) -> Result<LoadedAlignments> {
    let rtab_path = panaroo_dir.join("gene_presence_absence.Rtab");
    let sample_names = read_rtab_sample_names(&rtab_path)?;
    let genes = {
        let sample_indices = build_sample_indices(&sample_names, &rtab_path)?;
        let mut aln_paths = read_panaroo_dir(panaroo_dir)?;
        if let Some(threshold) = max_entropy {
            aln_paths = filter_alignments_by_entropy(panaroo_dir, aln_paths, threshold)?;
        }

        let pbar = get_progress_bar(aln_paths.len(), false, quiet);
        aln_paths
            .par_iter()
            .progress_with(pbar)
            .map(|aln| parse_gene_alignment(aln, &sample_indices, paralog_mode))
            .collect::<Result<Vec<_>>>()?
    };

    Ok(LoadedAlignments {
        sample_names,
        gene_sequences: genes,
    })
}

// Reads only the header row of Panaroo's Rtab and returns its sample columns.
fn read_rtab_sample_names(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path)
        .with_context(|| format!("failed to read Panaroo Rtab '{}'", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    let bytes = reader
        .read_line(&mut header)
        .with_context(|| format!("failed to read Panaroo Rtab '{}'", path.display()))?;

    if bytes == 0 {
        bail!("Panaroo Rtab '{}' contains no header row", path.display());
    }

    parse_rtab_header(path, &header)
}

// Parses the Rtab header and preserves Panaroo's sample column order.
fn parse_rtab_header(path: &Path, header: &str) -> Result<Vec<String>> {
    let header = header.trim_end_matches(['\r', '\n']);
    let columns: Vec<_> = header.split('\t').collect();

    if columns.first() != Some(&"Gene") {
        bail!(
            "Panaroo Rtab '{}' header must start with 'Gene'",
            path.display()
        );
    }

    if columns.len() < 2 {
        bail!(
            "Panaroo Rtab '{}' header must contain at least one sample column",
            path.display()
        );
    }

    let sample_names: Vec<_> = columns[1..]
        .iter()
        .map(|sample| (*sample).to_owned())
        .collect();
    build_sample_indices(&sample_names, path)?;

    Ok(sample_names)
}

// Builds the global sample index and rejects duplicate or empty Rtab samples.
fn build_sample_indices<'a>(
    sample_names: &'a [String],
    rtab_path: &Path,
) -> Result<HashMap<&'a str, usize>> {
    let mut sample_indices = HashMap::with_capacity(sample_names.len());
    for (index, sample) in sample_names.iter().enumerate() {
        if sample.is_empty() {
            bail!(
                "Panaroo Rtab '{}' header contains an empty sample name",
                rtab_path.display()
            );
        }

        if sample_indices.insert(sample.as_str(), index).is_some() {
            bail!(
                "Panaroo Rtab '{}' header contains duplicate sample name '{sample}'",
                rtab_path.display()
            );
        }
    }

    Ok(sample_indices)
}

// Parses one FASTA alignment and validates its records.
fn parse_gene_alignment(
    path: &Path,
    sample_indices: &HashMap<&str, usize>,
    paralog_mode: ParalogMode,
) -> Result<Gene> {
    let path = path.to_path_buf();
    let gene_name = gene_name_from_path(&path)?;
    let reader = open_alignment_reader(&path)?;
    let mut reader = Reader::new(reader);
    let mut parsed_sequences: HashMap<usize, SampleBases> = HashMap::new();
    let mut paralog_counts: HashMap<usize, usize> = HashMap::new();
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
        let Some(&global_sample_index) = sample_indices.get(sample) else {
            bail!(
                "sample '{sample}' in alignment '{}' does not appear in gene_presence_absence.Rtab",
                path.display()
            );
        };

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

        let new_sequence = SampleBases::from_sequence(sequence.as_ref());

        match parsed_sequences.entry(global_sample_index) {
            Entry::Occupied(mut entry) => {
                paralog_counts
                    .entry(global_sample_index)
                    .and_modify(|cnt| {
                        *cnt += 1;
                    })
                    .or_insert(1);
                match paralog_mode {
                    ParalogMode::First | ParalogMode::Skip => {
                        continue;
                    }
                    ParalogMode::Longest => {
                        if entry.get().non_gap_count(observed_len)
                            < new_sequence.non_gap_count(observed_len)
                        {
                            entry.insert(new_sequence);
                        }
                    }
                }
            }
            Entry::Vacant(entry) => {
                entry.insert(new_sequence);
            }
        }
    }

    let Some(alignment_len) = alignment_len else {
        bail!("alignment '{}' contains no FASTA records", path.display());
    };

    let paralog_count = paralog_counts.values().sum();
    if paralog_count > 0 && paralog_mode == ParalogMode::Skip {
        for sample_index in paralog_counts.keys() {
            parsed_sequences.remove(sample_index);
        }
    }

    Ok(Gene::new(
        gene_name,
        alignment_len,
        parsed_sequences,
        paralog_count,
    ))
}

// Drops alignments with entropy strictly greater than the requested threshold.
fn filter_alignments_by_entropy(
    panaroo_dir: &Path,
    aln_paths: Vec<PathBuf>,
    max_entropy: f64,
) -> Result<Vec<PathBuf>> {
    let entropy_path = panaroo_dir.join("alignment_entropy.csv");
    let entropies = read_alignment_entropy(&entropy_path)?;
    let mut retained = Vec::with_capacity(aln_paths.len());
    let mut removed_count = 0;
    let mut missing_count = 0;

    for aln_path in aln_paths {
        let gene_name = gene_name_from_path(&aln_path)?;
        match entropies.get(&gene_name) {
            Some(&entropy) if entropy > max_entropy => {
                removed_count += 1;
            }
            Some(_) => retained.push(aln_path),
            None => {
                missing_count += 1;
                retained.push(aln_path);
            }
        }
    }

    if missing_count > 0 {
        log::warn!(
            "{} alignments lacked entropy metadata in '{}'; keeping them",
            missing_count,
            entropy_path.display()
        );
    }
    log::info!(
        "Filtered {} alignments with entropy > {}",
        removed_count,
        max_entropy
    );

    Ok(retained)
}

// Reads Panaroo alignment entropy metadata keyed by normalized gene name.
fn read_alignment_entropy(path: &Path) -> Result<HashMap<String, f64>> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read alignment entropy CSV '{}'", path.display()))?;
    let mut lines = contents
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty());
    let Some((header_line_number, header_line)) = lines.next() else {
        bail!(
            "alignment entropy CSV '{}' contains no header row",
            path.display()
        );
    };
    let header = parse_csv_record(header_line).with_context(|| {
        format!(
            "malformed alignment entropy CSV '{}': line {}",
            path.display(),
            header_line_number + 1
        )
    })?;
    let (gene_column, entropy_column) = entropy_csv_columns(path, &header)?;
    let mut entropies = HashMap::new();

    for (line_number, line) in lines {
        let fields = parse_csv_record(line).with_context(|| {
            format!(
                "malformed alignment entropy CSV '{}': line {}",
                path.display(),
                line_number + 1
            )
        })?;
        if fields.len() <= gene_column.max(entropy_column) {
            bail!(
                "malformed alignment entropy CSV '{}': line {} has too few columns",
                path.display(),
                line_number + 1
            );
        }

        let raw_gene_name = fields[gene_column].trim();
        if raw_gene_name.is_empty() {
            bail!(
                "malformed alignment entropy CSV '{}': line {} has an empty gene name",
                path.display(),
                line_number + 1
            );
        }
        let gene_name = normalize_gene_name(raw_gene_name).with_context(|| {
            format!(
                "malformed alignment entropy CSV '{}': line {}",
                path.display(),
                line_number + 1
            )
        })?;

        let entropy = fields[entropy_column]
            .trim()
            .parse::<f64>()
            .with_context(|| {
                format!(
                    "malformed alignment entropy CSV '{}': line {} has invalid entropy",
                    path.display(),
                    line_number + 1
                )
            })?;
        if !entropy.is_finite() {
            bail!(
                "malformed alignment entropy CSV '{}': line {} has non-finite entropy",
                path.display(),
                line_number + 1
            );
        }

        if entropies.insert(gene_name.clone(), entropy).is_some() {
            bail!(
                "malformed alignment entropy CSV '{}' contains duplicate entropy rows for gene '{}'",
                path.display(),
                gene_name
            );
        }
    }

    Ok(entropies)
}

// Finds the gene and entropy columns in the CSV header.
fn entropy_csv_columns(path: &Path, header: &[String]) -> Result<(usize, usize)> {
    if header.len() < 2 {
        bail!(
            "malformed alignment entropy CSV '{}' header must contain at least two columns",
            path.display()
        );
    }

    let gene_column = header
        .iter()
        .position(|column| normalize_csv_header(column).contains("gene"))
        .unwrap_or(0);
    let entropy_column = header
        .iter()
        .position(|column| normalize_csv_header(column).contains("entropy"))
        .or_else(|| (header.len() == 2).then_some(1))
        .with_context(|| {
            format!(
                "malformed alignment entropy CSV '{}' header does not contain an entropy column",
                path.display()
            )
        })?;

    if gene_column == entropy_column {
        bail!(
            "malformed alignment entropy CSV '{}' uses the same column for gene and entropy",
            path.display()
        );
    }

    Ok((gene_column, entropy_column))
}

// Normalizes header names for loose matching.
fn normalize_csv_header(column: &str) -> String {
    column
        .chars()
        .filter(|character| character.is_ascii_alphanumeric())
        .flat_map(char::to_lowercase)
        .collect()
}

// Parses one simple CSV record, including quoted fields and escaped quotes.
fn parse_csv_record(line: &str) -> Result<Vec<String>> {
    let line = line.trim_end_matches('\r');
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;

    while let Some(character) = chars.next() {
        match character {
            '"' if in_quotes && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => {
                in_quotes = !in_quotes;
            }
            ',' if !in_quotes => {
                fields.push(field);
                field = String::new();
            }
            _ => field.push(character),
        }
    }

    if in_quotes {
        bail!("unterminated quoted CSV field");
    }

    fields.push(field);
    Ok(fields)
}

// Normalizes a FASTA record identifier to its sample name.
fn normalize_sample_id(record_id: &str) -> &str {
    record_id
        .split_once(';')
        .map_or(record_id, |(sample, _)| sample)
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

    normalize_gene_name(file_name)
}

// Derives a gene name by stripping known alignment/compression suffixes.
fn normalize_gene_name(name: &str) -> Result<String> {
    let Some(file_name) = Path::new(name).file_name().and_then(|name| name.to_str()) else {
        bail!("gene name is empty");
    };
    let mut name = file_name.to_owned();

    for suffix in ["gz", "bgz", "bz2", "xz", "zst"] {
        strip_extension(&mut name, suffix);
    }

    for suffix in ["aln.fas", "aln", "fa", "fasta", "fna", "fas"] {
        strip_extension(&mut name, suffix);
    }

    if name.is_empty() {
        bail!("gene name is empty");
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

    // Writes the required Panaroo Rtab header fixture.
    fn write_rtab(dir: &TempDir, sample_names: &[&str]) {
        fs::write(
            dir.path().join("gene_presence_absence.Rtab"),
            format!("Gene\t{}\n", sample_names.join("\t")),
        )
        .unwrap();
    }

    // Writes an alignment under Panaroo's standard alignment directory.
    fn write_panaroo_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let alignment_dir = dir.path().join("aligned_gene_sequences");
        fs::create_dir_all(&alignment_dir).unwrap();
        let path = alignment_dir.join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    // Parses a temporary alignment with explicit Rtab sample names.
    fn parse_with_samples(
        path: &Path,
        sample_names: &[&str],
        paralog_mode: ParalogMode,
    ) -> Result<Gene> {
        let sample_names: Vec<_> = sample_names
            .iter()
            .map(|sample| (*sample).to_owned())
            .collect();
        let sample_indices =
            build_sample_indices(&sample_names, Path::new("gene_presence_absence.Rtab")).unwrap();
        parse_gene_alignment(path, &sample_indices, paralog_mode)
    }

    // Parses a temporary alignment with the default paralog behavior.
    fn parse_default(path: &Path) -> Result<Gene> {
        parse_with_samples(path, &["sample1", "sample2"], ParalogMode::First)
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

        let gene = parse_default(&path).unwrap();
        assert_eq!(gene.name(), "gene");
        assert_eq!(gene.paralog_count(), None);
        assert!(gene.has_sample(0));
        assert!(gene.has_sample(1));
        assert_eq!(gene.snp_count(0, 1, false), (4, 4));
    }

    #[test]
    // Verifies Panaroo compound suffixes are stripped from gene names.
    fn parser_strips_panaroo_alignment_suffix_from_gene_name() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln.fas", ">sample1\nACGT\n>sample2\nTGCA\n");

        let gene = parse_default(&path).unwrap();

        assert_eq!(gene.name(), "gene");
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

        let gene =
            parse_with_samples(&path, &["cl7_bakta", "cl3_bakta"], ParalogMode::First).unwrap();

        assert!(gene.has_sample(0));
        assert!(gene.has_sample(1));
    }

    #[test]
    // Verifies duplicate normalized sample IDs keep the first sequence by default.
    fn parser_keeps_first_duplicate_normalized_sample_id_by_default() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nAAAA\n>sample2\nAAAA\n>sample1;second\nCCCC\n",
        );

        let gene = parse_default(&path).unwrap();

        assert_eq!(gene.paralog_count(), Some(1));
        assert_eq!(gene.snp_count(0, 1, false), (0, 4));
    }

    #[test]
    // Verifies skip mode marks duplicated samples for load-time filtering.
    fn parser_skip_mode_marks_duplicated_samples() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nACGT\n>sample2\nCCCC\n>sample1;second\nTGCA\n>sample3\nGGGG\n",
        );

        let gene = parse_with_samples(&path, &["sample1", "sample2", "sample3"], ParalogMode::Skip)
            .unwrap();

        assert_eq!(gene.paralog_count(), Some(1));
        assert!(gene.has_sample(0));
        assert!(gene.has_sample(1));
        assert!(gene.has_sample(2));
    }

    #[test]
    // Verifies longest mode keeps the duplicate with the most non-gap bytes.
    fn parser_longest_mode_keeps_duplicate_with_most_non_gap_bases() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;short\nAA--\n>sample2\nACGT\n>sample1;long\nACGT\n",
        );

        let gene =
            parse_with_samples(&path, &["sample1", "sample2"], ParalogMode::Longest).unwrap();

        assert_eq!(gene.paralog_count(), Some(1));
        assert_eq!(gene.snp_count(0, 1, true), (0, 4));
    }

    #[test]
    // Verifies longest mode keeps the first duplicate when non-gap counts tie.
    fn parser_longest_mode_keeps_first_duplicate_on_equal_non_gap_count() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(
            &dir,
            "gene.aln",
            ">sample1;first\nA--C\n>sample2\nA--C\n>sample1;second\nG--T\n",
        );

        let gene =
            parse_with_samples(&path, &["sample1", "sample2"], ParalogMode::Longest).unwrap();

        assert_eq!(gene.snp_count(0, 1, true), (0, 2));
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

        let error = parse_with_samples(&path, &["sample1"], ParalogMode::First).unwrap_err();
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
    // Verifies alignments cannot introduce samples absent from the Rtab header.
    fn parser_rejects_unknown_sample_names() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\nACGT\n>sample3\nACGT\n");

        let error = parse_default(&path).unwrap_err();

        assert!(
            error.to_string().contains("does not appear"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies unknown samples are rejected after Panaroo-style ID normalization.
    fn parser_rejects_unknown_normalized_sample_names() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\nACGT\n>sample3;copy\nACGT\n");

        let error = parse_default(&path).unwrap_err();
        let message = error.to_string();

        assert!(message.contains("sample 'sample3'"), "error: {message}");
        assert!(message.contains("does not appear"), "error: {message}");
    }

    #[test]
    // Verifies zero-length FASTA records are rejected.
    fn parser_rejects_zero_length_sequences() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", ">sample1\n\n");

        let error = parse_with_samples(&path, &["sample1"], ParalogMode::First).unwrap_err();

        assert!(
            error.to_string().contains("zero-length sequence"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies empty alignments are rejected.
    fn parser_rejects_empty_alignments() {
        let dir = TempDir::new().unwrap();
        let path = write_alignment(&dir, "gene.aln", "");

        let error = parse_with_samples(&path, &["sample1"], ParalogMode::First).unwrap_err();

        assert!(
            error.to_string().contains("contains no FASTA records"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies per-gene sample order is preserved independently.
    fn sample_order_can_differ_between_genes() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["s2", "s1"]);
        write_panaroo_alignment(&dir, "gene_a.aln.fas", ">s1\nAAAA\n>s2\nCCCC\n");
        write_panaroo_alignment(&dir, "gene_b.aln.fas", ">s2\nCCCC\n>s1\nAAAA\n");

        let loaded = load_genes(dir.path(), ParalogMode::First, None, true).unwrap();

        assert_eq!(loaded.sample_names, vec!["s2", "s1"]);
        assert_eq!(loaded.gene_sequences.len(), 2);
        assert!(loaded.gene_sequences[0].has_sample(1));
        assert!(loaded.gene_sequences[0].has_sample(0));
        assert!(loaded.gene_sequences[1].has_sample(0));
        assert!(loaded.gene_sequences[1].has_sample(1));
        assert_eq!(loaded.gene_sequences[0].snp_count(0, 1, false), (4, 4));
        assert_eq!(loaded.gene_sequences[1].snp_count(0, 1, false), (4, 4));
    }

    #[test]
    // Verifies genes may contain different subsets of samples.
    fn genes_may_have_different_sample_sets() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["s1", "s2", "s3"]);
        write_panaroo_alignment(
            &dir,
            "gene_a.aln.fas",
            ">s1;contig_a\nAAAA\n>s2;contig_a\nCCCC\n",
        );
        write_panaroo_alignment(
            &dir,
            "gene_b.aln.fas",
            ">s3;contig_b\nCCCC\n>s1;contig_b\nAAAA\n",
        );

        let loaded = load_genes(dir.path(), ParalogMode::First, None, true).unwrap();

        assert_eq!(loaded.sample_names, vec!["s1", "s2", "s3"]);
        assert_eq!(loaded.gene_sequences.len(), 2);
        assert!(loaded.gene_sequences[0].has_sample(0));
        assert!(loaded.gene_sequences[0].has_sample(1));
        assert!(!loaded.gene_sequences[0].has_sample(2));
        assert!(loaded.gene_sequences[1].has_sample(2));
        assert!(loaded.gene_sequences[1].has_sample(0));
        assert!(!loaded.gene_sequences[1].has_sample(1));
    }

    #[test]
    // Verifies loaded genes retain paralog metadata in sorted alignment order.
    fn load_genes_records_paralog_counts_in_alignment_order() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["s1", "s2"]);
        write_panaroo_alignment(
            &dir,
            "gene_a.aln.fas",
            ">s1;first\nAAAA\n>s2\nCCCC\n>s1;second\nTTTT\n",
        );
        write_panaroo_alignment(&dir, "gene_clean.aln.fas", ">s1\nAAAA\n>s2\nCCCC\n");
        write_panaroo_alignment(
            &dir,
            "gene_b.aln.fas",
            concat!(
                ">s1;first\nAAAA\n",
                ">s2;first\nCCCC\n",
                ">s1;second\nTTTT\n",
                ">s2;second\nGGGG\n",
            ),
        );

        let loaded = load_genes(dir.path(), ParalogMode::First, None, true).unwrap();
        let observed: Vec<_> = loaded
            .gene_sequences
            .iter()
            .filter_map(|gene| {
                gene.paralog_count()
                    .map(|paralog_count| (gene.name().to_owned(), paralog_count))
            })
            .collect();

        assert_eq!(
            observed,
            vec![("gene_a".to_string(), 1), ("gene_b".to_string(), 2)]
        );
    }

    #[test]
    // Verifies skip mode filters out genes with any paralog at load time.
    fn load_genes_skip_mode_filters_paralogous_genes() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["s1", "s2"]);
        write_panaroo_alignment(
            &dir,
            "gene_dup.aln.fas",
            ">s1;first\nAAAA\n>s2\nCCCC\n>s1;second\nTTTT\n",
        );
        write_panaroo_alignment(&dir, "gene_clean.aln.fas", ">s1\nAAAA\n>s2\nCCCC\n");

        let loaded = load_genes(dir.path(), ParalogMode::Skip, None, true).unwrap();

        assert_eq!(loaded.gene_sequences.len(), 1);
        assert_eq!(loaded.gene_sequences[0].name(), "gene_clean");
        assert_eq!(loaded.gene_sequences[0].paralog_count(), None);
    }

    #[test]
    // Verifies Rtab header parsing preserves sample column order.
    fn parse_rtab_header_preserves_sample_order() {
        let observed = parse_rtab_header(
            Path::new("gene_presence_absence.Rtab"),
            "Gene\tbeta\talpha\r\n",
        )
        .unwrap();

        assert_eq!(observed, vec!["beta", "alpha"]);
    }

    #[test]
    // Verifies Rtab headers must start with the Gene column.
    fn parse_rtab_header_requires_gene_first_column() {
        let error = parse_rtab_header(Path::new("gene_presence_absence.Rtab"), "gene\talpha\n")
            .unwrap_err();

        assert!(error.to_string().contains("must start with 'Gene'"));
    }

    #[test]
    // Verifies Rtab headers must contain at least one sample column.
    fn parse_rtab_header_requires_sample_columns() {
        let error =
            parse_rtab_header(Path::new("gene_presence_absence.Rtab"), "Gene\n").unwrap_err();

        assert!(error.to_string().contains("at least one sample"));
    }

    #[test]
    // Verifies duplicate Rtab sample names are rejected before loading.
    fn parse_rtab_header_rejects_duplicate_samples() {
        let error = parse_rtab_header(
            Path::new("gene_presence_absence.Rtab"),
            "Gene\talpha\talpha\n",
        )
        .unwrap_err();

        assert!(error.to_string().contains("duplicate sample name"));
    }

    #[test]
    // Verifies Panaroo discovery uses sorted aligned_gene_sequences .aln.fas files.
    fn read_panaroo_dir_uses_sorted_aligned_gene_sequences_files() {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("top_level.aln.fas"), "").unwrap();
        let alignment_dir = dir.path().join("aligned_gene_sequences");
        fs::create_dir(&alignment_dir).unwrap();
        let first = alignment_dir.join("alpha.aln.fas");
        let second = alignment_dir.join("zeta.aln.fas");
        fs::write(&second, "").unwrap();
        fs::write(alignment_dir.join("ignored.fas"), "").unwrap();
        fs::write(&first, "").unwrap();

        let observed = read_panaroo_dir(dir.path()).unwrap();

        assert_eq!(observed, vec![first, second]);
    }

    #[test]
    // Verifies Panaroo discovery rejects a missing aligned_gene_sequences directory.
    fn read_panaroo_dir_requires_aligned_gene_sequences() {
        let dir = TempDir::new().unwrap();
        fs::write(dir.path().join("gene.aln.fas"), "").unwrap();

        let error = read_panaroo_dir(dir.path()).unwrap_err();

        assert!(
            error.to_string().contains("aligned_gene_sequences"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies Panaroo discovery errors when the alignment directory is empty.
    fn read_panaroo_dir_rejects_empty_alignment_directory() {
        let dir = TempDir::new().unwrap();
        fs::create_dir(dir.path().join("aligned_gene_sequences")).unwrap();

        let error = read_panaroo_dir(dir.path()).unwrap_err();

        assert!(
            error
                .to_string()
                .contains("found no Panaroo alignment files ending in .aln.fas"),
            "error: {error}"
        );
    }

    #[test]
    // Verifies entropy filtering uses the same gene-name normalization as alignments.
    fn entropy_filter_matches_normalized_alignment_names() {
        let dir = TempDir::new().unwrap();
        let high = write_panaroo_alignment(&dir, "group_672.aln.fas", "");
        let equal = write_panaroo_alignment(&dir, "group_equal.aln.fas", "");
        fs::write(
            dir.path().join("alignment_entropy.csv"),
            "gene,entropy\ngroup_672.aln,0.51\ngroup_equal.aln,0.5\n",
        )
        .unwrap();

        let observed =
            filter_alignments_by_entropy(dir.path(), vec![high, equal.clone()], 0.5).unwrap();

        assert_eq!(observed, vec![equal]);
    }

    #[test]
    // Verifies alignments without entropy rows are retained.
    fn entropy_filter_keeps_alignments_without_metadata() {
        let dir = TempDir::new().unwrap();
        let known = write_panaroo_alignment(&dir, "known.aln.fas", "");
        let missing = write_panaroo_alignment(&dir, "missing.aln.fas", "");
        fs::write(
            dir.path().join("alignment_entropy.csv"),
            "gene,entropy\nknown,0.1\n",
        )
        .unwrap();

        let observed =
            filter_alignments_by_entropy(dir.path(), vec![known.clone(), missing.clone()], 0.5)
                .unwrap();

        assert_eq!(observed, vec![known, missing]);
    }

    #[test]
    // Verifies malformed entropy values are rejected.
    fn read_alignment_entropy_rejects_malformed_values() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("alignment_entropy.csv");
        fs::write(&path, "gene,entropy\ngene_a,not-a-number\n").unwrap();

        let error = read_alignment_entropy(&path).unwrap_err();

        assert!(error.to_string().contains("invalid entropy"));
    }

    #[test]
    // Verifies recombination tables are emitted as expected TSV.
    fn write_recombination_table_emits_tsv() {
        let sample_names = vec!["alpha".to_string(), "beta".to_string()];
        let genes = vec![
            Gene::new(
                "gene1".to_string(),
                1,
                HashMap::from([
                    (0, SampleBases::from_sequence(b"A")),
                    (1, SampleBases::from_sequence(b"A")),
                ]),
                0,
            ),
            Gene::new(
                "gene2".to_string(),
                1,
                HashMap::from([
                    (0, SampleBases::from_sequence(b"A")),
                    (1, SampleBases::from_sequence(b"A")),
                ]),
                0,
            ),
        ];
        let rows = vec![
            OutputRow {
                gene_index: 0,
                presence: vec![1, 0],
            },
            OutputRow {
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
