use crate::gene::GeneMetadata;
use anyhow::{Context, Result};
use hashbrown::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct OutputRow {
    pub(crate) gene_index: usize,
    pub(crate) presence: Vec<u8>,
}

// Orders rows by graph gene names and appends unmatched output genes by name.
pub(crate) fn order_rows_by_gene_names(
    rows: &[OutputRow],
    genes: &[GeneMetadata],
    ordered_gene_names: &[String],
) -> Result<Vec<OutputRow>> {
    let mut remaining = HashMap::with_capacity(rows.len());
    for (row_index, row) in rows.iter().enumerate() {
        let gene = genes.get(row.gene_index).with_context(|| {
            format!(
                "recombination table row references missing gene index {}",
                row.gene_index
            )
        })?;
        if remaining
            .insert(gene.name().to_owned(), row_index)
            .is_some()
        {
            anyhow::bail!(
                "multiple recombination rows have gene name '{}'",
                gene.name()
            );
        }
    }

    let mut ordered = Vec::with_capacity(rows.len());
    for gene_name in ordered_gene_names {
        if let Some(row_index) = remaining.remove(gene_name) {
            ordered.push(rows[row_index].clone());
        }
    }

    let mut unmatched: Vec<_> = remaining.into_iter().collect();
    unmatched.sort_unstable_by(|left, right| left.0.cmp(&right.0));
    ordered.extend(
        unmatched
            .into_iter()
            .map(|(_, row_index)| rows[row_index].clone()),
    );
    Ok(ordered)
}

// Writes a recombination presence table as tab-separated text.
pub fn write_recombination_table<W: Write>(
    sample_names: &[String],
    genes: &[GeneMetadata],
    rows: &[OutputRow],
    gene_sizes: Option<&HashMap<String, usize>>,
    mut writer: W,
) -> Result<()> {
    write!(writer, "gene")?;
    for sample in sample_names {
        write!(writer, "\t{sample}")?;
    }
    writeln!(writer)?;

    if let Some(gene_sizes) = gene_sizes {
        write!(writer, "#gene_frequency")?;
        for row in rows {
            let gene = genes.get(row.gene_index).with_context(|| {
                format!(
                    "recombination table row references missing gene index {}",
                    row.gene_index
                )
            })?;
            match gene_sizes.get(gene.name()) {
                Some(size) => write!(writer, "\t{size}/{}", sample_names.len())?,
                None => write!(writer, "\tNA")?,
            }
        }
        writeln!(writer)?;
    }

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
pub(crate) fn write_paralog_report(path: &Path, genes: &[GeneMetadata]) -> Result<usize> {
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

    writeln!(writer, "gene\tparalogs")
        .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    for (gene_name, paralog_count) in &paralog_rows {
        writeln!(writer, "{gene_name}\t{paralog_count}")
            .with_context(|| format!("failed to write paralog report '{}'", path.display()))?;
    }

    Ok(paralog_rows.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    // Verifies recombination tables are emitted as expected TSV.
    fn write_recombination_table_emits_tsv() {
        let sample_names = vec!["alpha".to_string(), "beta".to_string()];
        let genes = vec![
            GeneMetadata::new("gene1".to_string(), 0),
            GeneMetadata::new("gene2".to_string(), 0),
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

        write_recombination_table(&sample_names, &genes, &rows, None, &mut output).unwrap();

        assert_eq!(
            String::from_utf8(output).unwrap(),
            "gene\talpha\tbeta\ngene1\t1\t0\ngene2\t0\t1\n"
        );
    }

    #[test]
    fn orders_graph_matches_then_appends_unmatched_genes() {
        let genes = vec![
            GeneMetadata::new("gene_a".to_string(), 0),
            GeneMetadata::new("gene_b".to_string(), 0),
            GeneMetadata::new("gene_c".to_string(), 0),
        ];
        let rows = vec![
            OutputRow {
                gene_index: 0,
                presence: vec![1],
            },
            OutputRow {
                gene_index: 1,
                presence: vec![0],
            },
            OutputRow {
                gene_index: 2,
                presence: vec![1],
            },
        ];

        let ordered =
            order_rows_by_gene_names(&rows, &genes, &["gene_c".into(), "gene_a".into()]).unwrap();

        assert_eq!(
            ordered.iter().map(|row| row.gene_index).collect::<Vec<_>>(),
            vec![2, 0, 1]
        );
    }
}
