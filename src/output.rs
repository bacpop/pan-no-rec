use crate::gene::Gene;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::Write;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct OutputRow {
    pub(crate) gene_index: usize,
    pub(crate) presence: Vec<u8>,
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
    use crate::gene::SampleBases;
    use hashbrown::HashMap;

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
