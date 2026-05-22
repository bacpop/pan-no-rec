mod error;
mod gene;
mod io;
mod model;

use crate::gene::Gene;
use crate::io::load_genes;
use crate::model::{PairGeneStats, select_recombinant_gene_indices};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::ThreadPoolBuilder;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

pub use crate::error::CompareError;

pub type PairHits = HashMap<String, Vec<(String, String)>>;
type SamplePair = (String, String);

pub fn compare_alignments<P>(aln_paths: &[P], threads: usize) -> Result<PairHits, CompareError>
where
    P: AsRef<Path>,
{
    validate_threads(threads)?;

    let (sample_names, genes) = load_genes(aln_paths)?;

    if genes.is_empty() {
        return Ok(HashMap::new());
    }

    run_with_thread_pool(threads, || compare_loaded_alignments(&sample_names, &genes))
}

fn compare_loaded_alignments(sample_names: &[String], genes: &[Gene]) -> PairHits {
    let sample_pair_count = sample_pair_count(sample_names.len());
    let gene_pair_hits: Vec<_> = (0..sample_pair_count)
        .into_par_iter()
        .flat_map_iter(|pair_offset| {
            let (sample_a, sample_b) = sample_pair_indices(sample_names.len(), pair_offset);
            selected_pair_hits(genes, &sample_names[sample_a], &sample_names[sample_b])
        })
        .collect();

    let mut hits: PairHits = HashMap::new();
    for (gene, pair) in gene_pair_hits {
        hits.entry(gene).or_default().push(pair);
    }

    for pairs in hits.values_mut() {
        pairs.sort();
    }

    hits
}

fn validate_threads(threads: usize) -> Result<(), CompareError> {
    if threads == 0 {
        return Err(CompareError::InvalidThreadCount { threads });
    }

    Ok(())
}

fn run_with_thread_pool<T, F>(threads: usize, work: F) -> Result<T, CompareError>
where
    T: Send,
    F: FnOnce() -> T + Send,
{
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|source| CompareError::ThreadPoolBuild { threads, source })?;

    Ok(pool.install(work))
}

fn sample_pair_count(sample_count: usize) -> usize {
    sample_count * sample_count.saturating_sub(1) / 2
}

fn sample_pair_indices(sample_count: usize, pair_offset: usize) -> (usize, usize) {
    debug_assert!(sample_count >= 2);
    debug_assert!(pair_offset < sample_pair_count(sample_count));

    let mut low = 0;
    let mut high = sample_count - 2;
    while low < high {
        let midpoint = (low + high).div_ceil(2);
        if pairs_before_sample(sample_count, midpoint) <= pair_offset {
            low = midpoint;
        } else {
            high = midpoint - 1;
        }
    }

    let sample_a = low;
    let sample_b = sample_a + 1 + pair_offset - pairs_before_sample(sample_count, sample_a);
    (sample_a, sample_b)
}

fn pairs_before_sample(sample_count: usize, sample_index: usize) -> usize {
    sample_index * (2 * sample_count - sample_index - 1) / 2
}

fn selected_pair_hits(genes: &[Gene], sample_a: &str, sample_b: &str) -> Vec<(String, SamplePair)> {
    let pair_genes = collect_comparable_pair_gene_stats(genes, sample_a, sample_b);
    let recombinant_gene_indices = select_recombinant_gene_indices(pair_genes);
    let pair = (sample_a.to_owned(), sample_b.to_owned());

    recombinant_gene_indices
        .into_iter()
        .map(|gene_index| (genes[gene_index].name().to_owned(), pair.clone()))
        .collect()
}

fn collect_comparable_pair_gene_stats<'a>(
    genes: &'a [Gene],
    sample_a: &str,
    sample_b: &str,
) -> Vec<PairGeneStats<'a>> {
    genes
        .iter()
        .enumerate()
        .filter_map(|(gene_index, gene)| {
            let sample_a_index = gene.sample_index(sample_a)?;
            let sample_b_index = gene.sample_index(sample_b)?;

            Some(PairGeneStats {
                gene_index,
                gene_id: gene.name(),
                snps: gene.snp_count(sample_a_index, sample_b_index),
                length: gene.alignment_len(),
            })
        })
        .collect()
}

#[pyfunction(name = "compare_alignments", signature = (aln_paths, threads=1))]
fn compare_alignments_py(
    py: Python<'_>,
    aln_paths: Vec<String>,
    threads: usize,
) -> PyResult<PairHits> {
    let paths: Vec<_> = aln_paths.into_iter().map(PathBuf::from).collect();
    py.detach(|| {
        compare_alignments(&paths, threads)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    })
}

#[pymodule]
fn pangenome_recombination(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compare_alignments_py, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::{PyDict, PyModule};
    use std::fs;
    use tempfile::TempDir;

    fn write_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let path = dir.path().join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn sample_pair_offsets_iterate_global_pairs_in_lexicographic_order() {
        let observed: Vec<_> = (0..sample_pair_count(4))
            .map(|offset| sample_pair_indices(4, offset))
            .collect();

        assert_eq!(
            observed,
            vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        );
    }

    #[test]
    fn bayesian_threshold_returns_expected_gene_keys_and_pairs() {
        let dir = TempDir::new().unwrap();
        let mut paths = Vec::new();

        for mismatches in 0..12 {
            let s2 = format!("{}{}", "C".repeat(mismatches), "A".repeat(12 - mismatches));
            let contents = format!(">s1\n{}\n>s2\n{}\n", "A".repeat(12), s2);
            paths.push(write_alignment(
                &dir,
                &format!("gene{mismatches:02}.aln"),
                &contents,
            ));
        }

        let hits = compare_alignments(&paths, 1).unwrap();
        let mut observed: Vec<_> = hits.keys().cloned().collect();
        observed.sort();

        let expected: Vec<_> = (8..12).map(|index| format!("gene{index:02}")).collect();
        assert_eq!(observed, expected);

        for pairs in hits.values() {
            assert_eq!(pairs, &vec![("s1".to_string(), "s2".to_string())]);
        }
    }

    #[test]
    fn variable_sample_genes_accumulate_only_comparable_pairs() {
        let dir = TempDir::new().unwrap();
        let mut paths = Vec::new();

        for mismatches in 0..12 {
            let beta = format!("{}{}", "C".repeat(mismatches), "A".repeat(12 - mismatches));
            let contents = if mismatches == 0 {
                format!(
                    ">gamma;extra\n{}\n>beta;extra\n{}\n>alpha;extra\n{}\n",
                    "A".repeat(12),
                    beta,
                    "A".repeat(12)
                )
            } else if mismatches % 2 == 0 {
                format!(">beta;extra\n{}\n>alpha;extra\n{}\n", beta, "A".repeat(12))
            } else {
                format!(">alpha;extra\n{}\n>beta;extra\n{}\n", "A".repeat(12), beta)
            };

            paths.push(write_alignment(
                &dir,
                &format!("gene{mismatches:02}.aln"),
                &contents,
            ));
        }

        let hits = compare_alignments(&paths, 2).unwrap();
        let mut observed: Vec<_> = hits.keys().cloned().collect();
        observed.sort();

        let expected: Vec<_> = (8..12).map(|index| format!("gene{index:02}")).collect();
        assert_eq!(observed, expected);

        for pairs in hits.values() {
            assert_eq!(pairs, &vec![("alpha".to_string(), "beta".to_string())]);
        }
    }

    #[test]
    fn comparable_pair_stats_skip_genes_missing_either_sample() {
        let dir = TempDir::new().unwrap();
        let gene_ab = write_alignment(&dir, "gene_ab.aln", ">alpha\nAAAA\n>beta\nCCCC\n");
        let gene_ag = write_alignment(&dir, "gene_ag.aln", ">alpha\nAAAA\n>gamma\nTTTT\n");
        let gene_bg = write_alignment(&dir, "gene_bg.aln", ">beta\nCCCC\n>gamma\nTTTT\n");
        let (_sample_names, genes) = crate::io::load_genes(&[gene_ab, gene_ag, gene_bg]).unwrap();

        let observed: Vec<_> = collect_comparable_pair_gene_stats(&genes, "alpha", "beta")
            .into_iter()
            .map(|stats| stats.gene_id)
            .collect();

        assert_eq!(observed, vec!["gene_ab"]);
    }

    #[test]
    fn python_api_smoke_test_returns_dict() {
        let dir = TempDir::new().unwrap();
        let low = write_alignment(&dir, "gene_low.aln", ">s1\nAAAAAAAAAA\n>s2\nAAAAAAAAAA\n");
        let middle = write_alignment(
            &dir,
            "gene_middle.aln",
            ">s1\nAAAAAAAAAA\n>s2\nCAAAAAAAAA\n",
        );
        let high = write_alignment(&dir, "gene_high.aln", ">s1\nAAAAAAAAAA\n>s2\nCCCCCCCCAA\n");
        let paths: Vec<_> = [low, middle, high]
            .into_iter()
            .map(|path| path.to_string_lossy().to_string())
            .collect();

        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "pangenome_recombination").unwrap();
            pangenome_recombination(&module).unwrap();

            let result = module
                .getattr("compare_alignments")
                .unwrap()
                .call1((paths, 2usize))
                .unwrap();

            assert!(result.is_instance_of::<PyDict>());
            let dict = result.cast::<PyDict>().unwrap();
            assert!(dict.contains("gene_high").unwrap());
        });
    }

    #[test]
    fn compare_alignments_rejects_zero_threads() {
        let error = compare_alignments::<PathBuf>(&[], 0).unwrap_err();

        assert!(matches!(
            error,
            CompareError::InvalidThreadCount { threads: 0 }
        ));
    }

    #[test]
    fn rayon_pool_uses_requested_thread_count() {
        let threads = run_with_thread_pool(2, rayon::current_num_threads).unwrap();

        assert_eq!(threads, 2);
    }
}
