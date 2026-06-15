use crate::genome::Genome;
use crate::get_progress_bar;
use crate::model::select_recombinant_gene_indices;
use hashbrown::HashMap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

pub type PairHits = HashMap<usize, Vec<(usize, usize)>>;
type SamplePair = (usize, usize);

// Runs pairwise comparison across all sample pairs for loaded genes.
pub fn compare_loaded_alignments(
    sample_count: usize,
    genome: &Genome,
    gaps: bool,
    quiet: bool,
) -> PairHits {
    debug_assert_eq!(sample_count, genome.sample_count());

    // Flatten the upper triangle of the sample-by-sample matrix into one Rayon
    // range. That gives Rayon similarly sized units of work for each sample pair,
    // instead of parallelizing only by the outer sample index where early rows are
    // much larger than later rows.
    let sample_pair_count = sample_pair_count(sample_count);
    let progress_bar = get_progress_bar(sample_pair_count, true, quiet);
    let gene_pair_hits: Vec<_> = (0..sample_pair_count)
        .into_par_iter()
        .progress_with(progress_bar)
        .flat_map_iter(|pair_offset| {
            let (sample_a, sample_b) = sample_pair_indices(sample_count, pair_offset);
            selected_pair_hits(genome, sample_a, sample_b, gaps)
        })
        .collect();

    // Convert to a hash map keyed by global gene index.
    let mut hits: PairHits = HashMap::new();
    for (gene, pair) in gene_pair_hits {
        hits.entry(gene).or_default().push(pair);
    }

    hits
}

// Counts unique unordered sample pairs.
fn sample_pair_count(sample_count: usize) -> usize {
    sample_count * sample_count.saturating_sub(1) / 2
}

// Maps a flat upper-triangle offset to the corresponding sample indices.
fn sample_pair_indices(sample_count: usize, pair_offset: usize) -> (usize, usize) {
    debug_assert!(sample_count >= 2);
    debug_assert!(pair_offset < sample_pair_count(sample_count));

    // Find the row in the upper triangle that owns this flat offset. Row i
    // contains pairs (i, i + 1)..(i, sample_count - 1), so pairs_before_sample()
    // is the starting flat offset for row i.
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

// Counts pair offsets before a sample's row in the upper triangle.
fn pairs_before_sample(sample_count: usize, sample_index: usize) -> usize {
    // Number of upper-triangle pairs in rows before sample_index:
    // (sample_count - 1) + (sample_count - 2) + ... + (sample_count - sample_index).
    sample_index * (2 * sample_count - sample_index - 1) / 2
}

// Selects recombinant genes for one sample pair and tags them with that pair.
fn selected_pair_hits(
    genome: &Genome,
    sample_a: usize,
    sample_b: usize,
    gaps: bool,
) -> Vec<(usize, SamplePair)> {
    let pair_genes = genome.gene_snp_counts(sample_a, sample_b, gaps);
    let recombinant_gene_indices = select_recombinant_gene_indices(pair_genes);

    recombinant_gene_indices
        .into_iter()
        .map(|gene_index| (gene_index, (sample_a, sample_b)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::ParalogMode;
    use crate::gene::GeneMetadata;
    use crate::genome::load_genes;
    use crate::output::OutputRow;
    use crate::presence_table_from_pair_hits;
    use std::fs;
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    fn compare_alignments(panaroo_dir: &Path) -> (Vec<String>, Vec<GeneMetadata>, PairHits) {
        let (sample_names, genome) =
            load_genes(panaroo_dir, ParalogMode::First, None, true).expect("Test gene load failed");
        let gene_metadata = genome.gene_metadata().clone();
        let hits = compare_loaded_alignments(sample_names.len(), &genome, false, true);
        (sample_names, gene_metadata, hits)
    }

    // Normalizes Rayon-collected hit order for tests that compare pair vectors.
    fn sort_pair_hits(hits: &mut PairHits) {
        for pairs in hits.values_mut() {
            pairs.sort();
        }
    }

    fn infer_recombination_presence(panaroo_dir: &Path) -> Vec<OutputRow> {
        let (sample_names, genes, hits) = compare_alignments(panaroo_dir);
        presence_table_from_pair_hits(sample_names.len(), genes.len(), &hits, true)
    }

    // Writes the required Panaroo Rtab header fixture.
    fn write_rtab(dir: &TempDir, sample_names: &[&str]) {
        fs::write(
            dir.path().join("gene_presence_absence.Rtab"),
            format!("Gene\t{}\n", sample_names.join("\t")),
        )
        .unwrap();
    }

    // Writes a temporary FASTA alignment for tests.
    fn write_alignment(dir: &TempDir, name: &str, contents: &str) -> PathBuf {
        let alignment_dir = dir.path().join("aligned_gene_sequences");
        fs::create_dir_all(&alignment_dir).unwrap();
        let path = alignment_dir.join(name);
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    // Verifies flat pair offsets follow sample-pair lexicographic order.
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
    // Verifies thresholding returns the expected high-divergence gene hits.
    fn bayesian_threshold_returns_expected_gene_keys_and_pairs() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["s1", "s2"]);

        for mismatches in 0..12 {
            let s2 = format!("{}{}", "C".repeat(mismatches), "A".repeat(12 - mismatches));
            let contents = format!(">s1\n{}\n>s2\n{}\n", "A".repeat(12), s2);
            write_alignment(&dir, &format!("gene{mismatches:02}.aln.fas"), &contents);
        }

        let mut hits = compare_alignments(dir.path());
        sort_pair_hits(&mut hits.2);
        let mut observed: Vec<_> = hits.2.keys().copied().collect();
        observed.sort();

        assert_eq!(observed, vec![8, 9, 10, 11]);

        for pairs in hits.2.values() {
            assert_eq!(pairs, &vec![(0, 1)]);
        }
    }

    #[test]
    // Verifies missing samples limit comparisons to comparable gene pairs.
    fn variable_sample_genes_accumulate_only_comparable_pairs() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["alpha", "beta", "gamma"]);

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

            write_alignment(&dir, &format!("gene{mismatches:02}.aln.fas"), &contents);
        }

        let mut hits = compare_alignments(dir.path());
        sort_pair_hits(&mut hits.2);
        let mut observed: Vec<_> = hits.2.keys().copied().collect();
        observed.sort();

        assert_eq!(observed, vec![8, 9, 10, 11]);

        for pairs in hits.2.values() {
            assert_eq!(pairs, &vec![(0, 1)]);
        }
    }

    #[test]
    // Verifies pair statistics skip genes missing either requested sample.
    fn comparable_pair_stats_skip_genes_missing_either_sample() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["alpha", "beta", "gamma"]);
        write_alignment(&dir, "gene_ab.aln.fas", ">alpha\nAAAA\n>beta\nCCCC\n");
        write_alignment(&dir, "gene_ag.aln.fas", ">alpha\nAAAA\n>gamma\nTTTT\n");
        write_alignment(&dir, "gene_bg.aln.fas", ">beta\nCCCC\n>gamma\nTTTT\n");
        let (_sample_names, genome) =
            load_genes(dir.path(), ParalogMode::First, None, true).unwrap();

        let observed: Vec<_> = genome
            .gene_snp_counts(0, 1, false)
            .into_iter()
            .map(|stats| stats.gene_index)
            .collect();

        assert_eq!(observed, vec![0]);
    }

    #[test]
    // Verifies zero-length effective comparisons do not reach the model fit.
    fn comparable_pair_stats_skip_zero_length_alignments() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["alpha", "beta"]);
        write_alignment(&dir, "gene_zero.aln.fas", ">alpha\n----\n>beta\nAAAA\n");
        write_alignment(&dir, "gene_positive.aln.fas", ">alpha\nAAAA\n>beta\nAAAT\n");
        let (sample_names, genome) =
            load_genes(dir.path(), ParalogMode::First, None, true).unwrap();

        let observed: Vec<_> = genome
            .gene_snp_counts(0, 1, false)
            .into_iter()
            .map(|stats| (stats.gene_index, stats.snps, stats.length))
            .collect();

        assert_eq!(observed, vec![(0, 1, 4)]);

        let hits = compare_loaded_alignments(sample_names.len(), &genome, false, true);

        assert!(!hits.contains_key(&1));
    }

    #[test]
    // Verifies presence output keeps sorted alignment order.
    fn presence_table_keeps_all_loaded_genes_in_order() {
        let dir = TempDir::new().unwrap();
        write_rtab(&dir, &["alpha", "beta"]);
        write_alignment(
            &dir,
            "gene_low.aln.fas",
            ">beta\nAAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
        );
        write_alignment(
            &dir,
            "gene_middle.aln.fas",
            ">beta\nCAAAAAAAAA\n>alpha\nAAAAAAAAAA\n",
        );
        write_alignment(
            &dir,
            "gene_high.aln.fas",
            ">beta\nCCCCCCCCAA\n>alpha\nAAAAAAAAAA\n",
        );

        let rows = infer_recombination_presence(dir.path());

        assert_eq!(
            rows,
            vec![
                OutputRow {
                    gene_index: 0,
                    presence: vec![0, 0],
                },
                OutputRow {
                    gene_index: 1,
                    presence: vec![0, 0],
                },
                OutputRow {
                    gene_index: 2,
                    presence: vec![0, 0],
                },
            ]
        );
    }
}
