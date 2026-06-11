use crate::gene::{GeneMetadata, SampleBases, ParsedGeneAlignment};
use crate::model::PairGeneStats;
use crate::cli::ParalogMode;
use crate::panaroo_io::*;
use crate::get_progress_bar;

use std::path::Path;
use anyhow::{Result, bail};
use roaring::RoaringBitmap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
#[derive(Debug)]
pub(crate) struct Genome {
    sample_bases: Vec<SampleBases>,
    gene_masks: Vec<RoaringBitmap>,
    gene_lengths: Vec<usize>,
    sample_gene_presence: Vec<RoaringBitmap>,
    gene_sort_ranks: Vec<usize>,
    gene_metadata: Vec<GeneMetadata>,
}

impl Genome {
    pub(crate) fn new(sample_count: usize, gene_count: usize) -> Self {
        Genome {
            sample_bases: vec![SampleBases::default(); sample_count],
            gene_masks: vec![RoaringBitmap::new(); gene_count],
            gene_lengths: vec![0; gene_count],
            sample_gene_presence: vec![RoaringBitmap::new(); sample_count],
            gene_sort_ranks: vec![0; sample_count],
            gene_metadata: vec![GeneMetadata::default(); gene_count],
        }
    }

    pub(crate) fn add_alignment(&mut self, alignment: ParsedGeneAlignment) {
        let ParsedGeneAlignment {
            gene_index,
            metadata,
            alignment_len,
            offset,
            sequences,
        } = alignment;

        for (sample_index, bases) in sequences {
            let Some(target) = self.sample_bases.get_mut(sample_index) else {
                panic!(
                    "gene index {gene_index} references sample index {sample_index}, but only {} samples were loaded",
                    self.sample_bases.len()
                );
            };
            // Concats the bitvecs by ORing them
            target.union_assign(bases);
            self.sample_gene_presence[sample_index].insert(gene_index as u32);
        }

        let mut mask = RoaringBitmap::new();
        mask.insert_range(offset..(offset + alignment_len as u32));
        self.gene_masks[gene_index] = mask;
        self.gene_lengths[gene_index] = alignment_len;
        self.gene_metadata[gene_index] = metadata;
    }

    pub(crate) fn merge(mut self, other: Self) -> Self {
        debug_assert_eq!(self.sample_bases.len(), other.sample_bases.len());
        debug_assert_eq!(self.gene_lengths.len(), other.gene_lengths.len());

        // Combine all bitvecs
        for (left, right) in self.sample_bases.iter_mut().zip(other.sample_bases) {
            left.union_assign(right);
        }
        // Combine pres/abs
        for (left, right) in self
            .sample_gene_presence
            .iter_mut()
            .zip(other.sample_gene_presence)
        {
            *left |= &right;
        }

        // Combine the Vecs
        for (gene_index, ((other_length, other_mask), other_metadata)) in other
            .gene_lengths
            .into_iter()
            .zip(other.gene_masks)
            .zip(other.gene_metadata)
            .enumerate()
        {
            self.gene_lengths[gene_index] = other_length;
            self.gene_masks[gene_index] = other_mask;
            self.gene_metadata[gene_index] = other_metadata;
        }

        self
    }

    pub(crate) fn finish(&mut self) {
        self.gene_sort_ranks = gene_sort_ranks(&self.gene_metadata);
    }

    pub(crate) fn sample_count(&self) -> usize {
        self.sample_bases.len()
    }

    pub(crate) fn gene_snp_counts(
        &self,
        sample_a: usize,
        sample_b: usize,
        gaps: bool,
    ) -> Vec<PairGeneStats> {
        let left = &self.sample_bases[sample_a];
        let right = &self.sample_bases[sample_b];
        let matches = left.matching_sites(right);
        let both_gap = left.both_gap_sites(right);
        let either_gap = (!gaps).then(|| left.either_gap_sites(right));
        let comparable_genes =
            &self.sample_gene_presence[sample_a] & &self.sample_gene_presence[sample_b];

        // Iterate over distances using gene-wise maps
        comparable_genes
            .into_iter()
            .filter_map(|gene_index| {
                let gene_index = gene_index as usize;
                let mask = &self.gene_masks[gene_index];
                let excluded_gap_count = if gaps {
                    both_gap.intersection_len(mask)
                } else {
                    either_gap
                        .as_ref()
                        .expect("gap-excluding mode should build an either-gap bitmap")
                        .intersection_len(mask)
                } as usize;
                let length = self.gene_lengths[gene_index] - excluded_gap_count;
                if length == 0 {
                    return None;
                }

                let matching_count = matches.intersection_len(mask) as usize;
                Some(PairGeneStats {
                    gene_index,
                    gene_sort_rank: self.gene_sort_ranks[gene_index],
                    snps: length - matching_count,
                    length,
                })
            })
            .collect()
    }

    pub(crate) fn get_summary(&self) -> Result<(usize, usize)> {
        if self.gene_metadata.is_empty() {
            bail!("No valid genes loaded");
        } else if self.sample_bases.is_empty() {
            bail!("Alignments are empty");
        }
        Ok((self.gene_metadata.len(), self.sample_bases.len()))
    }

    pub(crate) fn gene_metadata(&self) -> &Vec<GeneMetadata> {
        &self.gene_metadata
    }

}

fn gene_sort_ranks(genes: &[GeneMetadata]) -> Vec<usize> {
    let mut ordered_indices: Vec<_> = (0..genes.len()).collect();
    ordered_indices.sort_by(|&left, &right| {
        genes[left]
            .name()
            .cmp(genes[right].name())
            .then_with(|| left.cmp(&right))
    });

    let mut ranks = vec![0; genes.len()];
    for (rank, gene_index) in ordered_indices.into_iter().enumerate() {
        ranks[gene_index] = rank;
    }
    ranks
}

// Loads all Panaroo alignments into a concatenated genome using the Rtab sample order.
pub(crate) fn load_genes(
    panaroo_dir: &Path,
    paralog_mode: ParalogMode,
    max_entropy: Option<f64>,
    quiet: bool,
) -> Result<(Vec<String>, Genome)> {
    let rtab_path = panaroo_dir.join("gene_presence_absence.Rtab");
    let sample_names = read_rtab_sample_names(&rtab_path)?;
    let sample_indices = build_sample_indices(&sample_names, &rtab_path)?;
    let mut aln_paths = read_panaroo_dir(panaroo_dir)?;
    if let Some(threshold) = max_entropy {
        aln_paths = filter_alignments_by_entropy(panaroo_dir, aln_paths, threshold)?;
    }

    let alignment_lengths = read_alignment_lengths(&aln_paths)?;
    let alignment_offsets = alignment_offsets(&alignment_lengths)?;
    let gene_count = aln_paths.len();
    let pbar = get_progress_bar(aln_paths.len(), false, quiet);
    let mut genome_accumulator = aln_paths
        .par_iter()
        .enumerate()
        .progress_with(pbar)
        .try_fold(
            || Genome::new(sample_names.len(), gene_count),
            |mut accumulator, (gene_index, aln)| {
                let parsed = parse_gene_alignment(
                    aln,
                    gene_index,
                    alignment_offsets[gene_index],
                    alignment_lengths[gene_index],
                    &sample_indices,
                    paralog_mode,
                )?;
                accumulator.add_alignment(parsed);
                Ok::<_, anyhow::Error>(accumulator)
            },
        )
        .try_reduce(
            || Genome::new(sample_names.len(), gene_count),
            |left, right| Ok(left.merge(right)),
        )?;

    genome_accumulator.finish();
    Ok((sample_names, genome_accumulator))
}

#[cfg(test)]
mod tests {
    use super::*;
    use hashbrown::HashMap;

    fn parsed_gene(
        gene_index: usize,
        name: &str,
        alignment_len: usize,
        offset: u32,
        sequences: &[(usize, &[u8])],
    ) -> ParsedGeneAlignment {
        ParsedGeneAlignment::new(
            gene_index,
            GeneMetadata::new(name.to_string(), 0),
            alignment_len,
            offset,
            sequences
                .iter()
                .map(|(sample, sequence)| {
                    (*sample, SampleBases::from_sequence_at(sequence, offset))
                })
                .collect(),
        )
    }

    fn build_genome(
        sample_count: usize,
        gene_count: usize,
        alignments: Vec<ParsedGeneAlignment>,
    ) -> Genome {
        let mut accumulator = GenomeAccumulator::new(sample_count, gene_count);
        for alignment in alignments {
            accumulator.add_alignment(alignment).unwrap();
        }
        accumulator.finish().unwrap().0
    }

    fn stats_tuples(stats: Vec<PairGeneStats>) -> Vec<(usize, usize, usize)> {
        stats
            .into_iter()
            .map(|stats| (stats.gene_index, stats.snps, stats.length))
            .collect()
    }

    #[test]
    fn concatenated_counts_match_per_gene_gap_excluding_mode() {
        let genome = build_genome(
            3,
            4,
            vec![
                parsed_gene(0, "zeta", 10, 0, &[(0, b"A--CGGTTT-"), (1, b"ACCCTG----")]),
                parsed_gene(
                    1,
                    "alpha",
                    10,
                    10,
                    &[(0, b"MRWSYKVHDB"), (1, b"AATCTGCATG")],
                ),
                parsed_gene(2, "missing", 4, 20, &[(0, b"AAAA"), (2, b"CCCC")]),
                parsed_gene(3, "zero", 4, 24, &[(0, b"----"), (1, b"AAAA")]),
            ],
        );

        let observed = stats_tuples(genome.gene_snp_counts(0, 1, false));

        assert_eq!(observed, vec![(0, 1, 4), (1, 0, 10)]);
    }

    #[test]
    fn concatenated_counts_match_per_gene_gap_including_mode() {
        let genome = build_genome(
            3,
            4,
            vec![
                parsed_gene(0, "zeta", 10, 0, &[(0, b"A--CGGTTT-"), (1, b"ACCCTG----")]),
                parsed_gene(
                    1,
                    "alpha",
                    10,
                    10,
                    &[(0, b"MRWSYKVHDB"), (1, b"AATCTGCATG")],
                ),
                parsed_gene(2, "missing", 4, 20, &[(0, b"AAAA"), (2, b"CCCC")]),
                parsed_gene(3, "zero", 4, 24, &[(0, b"----"), (1, b"AAAA")]),
            ],
        );

        let observed = stats_tuples(genome.gene_snp_counts(0, 1, true));

        assert_eq!(observed, vec![(0, 6, 9), (1, 0, 10), (3, 4, 4)]);
    }

    #[test]
    fn gene_sort_ranks_follow_gene_names() {
        let genome = build_genome(
            2,
            2,
            vec![
                parsed_gene(0, "zeta", 1, 0, &[(0, b"A"), (1, b"A")]),
                parsed_gene(1, "alpha", 1, 1, &[(0, b"A"), (1, b"A")]),
            ],
        );

        let observed: Vec<_> = genome
            .gene_snp_counts(0, 1, false)
            .into_iter()
            .map(|stats| (stats.gene_index, stats.gene_sort_rank))
            .collect();

        assert_eq!(observed, vec![(0, 1), (1, 0)]);
    }

    #[test]
    fn accumulator_merge_combines_parallel_chunks() {
        let mut left = GenomeAccumulator::new(2, 2);
        left.add_alignment(parsed_gene(0, "zeta", 3, 0, &[(0, b"AAA"), (1, b"AAA")]))
            .unwrap();
        let mut right = GenomeAccumulator::new(2, 2);
        right
            .add_alignment(parsed_gene(1, "alpha", 1, 3, &[(0, b"C"), (1, b"G")]))
            .unwrap();
        let (genome, metadata) = left.merge(right).unwrap().finish().unwrap();

        let observed = stats_tuples(genome.gene_snp_counts(0, 1, false));

        assert_eq!(metadata[0].name(), "zeta");
        assert_eq!(metadata[1].name(), "alpha");
        assert_eq!(observed, vec![(0, 0, 3), (1, 1, 1)]);
    }

    #[test]
    fn accumulator_rejects_out_of_range_sample_indices() {
        let mut accumulator = GenomeAccumulator::new(2, 1);
        let alignment = ParsedGeneAlignment::new(
            0,
            GeneMetadata::new("gene".to_string(), 0),
            1,
            0,
            HashMap::from([(2, SampleBases::from_sequence(b"A"))]),
        );

        let error = accumulator.add_alignment(alignment).unwrap_err();

        assert!(
            error.to_string().contains("sample index 2"),
            "error: {error}"
        );
    }

    #[test]
    fn accumulator_rejects_alignment_past_bitmap_limit() {
        let mut accumulator = GenomeAccumulator::new(0, 1);
        let alignment = ParsedGeneAlignment::new(
            0,
            GeneMetadata::new("large".to_string(), 0),
            1,
            u32::MAX,
            HashMap::new(),
        );

        let error = accumulator.add_alignment(alignment).unwrap_err();

        assert!(error.to_string().contains("exceeding"), "error: {error}");
    }
}
