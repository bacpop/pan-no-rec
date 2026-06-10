use crate::gene::{Gene, GeneMetadata, SampleBases};
use crate::model::PairGeneStats;
use anyhow::{Result, bail};
use roaring::RoaringBitmap;

#[derive(Debug)]
pub(crate) struct Genome {
    sample_bases: Vec<SampleBases>,
    gene_masks: Vec<RoaringBitmap>,
    gene_lengths: Vec<usize>,
    sample_gene_presence: Vec<RoaringBitmap>,
    gene_sort_ranks: Vec<usize>,
}

impl Genome {
    pub(crate) fn try_from_genes(
        sample_count: usize,
        genes: Vec<Gene>,
    ) -> Result<(Self, Vec<GeneMetadata>)> {
        validate_total_alignment_len(&genes)?;

        let mut sample_bases = vec![SampleBases::default(); sample_count];
        let mut sample_gene_presence = vec![RoaringBitmap::new(); sample_count];
        let mut gene_masks = Vec::with_capacity(genes.len());
        let mut gene_lengths = Vec::with_capacity(genes.len());
        let mut gene_metadata = Vec::with_capacity(genes.len());
        let mut offset = 0u32;

        for (gene_index, gene) in genes.into_iter().enumerate() {
            let gene_index = u32::try_from(gene_index)
                .map_err(|_| anyhow::anyhow!("loaded gene count exceeds the u32 index limit"))?;
            let (metadata, alignment_len, sequences) = gene.into_parts();
            let next_offset = offset as usize + alignment_len;
            let end = next_offset as u32;

            let mut mask = RoaringBitmap::new();
            mask.insert_range(offset..end);
            gene_masks.push(mask);
            gene_lengths.push(alignment_len);
            gene_metadata.push(metadata);

            for (sample_index, bases) in sequences {
                let Some(target) = sample_bases.get_mut(sample_index) else {
                    bail!(
                        "gene references sample index {sample_index}, but only {sample_count} samples were loaded"
                    );
                };
                target.append_shifted(bases, offset);
                sample_gene_presence[sample_index].insert(gene_index);
            }

            offset = end;
        }

        let gene_sort_ranks = gene_sort_ranks(&gene_metadata);
        let genome = Genome {
            sample_bases,
            gene_masks,
            gene_lengths,
            sample_gene_presence,
            gene_sort_ranks,
        };

        Ok((genome, gene_metadata))
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
}

fn validate_total_alignment_len(genes: &[Gene]) -> Result<()> {
    let mut total_len = 0usize;
    for gene in genes {
        total_len = total_len.checked_add(gene.alignment_len()).ok_or_else(|| {
            anyhow::anyhow!("total concatenated alignment length overflows usize")
        })?;
    }

    if total_len > u32::MAX as usize {
        bail!(
            "total concatenated alignment length {total_len} exceeds the {} column bitmap limit",
            u32::MAX
        );
    }

    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;
    use hashbrown::HashMap;

    fn gene(name: &str, alignment_len: usize, sequences: &[(usize, &[u8])]) -> Gene {
        Gene::new(
            name.to_string(),
            alignment_len,
            sequences
                .iter()
                .map(|(sample, sequence)| (*sample, SampleBases::from_sequence(sequence)))
                .collect(),
            0,
        )
    }

    fn stats_tuples(stats: Vec<PairGeneStats>) -> Vec<(usize, usize, usize)> {
        stats
            .into_iter()
            .map(|stats| (stats.gene_index, stats.snps, stats.length))
            .collect()
    }

    #[test]
    fn concatenated_counts_match_per_gene_gap_excluding_mode() {
        let genes = vec![
            gene("zeta", 10, &[(0, b"A--CGGTTT-"), (1, b"ACCCTG----")]),
            gene("alpha", 10, &[(0, b"MRWSYKVHDB"), (1, b"AATCTGCATG")]),
            gene("missing", 4, &[(0, b"AAAA"), (2, b"CCCC")]),
            gene("zero", 4, &[(0, b"----"), (1, b"AAAA")]),
        ];
        let (genome, _) = Genome::try_from_genes(3, genes).unwrap();

        let observed = stats_tuples(genome.gene_snp_counts(0, 1, false));

        assert_eq!(observed, vec![(0, 1, 4), (1, 0, 10)]);
    }

    #[test]
    fn concatenated_counts_match_per_gene_gap_including_mode() {
        let genes = vec![
            gene("zeta", 10, &[(0, b"A--CGGTTT-"), (1, b"ACCCTG----")]),
            gene("alpha", 10, &[(0, b"MRWSYKVHDB"), (1, b"AATCTGCATG")]),
            gene("missing", 4, &[(0, b"AAAA"), (2, b"CCCC")]),
            gene("zero", 4, &[(0, b"----"), (1, b"AAAA")]),
        ];
        let (genome, _) = Genome::try_from_genes(3, genes).unwrap();

        let observed = stats_tuples(genome.gene_snp_counts(0, 1, true));

        assert_eq!(observed, vec![(0, 6, 9), (1, 0, 10), (3, 4, 4)]);
    }

    #[test]
    fn gene_sort_ranks_follow_gene_names() {
        let genes = vec![
            gene("zeta", 1, &[(0, b"A"), (1, b"A")]),
            gene("alpha", 1, &[(0, b"A"), (1, b"A")]),
        ];
        let (genome, _) = Genome::try_from_genes(2, genes).unwrap();

        let observed: Vec<_> = genome
            .gene_snp_counts(0, 1, false)
            .into_iter()
            .map(|stats| (stats.gene_index, stats.gene_sort_rank))
            .collect();

        assert_eq!(observed, vec![(0, 1), (1, 0)]);
    }

    #[test]
    fn constructor_rejects_out_of_range_sample_indices() {
        let gene = Gene::new(
            "gene".to_string(),
            1,
            HashMap::from([(2, SampleBases::from_sequence(b"A"))]),
            0,
        );

        let error = Genome::try_from_genes(2, vec![gene]).unwrap_err();

        assert!(
            error.to_string().contains("sample index 2"),
            "error: {error}"
        );
    }

    #[test]
    fn constructor_rejects_total_alignment_len_past_bitmap_limit() {
        let genes = vec![
            Gene::new("large".to_string(), u32::MAX as usize, HashMap::new(), 0),
            Gene::new("extra".to_string(), 1, HashMap::new(), 0),
        ];

        let error = Genome::try_from_genes(0, genes).unwrap_err();

        assert!(error.to_string().contains("exceeds"), "error: {error}");
    }
}
