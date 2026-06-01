use statrs::function::gamma::ln_gamma;
use std::cmp::Ordering;

const DEFAULT_A0: f64 = 9.0;
const DEFAULT_A1: f64 = 1.0;

#[derive(Clone, Copy, Debug)]
pub(crate) struct PairGeneStats {
    pub(crate) gene_index: usize,
    pub(crate) gene_sort_rank: usize,
    pub(crate) snps: usize,
    pub(crate) length: usize,
}

// Selects genes above the inferred recombination threshold.
pub(crate) fn select_recombinant_gene_indices(
    genes: impl IntoIterator<Item = PairGeneStats>,
) -> Vec<usize> {
    let mut ordered_genes: Vec<_> = genes.into_iter().collect();
    order_pair_genes(&mut ordered_genes);

    let model_probs = threshold_model_probabilities(&ordered_genes);
    let threshold = find_threshold(&model_probs);

    ordered_genes
        .into_iter()
        .skip(threshold)
        .map(|gene| gene.gene_index)
        .collect()
}

// Orders genes by divergence and deterministic tie-breakers.
fn order_pair_genes(genes: &mut [PairGeneStats]) {
    genes.sort_by(|left, right| {
        snp_proportion(*left)
            .total_cmp(&snp_proportion(*right))
            .then_with(|| right.length.cmp(&left.length))
            .then_with(|| left.gene_sort_rank.cmp(&right.gene_sort_rank))
    });
}

// Computes the SNP proportion for one pairwise gene comparison.
fn snp_proportion(gene: PairGeneStats) -> f64 {
    gene.snps as f64 / gene.length as f64
}

// Computes normalized threshold model probabilities for ordered genes.
fn threshold_model_probabilities(ordered_genes: &[PairGeneStats]) -> Vec<f64> {
    if ordered_genes.is_empty() {
        return Vec::new();
    }

    let individual_logmls: Vec<_> = ordered_genes
        .iter()
        .map(|gene| {
            alt_log_likelihood(
                DEFAULT_A1,
                DEFAULT_A0,
                (gene.length - gene.snps) as f64,
                gene.snps as f64,
            )
        })
        .collect();

    let mut cumulative_lengths = 0usize;
    let mut cumulative_snps = 0usize;
    let joint_threshold_logmls: Vec<_> = ordered_genes
        .iter()
        .map(|gene| {
            cumulative_lengths += gene.length;
            cumulative_snps += gene.snps;
            alt_log_likelihood(
                DEFAULT_A1,
                DEFAULT_A0,
                (cumulative_lengths - cumulative_snps) as f64,
                cumulative_snps as f64,
            )
        })
        .collect();

    let mut suffix_individual_logmls = vec![0.0; ordered_genes.len()];
    let mut suffix_sum = 0.0;
    for index in (0..ordered_genes.len()).rev() {
        suffix_individual_logmls[index] = suffix_sum;
        suffix_sum += individual_logmls[index];
    }

    let threshold_logmls: Vec<_> = joint_threshold_logmls
        .into_iter()
        .zip(suffix_individual_logmls)
        .map(|(joint, suffix)| joint + suffix)
        .collect();

    let max_logml = threshold_logmls
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut normalised_weights: Vec<_> = threshold_logmls
        .iter()
        .map(|logml| (logml - max_logml).exp())
        .collect();
    let weight_sum: f64 = normalised_weights.iter().sum();

    for weight in &mut normalised_weights {
        *weight /= weight_sum;
    }

    normalised_weights
}

// Computes the beta-binomial alternative log likelihood.
fn alt_log_likelihood(a0: f64, a1: f64, n0: f64, n1: f64) -> f64 {
    ln_gamma(a0 + a1) - ln_gamma(a0 + a1 + n0 + n1) + ln_gamma(a0 + n0) + ln_gamma(a1 + n1)
        - ln_gamma(a0)
        - ln_gamma(a1)
}

// Converts model probabilities into an integer threshold.
fn find_threshold(model_probs: &[f64]) -> usize {
    let expected_threshold: f64 = model_probs
        .iter()
        .enumerate()
        .map(|(index, prob)| prob * (index + 1) as f64)
        .sum();

    round_ties_to_even(expected_threshold)
}

// Rounds a non-negative value using ties-to-even semantics.
fn round_ties_to_even(value: f64) -> usize {
    debug_assert!(value >= 0.0);

    let floor = value.floor();
    let fraction = value - floor;
    match fraction.total_cmp(&0.5) {
        Ordering::Less => floor as usize,
        Ordering::Greater => floor as usize + 1,
        Ordering::Equal => {
            let floor = floor as usize;
            if floor % 2 == 0 { floor } else { floor + 1 }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Builds pairwise gene stats for model tests.
    fn gene(gene_index: usize, gene_sort_rank: usize, snps: usize, length: usize) -> PairGeneStats {
        PairGeneStats {
            gene_index,
            gene_sort_rank,
            snps,
            length,
        }
    }

    // Returns gene indices after applying model ordering.
    fn ordered_indices(mut genes: Vec<PairGeneStats>) -> Vec<usize> {
        order_pair_genes(&mut genes);
        genes.into_iter().map(|gene| gene.gene_index).collect()
    }

    #[test]
    // Verifies genes sort by increasing SNP proportion.
    fn sorts_by_increasing_snp_proportion() {
        let observed = ordered_indices(vec![gene(0, 0, 2, 10), gene(1, 1, 1, 10)]);

        assert_eq!(observed, vec![1, 0]);
    }

    #[test]
    // Verifies equal SNP proportions prefer longer alignments.
    fn tie_breaks_by_longer_alignment_first() {
        let observed = ordered_indices(vec![gene(0, 0, 1, 10), gene(1, 1, 2, 20)]);

        assert_eq!(observed, vec![1, 0]);
    }

    #[test]
    // Verifies remaining ties are broken by precomputed gene-name sort rank.
    fn final_tie_breaker_is_gene_sort_rank() {
        let observed = ordered_indices(vec![gene(0, 1, 1, 10), gene(1, 0, 1, 10)]);

        assert_eq!(observed, vec![1, 0]);
    }

    #[test]
    // Verifies thresholding selects the reference recombinant tail.
    fn bayesian_threshold_selects_python_reference_tail() {
        let genes = vec![gene(0, 0, 8, 10), gene(1, 1, 0, 10), gene(2, 2, 1, 10)];

        let observed = select_recombinant_gene_indices(genes);

        assert_eq!(observed, vec![0]);
    }

    #[test]
    // Verifies threshold rounding matches Python's ties-to-even behavior.
    fn threshold_uses_python_rounding_ties_to_even() {
        assert_eq!(round_ties_to_even(2.5), 2);
        assert_eq!(round_ties_to_even(3.5), 4);
    }
}
