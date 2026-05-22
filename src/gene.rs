use roaring::RoaringBitmap;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
struct SampleBases {
    a: RoaringBitmap,
    c: RoaringBitmap,
    g: RoaringBitmap,
    t: RoaringBitmap,
    gap: RoaringBitmap,
}

impl SampleBases {
    fn from_sequence(sequence: &[u8]) -> Self {
        let mut bases = SampleBases::default();

        for (position, &base) in sequence.iter().enumerate() {
            bases.insert_base(position as u32, base);
        }

        bases
    }

    fn insert_base(&mut self, position: u32, base: u8) {
        match base.to_ascii_uppercase() {
            b'A' => {
                self.a.insert(position);
            }
            b'C' => {
                self.c.insert(position);
            }
            b'G' => {
                self.g.insert(position);
            }
            b'T' => {
                self.t.insert(position);
            }
            b'M' => {
                self.a.insert(position);
                self.c.insert(position);
            }
            b'R' => {
                self.a.insert(position);
                self.g.insert(position);
            }
            b'W' => {
                self.a.insert(position);
                self.t.insert(position);
            }
            b'S' => {
                self.c.insert(position);
                self.g.insert(position);
            }
            b'Y' => {
                self.c.insert(position);
                self.t.insert(position);
            }
            b'K' => {
                self.g.insert(position);
                self.t.insert(position);
            }
            b'V' => {
                self.a.insert(position);
                self.c.insert(position);
                self.g.insert(position);
            }
            b'H' => {
                self.a.insert(position);
                self.c.insert(position);
                self.t.insert(position);
            }
            b'D' => {
                self.a.insert(position);
                self.g.insert(position);
                self.t.insert(position);
            }
            b'B' => {
                self.c.insert(position);
                self.g.insert(position);
                self.t.insert(position);
            }
            b'-' => {
                self.gap.insert(position);
            }
            _ => {
                // includes N
                self.a.insert(position);
                self.c.insert(position);
                self.g.insert(position);
                self.t.insert(position);
            }
        }
    }
}

#[derive(Debug)]
pub(crate) struct Gene {
    name: String,
    alignment_len: usize,
    sample_indices: HashMap<String, usize>,
    samples: Vec<SampleBases>,
}

impl Gene {
    pub(crate) fn new(
        name: String,
        alignment_len: usize,
        sample_names: Vec<String>,
        ordered_sequences: Vec<Vec<u8>>,
    ) -> Self {
        debug_assert_eq!(sample_names.len(), ordered_sequences.len());

        let sample_indices = sample_names
            .iter()
            .enumerate()
            .map(|(index, sample)| (sample.clone(), index))
            .collect();
        let samples = ordered_sequences
            .iter()
            .map(|sequence| SampleBases::from_sequence(sequence))
            .collect();

        Gene {
            name,
            alignment_len,
            sample_indices,
            samples,
        }
    }

    pub(crate) fn snp_count(&self, sample_a: usize, sample_b: usize) -> usize {
        let left = &self.samples[sample_a];
        let right = &self.samples[sample_b];

        let mut matches = &left.a & &right.a;
        matches |= &left.c & &right.c;
        matches |= &left.g & &right.g;
        matches |= &left.t & &right.t;
        matches |= &left.gap & &right.gap;

        self.alignment_len - matches.len() as usize
    }

    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    pub(crate) fn alignment_len(&self) -> usize {
        self.alignment_len
    }

    pub(crate) fn sample_index(&self, sample: &str) -> Option<usize> {
        self.sample_indices.get(sample).copied()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gene(sequence_a: &[u8], sequence_b: &[u8]) -> Gene {
        Gene::new(
            "gene".to_string(),
            sequence_a.len(),
            vec!["sample_a".to_string(), "sample_b".to_string()],
            vec![sequence_a.to_vec(), sequence_b.to_vec()],
        )
    }

    fn one_base_snp_count(left: u8, right: u8) -> usize {
        gene(&[left], &[right]).snp_count(0, 1)
    }

    fn assert_membership(
        base: u8,
        expected_a: bool,
        expected_c: bool,
        expected_g: bool,
        expected_t: bool,
        expected_gap: bool,
    ) {
        let bases = SampleBases::from_sequence(&[base]);

        assert_eq!(bases.a.contains(0), expected_a, "{} A", base as char);
        assert_eq!(bases.c.contains(0), expected_c, "{} C", base as char);
        assert_eq!(bases.g.contains(0), expected_g, "{} G", base as char);
        assert_eq!(bases.t.contains(0), expected_t, "{} T", base as char);
        assert_eq!(bases.gap.contains(0), expected_gap, "{} gap", base as char);
    }

    #[test]
    fn exact_bases_match_only_themselves() {
        assert_eq!(gene(b"ACGT", b"ACGT").snp_count(0, 1), 0);
        assert_eq!(gene(b"ACGT", b"TGCA").snp_count(0, 1), 4);
    }

    #[test]
    fn iupac_bases_encode_expected_memberships() {
        let expectations = [
            (b'A', true, false, false, false, false),
            (b'C', false, true, false, false, false),
            (b'G', false, false, true, false, false),
            (b'T', false, false, false, true, false),
            (b'M', true, true, false, false, false),
            (b'R', true, false, true, false, false),
            (b'W', true, false, false, true, false),
            (b'S', false, true, true, false, false),
            (b'Y', false, true, false, true, false),
            (b'K', false, false, true, true, false),
            (b'V', true, true, true, false, false),
            (b'H', true, true, false, true, false),
            (b'D', true, false, true, true, false),
            (b'B', false, true, true, true, false),
            (b'N', true, true, true, true, false),
        ];

        for (base, expected_a, expected_c, expected_g, expected_t, expected_gap) in expectations {
            assert_membership(
                base,
                expected_a,
                expected_c,
                expected_g,
                expected_t,
                expected_gap,
            );
        }
    }

    #[test]
    fn iupac_bases_match_when_any_membership_overlaps() {
        assert_eq!(one_base_snp_count(b'R', b'A'), 0);
        assert_eq!(one_base_snp_count(b'R', b'G'), 0);
        assert_eq!(one_base_snp_count(b'R', b'C'), 1);
        assert_eq!(one_base_snp_count(b'B', b'A'), 1);
        assert_eq!(one_base_snp_count(b'B', b'T'), 0);
        assert_eq!(gene(b"MRWSYKVHDB", b"AATCTGCATG").snp_count(0, 1), 0);
    }

    #[test]
    fn n_and_unknown_non_gap_characters_match_any_ordinary_base() {
        assert_eq!(gene(b"NX?z", b"ACGT").snp_count(0, 1), 0);
        assert_eq!(one_base_snp_count(b'n', b't'), 0);
        assert_eq!(one_base_snp_count(b'?', b'A'), 0);
        assert_eq!(one_base_snp_count(b'X', b'C'), 0);
    }

    #[test]
    fn gap_matches_only_gap() {
        assert_membership(b'-', false, false, false, false, true);
        assert_eq!(one_base_snp_count(b'-', b'-'), 0);
        assert_eq!(one_base_snp_count(b'-', b'A'), 1);
        assert_eq!(one_base_snp_count(b'-', b'R'), 1);
        assert_eq!(one_base_snp_count(b'-', b'N'), 1);
        assert_eq!(one_base_snp_count(b'-', b'?'), 1);
    }

    #[test]
    fn two_sample_one_gene_has_expected_snp_count() {
        let gene = gene(b"AAAA", b"AACC");

        assert_eq!(gene.snp_count(0, 1), 2);
    }

    #[test]
    fn gap_mismatch_positions_are_counted() {
        let gene = gene(b"A-N?", b"AA-G");

        assert_eq!(gene.snp_count(0, 1), 2);
    }
}
