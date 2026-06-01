use hashbrown::HashMap;
use roaring::RoaringBitmap;

#[derive(Clone, Debug, Default)]
pub(crate) struct SampleBases {
    a: RoaringBitmap,
    c: RoaringBitmap,
    g: RoaringBitmap,
    t: RoaringBitmap,
    gap: RoaringBitmap,
}

impl SampleBases {
    // Encodes a sequence into per-base bitmap positions.
    pub(crate) fn from_sequence(sequence: &[u8]) -> Self {
        let mut bases = SampleBases::default();

        for (position, &base) in sequence.iter().enumerate() {
            bases.insert_base(position as u32, base);
        }

        bases
    }

    // Records one alignment position under its matching base set.
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

    // Counts positions encoded as any non-gap base.
    pub(crate) fn non_gap_count(&self, alignment_len: usize) -> usize {
        alignment_len - self.gap.len() as usize
    }
}

#[derive(Debug)]
pub(crate) struct Gene {
    name: String,
    alignment_len: usize,
    sequences: HashMap<usize, SampleBases>,
    paralog_count: usize,
}

impl Gene {
    // Builds a gene from ordered global sample indices and pre-encoded samples.
    pub(crate) fn new(
        name: String,
        alignment_len: usize,
        sequences: HashMap<usize, SampleBases>,
        paralog_count: usize,
    ) -> Self {
        Gene {
            name,
            alignment_len,
            sequences,
            paralog_count,
        }
    }

    // Counts non-matching alignment columns between two samples.
    // Ignores sites which are both gaps
    // If gaps == true, counts a gap vs base as a mismatch and includes in the length
    // If gaps = false, gap positions are ignored completely
    // Returns the number of mismatches, and the aligned length
    pub(crate) fn snp_count(&self, sample_a: usize, sample_b: usize, gaps: bool) -> (usize, usize) {
        let left = &self.sequences[&sample_a];
        let right = &self.sequences[&sample_b];

        let mut matches = &left.a & &right.a;
        matches |= &left.c & &right.c;
        matches |= &left.g & &right.g;
        matches |= &left.t & &right.t;

        let both_gap = &left.gap & &right.gap;
        let length = if gaps {
            self.alignment_len - both_gap.len() as usize
        } else {
            self.alignment_len - (left.gap.len() + right.gap.len() - both_gap.len()) as usize
        };

        (length - matches.len() as usize, length)
    }

    // Returns the gene identifier.
    pub(crate) fn name(&self) -> &str {
        &self.name
    }

    pub(crate) fn paralog_count(&self) -> Option<usize> {
        if self.paralog_count == 0 { None } else { Some(self.paralog_count) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Builds a two-sample test gene from byte sequences.
    fn gene(sequence_a: &[u8], sequence_b: &[u8]) -> Gene {
        Gene::new(
            "gene".to_string(),
            sequence_a.len(),
            vec![0, 1],
            vec![sequence_a.to_vec(), sequence_b.to_vec()],
        )
    }

    // Counts SNPs for a single-base two-sample gene using default gap handling.
    fn one_base_snp_count(left: u8, right: u8) -> (usize, usize) {
        gene(&[left], &[right]).snp_count(0, 1, false)
    }

    // Asserts bitmap membership for one encoded base.
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
    // Verifies sample lookup maps global sample indices to local sequence indices.
    fn sample_lookup_uses_global_sample_indices() {
        let gene = Gene::new(
            "gene".to_string(),
            4,
            vec![3, 1],
            vec![b"AAAA".to_vec(), b"CCCC".to_vec()],
        );

        assert_eq!(gene.sample_index(3), Some(0));
        assert_eq!(gene.sample_index(1), Some(1));
        assert_eq!(gene.sample_index(0), None);
    }

    #[test]
    // Verifies exact DNA bases only match identical bases.
    fn exact_bases_match_only_themselves() {
        assert_eq!(gene(b"ACGT", b"ACGT").snp_count(0, 1, false), (0, 4));
        assert_eq!(gene(b"ACGT", b"TGCA").snp_count(0, 1, false), (4, 4));
    }

    #[test]
    // Verifies IUPAC ambiguity codes map to the expected base sets.
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
    // Verifies ambiguous bases match when their base sets overlap.
    fn iupac_bases_match_when_any_membership_overlaps() {
        assert_eq!(one_base_snp_count(b'R', b'A'), (0, 1));
        assert_eq!(one_base_snp_count(b'R', b'G'), (0, 1));
        assert_eq!(one_base_snp_count(b'R', b'C'), (1, 1));
        assert_eq!(one_base_snp_count(b'B', b'A'), (1, 1));
        assert_eq!(one_base_snp_count(b'B', b'T'), (0, 1));
        assert_eq!(
            gene(b"MRWSYKVHDB", b"AATCTGCATG").snp_count(0, 1, false),
            (0, 10)
        );
    }

    #[test]
    // Verifies N and unknown non-gap bases match ordinary bases.
    fn n_and_unknown_non_gap_characters_match_any_ordinary_base() {
        assert_eq!(gene(b"NX?z", b"ACGT").snp_count(0, 1, false), (0, 4));
        assert_eq!(one_base_snp_count(b'n', b't'), (0, 1));
        assert_eq!(one_base_snp_count(b'?', b'A'), (0, 1));
        assert_eq!(one_base_snp_count(b'X', b'C'), (0, 1));
    }

    #[test]
    // Verifies gaps are ignored by default, and counted only in gap-inclusive mode.
    fn gap_handling_depends_on_gap_mode() {
        assert_membership(b'-', false, false, false, false, true);

        let gene = gene(b"A--CGGTTT-", b"ACCCTG----");

        assert_eq!(gene.snp_count(0, 1, false), (1, 4));
        assert_eq!(gene.snp_count(0, 1, true), (6, 9));
    }
}
