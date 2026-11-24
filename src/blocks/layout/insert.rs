//! Final placement of degree-two accessory chains between existing core flanks.
use super::*;
use std::collections::BTreeMap;

const MIN_MEAN_PERCENT: u128 = 3;

#[derive(Debug, Default, PartialEq, Eq)]
pub(super) struct Stats {
    pub eligible_nodes: usize,
    pub eligible_chains: usize,
    pub insertion_positions: usize,
    pub competing_positions: usize,
    pub inserted_nodes: usize,
    pub inserted_chains: usize,
}

pub(super) fn insert_accessory(
    paths: &mut Segments,
    adjacency: &[Vec<usize>],
    core: &[bool],
    sizes: &[usize],
    node_ids: &[i64],
    total_samples: usize,
) -> Stats {
    let mut positions = vec![None; adjacency.len()];
    for (segment, path) in paths.iter().enumerate() {
        for (offset, &node) in path.iter().enumerate() {
            positions[node] = Some((segment, offset));
        }
    }
    // Discover whole components before applying the frequency filter. A rare
    // interior gene must not split a candidate or cause a partial insertion.
    let available: Vec<_> = (0..adjacency.len())
        .map(|n| positions[n].is_none() && !core[n] && adjacency[n].len() == 2)
        .collect();
    let mut seen = vec![false; adjacency.len()];
    let mut stats = Stats::default();
    let mut winners: BTreeMap<(usize, usize), (Vec<usize>, u128, usize)> = BTreeMap::new();
    for node in 0..adjacency.len() {
        if !available[node] || seen[node] {
            continue;
        }
        let mut stack = vec![node];
        let mut component = Vec::new();
        let mut boundary = Vec::new();
        seen[node] = true;
        while let Some(n) = stack.pop() {
            component.push(n);
            for &next in &adjacency[n] {
                if !available[next] {
                    boundary.push((n, next));
                } else if !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        // Cycles have no external edges; dangling and branched attachments fail
        // the two distinct core flanks requirement.
        if boundary.len() != 2 {
            continue;
        }
        let (a, b) = (boundary[0].1, boundary[1].1);
        if a == b || !core[a] || !core[b] {
            continue;
        }
        let (Some((sa, ia)), Some((sb, ib))) = (positions[a], positions[b]) else {
            continue;
        };
        let left = ia.min(ib);
        let right = ia.max(ib);
        if sa != sb || right - left != 1 || left == 0 || right + 1 == paths[sa].len() {
            continue;
        }
        let weight: u128 = component
            .iter()
            .map(|&n| sizes[n].min(total_samples) as u128)
            .sum();
        if weight * 100 < MIN_MEAN_PERCENT * total_samples as u128 * component.len() as u128 {
            continue;
        }
        let (mut current, mut previous) = if ia < ib { boundary[0] } else { boundary[1] };
        let mut chain = Vec::with_capacity(component.len());
        loop {
            chain.push(current);
            let next = *adjacency[current].iter().find(|&&n| n != previous).unwrap();
            if !available[next] {
                break;
            }
            previous = current;
            current = next;
        }
        debug_assert_eq!(chain.len(), component.len());
        stats.eligible_nodes += chain.len();
        stats.eligible_chains += 1;
        match winners.entry((sa, left)) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert((chain, weight, 1));
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                let (incumbent, incumbent_weight, count) = entry.get_mut();
                *count += 1;
                if weight > *incumbent_weight
                    || (weight == *incumbent_weight
                        && chain
                            .iter()
                            .map(|&n| node_ids[n])
                            .cmp(incumbent.iter().map(|&n| node_ids[n]))
                            .is_lt())
                {
                    *incumbent = chain;
                    *incumbent_weight = weight;
                }
            }
        }
    }
    stats.insertion_positions = winners.len();
    stats.competing_positions = winners.values().filter(|(_, _, count)| *count > 1).count();
    let mut inserted_per_path = vec![0; paths.len()];
    for ((segment, offset), (chain, _, _)) in winners {
        let position = offset + inserted_per_path[segment] + 1;
        let chain_len = chain.len();
        paths[segment].splice(position..position, chain);
        inserted_per_path[segment] += chain_len;
        stats.inserted_nodes += chain_len;
        stats.inserted_chains += 1;
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    fn insert_accessory(
        paths: &mut Segments,
        adj: &[Vec<usize>],
        core: &[bool],
        sizes: &[usize],
        ids: &[i64],
    ) -> Stats {
        super::insert_accessory(paths, adj, core, sizes, ids, 20)
    }
    fn adjacency(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut result = vec![vec![]; n];
        for &(a, b) in edges {
            if a != b {
                result[a].push(b);
                result[b].push(a);
            }
        }
        for neighbors in &mut result {
            neighbors.sort_unstable();
            neighbors.dedup();
        }
        result
    }
    fn valid(paths: &Segments, adj: &[Vec<usize>]) {
        let flat: Vec<_> = paths.iter().flatten().collect();
        assert_eq!(
            flat.len(),
            flat.iter().collect::<std::collections::HashSet<_>>().len()
        );
        assert!(
            paths
                .iter()
                .all(|path| path.windows(2).all(|pair| adj[pair[0]].contains(&pair[1])))
        );
    }
    #[test]
    fn inserts_rare_genes_at_multiple_positions_and_preserves_flanks() {
        let adj = adjacency(
            8,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (1, 5),
                (5, 2),
                (2, 6),
                (6, 3),
                (3, 7),
                (7, 4),
            ],
        );
        let mut paths = vec![vec![0, 1, 2, 3, 4]];
        let core = vec![true, true, true, true, true, false, false, false];
        let stats = insert_accessory(
            &mut paths,
            &adj,
            &core,
            &[20, 20, 20, 20, 20, 1, 2, 1],
            &(0..8).collect::<Vec<_>>(),
        );
        assert_eq!(paths, vec![vec![0, 1, 5, 2, 6, 3, 4]]); // 7 has a terminal flank.
        assert_eq!(stats.inserted_nodes, 2);
        assert_eq!(stats.eligible_nodes, 2);
        valid(&paths, &adj);
        let mut reversed = vec![vec![4, 3, 2, 1, 0]];
        insert_accessory(
            &mut reversed,
            &adj,
            &core,
            &[20, 20, 20, 20, 20, 1, 2, 1],
            &(0..8).collect::<Vec<_>>(),
        );
        assert_eq!(
            reversed[0],
            paths[0].iter().rev().copied().collect::<Vec<_>>()
        );
        assert_eq!(
            insert_accessory(
                &mut paths,
                &adj,
                &core,
                &[20; 8],
                &(0..8).collect::<Vec<_>>()
            )
            .inserted_nodes,
            0
        );
    }
    #[test]
    fn competition_prefers_size_then_lowest_gml_id() {
        let adj = adjacency(
            7,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (1, 4),
                (4, 2),
                (1, 5),
                (5, 2),
                (1, 6),
                (6, 2),
            ],
        );
        let core = [true, true, true, true, false, false, false];
        let ids = [0, 1, 2, 3, 100, 50, 60];
        let mut paths = vec![vec![0, 1, 2, 3]];
        let stats = insert_accessory(&mut paths, &adj, &core, &[20, 20, 20, 20, 5, 8, 8], &ids);
        assert_eq!(paths[0], vec![0, 1, 5, 2, 3]);
        assert_eq!(
            stats,
            Stats {
                eligible_nodes: 3,
                eligible_chains: 3,
                insertion_positions: 1,
                competing_positions: 1,
                inserted_nodes: 1,
                inserted_chains: 1
            }
        );
    }
    #[test]
    fn excludes_nonconsecutive_cross_segment_noncore_and_wrong_degree() {
        let adj = adjacency(
            12,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (5, 6),
                (1, 7),
                (7, 3),
                (2, 8),
                (8, 5),
                (1, 9),
                (9, 2),
                (9, 4),
                (1, 10),
                (10, 11),
            ],
        );
        let core = [
            true, true, true, true, true, true, true, false, false, false, false, false,
        ];
        let mut paths = vec![vec![0, 1, 2, 3, 4], vec![5, 6]];
        let original = paths.clone();
        let stats = insert_accessory(
            &mut paths,
            &adj,
            &core,
            &[20; 12],
            &(0..12).collect::<Vec<_>>(),
        );
        assert_eq!(stats.inserted_nodes, 0);
        assert_eq!(paths, original);
    }
    #[test]
    fn does_not_insert_already_placed_or_core_nodes() {
        let adj = adjacency(6, &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2), (1, 5), (5, 2)]);
        let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
        let original = paths.clone();
        assert_eq!(
            insert_accessory(
                &mut paths,
                &adj,
                &[true, true, true, true, false, true],
                &[20; 6],
                &(0..6).collect::<Vec<_>>()
            )
            .inserted_nodes,
            0
        );
        assert_eq!(paths, original);
    }
    #[test]
    fn chain_mean_filter_is_inclusive_and_never_inserts_partial_chains() {
        let adj = adjacency(6, &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 2)]);
        let core = [true, true, true, true, false, false];
        let ids = [0, 1, 2, 3, 4, 5];
        for (frequencies, expected_count) in [([1, 5], 2), ([1, 4], 0), ([0, 0], 0), ([3, 3], 2)] {
            let mut paths = vec![vec![0, 1, 2, 3]];
            let sizes = [100, 100, 100, 100, frequencies[0], frequencies[1]];
            let stats = super::insert_accessory(&mut paths, &adj, &core, &sizes, &ids, 100);
            assert_eq!(stats.inserted_nodes, expected_count);
            assert_eq!(
                paths[0],
                if expected_count == 2 {
                    vec![0, 1, 4, 5, 2, 3]
                } else {
                    vec![0, 1, 2, 3]
                }
            );
            valid(&paths, &adj);
        }
        let mut reversed = vec![vec![3, 2, 1, 0]];
        super::insert_accessory(
            &mut reversed,
            &adj,
            &core,
            &[100, 100, 100, 100, 1, 5],
            &ids,
            100,
        );
        assert_eq!(reversed[0], vec![3, 2, 5, 4, 1, 0]);
    }
    #[test]
    fn singleton_threshold_and_sample_count_clamping() {
        let adj = adjacency(5, &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2)]);
        for (count, accepted) in [(2, false), (3, true), (4, true), (200, true)] {
            let mut paths = vec![vec![0, 1, 2, 3]];
            let stats = super::insert_accessory(
                &mut paths,
                &adj,
                &[true, true, true, true, false],
                &[100, 100, 100, 100, count],
                &[0, 1, 2, 3, 4],
                100,
            );
            assert_eq!(stats.inserted_nodes, usize::from(accepted));
        }
    }
    #[test]
    fn competing_chains_use_total_weight_then_oriented_ids() {
        // A singleton competes with two chains. Weight wins before mean;
        // equal weights are resolved using IDs in backbone orientation.
        let adj = adjacency(
            9,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (1, 4),
                (4, 2),
                (1, 5),
                (5, 6),
                (6, 2),
                (1, 7),
                (7, 8),
                (8, 2),
            ],
        );
        let core = [true, true, true, true, false, false, false, false, false];
        let ids = [0, 1, 2, 3, 4, 60, 10, 50, 70];
        let mut paths = vec![vec![0, 1, 2, 3]];
        let sizes = [100, 100, 100, 100, 9, 5, 5, 5, 5];
        let stats = super::insert_accessory(&mut paths, &adj, &core, &sizes, &ids, 100);
        assert_eq!(paths[0], vec![0, 1, 7, 8, 2, 3]);
        assert_eq!(stats.eligible_chains, 3);
        assert_eq!(stats.eligible_nodes, 5);
        assert_eq!(stats.competing_positions, 1);
        assert_eq!(stats.inserted_chains, 1);
        valid(&paths, &adj);
        let mut reversed = vec![vec![3, 2, 1, 0]];
        super::insert_accessory(&mut reversed, &adj, &core, &sizes, &ids, 100);
        assert_eq!(reversed[0], vec![3, 2, 6, 5, 1, 0]);
    }
    #[test]
    fn rejects_cycles_branches_and_chains_touching_placed_accessory() {
        let adj = adjacency(
            13,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (4, 5),
                (5, 6),
                (6, 4),
                (1, 7),
                (7, 8),
                (8, 2),
                (7, 9),
                (1, 10),
                (10, 11),
                (11, 2),
                (1, 12),
                (12, 1),
            ],
        );
        let core = [
            true, true, true, true, false, false, false, false, false, false, false, false, false,
        ];
        let mut paths = vec![vec![0, 1, 2, 3], vec![11]];
        let original = paths.clone();
        let stats = super::insert_accessory(
            &mut paths,
            &adj,
            &core,
            &[100; 13],
            &(0..13).collect::<Vec<_>>(),
            100,
        );
        assert_eq!(paths, original);
        assert_eq!(stats.inserted_nodes, 0);
    }
}
