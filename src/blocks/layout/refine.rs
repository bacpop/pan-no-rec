//! Deterministic, bounded rearrangements that preserve the backbone node set.
use super::*;
use std::collections::{BTreeMap, HashSet};
const VARIANTS: usize = 64;
const ROTATION_BUDGET: usize = 20_000;
const SPLICE_BUDGET: usize = 20_000;
pub(super) const DETOUR_BUDGET: usize = 20_000;

#[derive(Default, Debug)]
pub(super) struct Stats {
    pub rotations: usize,
    pub rotation_checks: usize,
    pub splices: usize,
    pub splice_checks: usize,
    pub rotation_exhausted: bool,
    pub splice_exhausted: bool,
}

fn endpoints(paths: &Segments) -> BTreeMap<usize, Vec<(usize, bool)>> {
    let mut map: BTreeMap<usize, Vec<(usize, bool)>> = BTreeMap::new();
    for (i, p) in paths.iter().enumerate() {
        map.entry(p[0]).or_default().push((i, true));
        if p.len() > 1 {
            map.entry(*p.last().unwrap()).or_default().push((i, false));
        }
    }
    map
}
// The orientation flag says that the endpoint is the first node of the path.
fn merge(paths: &mut Segments, a: usize, b: usize, first_a: bool, first_b: bool) {
    let mut second = paths.remove(b);
    if !first_b {
        second.reverse();
    }
    if first_a {
        paths[a].reverse();
    }
    paths[a].extend(second);
}
fn join_longest(paths: &mut Segments, adj: &[Vec<usize>]) {
    loop {
        let ends = endpoints(paths);
        let mut best = None;
        let mut length = 0;
        for (&x, a) in &ends {
            for &y in &adj[x] {
                if let Some(b) = ends.get(&y) {
                    for &(i, ra) in a {
                        for &(j, rb) in b {
                            if i < j && paths[i].len() + paths[j].len() > length {
                                length = paths[i].len() + paths[j].len();
                                best = Some((i, j, ra, rb));
                            }
                        }
                    }
                }
            }
        }
        let Some((a, b, ra, rb)) = best else { break };
        merge(paths, a, b, ra, rb);
    }
}
fn variants(path: &[usize], adj: &[Vec<usize>], budget: &mut usize) -> Vec<(Vec<usize>, usize)> {
    let key = |p: &[usize]| (p[0].min(*p.last().unwrap()), p[0].max(*p.last().unwrap()));
    let mut seen = HashSet::from([key(path)]);
    let mut result = vec![(path.to_vec(), 0)];
    let mut next = 0;
    while next < result.len() && result.len() < VARIANTS && *budget > 0 {
        let (path, depth) = result[next].clone();
        next += 1;
        for reverse in [false, true] {
            let mut oriented = path.clone();
            if reverse {
                oriented.reverse();
            }
            let indices: HashMap<_, _> =
                oriented.iter().enumerate().map(|(i, &x)| (x, i)).collect();
            for &neighbor in &adj[*oriented.last().unwrap()] {
                if *budget == 0 || result.len() == VARIANTS {
                    break;
                }
                *budget -= 1;
                let Some(&i) = indices.get(&neighbor) else {
                    continue;
                };
                if i + 2 >= path.len() {
                    continue;
                }
                let mut rotated = oriented.clone();
                rotated[i + 1..].reverse();
                if seen.insert(key(&rotated)) {
                    result.push((rotated, depth + 1));
                }
            }
        }
    }
    result
}
pub(super) fn refine(paths: &mut Segments, adj: &[Vec<usize>]) -> Stats {
    refine_with_limits(paths, adj, ROTATION_BUDGET, SPLICE_BUDGET)
}
fn refine_with_limits(
    paths: &mut Segments,
    adj: &[Vec<usize>],
    rotation_limit: usize,
    splice_limit: usize,
) -> Stats {
    let mut stats = Stats::default();
    let mut budget = rotation_limit;
    loop {
        let variants: Vec<_> = paths
            .iter()
            .map(|p| variants(p, adj, &mut budget))
            .collect();
        let mut ends: BTreeMap<usize, Vec<(usize, usize, bool)>> = BTreeMap::new();
        for (i, options) in variants.iter().enumerate() {
            for (v, (p, _)) in options.iter().enumerate() {
                ends.entry(p[0]).or_default().push((i, v, true));
                if p.len() > 1 {
                    ends.entry(*p.last().unwrap())
                        .or_default()
                        .push((i, v, false));
                }
            }
        }
        let mut best = None;
        let mut score = None;
        for (&x, a) in &ends {
            for &y in &adj[x] {
                if let Some(b) = ends.get(&y) {
                    for &(i, vi, ra) in a {
                        for &(j, vj, rb) in b {
                            if i >= j {
                                continue;
                            }
                            let (pa, da) = &variants[i][vi];
                            let (pb, db) = &variants[j][vj];
                            let value = (
                                pa.len() + pb.len(),
                                std::cmp::Reverse(da + db),
                                std::cmp::Reverse((i, j, vi, vj)),
                            );
                            if score.is_none_or(|s| value > s) {
                                score = Some(value);
                                best = Some((i, j, vi, vj, ra, rb));
                            }
                        }
                    }
                }
            }
        }
        let Some((i, j, vi, vj, ra, rb)) = best else {
            break;
        };
        paths[i] = variants[i][vi].0.clone();
        paths[j] = variants[j][vj].0.clone();
        merge(paths, i, j, ra, rb);
        stats.rotations += 1;
        if budget == 0 {
            break;
        }
    }
    stats.rotation_checks = rotation_limit - budget;
    stats.rotation_exhausted = budget == 0;
    'refinement: loop {
        let ends = endpoints(paths);
        let mut accepted = None;
        'candidates: for (j, path) in paths.iter().enumerate() {
            for reverse in [false, true] {
                let mut oriented = path.clone();
                if reverse {
                    oriented.reverse();
                }
                for k in 1..oriented.len().saturating_sub(1) {
                    for &x in &adj[oriented[k]] {
                        if let Some(options) = ends.get(&x) {
                            for &(i, ra) in options {
                                if i == j {
                                    continue;
                                }
                                if stats.splice_checks == splice_limit {
                                    stats.splice_exhausted = true;
                                    break 'refinement;
                                }
                                stats.splice_checks += 1;
                                let mut combined = paths[i].clone();
                                if ra {
                                    combined.reverse();
                                }
                                combined.extend(oriented[..=k].iter().rev());
                                let mut trial: Segments = paths
                                    .iter()
                                    .enumerate()
                                    .filter(|&(idx, _)| idx != i && idx != j)
                                    .map(|(_, p)| p.clone())
                                    .collect();
                                trial.push(combined);
                                trial.push(oriented[k + 1..].to_vec());
                                join_longest(&mut trial, adj);
                                if trial.len() < paths.len() {
                                    accepted = Some(trial);
                                    break 'candidates;
                                }
                            }
                        }
                    }
                }
            }
        }
        let Some(trial) = accepted else { break };
        *paths = trial;
        stats.splices += 1;
    }
    stats
}

// Enrich a connector without revisiting any placed or already routed node.
pub(super) fn detours(
    start: usize,
    route: Vec<usize>,
    used: &[bool],
    allowed: &[bool],
    adj: &[Vec<usize>],
    sizes: &[usize],
    budget: &mut usize,
) -> Vec<usize> {
    let mut occupied = used.to_vec();
    occupied[start] = true;
    for &x in &route {
        occupied[x] = true;
    }
    let mut result = Vec::new();
    let mut previous = start;
    for target in route {
        let mut stack = vec![(previous, Vec::new(), 0usize)];
        let mut best = Vec::new();
        let mut best_score = (0, 0);
        let mut checks = 0;
        while let Some((node, path, weight)) = stack.pop() {
            if *budget == 0 || checks == 256 {
                break;
            }
            checks += 1;
            *budget -= 1;
            for &next in &adj[node] {
                if next == target {
                    let score = (weight, path.len());
                    if score > best_score {
                        best = path.clone();
                        best_score = score;
                    }
                    continue;
                }
                if !allowed[next] || occupied[next] || path.contains(&next) || path.len() == 8 {
                    continue;
                }
                let mut extended = path.clone();
                extended.push(next);
                stack.push((next, extended, weight + sizes[next]));
            }
        }
        for x in best {
            occupied[x] = true;
            result.push(x);
        }
        result.push(target);
        previous = target;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    fn graph(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut adj = vec![vec![]; n];
        for &(a, b) in edges {
            adj[a].push(b);
            adj[b].push(a);
        }
        for neighbors in &mut adj {
            neighbors.sort_unstable();
        }
        adj
    }
    fn valid(paths: &Segments, adj: &[Vec<usize>], n: usize) {
        let flat: Vec<_> = paths.iter().flatten().copied().collect();
        assert_eq!(flat.len(), n);
        assert_eq!(flat.iter().collect::<HashSet<_>>().len(), n);
        assert!(
            paths
                .iter()
                .all(|p| !p.is_empty() && p.windows(2).all(|e| adj[e[0]].contains(&e[1])))
        );
    }
    #[test]
    fn rotation_exposes_an_internal_node_and_joins_without_repetition() {
        let adj = graph(6, &[(0, 1), (1, 2), (2, 3), (3, 1), (2, 4), (4, 5)]);
        let original = vec![vec![0, 1, 2, 3], vec![4, 5]];
        let mut paths = original.clone();
        let stats = refine(&mut paths, &adj);
        assert_eq!(paths.len(), 1);
        assert_eq!(stats.rotations, 1);
        valid(&paths, &adj, 6);
        let mut repeated = original;
        refine(&mut repeated, &adj);
        assert_eq!(paths, repeated);
    }
    #[test]
    fn splicing_frees_fragment_for_an_additional_join() {
        let adj = graph(
            9,
            &[
                (0, 1),
                (1, 2),
                (3, 4),
                (4, 5),
                (5, 6),
                (7, 8),
                (2, 4),
                (5, 7),
            ],
        );
        let mut paths = vec![vec![0, 1, 2], vec![3, 4, 5, 6], vec![7, 8]];
        let stats = refine(&mut paths, &adj);
        assert!(paths.len() < 3);
        assert!(stats.splices > 0);
        valid(&paths, &adj, 9);
    }
    #[test]
    fn non_improving_splice_is_rejected() {
        let adj = graph(6, &[(0, 1), (1, 2), (3, 4), (4, 5), (2, 4)]);
        let original = vec![vec![0, 1, 2], vec![3, 4, 5]];
        let mut paths = original.clone();
        let stats = refine(&mut paths, &adj);
        assert_eq!(paths, original);
        assert_eq!(stats.splices, 0);
        assert!(stats.splice_checks > 0);
    }
    #[test]
    fn longer_connectors_respect_masks_future_route_nodes_and_budget() {
        let adj = graph(4, &[(0, 2), (0, 1), (1, 2), (1, 3), (2, 3)]);
        let sizes = vec![20; 4];
        let allowed = vec![true; 4];
        let used = vec![true, false, false, false];
        let mut budget = 100;
        let route = detours(0, vec![2, 3], &used, &allowed, &adj, &sizes, &mut budget);
        assert_eq!(route, vec![1, 2, 3]);
        assert_eq!(used, vec![true, false, false, false]);
        let mut occupied = used.clone();
        occupied[1] = true;
        assert_eq!(
            detours(0, vec![2], &occupied, &allowed, &adj, &sizes, &mut 100),
            vec![2]
        );
        let mut eligible = allowed.clone();
        eligible[1] = false;
        eligible[3] = false;
        assert_eq!(
            detours(0, vec![2], &used, &eligible, &adj, &sizes, &mut 100),
            vec![2]
        );
        let mut budget = 1;
        assert_eq!(
            detours(0, vec![2], &used, &allowed, &adj, &sizes, &mut budget),
            vec![2]
        );
        assert_eq!(budget, 0);
        assert_eq!(
            detours(0, vec![2], &used, &allowed, &adj, &sizes, &mut budget),
            vec![2]
        );
    }
    #[test]
    fn rotation_budget_is_reproducible() {
        let adj = graph(5, &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 1), (4, 0)]);
        let mut budget = 1;
        let a = variants(&[0, 1, 2, 3, 4], &adj, &mut budget);
        assert_eq!(budget, 0);
        let mut budget = 1;
        let b = variants(&[0, 1, 2, 3, 4], &adj, &mut budget);
        assert_eq!(a, b);
        assert_eq!(variants(&[0], &adj, &mut 0), vec![(vec![0], 0)]);
    }
    #[test]
    fn exhausted_refinement_budgets_preserve_a_valid_incumbent() {
        let adj = graph(6, &[(0, 1), (1, 2), (3, 4), (4, 5), (2, 4)]);
        let original = vec![vec![0, 1, 2], vec![3, 4, 5]];
        let mut paths = original.clone();
        let stats = refine_with_limits(&mut paths, &adj, 0, 1);
        assert!(stats.rotation_exhausted && stats.splice_exhausted);
        assert_eq!(stats.splice_checks, 1);
        assert_eq!(paths, original);
        valid(&paths, &adj, 6);
    }
    #[test]
    fn layout_selection_prioritises_breaks_then_weight_and_continuity() {
        assert!(
            layout_score(&vec![vec![0, 1]], &[1, 1, 100])
                > layout_score(&vec![vec![0], vec![2]], &[1, 1, 100])
        );
        assert!(
            layout_score(&vec![vec![0, 2]], &[1, 1, 100])
                > layout_score(&vec![vec![0, 1]], &[1, 1, 100])
        );
        assert!(
            layout_score(&vec![vec![0, 1, 2], vec![3]], &[1; 4])
                > layout_score(&vec![vec![0, 1], vec![2, 3]], &[1; 4])
        );
    }
}
