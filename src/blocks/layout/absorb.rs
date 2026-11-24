//! Supported local absorption of short backbone segments, before accessory insertion.
use super::*;
use roaring::RoaringBitmap;
use std::collections::{BTreeMap, BTreeSet};

const MAX_SOURCE_GENES: usize = 25; // Exclusive.
const MIN_SUPPORT_PERCENT: u128 = 3;
#[derive(Clone, Copy)]
struct Limits {
    local_nodes: usize,
    connectors: usize,
    per_source: usize,
    global: usize,
}
const LIMITS: Limits = Limits {
    local_nodes: 64,
    connectors: 16,
    per_source: 20_000,
    global: 200_000,
};
#[cfg(test)]
use super::support::edge_members;
use super::support::{EdgeMembers, key};

#[derive(Debug, Default)]
pub(super) struct Stats {
    pub accepted: usize,
    pub added_nodes: usize,
    pub searched_sources: usize,
    pub multiple_hosts: usize,
    pub no_host: usize,
    pub no_slot: usize,
    pub oversized: usize,
    pub exhausted_sources: usize,
    pub expansions: usize,
    pub global_exhausted: bool,
    pub unavailable_edges: usize,
}
struct Input<'a> {
    adjacency: &'a [Vec<usize>],
    core: &'a [bool],
    sizes: &'a [usize],
    ids: Vec<i64>,
    members: &'a EdgeMembers,
    samples: usize,
}

pub(super) fn absorb(
    paths: &mut Segments,
    graph: &BlockGraph,
    adjacency: &[Vec<usize>],
    core: &[bool],
    sizes: &[usize],
    total_samples: usize,
    members: &EdgeMembers,
) -> Stats {
    let input = Input {
        adjacency,
        core,
        sizes,
        ids: graph.nodes.iter().map(|n| n.id).collect(),
        members,
        samples: total_samples,
    };
    run(paths, &input, LIMITS)
}

struct Local {
    allowed: Vec<bool>,
    required: Vec<bool>,
    host: usize,
    slots: Vec<usize>,
}
fn discover(
    source: usize,
    paths: &Segments,
    positions: &[Option<(usize, usize)>],
    input: &Input<'_>,
    stats: &mut Stats,
    limits: Limits,
) -> Option<Local> {
    let mut allowed = vec![false; positions.len()];
    let mut required = vec![false; positions.len()];
    let mut stack = paths[source].clone();
    for &n in &stack {
        allowed[n] = true;
        required[n] = true;
    }
    let mut count = stack.len();
    let mut contacts: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    while let Some(n) = stack.pop() {
        for &x in &input.adjacency[n] {
            if let Some((segment, offset)) = positions[x] {
                if segment != source {
                    contacts.entry(segment).or_default().insert(offset);
                }
            } else if !input.core[x] && !allowed[x] {
                allowed[x] = true;
                count += 1;
                stack.push(x);
            }
        }
    }
    // Connectivity uses all edges and all unused accessory nodes, including
    // rare/unsupported ones. Missing evidence cannot create a false sole host.
    if contacts.is_empty() {
        stats.no_host += 1;
        return None;
    }
    if contacts.len() != 1 {
        stats.multiple_hosts += 1;
        return None;
    }
    if count > limits.local_nodes {
        stats.oversized += 1;
        return None;
    }
    let (host, contacts) = contacts.into_iter().next().unwrap();
    let slots: Vec<_> = contacts
        .iter()
        .copied()
        .filter(|i| contacts.contains(&(i + 1)))
        .collect();
    if slots.is_empty() {
        stats.no_slot += 1;
        return None;
    }
    Some(Local {
        allowed,
        required,
        host,
        slots,
    })
}

#[derive(Debug)]
struct Candidate {
    source: usize,
    host: usize,
    slot: usize,
    route: Vec<usize>, // Interior only; host flanks are not inserted.
    support: u64,
    weight: u128,
    host_id: i64,
    route_ids: Vec<i64>,
}
impl Candidate {
    fn compare(&self, other: &Self) -> std::cmp::Ordering {
        self.support
            .cmp(&other.support)
            .then(self.weight.cmp(&other.weight))
            .then_with(|| other.host_id.cmp(&self.host_id))
            .then_with(|| other.slot.cmp(&self.slot))
            .then_with(|| other.route_ids.cmp(&self.route_ids))
    }
}
struct Search<'a> {
    input: &'a Input<'a>,
    local: &'a Local,
    paths: &'a Segments,
    source: usize,
    slot: usize,
    limits: Limits,
    remaining: usize,
    global: &'a mut usize,
    exhausted: bool,
    visited: Vec<bool>,
    route: Vec<usize>,
    best: Option<Candidate>,
    alternatives: usize,
    support_ties: usize,
}
impl Search<'_> {
    fn visit(
        &mut self,
        n: usize,
        target: usize,
        shared: Option<&RoaringBitmap>,
        required: usize,
        connectors: usize,
        weight: u128,
    ) {
        if self.remaining == 0 || *self.global == 0 {
            self.exhausted = true;
            return;
        }
        self.remaining -= 1;
        *self.global -= 1;
        if n == target {
            if required != 0 {
                return;
            }
            let support = shared.unwrap().len();
            let interior = self.route[1..self.route.len() - 1].to_vec();
            let candidate = Candidate {
                source: self.source,
                host: self.local.host,
                slot: self.slot,
                support,
                weight,
                host_id: self.input.ids[self.paths[self.local.host][0]],
                route_ids: interior.iter().map(|&n| self.input.ids[n]).collect(),
                route: interior,
            };
            self.alternatives += 1;
            match &self.best {
                Some(best) if support == best.support => self.support_ties += 1,
                Some(best) if support < best.support => (),
                _ => self.support_ties = 1,
            }
            log::debug!("Short segment insertion alternative: {:?}", candidate);
            if self
                .best
                .as_ref()
                .is_none_or(|best| candidate.compare(best).is_gt())
            {
                self.best = Some(candidate);
            }
            return;
        }
        for &next in &self.input.adjacency[n] {
            if self.visited[next] || (next != target && !self.local.allowed[next]) {
                continue;
            }
            let connector = next != target && !self.local.required[next];
            if connector && connectors == self.limits.connectors {
                continue;
            }
            let Some(Some(members)) = self.input.members.get(&key(n, next)) else {
                continue;
            };
            let mut intersection = members.clone();
            if let Some(shared) = shared {
                intersection &= shared;
            }
            if intersection.len() as u128 * 100 < MIN_SUPPORT_PERCENT * self.input.samples as u128 {
                continue;
            }
            self.visited[next] = true;
            self.route.push(next);
            let added_weight = if next == target {
                0
            } else {
                self.input.sizes[next].min(self.input.samples) as u128
            };
            self.visit(
                next,
                target,
                Some(&intersection),
                required - usize::from(self.local.required[next]),
                connectors + usize::from(connector),
                weight + added_weight,
            );
            self.route.pop();
            self.visited[next] = false;
            if self.exhausted {
                break;
            }
        }
    }
}

fn run(paths: &mut Segments, input: &Input<'_>, limits: Limits) -> Stats {
    let mut stats = Stats {
        unavailable_edges: input.members.values().filter(|m| m.is_none()).count(),
        ..Stats::default()
    };
    if stats.unavailable_edges > 0 {
        log::warn!(
            "Short segment insertion: {} normalised edges have missing, empty or malformed members and cannot support replacement routes",
            stats.unavailable_edges
        );
    }
    if input.members.values().all(Option::is_none) {
        return stats;
    }
    let mut global = limits.global;
    loop {
        let mut positions = vec![None; input.adjacency.len()];
        for (i, path) in paths.iter().enumerate() {
            for (j, &n) in path.iter().enumerate() {
                positions[n] = Some((i, j));
            }
        }
        let mut best: Option<Candidate> = None;
        for source in 0..paths.len() {
            if paths[source].len() >= MAX_SOURCE_GENES {
                continue;
            }
            let Some(local) = discover(source, paths, &positions, input, &mut stats, limits) else {
                continue;
            };
            stats.searched_sources += 1;
            let mut search = Search {
                input,
                local: &local,
                paths,
                source,
                slot: 0,
                limits,
                remaining: limits.per_source,
                global: &mut global,
                exhausted: false,
                visited: vec![false; input.adjacency.len()],
                route: Vec::new(),
                best: None,
                alternatives: 0,
                support_ties: 0,
            };
            for &slot in &local.slots {
                search.slot = slot;
                let start = paths[local.host][slot];
                let target = paths[local.host][slot + 1];
                search.visited[start] = true;
                search.route.push(start);
                search.visit(start, target, None, paths[source].len(), 0, 0);
                search.route.pop();
                search.visited[start] = false;
                if search.exhausted {
                    break;
                }
            }
            log::debug!(
                "Short source {}: {} supported alternatives, {} sharing best support, exhausted {}",
                input.ids[paths[source][0]],
                search.alternatives,
                search.support_ties,
                search.exhausted
            );
            if search.exhausted {
                stats.exhausted_sources += 1;
                if *search.global == 0 {
                    stats.global_exhausted = true;
                    break;
                }
                continue;
            }
            if let Some(candidate) = search.best {
                log::info!(
                    "Short source {}: {} supported alternatives, {} sharing best support",
                    input.ids[paths[source][0]],
                    search.alternatives,
                    search.support_ties
                );
                if best
                    .as_ref()
                    .is_none_or(|best| candidate.compare(best).is_gt())
                {
                    best = Some(candidate);
                }
            }
        }
        if stats.global_exhausted {
            break;
        }
        let Some(best) = best else {
            break;
        };
        let host = &paths[best.host];
        let flanks = (host[best.slot], host[best.slot + 1]);
        let old_support = input
            .members
            .get(&key(flanks.0, flanks.1))
            .and_then(Option::as_ref)
            .map(RoaringBitmap::len);
        let added = best.route.len() - paths[best.source].len();
        log::info!(
            "Absorb short segment starting at GML {} ({} genes) into host starting at GML {}, between GML {} and {}: route {:?}, shared support {}/{}, original edge support {:?}, added accessory genes {}",
            input.ids[paths[best.source][0]],
            paths[best.source].len(),
            best.host_id,
            input.ids[flanks.0],
            input.ids[flanks.1],
            best.route_ids,
            best.support,
            input.samples,
            old_support,
            added
        );
        paths[best.host].splice(best.slot + 1..best.slot + 1, best.route);
        paths.remove(best.source);
        stats.accepted += 1;
        stats.added_nodes += added;
    }
    stats.expansions = limits.global - global;
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Fixture {
        adjacency: Vec<Vec<usize>>,
        members: EdgeMembers,
        sizes: Vec<usize>,
        ids: Vec<i64>,
    }
    impl Fixture {
        fn new(n: usize, edges: &[(usize, usize)], support: &[u32]) -> Self {
            let mut adjacency = vec![Vec::new(); n];
            let mut members = BTreeMap::new();
            for &(a, b) in edges {
                adjacency[a].push(b);
                adjacency[b].push(a);
                members.insert(key(a, b), Some(support.iter().copied().collect()));
            }
            for neighbors in &mut adjacency {
                neighbors.sort_unstable();
                neighbors.dedup();
            }
            Self {
                adjacency,
                members,
                sizes: vec![100; n],
                ids: (0..n as i64).collect(),
            }
        }
        fn run(&self, paths: &mut Segments, limits: Limits) -> Stats {
            let mut core = vec![false; self.adjacency.len()];
            for &n in paths.iter().flatten() {
                core[n] = true;
            }
            let input = Input {
                adjacency: &self.adjacency,
                core: &core,
                sizes: &self.sizes,
                ids: self.ids.clone(),
                members: &self.members,
                samples: 100,
            };
            let before: BTreeSet<_> = paths.iter().flatten().copied().collect();
            let old_count = paths.len();
            let stats = run(paths, &input, limits);
            let after: BTreeSet<_> = paths.iter().flatten().copied().collect();
            assert!(before.is_subset(&after));
            assert_eq!(after.len(), paths.iter().map(Vec::len).sum());
            assert_eq!(old_count - paths.len(), stats.accepted);
            assert_eq!(after.len() - before.len(), stats.added_nodes);
            assert!(paths.iter().all(|p| {
                p.windows(2)
                    .all(|pair| self.adjacency[pair[0]].contains(&pair[1]))
            }));
            stats
        }
    }
    fn simple(support: &[u32]) -> Fixture {
        Fixture::new(5, &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2)], support)
    }
    #[test]
    fn exact_threshold_passes_even_with_stronger_original_edge() {
        let mut f = simple(&[0, 1, 2]);
        f.members.insert(key(1, 2), Some((0..100).collect()));
        let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
        assert_eq!(f.run(&mut paths, LIMITS).accepted, 1);
        assert_eq!(paths, vec![vec![0, 1, 4, 2, 3]]);
    }
    #[test]
    fn below_threshold_or_disjoint_edge_members_do_not_pass() {
        for disjoint in [false, true] {
            let mut f = simple(&[0, 1, 2]);
            f.members.insert(
                key(4, 2),
                Some(if disjoint {
                    [3, 4, 5].into_iter().collect()
                } else {
                    [0, 1].into_iter().collect()
                }),
            );
            let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
            assert_eq!(f.run(&mut paths, LIMITS).accepted, 0);
        }
    }
    #[test]
    fn inserts_all_source_genes_in_a_new_order_and_adds_connectors() {
        let f = Fixture::new(
            9,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (4, 5),
                (5, 6),
                (4, 6),
                (1, 5),
                (6, 2),
                (6, 7),
                (7, 8),
                (8, 2),
            ],
            &[0, 1, 2],
        );
        let mut paths = vec![vec![0, 1, 2, 3], vec![4, 5, 6]];
        let stats = f.run(&mut paths, LIMITS);
        assert_eq!(stats.accepted, 1);
        assert_eq!(stats.added_nodes, 2);
        assert_eq!(paths, vec![vec![0, 1, 5, 4, 6, 7, 8, 2, 3]]);
    }
    #[test]
    fn complete_topology_excludes_a_second_host_even_through_unsupported_nodes() {
        let mut f = Fixture::new(
            7,
            &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2), (4, 6), (6, 5)],
            &[0, 1, 2],
        );
        f.sizes[6] = 1;
        f.members.insert(key(6, 5), None);
        let mut paths = vec![vec![0, 1, 2, 3], vec![4], vec![5]];
        let stats = f.run(&mut paths, LIMITS);
        assert_eq!(stats.accepted, 0);
        assert!(stats.multiple_hosts > 0);
    }
    #[test]
    fn dead_end_branch_stays_separate() {
        let f = Fixture::new(6, &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 5)], &[0, 1, 2]);
        let mut paths = vec![vec![0, 1, 2, 3], vec![4, 5]];
        assert_eq!(f.run(&mut paths, LIMITS).accepted, 0);
    }
    #[test]
    fn conflicts_are_revalidated_after_highest_support_insertion() {
        let mut f = Fixture::new(
            6,
            &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2), (1, 5), (5, 2)],
            &[0, 1, 2],
        );
        for pair in [(1, 5), (5, 2)] {
            f.members
                .insert(key(pair.0, pair.1), Some((0..10).collect()));
        }
        let mut paths = vec![vec![0, 1, 2, 3], vec![4], vec![5]];
        assert_eq!(f.run(&mut paths, LIMITS).accepted, 1);
        assert_eq!(paths, vec![vec![0, 1, 5, 2, 3], vec![4]]);
    }
    #[test]
    fn support_ties_prefer_weight_then_route_ids_deterministically() {
        let mut f = Fixture::new(
            7,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (1, 4),
                (4, 5),
                (5, 2),
                (4, 6),
                (6, 2),
            ],
            &[0, 1, 2],
        );
        f.sizes[5] = 50;
        f.sizes[6] = 40;
        let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
        f.run(&mut paths, LIMITS);
        assert_eq!(paths, vec![vec![0, 1, 4, 5, 2, 3]]);
        f.sizes[6] = 50;
        f.ids[5] = 60;
        f.ids[6] = 50;
        for _ in 0..2 {
            let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
            f.run(&mut paths, LIMITS);
            assert_eq!(paths, vec![vec![0, 1, 4, 6, 2, 3]]);
        }
    }
    #[test]
    fn missing_members_do_not_fall_back_to_sizes() {
        let mut f = simple(&[0, 1, 2]);
        for m in f.members.values_mut() {
            *m = None;
        }
        let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
        let stats = f.run(&mut paths, LIMITS);
        assert_eq!(stats.accepted, 0);
        assert_eq!(stats.unavailable_edges, 5);
    }
    #[test]
    fn budgets_never_publish_a_partially_searched_candidate() {
        let f = Fixture::new(
            6,
            &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 2), (4, 5)],
            &[0, 1, 2],
        );
        // The source already has a valid candidate before its final prefix
        // expansion, but an incomplete search must discard it.
        for limits in [
            Limits {
                per_source: 4,
                ..LIMITS
            },
            Limits {
                global: 1,
                ..LIMITS
            },
            Limits {
                local_nodes: 0,
                ..LIMITS
            },
        ] {
            let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
            let stats = f.run(&mut paths, limits);
            assert_eq!(stats.accepted, 0);
            assert!(stats.exhausted_sources > 0 || stats.oversized > 0);
            if limits.global == 1 {
                assert!(stats.global_exhausted);
            }
        }
    }
    #[test]
    fn connector_limit_and_exclusive_source_length_cutoff() {
        let f = Fixture::new(
            6,
            &[(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 2)],
            &[0, 1, 2],
        );
        let mut paths = vec![vec![0, 1, 2, 3], vec![4]];
        assert_eq!(
            f.run(
                &mut paths,
                Limits {
                    connectors: 0,
                    ..LIMITS
                }
            )
            .accepted,
            0
        );
        for length in [24, 25] {
            let mut edges: Vec<_> = (2..length + 1).map(|n| (n, n + 1)).collect();
            edges.extend([(0, 1), (0, 2), (length + 1, 1)]);
            let f = Fixture::new(length + 2, &edges, &[0, 1, 2]);
            let mut paths = vec![vec![0, 1], (2..length + 2).collect()];
            assert_eq!(
                f.run(&mut paths, LIMITS).accepted,
                usize::from(length == 24)
            );
        }
    }
    #[test]
    fn parses_deduplicates_unions_and_rejects_invalid_edge_members() {
        let graph=BlockGraph::from_gml_str(r#"graph [
            node [id 10 size 100] node [id 20 size 100] node [id 30 size 100] node [id 40 size 100] node [id 50 size 100] node [id 60 size 100]
            edge [source 10 target 20 members 0 members 0 members 1]
            edge [source 20 target 10 members 2]
            edge [source 20 target 30 size 100]
            edge [source 30 target 40 members -1]
            edge [source 40 target 50 members "bad"]
            edge [source 50 target 60 members 0]
            edge [source 60 target 50 members 0.5]
            edge [source 10 target 60 members "_networkx_list_start" members 7]
        ]"#).unwrap();
        let members = edge_members(&graph);
        assert_eq!(
            members[&(0, 1)]
                .as_ref()
                .unwrap()
                .iter()
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(
            members[&(0, 5)]
                .as_ref()
                .unwrap()
                .iter()
                .collect::<Vec<_>>(),
            vec![7]
        );
        for edge in [(1, 2), (2, 3), (3, 4), (4, 5)] {
            assert!(members[&edge].is_none());
        }
    }
    #[test]
    fn graph_members_flow_through_public_pass_and_original_host_endpoints_survive() {
        let mut text = String::from("graph [");
        for n in 0..5 {
            text.push_str(&format!("node [id {n} size 100]"));
        }
        for (a, b) in [(0, 1), (1, 2), (2, 3), (1, 4), (4, 2)] {
            text.push_str(&format!(
                "edge [source {a} target {b} members 0 members 1 members 2]"
            ));
        }
        text.push(']');
        let graph = BlockGraph::from_gml_str(&text).unwrap();
        let mut paths = vec![vec![4], vec![0, 1, 2, 3]]; // source before host exercises removal indices.
        let stats = absorb(
            &mut paths,
            &graph,
            &graph.adjacency(),
            &[true; 5],
            &[100; 5],
            100,
            &edge_members(&graph),
        );
        assert_eq!(stats.accepted, 1);
        assert_eq!(paths, vec![vec![0, 1, 4, 2, 3]]);
    }
}
