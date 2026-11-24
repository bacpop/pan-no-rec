//! Gap-preserving proximity ordering of complete backbone segments.
use super::support::EdgeMembers;
use super::*;
use rayon::prelude::*;

const UNREACHABLE: usize = usize::MAX;
#[derive(Debug, Default)]
pub(super) struct GapStats {
    pub total: usize,
    pub median: f64,
    pub maximum: usize,
    pub disconnected: usize,
}
#[derive(Debug, Default)]
pub(super) struct Stats {
    pub before: GapStats,
    pub after: GapStats,
    pub mutual_unique: usize,
    pub support_ties: usize,
    pub deterministic_ties: usize,
    pub rounds: usize,
    pub components: usize,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Oriented {
    segment: usize,
    reversed: bool,
}
type Chain = Vec<Oriented>;
struct Row {
    distances: Vec<usize>,
    supports: Vec<u64>,
}
struct Metrics {
    nodes: Vec<usize>,
    ends: Vec<[usize; 2]>,
    rows: Vec<Row>,
}
// FIFO layers ensure all equal-distance predecessors update support before a
// node is expanded. This finds the widest bottleneck among ALL shortest paths.
fn endpoint_row(
    source: usize,
    nodes: &[usize],
    adjacency: &[Vec<usize>],
    weights: &[Vec<u64>],
) -> Row {
    let mut distance = vec![UNREACHABLE; adjacency.len()];
    let mut support = vec![0; adjacency.len()];
    let mut queue = VecDeque::from([source]);
    distance[source] = 0;
    support[source] = u64::MAX;
    while let Some(n) = queue.pop_front() {
        for (edge, &other) in adjacency[n].iter().enumerate() {
            let next = distance[n] + 1;
            if distance[other] == UNREACHABLE {
                distance[other] = next;
                queue.push_back(other);
            }
            if distance[other] == next {
                support[other] = support[other].max(support[n].min(weights[n][edge]));
            }
        }
    }
    Row {
        distances: nodes.iter().map(|&n| distance[n]).collect(),
        supports: nodes.iter().map(|&n| support[n]).collect(),
    }
}
impl Metrics {
    fn new(paths: &Segments, adjacency: &[Vec<usize>], members: &EdgeMembers) -> Self {
        let mut nodes = Vec::new();
        let mut ends = Vec::new();
        for p in paths {
            let left = nodes.len();
            nodes.push(p[0]);
            let right = if p.len() == 1 {
                left
            } else {
                nodes.push(*p.last().unwrap());
                nodes.len() - 1
            };
            ends.push([left, right]);
        }
        let weights: Vec<Vec<_>> = adjacency
            .iter()
            .enumerate()
            .map(|(a, neighbors)| {
                neighbors
                    .iter()
                    .map(|&b| {
                        members
                            .get(&(a.min(b), a.max(b)))
                            .and_then(Option::as_ref)
                            .map_or(0, |s| s.len())
                    })
                    .collect()
            })
            .collect();
        let rows = nodes
            .par_iter()
            .map(|&n| endpoint_row(n, &nodes, adjacency, &weights))
            .collect();
        Self { nodes, ends, rows }
    }
    fn chain_ends(&self, chain: &[Oriented]) -> [usize; 2] {
        let first = chain[0];
        let last = *chain.last().unwrap();
        [
            self.ends[first.segment][usize::from(first.reversed)],
            self.ends[last.segment][usize::from(!last.reversed)],
        ]
    }
    fn gaps(&self, order: &[Oriented]) -> GapStats {
        let mut finite = Vec::new();
        let mut disconnected = 0;
        for pair in order.windows(2) {
            let a = self.chain_ends(&pair[..1])[1];
            let b = self.chain_ends(&pair[1..])[0];
            let d = self.rows[a].distances[b];
            if d == UNREACHABLE {
                disconnected += 1;
            } else {
                finite.push(d);
            }
        }
        finite.sort_unstable();
        let count = finite.len();
        GapStats {
            total: finite.iter().sum(),
            maximum: finite.last().copied().unwrap_or(0),
            median: if count == 0 {
                0.0
            } else {
                (finite[(count - 1) / 2] as f64 + finite[count / 2] as f64) / 2.0
            },
            disconnected,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Choice {
    a: usize,
    b: usize,
    sa: usize,
    sb: usize,
    ea: usize,
    eb: usize,
    distance: usize,
}
impl Choice {
    fn key(&self) -> (usize, usize) {
        (self.ea.min(self.eb), self.ea.max(self.eb))
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Rule {
    MutualUnique,
    Support,
    Deterministic,
}
impl Rule {
    fn label(self) -> &'static str {
        match self {
            Self::MutualUnique => "mutual unique",
            Self::Support => "bottleneck support",
            Self::Deterministic => "deterministic IDs",
        }
    }
}
fn reverse(chain: &mut Chain) {
    chain.reverse();
    for part in chain {
        part.reversed = !part.reversed;
    }
}
fn candidates(chains: &[Chain], metrics: &Metrics) -> Vec<Choice> {
    let mut result = Vec::new();
    for (a, first) in chains.iter().enumerate() {
        let aa = metrics.chain_ends(first);
        for (b, second) in chains.iter().enumerate().skip(a + 1) {
            let bb = metrics.chain_ends(second);
            for sa in 0..2 {
                if sa == 1 && aa[0] == aa[1] {
                    continue;
                }
                for sb in 0..2 {
                    if sb == 1 && bb[0] == bb[1] {
                        continue;
                    }
                    let distance = metrics.rows[aa[sa]].distances[bb[sb]];
                    if distance != UNREACHABLE {
                        result.push(Choice {
                            a,
                            b,
                            sa,
                            sb,
                            ea: aa[sa],
                            eb: bb[sb],
                            distance,
                        });
                    }
                }
            }
        }
    }
    result
}
fn arrange(metrics: &Metrics, stats: &mut Stats) -> Vec<Oriented> {
    let mut chains: Vec<Chain> = (0..metrics.ends.len())
        .map(|segment| {
            vec![Oriented {
                segment,
                reversed: false,
            }]
        })
        .collect();
    loop {
        let choices = candidates(&chains, metrics);
        if choices.is_empty() {
            break;
        }
        let mut nearest = vec![None; chains.len()];
        let mut counts = vec![0; chains.len()];
        for (index, c) in choices.iter().enumerate() {
            for n in [c.a, c.b] {
                match nearest[n] {
                    None => {
                        nearest[n] = Some(index);
                        counts[n] = 1;
                    }
                    Some(old) if c.distance < choices[old].distance => {
                        nearest[n] = Some(index);
                        counts[n] = 1;
                    }
                    Some(old) if c.distance == choices[old].distance => counts[n] += 1,
                    _ => (),
                }
            }
        }
        let mut batch: Vec<_> = choices
            .iter()
            .enumerate()
            .filter_map(|(i, c)| {
                (counts[c.a] == 1
                    && counts[c.b] == 1
                    && nearest[c.a] == Some(i)
                    && nearest[c.b] == Some(i))
                .then_some(*c)
            })
            .collect();
        let rule;
        if batch.is_empty() {
            let minimum = choices.iter().map(|c| c.distance).min().unwrap();
            let tied: Vec<_> = choices
                .iter()
                .filter(|c| c.distance == minimum)
                .copied()
                .collect();
            let support = tied
                .iter()
                .map(|c| metrics.rows[c.ea].supports[c.eb])
                .max()
                .unwrap();
            let mut finalists: Vec<_> = tied
                .into_iter()
                .filter(|c| metrics.rows[c.ea].supports[c.eb] == support)
                .collect();
            finalists.sort_by_key(Choice::key);
            rule = if finalists.len() == 1 {
                Rule::Support
            } else {
                Rule::Deterministic
            };
            batch.push(finalists[0]);
        } else {
            batch.sort_by_key(|c| (c.distance, c.key()));
            rule = Rule::MutualUnique;
        }
        stats.rounds += 1;
        // Mutual pairs are chain-disjoint; stable segment anchors survive index
        // shifts while each accepted pairing is applied separately.
        let pending: Vec<_> = batch
            .iter()
            .map(|c| (*c, chains[c.a][0].segment, chains[c.b][0].segment))
            .collect();
        for (choice, anchor_a, anchor_b) in pending {
            let a = chains
                .iter()
                .position(|c| c.iter().any(|p| p.segment == anchor_a))
                .unwrap();
            let b = chains
                .iter()
                .position(|c| c.iter().any(|p| p.segment == anchor_b))
                .unwrap();
            let mut first = chains[a].clone();
            let mut second = chains[b].clone();
            if choice.sa == 0 {
                reverse(&mut first);
            }
            if choice.sb == 1 {
                reverse(&mut second);
            }
            first.extend(second);
            chains.remove(a.max(b));
            chains.remove(a.min(b));
            chains.push(first);
            chains.sort_by_key(|c| c.iter().map(|p| p.segment).min().unwrap());
            match rule {
                Rule::MutualUnique => stats.mutual_unique += 1,
                Rule::Support => stats.support_ties += 1,
                Rule::Deterministic => stats.deterministic_ties += 1,
            }
            log::debug!(
                "Segment ordering round {}: {} between GML {} and {}, distance {}, bottleneck {}",
                stats.rounds,
                rule.label(),
                metrics.nodes[choice.ea],
                metrics.nodes[choice.eb],
                choice.distance,
                metrics.rows[choice.ea].supports[choice.eb]
            );
        }
    }
    stats.components = chains.len();
    for chain in &mut chains {
        let ends = metrics.chain_ends(chain);
        if ends[0] > ends[1] {
            reverse(chain);
        }
        for part in chain {
            if metrics.ends[part.segment][0] == metrics.ends[part.segment][1] {
                part.reversed = false;
            }
        }
    }
    chains.into_iter().flatten().collect()
}
pub(super) fn order_segments(
    paths: &mut Segments,
    adjacency: &[Vec<usize>],
    members: &EdgeMembers,
) -> Stats {
    let mut stats = Stats::default();
    if paths.len() < 2 {
        stats.components = paths.len();
        return stats;
    }
    let metrics = Metrics::new(paths, adjacency, members);
    let original: Vec<_> = (0..paths.len())
        .map(|segment| Oriented {
            segment,
            reversed: false,
        })
        .collect();
    stats.before = metrics.gaps(&original);
    let ordered = arrange(&metrics, &mut stats);
    stats.after = metrics.gaps(&ordered);
    let mut storage: Vec<_> = std::mem::take(paths).into_iter().map(Some).collect();
    for part in ordered {
        let mut path = storage[part.segment].take().unwrap();
        if part.reversed {
            path.reverse();
        }
        paths.push(path);
    }
    stats
}

#[cfg(test)]
mod tests {
    use super::*;
    fn graph(n: usize, edges: &[(usize, usize, u32)]) -> (Vec<Vec<usize>>, EdgeMembers) {
        let mut adj = vec![Vec::new(); n];
        let mut members = std::collections::BTreeMap::new();
        for &(a, b, support) in edges {
            adj[a].push(b);
            adj[b].push(a);
            members.insert((a.min(b), a.max(b)), Some((0..support).collect()));
        }
        for ns in &mut adj {
            ns.sort_unstable();
            ns.dedup();
        }
        (adj, members)
    }
    fn matrix(n: usize, edges: &[(usize, usize, u64)]) -> Metrics {
        let mut rows: Vec<_> = (0..n)
            .map(|_| Row {
                distances: vec![10; n],
                supports: vec![0; n],
            })
            .collect();
        for (i, row) in rows.iter_mut().enumerate() {
            row.distances[i] = 0;
        }
        for &(a, b, support) in edges {
            rows[a].distances[b] = 1;
            rows[b].distances[a] = 1;
            rows[a].supports[b] = support;
            rows[b].supports[a] = support;
        }
        Metrics {
            nodes: (0..n).collect(),
            ends: (0..n).map(|n| [n, n]).collect(),
            rows,
        }
    }
    #[test]
    fn breadth_first_support_considers_all_shortest_paths_and_later_predecessors() {
        let (adj, members) = graph(5, &[(0, 1, 2), (1, 3, 2), (0, 2, 8), (2, 3, 7), (3, 4, 6)]);
        let paths = vec![vec![0], vec![3], vec![4]];
        let metrics = Metrics::new(&paths, &adj, &members);
        assert_eq!(metrics.rows[0].distances, vec![0, 2, 3]);
        assert_eq!(metrics.rows[0].supports[1..], [7, 6]);
    }
    #[test]
    fn stronger_longer_route_does_not_override_shortest_distance() {
        let (adj, members) = graph(4, &[(0, 3, 1), (0, 1, 100), (1, 2, 100), (2, 3, 100)]);
        let metrics = Metrics::new(&vec![vec![0], vec![3]], &adj, &members);
        assert_eq!(metrics.rows[0].distances[1], 1);
        assert_eq!(metrics.rows[0].supports[1], 1);
    }
    #[test]
    fn tied_distances_use_support_then_stable_endpoint_indices() {
        let metrics = matrix(3, &[(0, 1, 3), (0, 2, 8), (1, 2, 4)]);
        let mut stats = Stats::default();
        let order = arrange(&metrics, &mut stats);
        assert_eq!(
            order.iter().map(|p| p.segment).collect::<Vec<_>>(),
            vec![0, 2, 1]
        );
        assert_eq!(stats.support_ties, 2);
        assert_eq!(stats.deterministic_ties, 0);
        let metrics = matrix(3, &[(0, 1, 0), (0, 2, 0), (1, 2, 0)]);
        let mut stats = Stats::default();
        let order = arrange(&metrics, &mut stats);
        assert_eq!(
            order.iter().map(|p| p.segment).collect::<Vec<_>>(),
            vec![1, 0, 2]
        );
        assert_eq!(stats.deterministic_ties, 2);
    }
    #[test]
    fn mutual_unique_pairs_precede_ambiguous_shorter_connections() {
        let mut metrics = matrix(5, &[(0, 1, 0), (0, 2, 0), (1, 2, 0)]);
        metrics.rows[3].distances[4] = 5;
        metrics.rows[4].distances[3] = 5;
        let mut stats = Stats::default();
        let order = arrange(&metrics, &mut stats);
        assert_eq!(stats.mutual_unique, 1);
        assert_eq!(stats.deterministic_ties, 3);
        assert_eq!(order.len(), 5);
    }
    #[test]
    fn whole_segments_can_reverse_without_changing_membership_or_boundaries() {
        let (adj, members) = graph(
            6,
            &[(0, 1, 10), (2, 3, 10), (4, 5, 10), (1, 3, 10), (2, 4, 5)],
        );
        let original = vec![vec![0, 1], vec![2, 3], vec![4, 5]];
        let mut paths = original.clone();
        let stats = order_segments(&mut paths, &adj, &members);
        assert_eq!(paths, vec![vec![0, 1], vec![3, 2], vec![4, 5]]);
        assert_eq!(stats.after.total, 2);
        assert_eq!(paths.len(), 3);
        let mut sorted: Vec<_> = paths.iter().flatten().copied().collect();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..6).collect::<Vec<_>>());
        for p in &paths {
            assert!(
                original
                    .iter()
                    .any(|q| p == q || p.iter().eq(q.iter().rev()))
            );
        }
    }
    #[test]
    fn disconnected_components_stay_grouped_with_explicit_infinite_gaps() {
        let (adj, members) = graph(6, &[(0, 1, 2), (2, 3, 2), (4, 5, 2)]);
        let mut paths = vec![vec![0, 1], vec![2], vec![4, 5], vec![3]];
        let stats = order_segments(&mut paths, &adj, &members);
        assert_eq!(paths, vec![vec![0, 1], vec![2], vec![3], vec![4, 5]]);
        assert_eq!(stats.components, 3);
        assert_eq!(stats.after.disconnected, 2);
        assert_eq!(stats.after.total, 1);
    }
    #[test]
    fn missing_support_does_not_remove_topological_connections() {
        let (adj, mut members) = graph(3, &[(0, 1, 2), (1, 2, 2)]);
        for value in members.values_mut() {
            *value = None;
        }
        let metrics = Metrics::new(&vec![vec![0], vec![2]], &adj, &members);
        assert_eq!(metrics.rows[0].distances[1], 2);
        assert_eq!(metrics.rows[0].supports[1], 0);
    }
    #[test]
    fn empty_singleton_and_thread_counts_are_deterministic() {
        let (adj, members) = graph(
            6,
            &[(0, 1, 10), (2, 3, 10), (4, 5, 10), (1, 3, 10), (2, 4, 5)],
        );
        let mut empty = vec![];
        assert_eq!(order_segments(&mut empty, &adj, &members).components, 0);
        let mut single = vec![vec![0, 1]];
        assert_eq!(order_segments(&mut single, &adj, &members).components, 1);
        assert_eq!(single, vec![vec![0, 1]]);
        let run = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let mut paths = vec![vec![0, 1], vec![2, 3], vec![4, 5]];
                    order_segments(&mut paths, &adj, &members);
                    paths
                })
        };
        assert_eq!(run(1), run(4));
    }
}
