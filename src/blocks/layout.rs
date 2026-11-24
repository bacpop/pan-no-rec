//! Bounded continuity improvements around the original layout incumbent.
use super::*;
mod absorb;
mod insert;
mod order;
mod refine;
mod support;

const BACKTRACK_DEPTH: usize = 2;
const ALTERNATIVES_PER_END: usize = 3;
const BACKTRACK_SEARCH_BUDGET: usize = 256;

type Segments = Vec<Vec<usize>>;

impl BlockGraph {
    /// Lay out every node once, growing backbone paths from both ends, retaining
    /// the original layout if the revised backbone has more graph breaks.
    pub fn layout_core_genes(&mut self, total_samples: usize) -> Result<Vec<i64>> {
        Ok(self.layout_with_metadata(total_samples)?.node_order)
    }
    pub(crate) fn layout_with_metadata(&mut self, total_samples: usize) -> Result<LayoutResult> {
        // Preserve the original RNG stream and suffix as well as its backbone.
        // This also performs all existing size/sample-count validation.
        let baseline: Vec<_> = self
            .layout_baseline(total_samples)?
            .into_iter()
            .map(|id| self.node_indices[&id])
            .collect();
        let sizes: Vec<_> = self
            .nodes
            .iter()
            .map(|n| n.size().map(|s| s.min(total_samples)))
            .collect::<Result<_>>()?;
        let core: Vec<_> = sizes
            .iter()
            .map(|&s| is_at_least_percent(s, total_samples, CORE_PERCENT))
            .collect();
        let eligible: Vec<_> = sizes
            .iter()
            .map(|&s| is_above_percent(s, total_samples, ELIGIBLE_PERCENT))
            .collect();
        let adjacency = self.adjacency();
        let end = baseline.iter().rposition(|&i| core[i]).map_or(0, |i| i + 1);
        let original = split_paths(&baseline[..end], &adjacency);
        let original_breaks = original.len().saturating_sub(1);
        let mut rng = LayoutRng::default();
        let mut search = Search {
            adjacency: &adjacency,
            sizes: &sizes,
            core: &core,
            eligible: &eligible,
            rng: &mut rng,
            retries: BACKTRACK_SEARCH_BUDGET,
            detour_budget: 0,
        };
        let revised = search.backbone();
        let retry_searches = BACKTRACK_SEARCH_BUDGET - search.retries;
        let revised_breaks = revised.len().saturating_sub(1);
        let improved = revised_breaks < original_breaks
            || (revised_breaks == original_breaks
                && weight(&revised, &sizes) > weight(&original, &sizes));
        let mut incumbent = if improved {
            revised.clone()
        } else {
            original.clone()
        };
        join_endpoints(&mut incumbent, &adjacency, &sizes);
        let before_refinement = incumbent.len().saturating_sub(1);
        let mut candidates = vec![("current", incumbent.clone()), ("revised", revised)];
        let mut long_rng = LayoutRng::default();
        let mut long_search = Search {
            adjacency: &adjacency,
            sizes: &sizes,
            core: &core,
            eligible: &eligible,
            rng: &mut long_rng,
            retries: BACKTRACK_SEARCH_BUDGET,
            detour_budget: refine::DETOUR_BUDGET,
        };
        let longer = long_search.backbone();
        log::info!(
            "Longer connectors: {} expansions, exhausted: {}",
            refine::DETOUR_BUDGET - long_search.detour_budget,
            long_search.detour_budget == 0
        );
        candidates.push(("longer", longer));
        let mut selected = "current (unmodified)";
        let mut seen = Vec::new();
        for (name, mut candidate) in candidates {
            if seen.contains(&candidate) {
                continue;
            }
            seen.push(candidate.clone());
            join_endpoints(&mut candidate, &adjacency, &sizes);
            let before = candidate.len().saturating_sub(1);
            let stats = refine::refine(&mut candidate, &adjacency);
            log::info!(
                "Layout refinement {}: {} -> {} breaks; {:?}",
                name,
                before,
                candidate.len().saturating_sub(1),
                stats
            );
            if layout_score(&candidate, &sizes) > layout_score(&incumbent, &sizes) {
                incumbent = candidate;
                selected = name;
            }
        }
        let mut segments = incumbent;
        let members = support::edge_members(self);
        let absorptions = absorb::absorb(
            &mut segments,
            self,
            &adjacency,
            &core,
            &sizes,
            total_samples,
            &members,
        );
        log::info!("Supported short segment insertions: {:?}", absorptions);
        let before_insertions: usize = segments.iter().map(Vec::len).sum();
        let node_ids: Vec<_> = self.nodes.iter().map(|node| node.id).collect();
        let insertions = insert::insert_accessory(
            &mut segments,
            &adjacency,
            &core,
            &sizes,
            &node_ids,
            total_samples,
        );
        log::info!(
            "Final accessory insertions: {:?}; backbone {} -> {} nodes",
            insertions,
            before_insertions,
            before_insertions + insertions.inserted_nodes
        );
        let ordering = order::order_segments(&mut segments, &adjacency, &members);
        log::info!("Backbone segment ordering: {:?}", ordering);
        log::info!(
            "Backbone gap distances: total {} -> {}, median {} -> {}, maximum {} -> {}, disconnected {} -> {}",
            ordering.before.total,
            ordering.after.total,
            ordering.before.median,
            ordering.after.median,
            ordering.before.maximum,
            ordering.after.maximum,
            ordering.before.disconnected,
            ordering.after.disconnected
        );
        let joined = before_refinement.saturating_sub(segments.len().saturating_sub(1));
        let backbone_len: usize = segments.iter().map(Vec::len).sum();
        let backbone_breaks = segments.len().saturating_sub(1);
        let mut offset = 0;
        let backbone_segment_starts = segments
            .iter()
            .map(|path| {
                let start = offset;
                offset += path.len();
                start
            })
            .collect();
        let mut order: Vec<_> = segments.into_iter().flatten().collect();
        let mut used = node_mask(self.nodes.len(), &order);
        let same_nodes = used == node_mask(self.nodes.len(), &baseline[..end]);
        if same_nodes {
            order.extend_from_slice(&baseline[end..]);
        } else {
            append_accessory(&mut order, &mut used, &adjacency, &sizes, &mut rng);
        }
        log::info!(
            "Layout backbone: {} nodes, {} continuous segments, {} breaks (baseline {}, revised {}); selected {}; break reduction vs incumbent {}; preserved accessory suffix: {}; extra backtracking searches: {}",
            backbone_len,
            backbone_breaks + usize::from(backbone_len > 0),
            backbone_breaks,
            original_breaks,
            revised_breaks,
            selected,
            joined,
            same_nodes,
            retry_searches
        );
        debug_assert_eq!(order.len(), self.nodes.len());
        debug_assert!(node_mask(self.nodes.len(), &order).iter().all(|&v| v));
        self.clear_coordinates();
        for (coordinate, &index) in order.iter().enumerate() {
            self.nodes[index].coordinate = Some(coordinate);
        }
        Ok(LayoutResult {
            node_order: order
                .into_iter()
                .map(|index| self.nodes[index].id)
                .collect(),
            backbone_len,
            backbone_segment_starts,
        })
    }
}

fn layout_score(paths: &Segments, sizes: &[usize]) -> (std::cmp::Reverse<usize>, usize, u128) {
    (
        std::cmp::Reverse(paths.len().saturating_sub(1)),
        weight(paths, sizes),
        paths.iter().map(|p| (p.len() as u128).pow(2)).sum(),
    )
}

fn weight(paths: &[Vec<usize>], sizes: &[usize]) -> usize {
    paths.iter().flatten().map(|&i| sizes[i]).sum()
}
fn split_paths(order: &[usize], adjacency: &[Vec<usize>]) -> Segments {
    let mut paths: Segments = Vec::new();
    for &node in order {
        if paths
            .last()
            .and_then(|p| p.last())
            .is_none_or(|&last| adjacency[last].binary_search(&node).is_err())
        {
            paths.push(Vec::new());
        }
        paths.last_mut().unwrap().push(node);
    }
    paths
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Connector {
    nodes: Vec<usize>,
    prepend: bool,
}
impl Connector {
    fn extend(&self, path: &mut Vec<usize>) {
        if self.prepend {
            let mut prefix = self.nodes.clone();
            prefix.reverse();
            prefix.extend_from_slice(path);
            *path = prefix;
        } else {
            path.extend_from_slice(&self.nodes);
        }
    }
}
#[derive(Clone)]
struct Decision {
    before: Vec<usize>,
    chosen: Connector,
}

struct Search<'a> {
    adjacency: &'a [Vec<usize>],
    sizes: &'a [usize],
    core: &'a [bool],
    eligible: &'a [bool],
    rng: &'a mut LayoutRng,
    retries: usize,
    detour_budget: usize,
}
impl Search<'_> {
    fn used(&self, locked: &[bool], path: &[usize]) -> Vec<bool> {
        let mut used = locked.to_vec();
        for &node in path {
            used[node] = true;
        }
        used
    }
    fn route(
        &mut self,
        start: usize,
        targets: &[bool],
        allowed: &[bool],
        used: &[bool],
        retry: bool,
    ) -> Option<Vec<usize>> {
        if retry {
            if self.retries == 0 {
                return None;
            }
            self.retries -= 1;
        }
        let route = best_path_to_core(
            start,
            targets,
            allowed,
            self.adjacency,
            self.sizes,
            self.eligible,
            used,
            self.rng,
        )?;
        if self.detour_budget == 0 {
            Some(route)
        } else {
            Some(refine::detours(
                start,
                route,
                used,
                allowed,
                self.adjacency,
                self.sizes,
                &mut self.detour_budget,
            ))
        }
    }
    // Check both endpoints with the >5% mask before permitting rare connectors.
    fn choices(
        &mut self,
        path: &[usize],
        used: &[bool],
        retry: bool,
        number: usize,
    ) -> Vec<Connector> {
        let targets: Vec<_> = self
            .core
            .iter()
            .enumerate()
            .map(|(i, &c)| c && !used[i])
            .collect();
        if !targets.iter().any(|&v| v) {
            return Vec::new();
        }
        let all = vec![true; used.len()];
        let eligible = self.eligible.to_vec();
        for allowed in [&eligible, &all] {
            let mut choices = Vec::new();
            for prepend in [false, true] {
                let start = if prepend {
                    path[0]
                } else {
                    *path.last().unwrap()
                };
                let mut remaining = targets.clone();
                for _ in 0..number {
                    let Some(nodes) = self.route(start, &remaining, allowed, used, retry) else {
                        break;
                    };
                    remaining[*nodes.last().unwrap()] = false;
                    let choice = Connector { nodes, prepend };
                    if !choices.contains(&choice) {
                        choices.push(choice);
                    }
                }
                // The main pass preserves successful append choices; the other
                // whole-path endpoint is tried when this endpoint cannot extend.
                if !retry && !choices.is_empty() {
                    return choices;
                }
            }
            if !choices.is_empty() {
                return choices;
            }
        }
        Vec::new()
    }
    fn grow(
        &mut self,
        path: &mut Vec<usize>,
        locked: &[bool],
        history: &mut Vec<Decision>,
        retry: bool,
    ) {
        loop {
            let used = self.used(locked, path);
            let Some(choice) = self.choices(path, &used, retry, 1).into_iter().next() else {
                break;
            };
            history.push(Decision {
                before: path.clone(),
                chosen: choice.clone(),
            });
            if history.len() > BACKTRACK_DEPTH {
                history.remove(0);
            }
            choice.extend(path);
        }
    }
    fn repair(
        &mut self,
        path: &[usize],
        locked: &[bool],
        history: &[Decision],
    ) -> Option<(Vec<usize>, Vec<Decision>)> {
        let mut best = path.to_vec();
        let mut best_history = Vec::new();
        for decision in history.iter().rev() {
            let used = self.used(locked, &decision.before);
            let alternatives = self.choices(&decision.before, &used, true, ALTERNATIVES_PER_END);
            for choice in alternatives {
                if choice == decision.chosen || self.retries == 0 {
                    continue;
                }
                let mut trial = decision.before.clone();
                let mut trial_history = vec![Decision {
                    before: trial.clone(),
                    chosen: choice.clone(),
                }];
                choice.extend(&mut trial);
                self.grow(&mut trial, locked, &mut trial_history, true);
                if self.score(&trial) > self.score(&best) {
                    best = trial;
                    best_history = trial_history;
                }
            }
        }
        // Extra accessory weight alone is not enough to accept a repair.
        if self.score(&best).0 > self.score(path).0 {
            Some((best, best_history))
        } else {
            None
        }
    }
    fn score(&self, path: &[usize]) -> (usize, usize) {
        (
            path.iter().filter(|&&i| self.core[i]).count(),
            path.iter().map(|&i| self.sizes[i]).sum(),
        )
    }
    fn backbone(&mut self) -> Segments {
        let mut used = vec![false; self.core.len()];
        let mut flat = Vec::new();
        loop {
            let remaining: Vec<_> = self
                .core
                .iter()
                .enumerate()
                .map(|(i, &c)| c && !used[i])
                .collect();
            let mut components = connected_components(self.adjacency, &remaining);
            if components.is_empty() {
                break;
            }
            components.sort_by_key(|c| std::cmp::Reverse(c.len()));
            let largest = components[0].len();
            let tied: Vec<_> = components
                .iter()
                .take_while(|c| c.len() == largest)
                .collect();
            let component = tied[self.rng.choose_index(&tied)];
            let mut path = greedy_path(
                self.adjacency,
                &node_mask(used.len(), component),
                self.sizes,
                self.rng,
            );
            let mut history = Vec::new();
            loop {
                self.grow(&mut path, &used, &mut history, false);
                if self.retries == 0 {
                    break;
                }
                if let Some((repaired, decisions)) = self.repair(&path, &used, &history) {
                    path = repaired;
                    history = decisions;
                } else {
                    break;
                }
            }
            append_nodes(path, &mut flat, &mut used);
        }
        split_paths(&flat, self.adjacency)
    }
}

// Reorient whole segments only; an interior edge is never a valid join.
fn join_endpoints(paths: &mut Segments, adjacency: &[Vec<usize>], sizes: &[usize]) {
    loop {
        let weights: Vec<usize> = paths
            .iter()
            .map(|p| p.iter().map(|&i| sizes[i]).sum())
            .collect();
        let mut best = None;
        let mut best_score = (0, 0);
        for a in 0..paths.len() {
            for b in a + 1..paths.len() {
                for reverse_a in [false, true] {
                    for reverse_b in [false, true] {
                        let x = if reverse_a {
                            paths[a][0]
                        } else {
                            *paths[a].last().unwrap()
                        };
                        let y = if reverse_b {
                            *paths[b].last().unwrap()
                        } else {
                            paths[b][0]
                        };
                        let score = (weights[a] + weights[b], paths[a].len() + paths[b].len());
                        if adjacency[x].binary_search(&y).is_ok()
                            && (best.is_none() || score > best_score)
                        {
                            best = Some((a, b, reverse_a, reverse_b));
                            best_score = score;
                        }
                    }
                }
            }
        }
        let Some((a, b, reverse_a, reverse_b)) = best else {
            break;
        };
        let mut second = paths.remove(b);
        if reverse_b {
            second.reverse();
        }
        if reverse_a {
            paths[a].reverse();
        }
        paths[a].extend(second);
    }
}

pub(super) fn append_accessory(
    order: &mut Vec<usize>,
    used: &mut [bool],
    adjacency: &[Vec<usize>],
    sizes: &[usize],
    rng: &mut LayoutRng,
) {
    loop {
        let remaining: Vec<_> = used.iter().map(|&u| !u).collect();
        let mut components = connected_components(adjacency, &remaining);
        if components.is_empty() {
            break;
        }
        components.sort_by(|a, b| {
            let aw: usize = a.iter().map(|&i| sizes[i]).sum();
            let bw: usize = b.iter().map(|&i| sizes[i]).sum();
            bw.cmp(&aw).then(b.len().cmp(&a.len()))
        });
        let best_weight: usize = components[0].iter().map(|&i| sizes[i]).sum();
        let best_count = components[0].len();
        let tied: Vec<_> = components
            .iter()
            .take_while(|c| {
                c.len() == best_count && c.iter().map(|&i| sizes[i]).sum::<usize>() == best_weight
            })
            .collect();
        let component = tied[rng.choose_index(&tied)];
        let path = greedy_path(adjacency, &node_mask(used.len(), component), sizes, rng);
        append_nodes(path, order, used);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn adjacency(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
        let mut result = vec![vec![]; n];
        for &(a, b) in edges {
            result[a].push(b);
            result[b].push(a);
        }
        for list in &mut result {
            list.sort_unstable();
        }
        result
    }
    fn breaks(order: &[i64], graph: &BlockGraph) -> usize {
        let adj = graph.adjacency();
        order
            .windows(2)
            .filter(|pair| {
                adj[graph.node_indices[&pair[0]]]
                    .binary_search(&graph.node_indices[&pair[1]])
                    .is_err()
            })
            .count()
    }
    #[test]
    fn extends_opposite_end_of_whole_path() {
        let adj = adjacency(4, &[(0, 1), (1, 2), (0, 3)]);
        let core = vec![true; 4];
        let sizes = vec![20; 4];
        let mut rng = LayoutRng::default();
        let mut search = Search {
            adjacency: &adj,
            sizes: &sizes,
            core: &core,
            eligible: &core,
            rng: &mut rng,
            retries: 256,
            detour_budget: 0,
        };
        let mut path = vec![0, 1, 2];
        search.grow(&mut path, &[false; 4], &mut vec![], false);
        assert_eq!(path, vec![3, 0, 1, 2]);
        assert_eq!(search.retries, 256);
    }
    #[test]
    fn eligible_route_at_other_end_precedes_rare_bridge() {
        let adj = adjacency(4, &[(0, 1), (1, 3), (3, 2), (0, 2)]);
        let core = vec![true, true, true, false];
        let eligible = core.clone();
        let sizes = vec![20, 20, 20, 1];
        let mut rng = LayoutRng::default();
        let mut search = Search {
            adjacency: &adj,
            sizes: &sizes,
            core: &core,
            eligible: &eligible,
            rng: &mut rng,
            retries: 256,
            detour_budget: 0,
        };
        let mut path = vec![0, 1];
        search.grow(&mut path, &[false; 4], &mut vec![], false);
        assert_eq!(path, vec![2, 0, 1]);
    }
    #[test]
    fn no_core_nodes_preserves_original_accessory_order() {
        let mut graph = BlockGraph::from_gml_str(
            "graph [node [id 1 size 2] node [id 2 size 1] edge [source 1 target 2]]",
        )
        .unwrap();
        let original = graph.layout_baseline(20).unwrap();
        assert_eq!(graph.layout_core_genes(20).unwrap(), original);
    }
    #[test]
    fn retries_recent_connection_without_leaking_used_nodes() {
        let adj = adjacency(5, &[(0, 1), (1, 2), (1, 3), (3, 4)]);
        let core = vec![true; 5];
        let sizes = vec![20; 5];
        let mut rng = LayoutRng::default();
        let mut search = Search {
            adjacency: &adj,
            sizes: &sizes,
            core: &core,
            eligible: &core,
            rng: &mut rng,
            retries: 256,
            detour_budget: 0,
        };
        let history = vec![Decision {
            before: vec![0, 1],
            chosen: Connector {
                nodes: vec![2],
                prepend: false,
            },
        }];
        let original = vec![0, 1, 2];
        let locked = vec![false; 5];
        let (repaired, _) = search.repair(&original, &locked, &history).unwrap();
        assert_eq!(repaired, vec![0, 1, 3, 4]);
        assert!(!search.used(&locked, &repaired)[2]);
        assert_eq!(original, vec![0, 1, 2]);
        assert!(locked.iter().all(|&v| !v));
        search.retries = 0;
        assert!(
            search
                .choices(&original, &search.used(&locked, &original), true, 1)
                .is_empty()
        );
    }
    #[test]
    fn joins_reversed_and_singleton_segments_but_not_interior_edges() {
        let sizes = vec![1; 5];
        let mut paths = vec![vec![0, 1, 2], vec![3, 4]];
        join_endpoints(
            &mut paths,
            &adjacency(5, &[(0, 1), (1, 2), (1, 3), (3, 4)]),
            &sizes,
        );
        assert_eq!(paths.len(), 2);
        join_endpoints(
            &mut paths,
            &adjacency(5, &[(0, 1), (1, 2), (2, 4), (3, 4)]),
            &sizes,
        );
        assert_eq!(paths, vec![vec![0, 1, 2, 4, 3]]);
        let mut singletons = vec![vec![0], vec![1], vec![2]];
        join_endpoints(&mut singletons, &adjacency(3, &[(0, 1)]), &[1; 3]);
        assert_eq!(singletons, vec![vec![0, 1], vec![2]]);
    }
    #[test]
    fn failed_repair_is_unchanged_and_budget_is_bounded() {
        let adj = adjacency(3, &[(0, 1), (1, 2)]);
        let core = vec![true; 3];
        let sizes = vec![20; 3];
        let mut rng = LayoutRng::default();
        let mut search = Search {
            adjacency: &adj,
            sizes: &sizes,
            core: &core,
            eligible: &core,
            rng: &mut rng,
            retries: 1,
            detour_budget: 0,
        };
        let history = vec![Decision {
            before: vec![0, 1],
            chosen: Connector {
                nodes: vec![2],
                prepend: false,
            },
        }];
        assert!(search.repair(&[0, 1, 2], &[false; 3], &history).is_none());
        assert_eq!(search.retries, 0);
    }
    #[test]
    fn layouts_cover_nodes_and_never_increase_backbone_breaks() {
        let mut rng = LayoutRng::default();
        for _ in 0..30 {
            let mut text = String::from("graph [");
            for i in 0..18 {
                let size = if i % 3 == 0 { 3 } else { 20 };
                text.push_str(&format!(" node [ id {i} size {size} ]"));
            }
            for i in 0..18 {
                for j in i + 1..18 {
                    if rng.next() % 9 == 0 {
                        text.push_str(&format!(" edge [ source {i} target {j} ]"));
                    }
                }
            }
            text.push(']');
            let mut graph = BlockGraph::from_gml_str(&text).unwrap();
            let baseline = graph.layout_baseline(20).unwrap();
            let core_end = |order: &[i64]| {
                order
                    .iter()
                    .rposition(|&id| graph.node_size(id).unwrap().unwrap() >= 19)
                    .map_or(0, |i| i + 1)
            };
            let base_end = core_end(&baseline);
            let before = breaks(&baseline[..base_end], &graph);
            let order = graph.layout_core_genes(20).unwrap();
            let end = order
                .iter()
                .rposition(|&id| graph.node_size(id).unwrap().unwrap() >= 19)
                .map_or(0, |i| i + 1);
            assert!(breaks(&order[..end], &graph) <= before);
            assert_eq!(order.len(), 18);
            assert_eq!(
                order.iter().collect::<std::collections::HashSet<_>>().len(),
                18
            );
            for (i, &id) in order.iter().enumerate() {
                assert_eq!(graph.node_coordinate(id), Some(i));
            }
            assert_eq!(order, graph.layout_core_genes(20).unwrap());
        }
    }
    #[test]
    #[ignore = "requires PMEN1_LAYOUT_GRAPH; writes benchmark orders only under /tmp"]
    fn pmen1_layout_benchmark() {
        let _ = simple_logger::init_with_level(log::Level::Info);
        let path = std::env::var("PMEN1_LAYOUT_GRAPH").expect("set PMEN1_LAYOUT_GRAPH");
        let mut graph = load_gml(Path::new(&path)).unwrap();
        let start = std::time::Instant::now();
        let baseline = graph.layout_baseline(221).unwrap();
        let baseline_time = start.elapsed();
        let start = std::time::Instant::now();
        let metadata = graph.layout_with_metadata(221).unwrap();
        let order = metadata.node_order.clone();
        let improved_time = start.elapsed();
        let end = |order: &[i64]| {
            order
                .iter()
                .rposition(|&id| graph.node_size(id).unwrap().unwrap().min(221) * 100 >= 221 * 95)
                .map_or(0, |i| i + 1)
        };
        let base_end = end(&baseline);
        let improved_end = metadata.backbone_len;
        let names: Vec<_> = order
            .iter()
            .map(|&id| graph.node_name(id).unwrap().to_owned())
            .collect();
        let segment_file = Path::new("/tmp/pmen1_layout_segments.tsv");
        crate::layout_segments::write(segment_file, &graph, &metadata, &names).unwrap();
        let segments = crate::layout_segments::read(segment_file, &names).unwrap();
        assert_eq!(segments.backbone_breaks().len(), 23);
        assert_eq!(metadata.backbone_len, 1927);
        // Compare against the independently cached Python experiment, including
        // whole-segment orientations, then replay from that experiment's input.
        let reference =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("dev/layout_viewer/endpoint-distances");
        let endpoints =
            std::fs::read_to_string(reference.join("pmen1-final-segments.endpoints.tsv")).unwrap();
        let endpoint_ids: HashMap<_, _> = endpoints
            .lines()
            .skip(1)
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                (
                    format!("S{}{}", f[0], if f[3] == "only" { "*" } else { f[3] }),
                    f[4].parse::<i64>().unwrap(),
                )
            })
            .collect();
        let expected =
            std::fs::read_to_string(reference.join("pmen1-supported-tie-segment-order.order.tsv"))
                .unwrap();
        let rows: Vec<_> = expected
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect();
        assert_eq!(rows.len(), metadata.backbone_segment_starts.len());
        let mut restored = vec![Vec::new(); rows.len()];
        for (i, row) in rows.iter().enumerate() {
            let start = metadata.backbone_segment_starts[i];
            let end = metadata
                .backbone_segment_starts
                .get(i + 1)
                .copied()
                .unwrap_or(metadata.backbone_len);
            assert_eq!(order[start], endpoint_ids[row[5]]);
            assert_eq!(order[end - 1], endpoint_ids[row[6]]);
            assert_eq!(end - start, row[3].parse::<usize>().unwrap());
            let mut path: Vec<_> = order[start..end]
                .iter()
                .map(|id| graph.node_index(*id).unwrap())
                .collect();
            if row[4] == "reverse" {
                path.reverse();
            }
            restored[row[1].parse::<usize>().unwrap() - 1] = path;
        }
        let stats = order::order_segments(
            &mut restored,
            &graph.adjacency(),
            &support::edge_members(&graph),
        );
        assert_eq!(
            (
                stats.mutual_unique,
                stats.support_ties,
                stats.deterministic_ties
            ),
            (18, 4, 1)
        );
        assert_eq!((stats.before.total, stats.after.total), (296, 178));
        assert_eq!((stats.after.median, stats.after.maximum), (5.0, 31));
        assert_eq!(
            restored
                .iter()
                .flatten()
                .map(|&i| graph.nodes[i].id)
                .collect::<Vec<_>>(),
            order[..metadata.backbone_len]
        );
        if baseline[..base_end]
            .iter()
            .collect::<std::collections::HashSet<_>>()
            == order[..improved_end]
                .iter()
                .collect::<std::collections::HashSet<_>>()
        {
            assert_eq!(baseline[base_end..], order[improved_end..]);
        }
        let old_breaks = breaks(&baseline[..base_end], &graph);
        let new_breaks = breaks(&order[..improved_end], &graph);
        println!(
            "Measured breaks: baseline {old_breaks}, refined {new_breaks}; backbone {improved_end}"
        );
        assert_eq!(old_breaks, 86);
        assert!(new_breaks <= 30, "PMEN1 regression: {new_breaks}");
        assert_eq!(order.len(), baseline.len());
        assert_eq!(
            order.iter().collect::<std::collections::HashSet<_>>().len(),
            order.len()
        );
        println!(
            "PMEN1 baseline: nodes={} backbone={} breaks={} time={:?}; revised: backbone={} breaks={} time={:?}; suffix unchanged={}",
            baseline.len(),
            base_end,
            old_breaks,
            baseline_time,
            improved_end,
            new_breaks,
            improved_time,
            baseline[base_end..] == order[improved_end..]
        );
        for (name, values) in [("baseline", &baseline), ("improved", &order)] {
            let data = values
                .iter()
                .enumerate()
                .map(|(i, &id)| format!("{}\t{}\t{}\n", i + 1, id, graph.node_name(id).unwrap()))
                .collect::<String>();
            std::fs::write(format!("/tmp/pmen1_layout_{name}.tsv"), data).unwrap();
        }
    }
}
