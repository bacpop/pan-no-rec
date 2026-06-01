use crate::get_progress_bar;
use hashbrown::HashMap;
use indicatif::ParallelProgressIterator;
use petgraph::graph::{NodeIndex, UnGraph};
use rayon::prelude::*;

use crate::dists::PairHits;
use crate::gene::Gene;
use crate::output::OutputRow;

type RecombinationGraph = UnGraph<usize, ()>;

// Builds the public presence table from pairwise recombinant hits.
pub fn presence_table_from_pair_hits(
    sample_count: usize,
    genes: &[Gene],
    hits: &PairHits,
    quiet: bool,
) -> Vec<OutputRow> {
    let progress_bar = get_progress_bar(genes.len(), false, quiet);
    genes
        .par_iter()
        .enumerate()
        .progress_with(progress_bar)
        .map(|(gene_index, _)| {
            let pairs = hits.get(&gene_index).map(Vec::as_slice).unwrap_or(&[]);

            OutputRow {
                gene_index,
                presence: infer_gene_presence(sample_count, pairs),
            }
        })
        .collect()
}

// Infers sample-level recombination presence from pairwise gene hits.
pub(crate) fn infer_gene_presence(sample_count: usize, pairs: &[(usize, usize)]) -> Vec<u8> {
    let graph = recombinant_pair_graph(sample_count, pairs);
    let mut presence = vec![0; sample_count];

    for sample_index in prune_recombinant_samples(&graph) {
        presence[sample_index] = 1;
    }

    presence
}

// Builds an undirected graph whose edges are recombinant sample pairs.
fn recombinant_pair_graph(sample_count: usize, pairs: &[(usize, usize)]) -> RecombinationGraph {
    let mut graph = RecombinationGraph::default();
    let mut sample_nodes = HashMap::new();

    for &(left_index, right_index) in pairs {
        if left_index >= sample_count || right_index >= sample_count || left_index == right_index {
            continue;
        }

        let left_node = *sample_nodes
            .entry(left_index)
            .or_insert_with(|| graph.add_node(left_index));
        let right_node = *sample_nodes
            .entry(right_index)
            .or_insert_with(|| graph.add_node(right_index));

        if graph.find_edge(left_node, right_node).is_none() {
            graph.add_edge(left_node, right_node, ());
        }
    }

    graph
}

// Prunes the graph to identify samples implicated by dense hit structure.
fn prune_recombinant_samples(graph: &RecombinationGraph) -> Vec<usize> {
    let mut active = vec![true; graph.node_count()];
    let mut active_count = graph.node_count();
    let mut recombinant_samples = Vec::new();

    while active_count > 3 {
        let degrees = active_degrees(graph, &active);
        let mut active_nodes: Vec<_> = graph
            .node_indices()
            .filter(|&node| is_active(&active, node))
            .map(|node| (node, degrees[node.index()]))
            .collect();

        active_nodes.sort_by(|(left_node, left_degree), (right_node, right_degree)| {
            right_degree
                .cmp(left_degree)
                .then_with(|| graph[*left_node].cmp(&graph[*right_node]))
        });

        let highest_degree = active_nodes[0].1;
        if highest_degree > 0
            && active_nodes
                .iter()
                .all(|(_, degree)| *degree == highest_degree)
        {
            remove_nodes(
                graph,
                &mut active,
                &mut active_count,
                active_nodes.into_iter().map(|(node, _)| node),
                &mut recombinant_samples,
            );
            break;
        }

        let second_highest_degree = active_nodes[1].1;
        if highest_degree > 3 * second_highest_degree {
            remove_nodes(
                graph,
                &mut active,
                &mut active_count,
                [active_nodes[0].0],
                &mut recombinant_samples,
            );
            continue;
        }

        if highest_degree == 0 {
            break;
        }

        let core_nodes = highest_non_empty_k_core(graph, &active, highest_degree);
        if core_nodes.is_empty() {
            break;
        }

        remove_nodes(
            graph,
            &mut active,
            &mut active_count,
            core_nodes,
            &mut recombinant_samples,
        );
    }

    recombinant_samples.sort_unstable();
    recombinant_samples
}

// Counts each active node's active neighbors.
fn active_degrees(graph: &RecombinationGraph, active: &[bool]) -> Vec<usize> {
    let mut degrees = vec![0; active.len()];

    for node in graph.node_indices().filter(|&node| is_active(active, node)) {
        degrees[node.index()] = graph
            .neighbors(node)
            .filter(|&neighbor| is_active(active, neighbor))
            .count();
    }

    degrees
}

// Finds the highest non-empty k-core among active nodes.
fn highest_non_empty_k_core(
    graph: &RecombinationGraph,
    active: &[bool],
    max_k: usize,
) -> Vec<NodeIndex> {
    for k in (1..=max_k).rev() {
        let core_nodes = k_core_nodes(graph, active, k);
        if !core_nodes.is_empty() {
            return core_nodes;
        }
    }

    Vec::new()
}

// Computes the active nodes that remain in a specific k-core.
fn k_core_nodes(graph: &RecombinationGraph, active: &[bool], k: usize) -> Vec<NodeIndex> {
    if k == 0 {
        return graph
            .node_indices()
            .filter(|&node| is_active(active, node))
            .collect();
    }

    let mut in_core = active.to_vec();
    loop {
        let mut changed = false;

        for node in graph.node_indices() {
            if !is_active(&in_core, node) {
                continue;
            }

            let active_degree = graph
                .neighbors(node)
                .filter(|&neighbor| is_active(&in_core, neighbor))
                .count();
            if active_degree < k {
                in_core[node.index()] = false;
                changed = true;
            }
        }

        if !changed {
            break;
        }
    }

    graph
        .node_indices()
        .filter(|&node| is_active(&in_core, node))
        .collect()
}

// Removes active graph nodes and records their sample indices.
fn remove_nodes(
    graph: &RecombinationGraph,
    active: &mut [bool],
    active_count: &mut usize,
    nodes: impl IntoIterator<Item = NodeIndex>,
    recombinant_samples: &mut Vec<usize>,
) {
    for node in nodes {
        if !is_active(active, node) {
            continue;
        }

        active[node.index()] = false;
        *active_count -= 1;
        recombinant_samples.push(graph[node]);
    }
}

// Checks whether a graph node is currently active.
fn is_active(active: &[bool], node: NodeIndex) -> bool {
    active.get(node.index()).copied().unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Builds a test graph from sample-indexed edges.
    fn graph_from_edges(node_count: usize, edges: &[(usize, usize)]) -> RecombinationGraph {
        let mut graph = RecombinationGraph::default();
        let nodes: Vec<_> = (0..node_count)
            .map(|sample| graph.add_node(sample))
            .collect();

        for &(left, right) in edges {
            graph.add_edge(nodes[left], nodes[right], ());
        }

        graph
    }

    // Builds a complete test graph with the requested node count.
    fn complete_graph(node_count: usize) -> RecombinationGraph {
        let mut edges = Vec::new();
        for left in 0..node_count {
            for right in left + 1..node_count {
                edges.push((left, right));
            }
        }

        graph_from_edges(node_count, &edges)
    }

    // Marks every node in a test graph as active.
    fn active_nodes(graph: &RecombinationGraph) -> Vec<bool> {
        vec![true; graph.node_count()]
    }

    // Converts graph nodes into sorted sample indices for assertions.
    fn node_samples(graph: &RecombinationGraph, nodes: Vec<NodeIndex>) -> Vec<usize> {
        let mut samples: Vec<_> = nodes.into_iter().map(|node| graph[node]).collect();
        samples.sort_unstable();
        samples
    }

    // Runs presence inference from sample-index pair fixtures.
    fn presence_for_pairs(sample_count: usize, pairs: &[(usize, usize)]) -> Vec<u8> {
        infer_gene_presence(sample_count, pairs)
    }

    #[test]
    // Verifies a complete graph remains intact at its maximum k-core.
    fn k_core_of_complete_graph_contains_all_nodes_at_max_k() {
        let graph = complete_graph(5);

        let observed = node_samples(&graph, k_core_nodes(&graph, &active_nodes(&graph), 4));

        assert_eq!(observed, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    // Verifies k-core search lowers k when a star has no max core.
    fn highest_non_empty_k_core_lowers_past_empty_star_max_core() {
        let graph = graph_from_edges(6, &[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]);
        let active = active_nodes(&graph);

        assert!(k_core_nodes(&graph, &active, 5).is_empty());
        let observed = node_samples(&graph, highest_non_empty_k_core(&graph, &active, 5));

        assert_eq!(observed, vec![0, 1, 2, 3, 4, 5]);
    }

    #[test]
    // Verifies k-core search lowers k for sparse path graphs.
    fn highest_non_empty_k_core_lowers_for_sparse_path() {
        let graph = graph_from_edges(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
        let active = active_nodes(&graph);

        assert!(k_core_nodes(&graph, &active, 2).is_empty());
        let observed = node_samples(&graph, highest_non_empty_k_core(&graph, &active, 2));

        assert_eq!(observed, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    // Verifies isolated residual nodes stop pruning cleanly.
    fn isolated_residual_nodes_do_not_loop_forever() {
        let graph = graph_from_edges(4, &[]);

        let observed = prune_recombinant_samples(&graph);

        assert!(observed.is_empty());
    }

    #[test]
    // Verifies small initial graphs do not mark recombination presence.
    fn pruning_skips_initial_graphs_with_at_most_three_nodes() {
        let observed = presence_for_pairs(4, &[(0, 1), (1, 2)]);

        assert_eq!(observed, vec![0, 0, 0, 0]);
    }

    #[test]
    // Verifies a dominant hub is marked without marking its leaves.
    fn hub_degree_more_than_three_times_second_highest_marks_only_hub() {
        let observed = presence_for_pairs(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]);

        assert_eq!(observed, vec![1, 0, 0, 0, 0]);
    }

    #[test]
    // Verifies regular dense structure marks all remaining nodes.
    fn regular_graph_marks_all_remaining_nodes() {
        let observed = presence_for_pairs(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);

        assert_eq!(observed, vec![1, 1, 1, 1]);
    }

    #[test]
    // Verifies pruning stops once three or fewer nodes remain.
    fn pruning_stops_when_removals_leave_at_most_three_nodes() {
        let observed = presence_for_pairs(5, &[(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)]);

        assert_eq!(observed, vec![1, 1, 1, 1, 0]);
    }

    #[test]
    // Verifies invalid and self-pair edges are skipped defensively.
    fn invalid_and_self_pairs_are_ignored() {
        let observed = presence_for_pairs(4, &[(0, 0), (0, 4), (4, 1), (0, 1), (1, 2)]);

        assert_eq!(observed, vec![0, 0, 0, 0]);
    }
}
