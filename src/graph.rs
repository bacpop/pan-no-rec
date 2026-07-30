use crate::gene::GeneMetadata;
use crate::get_progress_bar;
use anyhow::{Context, Result, bail};
use hashbrown::{HashMap, HashSet};
use indicatif::ParallelProgressIterator;
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use rayon::prelude::*;
use roaring::RoaringBitmap;

use crate::dists::PairHits;
use crate::output::OutputRow;

use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

type RecombinationGraph = UnGraph<usize, ()>;

#[derive(Clone, Debug)]
struct ComponentWork {
    id: usize,
    nodes: Vec<NodeIndex>,
    iteration: usize,
}

#[derive(Debug)]
pub(crate) struct GraphDebugConfig {
    output_dir: PathBuf,
    gene_names: HashSet<String>,
}

impl GraphDebugConfig {
    pub(crate) fn new(output_dir: PathBuf, gene_names: Vec<String>) -> Result<Self> {
        if gene_names.is_empty() {
            bail!("--debug-graph-dir requires at least one --debug-graph-genes value");
        }

        let mut selected_genes = HashSet::with_capacity(gene_names.len());
        for gene_name in gene_names {
            let gene_name = gene_name.trim();
            if gene_name.is_empty() {
                bail!("debug graph gene name must not be empty");
            }
            selected_genes.insert(gene_name.to_string());
        }

        fs::create_dir_all(&output_dir).with_context(|| {
            format!(
                "failed to create debug graph directory '{}'",
                output_dir.display()
            )
        })?;

        Ok(GraphDebugConfig {
            output_dir,
            gene_names: selected_genes,
        })
    }

    fn for_gene<'a>(
        &'a self,
        gene_index: usize,
        gene_name: &'a str,
        sample_names: &'a [String],
    ) -> Option<GeneGraphDebug<'a>> {
        self.gene_names
            .contains(gene_name)
            .then_some(GeneGraphDebug {
                output_dir: &self.output_dir,
                gene_index,
                gene_name,
                sample_names,
            })
    }
}

struct GeneGraphDebug<'a> {
    output_dir: &'a Path,
    gene_index: usize,
    gene_name: &'a str,
    sample_names: &'a [String],
}

// Builds the public presence table from pairwise recombinant hits.
pub fn presence_table_from_pair_hits(
    sample_names: &[String],
    genes: &[GeneMetadata],
    hits: &PairHits,
    debug_graphs: Option<&GraphDebugConfig>,
    quiet: bool,
) -> Result<Vec<OutputRow>> {
    let sample_count = sample_names.len();
    let progress_bar = get_progress_bar(genes.len(), false, quiet);
    (0..genes.len())
        .into_par_iter()
        .progress_with(progress_bar)
        .map(|gene_index| -> Result<OutputRow> {
            let gene = &genes[gene_index];
            let debug = debug_graphs
                .and_then(|config| config.for_gene(gene_index, gene.name(), sample_names));
            let presence = match hits.get(&gene_index) {
                Some(pair_offsets) => {
                    infer_gene_presence(sample_count, pair_offsets, debug.as_ref())?
                }
                None => {
                    if let Some(debug) = debug.as_ref() {
                        let empty_pair_offsets = RoaringBitmap::new();
                        infer_gene_presence(sample_count, &empty_pair_offsets, Some(debug))?
                    } else {
                        vec![0; sample_count]
                    }
                }
            };

            Ok(OutputRow {
                gene_index,
                presence,
            })
        })
        .collect()
}

// Infers sample-level recombination presence from pairwise gene hits.
fn infer_gene_presence(
    sample_count: usize,
    pair_offsets: &RoaringBitmap,
    debug: Option<&GeneGraphDebug<'_>>,
) -> Result<Vec<u8>> {
    let graph = recombinant_pair_graph(sample_count, pair_offsets);
    let mut presence = vec![0; sample_count];

    for sample_index in prune_recombinant_samples(&graph, debug)? {
        presence[sample_index] = 1;
    }

    Ok(presence)
}

// Builds an undirected graph whose edges are recombinant sample pairs.
fn recombinant_pair_graph(sample_count: usize, pair_offsets: &RoaringBitmap) -> RecombinationGraph {
    let mut graph = RecombinationGraph::default();
    let mut sample_nodes = HashMap::new();

    for pair_offset in pair_offsets {
        let (left_index, right_index) = sample_pair_indices(sample_count, pair_offset as usize);
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

// Maps a flat upper-triangle offset to sample indices.
fn sample_pair_indices(sample_count: usize, pair_offset: usize) -> (usize, usize) {
    debug_assert!(sample_count >= 2);

    let mut low = 0;
    let mut high = sample_count.saturating_sub(2);
    while low < high {
        let midpoint = (low + high).div_ceil(2);
        if pairs_before_sample(sample_count, midpoint) <= pair_offset {
            low = midpoint;
        } else {
            high = midpoint - 1;
        }
    }

    let sample_a = low;
    let sample_b = sample_a + 1 + pair_offset - pairs_before_sample(sample_count, sample_a);
    (sample_a, sample_b)
}

// Counts pair offsets before a sample's row in the upper triangle.
fn pairs_before_sample(sample_count: usize, sample_index: usize) -> usize {
    sample_index * (2 * sample_count - sample_index - 1) / 2
}

// Prunes the graph to identify samples implicated by dense hit structure.
fn prune_recombinant_samples(
    graph: &RecombinationGraph,
    debug: Option<&GeneGraphDebug<'_>>,
) -> Result<Vec<usize>> {
    let mut active = vec![true; graph.node_count()];
    let mut recombinant_samples = Vec::new();
    let initial_components = active_components(graph, &active, graph.node_indices());

    if initial_components.is_empty() {
        write_debug_snapshot(debug, graph, None, 0, "initial", &active, &[], &[])?;
        return Ok(recombinant_samples);
    }

    for (component_id, nodes) in initial_components.into_iter().enumerate() {
        prune_component(
            graph,
            &mut active,
            ComponentWork {
                id: component_id,
                nodes,
                iteration: 0,
            },
            debug,
            &mut recombinant_samples,
        )?;
    }

    recombinant_samples.sort_unstable();
    Ok(recombinant_samples)
}

fn prune_component(
    graph: &RecombinationGraph,
    active: &mut [bool],
    mut component: ComponentWork,
    debug: Option<&GeneGraphDebug<'_>>,
    recombinant_samples: &mut Vec<usize>,
) -> Result<()> {
    sort_nodes_by_sample(graph, &mut component.nodes);
    write_debug_snapshot(
        debug,
        graph,
        Some(component.id),
        component.iteration,
        "initial",
        active,
        &component.nodes,
        &[],
    )?;

    loop {
        let active_subcomponents =
            active_components(graph, active, component.nodes.iter().copied());
        let mut removals = Vec::new();

        for subcomponent in active_subcomponents {
            if subcomponent.len() < 3 {
                continue;
            }

            let (max_degree, nodes_to_remove) = max_degree_nodes(graph, active, &subcomponent);
            removals.push((subcomponent, max_degree, nodes_to_remove));
        }

        if removals.is_empty() {
            component.iteration += 1;
            let mut remaining_nodes: Vec<_> = component
                .nodes
                .iter()
                .copied()
                .filter(|&node| is_active(active, node))
                .collect();
            sort_nodes_by_sample(graph, &mut remaining_nodes);
            write_debug_snapshot(
                debug,
                graph,
                Some(component.id),
                component.iteration,
                "stop_too_small",
                active,
                &remaining_nodes,
                &[],
            )?;
            break;
        }

        for (subcomponent, max_degree, nodes_to_remove) in &removals {
            component.iteration += 1;
            let step = format!("remove_max_degree_d{max_degree}");
            write_debug_snapshot(
                debug,
                graph,
                Some(component.id),
                component.iteration,
                &step,
                active,
                subcomponent,
                nodes_to_remove,
            )?;
        }

        let mut removed_nodes = Vec::new();
        for (_, _, nodes_to_remove) in removals {
            removed_nodes.extend(nodes_to_remove);
        }
        sort_nodes_by_sample(graph, &mut removed_nodes);
        remove_nodes(graph, active, removed_nodes, recombinant_samples);
    }

    Ok(())
}

fn max_degree_nodes(
    graph: &RecombinationGraph,
    active: &[bool],
    component_nodes: &[NodeIndex],
) -> (usize, Vec<NodeIndex>) {
    let degrees = component_degrees(graph, active, component_nodes);
    let max_degree = component_nodes
        .iter()
        .copied()
        .filter(|&node| is_active(active, node))
        .map(|node| degrees[node.index()])
        .max()
        .unwrap_or(0);
    let mut nodes: Vec<_> = component_nodes
        .iter()
        .copied()
        .filter(|&node| is_active(active, node) && degrees[node.index()] == max_degree)
        .collect();
    sort_nodes_by_sample(graph, &mut nodes);
    (max_degree, nodes)
}

fn write_debug_snapshot(
    debug: Option<&GeneGraphDebug<'_>>,
    graph: &RecombinationGraph,
    component_id: Option<usize>,
    iteration: usize,
    step: &str,
    active: &[bool],
    component_nodes: &[NodeIndex],
    removed_this_step: &[NodeIndex],
) -> Result<()> {
    let Some(debug) = debug else {
        return Ok(());
    };

    let component_label = component_label(component_id);
    let path = debug.output_dir.join(format!(
        "{:06}_{}_component_{}_iter_{:03}_{}.graphml",
        debug.gene_index,
        sanitize_file_stem(debug.gene_name),
        component_label,
        iteration,
        sanitize_file_stem(step)
    ));

    write_graphml(
        graph,
        active,
        component_nodes,
        removed_this_step,
        debug,
        component_id,
        iteration,
        step,
        &path,
    )
}

fn write_graphml(
    graph: &RecombinationGraph,
    active: &[bool],
    component_nodes: &[NodeIndex],
    removed_this_step: &[NodeIndex],
    debug: &GeneGraphDebug<'_>,
    component_id: Option<usize>,
    iteration: usize,
    step: &str,
    path: &Path,
) -> Result<()> {
    let mut file = File::create(path)
        .with_context(|| format!("failed to create GraphML file '{}'", path.display()))?;
    let degrees = component_degrees(graph, active, component_nodes);
    let component_set = node_set(component_nodes);
    let component_label = component_label(component_id);

    writeln!(file, r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
    writeln!(
        file,
        r#"<graphml xmlns="http://graphml.graphdrawing.org/xmlns">"#
    )?;
    writeln!(
        file,
        r#"  <key id="gene_name" for="graph" attr.name="gene_name" attr.type="string"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="gene_index" for="graph" attr.name="gene_index" attr.type="int"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="component" for="graph" attr.name="component" attr.type="string"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="iteration" for="graph" attr.name="iteration" attr.type="int"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="step" for="graph" attr.name="step" attr.type="string"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="name" for="node" attr.name="name" attr.type="string"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="sample_index" for="node" attr.name="sample_index" attr.type="int"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="sample_name" for="node" attr.name="sample_name" attr.type="string"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="active" for="node" attr.name="active" attr.type="boolean"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="removed_this_step" for="node" attr.name="removed_this_step" attr.type="boolean"/>"#
    )?;
    writeln!(
        file,
        r#"  <key id="active_degree" for="node" attr.name="active_degree" attr.type="int"/>"#
    )?;
    writeln!(
        file,
        r#"  <graph id="g{}_c{}_i{}" edgedefault="undirected">"#,
        debug.gene_index, component_label, iteration
    )?;
    writeln!(
        file,
        r#"    <data key="gene_name">{}</data>"#,
        xml_escape(debug.gene_name)
    )?;
    writeln!(
        file,
        r#"    <data key="gene_index">{}</data>"#,
        debug.gene_index
    )?;
    writeln!(
        file,
        r#"    <data key="component">{}</data>"#,
        xml_escape(&component_label)
    )?;
    writeln!(file, r#"    <data key="iteration">{iteration}</data>"#)?;
    writeln!(file, r#"    <data key="step">{}</data>"#, xml_escape(step))?;

    for &node in component_nodes {
        let sample_index = graph[node];
        let fallback_name = sample_index.to_string();
        let sample_name = debug
            .sample_names
            .get(sample_index)
            .map(String::as_str)
            .unwrap_or(&fallback_name);
        let is_active = is_active(active, node);
        let removed_now = removed_this_step.contains(&node);

        writeln!(file, r#"    <node id="n{}">"#, node.index())?;
        writeln!(
            file,
            r#"      <data key="name">{}</data>"#,
            xml_escape(sample_name)
        )?;
        writeln!(
            file,
            r#"      <data key="sample_index">{sample_index}</data>"#
        )?;
        writeln!(
            file,
            r#"      <data key="sample_name">{}</data>"#,
            xml_escape(sample_name)
        )?;
        writeln!(file, r#"      <data key="active">{is_active}</data>"#)?;
        writeln!(
            file,
            r#"      <data key="removed_this_step">{removed_now}</data>"#
        )?;
        writeln!(
            file,
            r#"      <data key="active_degree">{}</data>"#,
            degrees[node.index()]
        )?;
        writeln!(file, r#"    </node>"#)?;
    }

    let mut edge_index = 0usize;
    for edge in graph.edge_references() {
        if !(component_set.contains(&edge.source()) && component_set.contains(&edge.target())) {
            continue;
        }
        writeln!(
            file,
            r#"    <edge id="e{}" source="n{}" target="n{}"/>"#,
            edge_index,
            edge.source().index(),
            edge.target().index()
        )?;
        edge_index += 1;
    }

    writeln!(file, r#"  </graph>"#)?;
    writeln!(file, r#"</graphml>"#)?;

    Ok(())
}

fn component_label(component_id: Option<usize>) -> String {
    component_id.map_or_else(|| "none".to_string(), |id| format!("{id:03}"))
}

fn sanitize_file_stem(value: &str) -> String {
    let sanitized: String = value
        .chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '-' | '_' | '.') {
                character
            } else {
                '_'
            }
        })
        .collect();

    if sanitized.is_empty() {
        "unnamed".to_string()
    } else {
        sanitized
    }
}

fn xml_escape(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&apos;"),
            _ => escaped.push(character),
        }
    }
    escaped
}

fn active_components(
    graph: &RecombinationGraph,
    active: &[bool],
    nodes: impl IntoIterator<Item = NodeIndex>,
) -> Vec<Vec<NodeIndex>> {
    let candidates: HashSet<_> = nodes
        .into_iter()
        .filter(|&node| is_active(active, node))
        .collect();
    let mut unseen = candidates.clone();
    let mut starts: Vec<_> = candidates.into_iter().collect();
    sort_nodes_by_sample(graph, &mut starts);
    let mut components = Vec::new();

    for start in starts {
        if !unseen.remove(&start) {
            continue;
        }

        let mut component = vec![start];
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            let mut neighbors: Vec<_> = graph
                .neighbors(node)
                .filter(|neighbor| unseen.contains(neighbor))
                .collect();
            sort_nodes_by_sample(graph, &mut neighbors);
            for neighbor in neighbors {
                if unseen.remove(&neighbor) {
                    component.push(neighbor);
                    stack.push(neighbor);
                }
            }
        }

        sort_nodes_by_sample(graph, &mut component);
        components.push(component);
    }

    components.sort_by_key(|component| {
        component
            .iter()
            .map(|&node| graph[node])
            .min()
            .unwrap_or(usize::MAX)
    });
    components
}

// Counts each component node's active neighbors within that same component.
fn component_degrees(
    graph: &RecombinationGraph,
    active: &[bool],
    component_nodes: &[NodeIndex],
) -> Vec<usize> {
    let component_set = node_set(component_nodes);
    let mut degrees = vec![0; active.len()];

    for &node in component_nodes {
        if !is_active(active, node) {
            continue;
        }

        degrees[node.index()] = graph
            .neighbors(node)
            .filter(|&neighbor| is_active(active, neighbor) && component_set.contains(&neighbor))
            .count();
    }

    degrees
}

fn node_set(nodes: &[NodeIndex]) -> HashSet<NodeIndex> {
    nodes.iter().copied().collect()
}

fn sort_nodes_by_sample(graph: &RecombinationGraph, nodes: &mut [NodeIndex]) {
    nodes.sort_by_key(|node| graph[*node]);
}

// Removes active graph nodes and records their sample indices.
fn remove_nodes(
    graph: &RecombinationGraph,
    active: &mut [bool],
    nodes: impl IntoIterator<Item = NodeIndex>,
    recombinant_samples: &mut Vec<usize>,
) {
    for node in nodes {
        if !is_active(active, node) {
            continue;
        }

        active[node.index()] = false;
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
    use std::fs;
    use tempfile::TempDir;

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

    // Runs presence inference from sample-index pair fixtures.
    fn presence_for_pairs(sample_count: usize, pairs: &[(usize, usize)]) -> Vec<u8> {
        let pair_offsets = pairs
            .iter()
            .map(|&(left, right)| pair_offset(sample_count, left, right))
            .collect();
        infer_gene_presence(sample_count, &pair_offsets, None).unwrap()
    }

    fn pair_offset(sample_count: usize, left: usize, right: usize) -> u32 {
        if left >= sample_count || right >= sample_count || left == right {
            return (sample_count * sample_count) as u32;
        }

        let (left, right) = if left < right {
            (left, right)
        } else {
            (right, left)
        };
        (pairs_before_sample(sample_count, left) + right - left - 1) as u32
    }

    #[test]
    // Verifies isolated residual nodes stop pruning cleanly.
    fn isolated_residual_nodes_do_not_loop_forever() {
        let graph = graph_from_edges(4, &[]);

        let observed = prune_recombinant_samples(&graph, None).unwrap();

        assert!(observed.is_empty());
    }

    #[test]
    // Verifies components with fewer than three nodes do not mark recombination presence.
    fn pruning_skips_components_with_fewer_than_three_nodes() {
        let observed = presence_for_pairs(4, &[(0, 1)]);

        assert_eq!(observed, vec![0, 0, 0, 0]);
    }

    #[test]
    // Verifies a three-node component removes its unique highest-degree node.
    fn three_node_path_marks_only_center_node() {
        let observed = presence_for_pairs(3, &[(0, 1), (1, 2)]);

        assert_eq!(observed, vec![0, 1, 0]);
    }

    #[test]
    // Verifies a three-node equal-degree component marks every node.
    fn three_node_triangle_marks_all_nodes() {
        let observed = presence_for_pairs(3, &[(0, 1), (1, 2), (0, 2)]);

        assert_eq!(observed, vec![1, 1, 1]);
    }

    #[test]
    // Verifies a star removes its max-degree center and leaves tiny components unmarked.
    fn star_marks_only_center_node() {
        let observed = presence_for_pairs(5, &[(0, 1), (0, 2), (0, 3), (0, 4)]);

        assert_eq!(observed, vec![1, 0, 0, 0, 0]);
    }

    #[test]
    // Verifies regular structures remove all nodes tied for maximum degree.
    fn regular_graph_marks_all_max_degree_nodes() {
        let observed = presence_for_pairs(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);

        assert_eq!(observed, vec![1, 1, 1, 1]);
    }

    #[test]
    // Verifies pruning repeats until only singleton or pair components remain.
    fn pruning_repeats_until_only_tiny_components_remain() {
        let observed = presence_for_pairs(5, &[(0, 1), (1, 2), (2, 3), (3, 0), (3, 4)]);

        assert_eq!(observed, vec![0, 1, 0, 1, 0]);
    }

    #[test]
    // Verifies disconnected components are pruned independently.
    fn disconnected_components_do_not_share_degree_decisions() {
        let observed = presence_for_pairs(8, &[(0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (6, 7)]);

        assert_eq!(observed, vec![1, 0, 0, 0, 0, 0, 1, 0]);
    }

    #[test]
    // Verifies pruning a component can split it into independently handled subcomponents.
    fn component_splits_are_reprocessed_independently() {
        let observed = presence_for_pairs(
            11,
            &[
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 5),
                (0, 6),
                (0, 7),
                (0, 8),
                (0, 9),
                (0, 10),
                (1, 2),
                (2, 3),
                (1, 3),
                (4, 5),
                (5, 6),
            ],
        );

        assert_eq!(observed, vec![1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0]);
    }

    #[test]
    // Verifies selected genes emit GraphML snapshots for pruning decisions.
    fn graph_debug_writes_graphml_snapshots_for_selected_gene() {
        let dir = TempDir::new().unwrap();
        let output_dir = dir.path().join("graphs");
        let config = GraphDebugConfig::new(output_dir.clone(), vec!["geneA".to_string()]).unwrap();
        let sample_names = vec![
            "alpha".to_string(),
            "beta".to_string(),
            "gamma".to_string(),
            "delta".to_string(),
            "epsilon".to_string(),
        ];
        let debug = config.for_gene(7, "geneA", &sample_names).unwrap();
        let pair_offsets = [(0, 1), (0, 2), (0, 3), (0, 4)]
            .into_iter()
            .map(|(left, right)| pair_offset(5, left, right))
            .collect();

        let observed = infer_gene_presence(5, &pair_offsets, Some(&debug)).unwrap();

        assert_eq!(observed, vec![1, 0, 0, 0, 0]);
        let mut files: Vec<_> = fs::read_dir(&output_dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        files.sort();
        assert_eq!(
            files,
            vec![
                "000007_geneA_component_000_iter_000_initial.graphml",
                "000007_geneA_component_000_iter_001_remove_max_degree_d4.graphml",
                "000007_geneA_component_000_iter_002_stop_too_small.graphml",
            ]
        );

        let hub_snapshot = fs::read_to_string(
            output_dir.join("000007_geneA_component_000_iter_001_remove_max_degree_d4.graphml"),
        )
        .unwrap();
        assert!(hub_snapshot.contains(r#"<data key="gene_name">geneA</data>"#));
        assert!(hub_snapshot.contains(r#"<data key="sample_name">alpha</data>"#));
        assert!(hub_snapshot.contains(r#"<data key="removed_this_step">true</data>"#));
    }

    #[test]
    // Verifies selected genes with no graph nodes emit a diagnostic empty component snapshot.
    fn graph_debug_writes_component_none_snapshot_for_empty_graph() {
        let dir = TempDir::new().unwrap();
        let output_dir = dir.path().join("graphs");
        let config = GraphDebugConfig::new(output_dir.clone(), vec!["geneA".to_string()]).unwrap();
        let sample_names = vec!["alpha".to_string(), "beta".to_string()];
        let debug = config.for_gene(3, "geneA", &sample_names).unwrap();
        let pair_offsets = RoaringBitmap::new();

        let observed = infer_gene_presence(2, &pair_offsets, Some(&debug)).unwrap();

        assert_eq!(observed, vec![0, 0]);
        let mut files: Vec<_> = fs::read_dir(&output_dir)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        files.sort();
        assert_eq!(
            files,
            vec!["000003_geneA_component_none_iter_000_initial.graphml"]
        );

        let snapshot = fs::read_to_string(
            output_dir.join("000003_geneA_component_none_iter_000_initial.graphml"),
        )
        .unwrap();
        assert!(snapshot.contains(r#"<data key="component">none</data>"#));
    }

    #[test]
    // Verifies invalid and self-pair edges are skipped defensively.
    fn invalid_and_self_pairs_are_ignored() {
        let observed = presence_for_pairs(4, &[(0, 0), (0, 4), (4, 1), (0, 1)]);

        assert_eq!(observed, vec![0, 0, 0, 0]);
    }
}
