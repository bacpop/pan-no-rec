use anyhow::{Context, Result, bail};
use std::collections::{HashMap, VecDeque};
use std::fs;
use std::path::Path;

#[derive(Clone, Debug, PartialEq)]
pub(crate) enum GmlValue {
    Integer(i64),
    Float(f64),
    String(String),
    List(Vec<GmlAttribute>),
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct GmlAttribute {
    pub(crate) key: String,
    pub(crate) value: GmlValue,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BlockNode {
    pub(crate) id: i64,
    pub(crate) label: Option<String>,
    pub(crate) coordinate: Option<usize>,
    pub(crate) attributes: Vec<GmlAttribute>,
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct BlockEdge {
    pub(crate) source: i64,
    pub(crate) target: i64,
    pub(crate) attributes: Vec<GmlAttribute>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BlockGraph {
    pub(crate) nodes: Vec<BlockNode>,
    pub(crate) edges: Vec<BlockEdge>,
    node_indices: HashMap<i64, usize>,
}

impl BlockGraph {
    pub fn from_gml_path(path: &Path) -> Result<Self> {
        let contents = fs::read_to_string(path)
            .with_context(|| format!("failed to read GML graph '{}'", path.display()))?;
        Self::from_gml_str(&contents)
            .with_context(|| format!("failed to parse GML graph '{}'", path.display()))
    }

    pub fn from_gml_str(contents: &str) -> Result<Self> {
        let tokens = tokenize(contents)?;
        let mut parser = Parser::new(tokens);
        let document = parser.parse_document()?;
        let graph_attrs = document
            .iter()
            .find_map(|attribute| match (&attribute.key[..], &attribute.value) {
                ("graph", GmlValue::List(attributes)) => Some(attributes.as_slice()),
                _ => None,
            })
            .with_context(|| "GML document does not contain a graph list")?;

        Self::from_graph_attributes(graph_attrs)
    }

    pub fn node_index(&self, node_id: i64) -> Option<usize> {
        self.node_indices.get(&node_id).copied()
    }

    pub fn node_coordinate(&self, node_id: i64) -> Option<usize> {
        self.node_index(node_id)
            .and_then(|index| self.nodes[index].coordinate)
    }

    pub fn node_name(&self, node_id: i64) -> Option<&str> {
        self.node_index(node_id).and_then(|index| {
            self.nodes[index]
                .string_attribute("name")
                .or(self.nodes[index].label.as_deref())
        })
    }

    pub fn node_size(&self, node_id: i64) -> Result<Option<usize>> {
        self.node_index(node_id)
            .map(|index| self.nodes[index].size())
            .transpose()
    }

    pub fn assign_named_node_coordinate(&mut self, name: &str, coordinate: usize) -> Result<i64> {
        let node_index = self.unique_node_with_string_attribute("name", name)?;
        for (index, node) in self.nodes.iter_mut().enumerate() {
            if index != node_index && node.coordinate == Some(coordinate) {
                node.coordinate = None;
            }
        }
        self.nodes[node_index].coordinate = Some(coordinate);
        Ok(self.nodes[node_index].id)
    }

    pub fn assign_dnaa_coordinate_zero(&mut self) -> Result<i64> {
        self.assign_named_node_coordinate("dnaA", 0)
    }

    pub fn clear_coordinates(&mut self) {
        for node in &mut self.nodes {
            node.coordinate = None;
        }
    }

    /// Assigns a sample-weighted, simple-path layout to the graph's genes.
    ///
    /// Core nodes (at least 95% of samples) are laid out first. Remaining nodes are
    /// appended in sample-weighted simple path segments so every graph node
    /// receives a coordinate.
    pub fn layout_core_genes(&mut self, total_samples: usize) -> Result<Vec<i64>> {
        if total_samples == 0 {
            bail!("total sample count must be greater than zero");
        }

        let sizes: Vec<_> = self
            .nodes
            .iter()
            .map(|node| node.size())
            .collect::<Result<_>>()?;
        for (node, &size) in self.nodes.iter().zip(&sizes) {
            if size > total_samples {
                bail!(
                    "GML node {} has size {size}, greater than total sample count {total_samples}",
                    node.id
                );
            }
        }

        self.clear_coordinates();
        let adjacency = self.adjacency();
        let core: Vec<_> = sizes
            .iter()
            .map(|&size| is_at_least_percent(size, total_samples, CORE_PERCENT))
            .collect();
        let eligible: Vec<_> = sizes
            .iter()
            .map(|&size| is_above_percent(size, total_samples, ELIGIBLE_PERCENT))
            .collect();
        let mut rng = LayoutRng::default();
        let mut ordered = Vec::with_capacity(self.nodes.len());
        let mut used = vec![false; self.nodes.len()];

        let mut core_components = connected_components(&adjacency, &core);
        if !core_components.is_empty() {
            let largest_size = core_components
                .iter()
                .map(Vec::len)
                .max()
                .unwrap_or_default();
            let anchor_candidates: Vec<_> = core_components
                .iter()
                .enumerate()
                .filter_map(|(index, component)| (component.len() == largest_size).then_some(index))
                .collect();
            let anchor_index = rng.choose_index(&anchor_candidates);
            let anchor = core_components.swap_remove(anchor_index);
            let anchor_set = node_mask(self.nodes.len(), &anchor);
            let anchor_path = greedy_path(&adjacency, &anchor_set, &sizes, &mut rng);
            append_nodes(anchor_path, &mut ordered, &mut used);

            // Grow from the free end of the anchor path. Each connector is a
            // shortest path; ties prefer greater new sample weight.
            loop {
                let Some(&current) = ordered.last() else {
                    break;
                };
                let remaining_core: Vec<_> = core
                    .iter()
                    .enumerate()
                    .map(|(index, &is_core)| is_core && !used[index])
                    .collect();
                if !remaining_core.iter().any(|&present| present) {
                    break;
                }

                let eligible_path = best_path_to_core(
                    current,
                    &remaining_core,
                    &eligible,
                    &adjacency,
                    &sizes,
                    &eligible,
                    &used,
                    &mut rng,
                );
                let path = if eligible_path.is_some() {
                    eligible_path
                } else {
                    best_path_to_core(
                        current,
                        &remaining_core,
                        &vec![true; self.nodes.len()],
                        &adjacency,
                        &sizes,
                        &eligible,
                        &used,
                        &mut rng,
                    )
                };

                if let Some(path) = path {
                    append_nodes(path, &mut ordered, &mut used);
                } else {
                    // Disconnected core components become separate segments,
                    // largest first. Their internal path is greedily selected.
                    let mut disconnected = connected_components(&adjacency, &remaining_core);
                    if disconnected.is_empty() {
                        break;
                    }
                    disconnected.sort_by(|left, right| right.len().cmp(&left.len()));
                    let largest = disconnected[0].len();
                    let tied: Vec<_> = disconnected
                        .iter()
                        .take_while(|component| component.len() == largest)
                        .collect();
                    let chosen = rng.choose_index(&tied);
                    let component = tied[chosen];
                    let component_mask = node_mask(self.nodes.len(), component);
                    let segment = greedy_path(&adjacency, &component_mask, &sizes, &mut rng);
                    if segment.is_empty() {
                        break;
                    }
                    append_nodes(segment, &mut ordered, &mut used);
                }
            }
        }

        // Lay out every remaining node in additional segments. This includes
        // branch nodes and rare genes that could not be used as connectors.
        loop {
            let remaining_nodes: Vec<_> = used
                .iter()
                .enumerate()
                .map(|(_, &is_used)| !is_used)
                .collect();
            if !remaining_nodes.iter().any(|&present| present) {
                break;
            }

            let mut components = connected_components(&adjacency, &remaining_nodes);
            components.sort_by(|left, right| {
                let left_weight: usize = left.iter().map(|&index| sizes[index]).sum();
                let right_weight: usize = right.iter().map(|&index| sizes[index]).sum();
                right_weight
                    .cmp(&left_weight)
                    .then_with(|| right.len().cmp(&left.len()))
            });
            let best_weight: usize = components[0].iter().map(|&index| sizes[index]).sum();
            let best_count = components[0].len();
            let tied: Vec<_> = components
                .iter()
                .take_while(|component| {
                    component.iter().map(|&index| sizes[index]).sum::<usize>() == best_weight
                        && component.len() == best_count
                })
                .collect();
            let component = tied[rng.choose_index(&tied)];
            let component_mask = node_mask(self.nodes.len(), component);
            let segment = greedy_path(&adjacency, &component_mask, &sizes, &mut rng);
            if segment.is_empty() {
                break;
            }
            append_nodes(segment, &mut ordered, &mut used);
        }

        for (coordinate, &node_index) in ordered.iter().enumerate() {
            self.nodes[node_index].coordinate = Some(coordinate);
        }

        Ok(ordered
            .into_iter()
            .map(|node_index| self.nodes[node_index].id)
            .collect())
    }

    fn adjacency(&self) -> Vec<Vec<usize>> {
        let mut adjacency = vec![Vec::new(); self.nodes.len()];
        for edge in &self.edges {
            let (Some(&source), Some(&target)) = (
                self.node_indices.get(&edge.source),
                self.node_indices.get(&edge.target),
            ) else {
                continue;
            };
            if source != target {
                adjacency[source].push(target);
                adjacency[target].push(source);
            }
        }
        for neighbors in &mut adjacency {
            neighbors.sort_unstable();
            neighbors.dedup();
        }
        adjacency
    }

    fn from_graph_attributes(attributes: &[GmlAttribute]) -> Result<Self> {
        let mut nodes = Vec::new();
        let mut edges = Vec::new();
        let mut node_indices = HashMap::new();

        for attribute in attributes {
            match (&attribute.key[..], &attribute.value) {
                ("node", GmlValue::List(node_attributes)) => {
                    let node = parse_node(node_attributes)?;
                    let node_index = nodes.len();
                    if node_indices.insert(node.id, node_index).is_some() {
                        bail!("GML graph contains duplicate node id {}", node.id);
                    }
                    nodes.push(node);
                }
                ("edge", GmlValue::List(edge_attributes)) => {
                    edges.push(parse_edge(edge_attributes)?);
                }
                _ => {}
            }
        }

        for edge in &edges {
            if !node_indices.contains_key(&edge.source) {
                bail!("GML edge references missing source node {}", edge.source);
            }
            if !node_indices.contains_key(&edge.target) {
                bail!("GML edge references missing target node {}", edge.target);
            }
        }

        Ok(BlockGraph {
            nodes,
            edges,
            node_indices,
        })
    }

    fn unique_node_with_string_attribute(&self, key: &str, value: &str) -> Result<usize> {
        let mut matching_indices = self.nodes.iter().enumerate().filter_map(|(index, node)| {
            (node.string_attribute(key) == Some(value)).then_some(index)
        });
        let Some(node_index) = matching_indices.next() else {
            bail!("GML graph does not contain a node with {key} {value}");
        };
        if matching_indices.next().is_some() {
            bail!("GML graph contains multiple nodes with {key} {value}");
        }

        Ok(node_index)
    }
}

impl BlockNode {
    pub(crate) fn string_attribute(&self, key: &str) -> Option<&str> {
        optional_string(&self.attributes, key)
    }

    fn size(&self) -> Result<usize> {
        self.attributes
            .iter()
            .find_map(|attribute| {
                if attribute.key == "size" {
                    match &attribute.value {
                        GmlValue::Integer(value) if *value >= 0 => Some(
                            usize::try_from(*value).context("node size does not fit this platform"),
                        ),
                        GmlValue::Integer(value) => Some(Err(anyhow::anyhow!(
                            "GML node {} has negative size {value}",
                            self.id
                        ))),
                        _ => Some(Err(anyhow::anyhow!(
                            "GML node {} has a non-integer size attribute",
                            self.id
                        ))),
                    }
                } else {
                    None
                }
            })
            .with_context(|| format!("GML node {} is missing its size attribute", self.id))?
    }
}

const CORE_PERCENT: u128 = 95;
const ELIGIBLE_PERCENT: u128 = 5;

fn is_above_percent(size: usize, total: usize, percent: u128) -> bool {
    (size as u128) * 100 > (total as u128) * percent
}

fn is_at_least_percent(size: usize, total: usize, percent: u128) -> bool {
    (size as u128) * 100 >= (total as u128) * percent
}

fn node_mask(node_count: usize, nodes: &[usize]) -> Vec<bool> {
    let mut mask = vec![false; node_count];
    for &node in nodes {
        mask[node] = true;
    }
    mask
}

fn connected_components(adjacency: &[Vec<usize>], included: &[bool]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; adjacency.len()];
    let mut components = Vec::new();

    for start in 0..adjacency.len() {
        if !included[start] || seen[start] {
            continue;
        }
        let mut component = Vec::new();
        let mut queue = VecDeque::from([start]);
        seen[start] = true;
        while let Some(node) = queue.pop_front() {
            component.push(node);
            for &neighbor in &adjacency[node] {
                if included[neighbor] && !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push_back(neighbor);
                }
            }
        }
        components.push(component);
    }

    components
}

fn append_nodes(nodes: Vec<usize>, ordered: &mut Vec<usize>, used: &mut [bool]) {
    for node in nodes {
        if !used[node] {
            used[node] = true;
            ordered.push(node);
        }
    }
}

// Builds a greedy simple path, trying a bounded number of randomized starts.
fn greedy_path(
    adjacency: &[Vec<usize>],
    included: &[bool],
    sizes: &[usize],
    rng: &mut LayoutRng,
) -> Vec<usize> {
    let mut starts: Vec<_> = included
        .iter()
        .enumerate()
        .filter_map(|(index, &include)| include.then_some(index))
        .collect();
    rng.shuffle(&mut starts);
    starts.truncate(16);

    let mut best = Vec::new();
    let mut best_weight = 0usize;
    for start in starts {
        let mut path = vec![start];
        let mut in_path = vec![false; adjacency.len()];
        in_path[start] = true;

        loop {
            let mut candidates: Vec<_> = adjacency[*path.last().unwrap()]
                .iter()
                .copied()
                .filter(|&neighbor| included[neighbor] && !in_path[neighbor])
                .collect();
            if candidates.is_empty() {
                break;
            }
            rng.shuffle(&mut candidates);
            candidates.sort_by(|&left, &right| {
                let left_options = adjacency[left]
                    .iter()
                    .filter(|&&neighbor| included[neighbor] && !in_path[neighbor])
                    .count();
                let right_options = adjacency[right]
                    .iter()
                    .filter(|&&neighbor| included[neighbor] && !in_path[neighbor])
                    .count();
                left_options
                    .cmp(&right_options)
                    .then_with(|| sizes[right].cmp(&sizes[left]))
            });
            let next = candidates[0];
            path.push(next);
            in_path[next] = true;
        }

        let weight: usize = path.iter().map(|&node| sizes[node]).sum();
        if weight > best_weight || (weight == best_weight && path.len() > best.len()) {
            best_weight = weight;
            best = path;
        }
    }

    best
}

// Finds shortest paths to unplaced core nodes, maximizing eligible sample
// weight among equal-length routes and then choosing the best target route.
fn best_path_to_core(
    start: usize,
    target: &[bool],
    allowed: &[bool],
    adjacency: &[Vec<usize>],
    sizes: &[usize],
    eligible: &[bool],
    used: &[bool],
    rng: &mut LayoutRng,
) -> Option<Vec<usize>> {
    let mut previous = vec![None; adjacency.len()];
    let mut distance = vec![usize::MAX; adjacency.len()];
    let mut path_weight = vec![0usize; adjacency.len()];
    let mut queue = VecDeque::from([start]);
    distance[start] = 0;

    while let Some(node) = queue.pop_front() {
        let mut neighbors = adjacency[node].clone();
        rng.shuffle(&mut neighbors);
        for neighbor in neighbors {
            if !allowed[neighbor] || used[neighbor] {
                continue;
            }
            let candidate_distance = distance[node] + 1;
            let candidate_weight = path_weight[node]
                + if eligible[neighbor] {
                    sizes[neighbor]
                } else {
                    0
                };
            if distance[neighbor] == usize::MAX {
                distance[neighbor] = candidate_distance;
                previous[neighbor] = Some(node);
                path_weight[neighbor] = candidate_weight;
                queue.push_back(neighbor);
            } else if distance[neighbor] == candidate_distance
                && (candidate_weight > path_weight[neighbor]
                    || (candidate_weight == path_weight[neighbor]
                        && rng.choose_index(&[false, true]) == 0))
            {
                previous[neighbor] = Some(node);
                path_weight[neighbor] = candidate_weight;
            }
        }
    }

    let mut candidates: Vec<_> = target
        .iter()
        .enumerate()
        .filter_map(|(index, &is_target)| {
            (is_target && distance[index] != usize::MAX).then_some(index)
        })
        .collect();
    if candidates.is_empty() {
        return None;
    }

    candidates.sort_by(|&left, &right| {
        path_weight[right]
            .cmp(&path_weight[left])
            .then_with(|| distance[left].cmp(&distance[right]))
    });
    let best_weight = path_weight[candidates[0]];
    let best_distance = distance[candidates[0]];
    let tied: Vec<_> = candidates
        .into_iter()
        .take_while(|&node| path_weight[node] == best_weight && distance[node] == best_distance)
        .collect();
    let target_node = tied[rng.choose_index(&tied)];

    let mut path = Vec::new();
    let mut current = target_node;
    while current != start {
        path.push(current);
        current = previous[current]?;
    }
    path.reverse();
    Some(path)
}

#[derive(Debug)]
struct LayoutRng(u64);

impl Default for LayoutRng {
    fn default() -> Self {
        // Fixed seed keeps layouts reproducible while distributing ties.
        Self(0x9e37_79b9_7f4a_7c15)
    }
}

impl LayoutRng {
    fn next(&mut self) -> u64 {
        let mut value = self.0;
        value ^= value << 13;
        value ^= value >> 7;
        value ^= value << 17;
        self.0 = value;
        value
    }

    fn choose_index<T>(&mut self, values: &[T]) -> usize {
        (self.next() as usize) % values.len()
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = self.choose_index(&values[..=index]);
            values.swap(index, other);
        }
    }
}

pub fn load_gml(path: &Path) -> Result<BlockGraph> {
    BlockGraph::from_gml_path(path)
}

fn parse_node(attributes: &[GmlAttribute]) -> Result<BlockNode> {
    let id = required_integer(attributes, "id")?;
    let label = optional_string(attributes, "label")
        .or_else(|| optional_string(attributes, "name"))
        .map(str::to_owned);

    Ok(BlockNode {
        id,
        label,
        coordinate: None,
        attributes: attributes.to_vec(),
    })
}

fn parse_edge(attributes: &[GmlAttribute]) -> Result<BlockEdge> {
    Ok(BlockEdge {
        source: required_integer(attributes, "source")?,
        target: required_integer(attributes, "target")?,
        attributes: attributes.to_vec(),
    })
}

fn required_integer(attributes: &[GmlAttribute], key: &str) -> Result<i64> {
    attributes
        .iter()
        .find_map(|attribute| {
            if attribute.key == key {
                attribute.value.as_integer()
            } else {
                None
            }
        })
        .with_context(|| format!("GML list is missing required integer key '{key}'"))
}

fn optional_string<'a>(attributes: &'a [GmlAttribute], key: &str) -> Option<&'a str> {
    attributes.iter().find_map(|attribute| {
        if attribute.key == key {
            attribute.value.as_str()
        } else {
            None
        }
    })
}

impl GmlValue {
    fn as_integer(&self) -> Option<i64> {
        match self {
            GmlValue::Integer(value) => Some(*value),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            GmlValue::String(value) => Some(value),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum Token {
    Atom(String),
    StringLiteral(String),
    OpenBracket,
    CloseBracket,
}

fn tokenize(contents: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let mut chars = contents.char_indices().peekable();

    while let Some((start, character)) = chars.next() {
        match character {
            '[' => tokens.push(Token::OpenBracket),
            ']' => tokens.push(Token::CloseBracket),
            '"' => tokens.push(Token::StringLiteral(read_string_literal(&mut chars)?)),
            '#' => skip_comment(&mut chars),
            character if character.is_whitespace() => {}
            _ => {
                let mut end = start + character.len_utf8();
                while let Some(&(next_index, next_character)) = chars.peek() {
                    if next_character.is_whitespace()
                        || matches!(next_character, '[' | ']' | '"' | '#')
                    {
                        break;
                    }
                    chars.next();
                    end = next_index + next_character.len_utf8();
                }
                tokens.push(Token::Atom(contents[start..end].to_string()));
            }
        }
    }

    Ok(tokens)
}

fn read_string_literal(
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
) -> Result<String> {
    let mut value = String::new();

    while let Some((_, character)) = chars.next() {
        match character {
            '"' => return Ok(value),
            '\\' => {
                let Some((_, escaped)) = chars.next() else {
                    bail!("unterminated escape sequence in GML string literal");
                };
                match escaped {
                    '"' => value.push('"'),
                    '\\' => value.push('\\'),
                    'n' => value.push('\n'),
                    'r' => value.push('\r'),
                    't' => value.push('\t'),
                    other => value.push(other),
                }
            }
            other => value.push(other),
        }
    }

    bail!("unterminated GML string literal")
}

fn skip_comment(chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>) {
    while let Some((_, character)) = chars.next() {
        if character == '\n' {
            break;
        }
    }
}

struct Parser {
    tokens: Vec<Token>,
    position: usize,
}

impl Parser {
    fn new(tokens: Vec<Token>) -> Self {
        Parser {
            tokens,
            position: 0,
        }
    }

    fn parse_document(&mut self) -> Result<Vec<GmlAttribute>> {
        let mut attributes = Vec::new();
        while self.position < self.tokens.len() {
            attributes.push(self.parse_attribute()?);
        }
        Ok(attributes)
    }

    fn parse_attribute(&mut self) -> Result<GmlAttribute> {
        let key = self.expect_atom("GML key")?;
        let value = self
            .parse_value()
            .with_context(|| format!("failed to parse value for GML key '{key}'"))?;
        Ok(GmlAttribute { key, value })
    }

    fn parse_value(&mut self) -> Result<GmlValue> {
        match self.next_token() {
            Some(Token::OpenBracket) => self.parse_list().map(GmlValue::List),
            Some(Token::StringLiteral(value)) => Ok(GmlValue::String(value)),
            Some(Token::Atom(value)) => Ok(parse_atom_value(value)),
            Some(Token::CloseBracket) => bail!("unexpected ']' while parsing GML value"),
            None => bail!("unexpected end of GML document while parsing value"),
        }
    }

    fn parse_list(&mut self) -> Result<Vec<GmlAttribute>> {
        let mut attributes = Vec::new();

        loop {
            match self.peek_token() {
                Some(Token::CloseBracket) => {
                    self.position += 1;
                    return Ok(attributes);
                }
                Some(_) => attributes.push(self.parse_attribute()?),
                None => bail!("unterminated GML list"),
            }
        }
    }

    fn expect_atom(&mut self, context: &str) -> Result<String> {
        match self.next_token() {
            Some(Token::Atom(value)) => Ok(value),
            Some(Token::StringLiteral(_)) => bail!("expected {context}, found string literal"),
            Some(Token::OpenBracket) => bail!("expected {context}, found '['"),
            Some(Token::CloseBracket) => bail!("expected {context}, found ']'"),
            None => bail!("expected {context}, found end of GML document"),
        }
    }

    fn peek_token(&self) -> Option<&Token> {
        self.tokens.get(self.position)
    }

    fn next_token(&mut self) -> Option<Token> {
        let token = self.tokens.get(self.position).cloned();
        if token.is_some() {
            self.position += 1;
        }
        token
    }
}

fn parse_atom_value(value: String) -> GmlValue {
    if let Ok(integer) = value.parse::<i64>() {
        return GmlValue::Integer(integer);
    }

    if let Ok(float) = value.parse::<f64>() {
        return GmlValue::Float(float);
    }

    GmlValue::String(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_networkx_style_gml_with_unassigned_coordinates() {
        let graph = BlockGraph::from_gml_str(
            r#"
            graph [
              directed 0
              node [
                id 10
                label "gene_a"
              ]
              node [
                id 4
                label "gene_b"
              ]
              edge [
                source 10
                target 4
              ]
            ]
            "#,
        )
        .unwrap();

        assert_eq!(graph.nodes.len(), 2);
        assert_eq!(graph.edges.len(), 1);
        assert_eq!(graph.nodes[0].label.as_deref(), Some("gene_a"));
        assert_eq!(graph.node_coordinate(10), None);
        assert_eq!(graph.node_coordinate(4), None);
    }

    #[test]
    fn assigns_dnaa_node_coordinate_zero_from_name_attribute() {
        let mut graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [
                id 1
                name "abc"
              ]
              node [
                id 2
                name "dnaA"
              ]
            ]
            "#,
        )
        .unwrap();

        let node_id = graph.assign_dnaa_coordinate_zero().unwrap();

        assert_eq!(node_id, 2);
        assert_eq!(graph.node_coordinate(1), None);
        assert_eq!(graph.node_coordinate(2), Some(0));
    }

    #[test]
    fn named_coordinate_assignment_rejects_ambiguous_names() {
        let mut graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [ id 1 name "dnaA" ]
              node [ id 2 name "dnaA" ]
            ]
            "#,
        )
        .unwrap();

        let error = graph.assign_dnaa_coordinate_zero().unwrap_err();

        assert!(
            error.to_string().contains("multiple nodes"),
            "error: {error}"
        );
    }

    #[test]
    fn preserves_nested_attributes() {
        let graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [
                id 0
                label "gene"
                metadata [
                  genome_count 3
                ]
              ]
            ]
            "#,
        )
        .unwrap();

        assert!(matches!(
            &graph.nodes[0].attributes[2].value,
            GmlValue::List(_)
        ));
    }

    #[test]
    fn rejects_edges_with_missing_nodes() {
        let error = BlockGraph::from_gml_str(
            r#"
            graph [
              node [ id 0 ]
              edge [ source 0 target 1 ]
            ]
            "#,
        )
        .unwrap_err();

        assert!(
            error.to_string().contains("missing target node"),
            "error: {error}"
        );
    }

    fn layout_fixture() -> &'static str {
        r#"
        graph [
          node [ id 1 size 20 name "core_a" ]
          node [ id 2 size 1 name "rare_bridge" ]
          node [ id 3 size 20 name "core_b" ]
          node [ id 4 size 2 name "common_gene" ]
          edge [ source 1 target 2 ]
          edge [ source 2 target 3 ]
          edge [ source 3 target 4 ]
        ]
        "#
    }

    #[test]
    fn layout_connects_core_with_low_frequency_bridge_and_assigns_coordinates() {
        let mut graph = BlockGraph::from_gml_str(layout_fixture()).unwrap();

        let layout = graph.layout_core_genes(20).unwrap();

        assert_eq!(layout.len(), 4);
        assert_eq!(
            layout
                .iter()
                .copied()
                .collect::<std::collections::HashSet<_>>()
                .len(),
            4
        );
        let positions: HashMap<_, _> = graph
            .nodes
            .iter()
            .map(|node| (node.id, node.coordinate.unwrap()))
            .collect();
        assert_eq!(positions.len(), 4);
        assert_eq!(positions[&layout[0]], 0);
        assert_eq!(positions[&layout[3]], 3);
        assert!(layout.contains(&2));
    }

    #[test]
    fn layout_prefers_greater_sample_weight_between_equal_length_routes() {
        let mut graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [ id 1 size 20 ]
              node [ id 2 size 2 ]
              node [ id 3 size 15 ]
              node [ id 4 size 20 ]
              edge [ source 1 target 2 ]
              edge [ source 2 target 4 ]
              edge [ source 1 target 3 ]
              edge [ source 3 target 4 ]
            ]
            "#,
        )
        .unwrap();

        let layout = graph.layout_core_genes(20).unwrap();

        assert_eq!(layout[1], 3);
    }

    #[test]
    fn layout_is_reproducible_and_uses_strict_frequency_thresholds() {
        let mut first = BlockGraph::from_gml_str(layout_fixture()).unwrap();
        let mut second = BlockGraph::from_gml_str(layout_fixture()).unwrap();

        assert_eq!(
            first.layout_core_genes(20).unwrap(),
            second.layout_core_genes(20).unwrap()
        );
        assert!(is_at_least_percent(19, 20, CORE_PERCENT));
        assert!(is_above_percent(20, 20, CORE_PERCENT));
        assert!(!is_above_percent(1, 20, ELIGIBLE_PERCENT));
        assert!(is_above_percent(2, 20, ELIGIBLE_PERCENT));
    }

    #[test]
    fn disconnected_core_components_are_kept_as_separate_segments() {
        let mut graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [ id 1 size 20 ]
              node [ id 2 size 20 ]
              node [ id 3 size 20 ]
              node [ id 4 size 20 ]
              node [ id 5 size 20 ]
              edge [ source 1 target 2 ]
              edge [ source 3 target 4 ]
            ]
            "#,
        )
        .unwrap();

        let layout = graph.layout_core_genes(20).unwrap();

        assert_eq!(layout.len(), 5);
        assert_eq!(
            layout
                .iter()
                .copied()
                .collect::<std::collections::HashSet<_>>()
                .len(),
            5
        );
        let component_pairs = [
            std::collections::HashSet::from([1, 2]),
            std::collections::HashSet::from([3, 4]),
        ];
        assert!(
            component_pairs
                .iter()
                .any(|pair| { pair.contains(&layout[0]) && pair.contains(&layout[1]) })
        );
        assert!(
            component_pairs
                .iter()
                .any(|pair| { pair.contains(&layout[2]) && pair.contains(&layout[3]) })
        );
        assert_eq!(layout[4], 5);
    }

    #[test]
    fn branched_core_graph_layout_never_repeats_nodes() {
        let mut graph = BlockGraph::from_gml_str(
            r#"
            graph [
              node [ id 1 size 20 ]
              node [ id 2 size 20 ]
              node [ id 3 size 20 ]
              edge [ source 1 target 2 ]
              edge [ source 1 target 3 ]
            ]
            "#,
        )
        .unwrap();

        let layout = graph.layout_core_genes(20).unwrap();

        assert_eq!(layout.len(), 3);
        assert_eq!(
            layout
                .iter()
                .copied()
                .collect::<std::collections::HashSet<_>>()
                .len(),
            3
        );
        assert_eq!(
            graph
                .nodes
                .iter()
                .filter(|node| node.coordinate.is_some())
                .count(),
            3
        );
    }

    #[test]
    fn layout_rejects_invalid_sample_counts_and_missing_sizes() {
        let mut graph = BlockGraph::from_gml_str("graph [ node [ id 1 size 1 ] ]").unwrap();
        assert!(graph.layout_core_genes(0).is_err());

        let mut missing_size = BlockGraph::from_gml_str("graph [ node [ id 1 ] ]").unwrap();
        assert!(missing_size.layout_core_genes(10).is_err());
    }
}
