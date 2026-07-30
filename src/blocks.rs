use anyhow::{bail, Context, Result};
use std::collections::HashMap;
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
pub(crate) struct BlockGraph {
    pub(crate) nodes: Vec<BlockNode>,
    pub(crate) edges: Vec<BlockEdge>,
    node_indices: HashMap<i64, usize>,
}

impl BlockGraph {
    pub(crate) fn from_gml_path(path: &Path) -> Result<Self> {
        let contents = fs::read_to_string(path)
            .with_context(|| format!("failed to read GML graph '{}'", path.display()))?;
        Self::from_gml_str(&contents)
            .with_context(|| format!("failed to parse GML graph '{}'", path.display()))
    }

    pub(crate) fn from_gml_str(contents: &str) -> Result<Self> {
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

    pub(crate) fn node_index(&self, node_id: i64) -> Option<usize> {
        self.node_indices.get(&node_id).copied()
    }

    pub(crate) fn node_coordinate(&self, node_id: i64) -> Option<usize> {
        self.node_index(node_id)
            .and_then(|index| self.nodes[index].coordinate)
    }

    pub(crate) fn assign_named_node_coordinate(
        &mut self,
        name: &str,
        coordinate: usize,
    ) -> Result<i64> {
        let node_index = self.unique_node_with_string_attribute("name", name)?;
        for (index, node) in self.nodes.iter_mut().enumerate() {
            if index != node_index && node.coordinate == Some(coordinate) {
                node.coordinate = None;
            }
        }
        self.nodes[node_index].coordinate = Some(coordinate);
        Ok(self.nodes[node_index].id)
    }

    pub(crate) fn assign_dnaa_coordinate_zero(&mut self) -> Result<i64> {
        self.assign_named_node_coordinate("dnaA", 0)
    }

    pub(crate) fn clear_coordinates(&mut self) {
        for node in &mut self.nodes {
            node.coordinate = None;
        }
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
}

pub(crate) fn load_gml(path: &Path) -> Result<BlockGraph> {
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
}
