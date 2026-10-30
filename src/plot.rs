use anyhow::{Context, Result, bail};
use clap::Parser;
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(about = "Plot a pan-no-rec TSV beside a Newick phylogeny")]
struct Args {
    /// pan-no-rec recombination table
    #[arg(long, value_name = "TSV")]
    tsv: PathBuf,

    /// Phylogeny in Newick format
    #[arg(long, value_name = "NEWICK")]
    tree: PathBuf,

    /// Output filename prefix; writes PREFIX.svg and PREFIX.png
    #[arg(long, value_name = "PREFIX", default_value = "recombinant_heatmap")]
    output: PathBuf,
}

struct Matrix {
    genes: Vec<String>,
    isolates: Vec<String>,
    values: Vec<Vec<bool>>,
    frequencies: Vec<Option<(usize, usize)>>,
}

#[derive(Debug)]
struct TreeNode {
    label: Option<String>,
    branch_length: Option<f64>,
    children: Vec<TreeNode>,
}

struct PositionedNode {
    x: f64,
    y: f64,
    children: Vec<PositionedNode>,
}

pub fn run() -> Result<()> {
    let args = Args::parse();
    let matrix = read_matrix(&args.tsv)?;
    let newick = fs::read_to_string(&args.tree)
        .with_context(|| format!("failed to read Newick tree '{}'", args.tree.display()))?;
    let tree = NewickParser::parse(&newick)
        .with_context(|| format!("failed to parse Newick tree '{}'", args.tree.display()))?;
    let svg = render_heatmap(&matrix, &tree)?;
    let png = rasterize_svg(&svg)?;
    let svg_path = svg_path_for(&args.output);
    let png_path = png_path_for(&args.output);
    fs::write(&svg_path, &svg)
        .with_context(|| format!("failed to write SVG '{}'", svg_path.display()))?;
    fs::write(&png_path, png)
        .with_context(|| format!("failed to write PNG '{}'", png_path.display()))?;
    eprintln!("Wrote recombinant heatmap to '{}'", svg_path.display());
    eprintln!("Wrote bitmap heatmap to '{}'", png_path.display());
    Ok(())
}

fn path_with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

fn svg_path_for(prefix: &Path) -> PathBuf {
    path_with_suffix(prefix, ".svg")
}

fn png_path_for(prefix: &Path) -> PathBuf {
    path_with_suffix(prefix, ".png")
}

fn rasterize_svg(svg: &str) -> Result<Vec<u8>> {
    let mut options = resvg::usvg::Options::default();
    options.fontdb_mut().load_system_fonts();
    let tree = resvg::usvg::Tree::from_str(svg, &options)
        .context("failed to parse generated heatmap SVG for PNG output")?;
    let size = tree.size().to_int_size();
    let mut pixmap = resvg::tiny_skia::Pixmap::new(size.width(), size.height())
        .context("failed to allocate PNG canvas for heatmap")?;
    resvg::render(
        &tree,
        resvg::tiny_skia::Transform::default(),
        &mut pixmap.as_mut(),
    );
    pixmap.encode_png().context("failed to encode heatmap PNG")
}

fn read_matrix(path: &PathBuf) -> Result<Matrix> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read pan-no-rec recombination TSV '{}'", path.display()))?;
    let mut lines = contents.lines().filter(|line| !line.is_empty());
    let header = lines
        .next()
        .with_context(|| format!("recombination TSV '{}' is empty", path.display()))?
        .trim_end_matches('\r');
    let header_fields: Vec<_> = header.split('\t').collect();
    if header_fields.first() != Some(&"gene") || header_fields.len() < 2 {
        bail!(
            "Malformed pan-no-rec TSV!'{}' must start with 'gene' and contain sample columns",
            path.display()
        );
    }

    let isolates: Vec<_> = header_fields[1..]
        .iter()
        .map(|name| (*name).to_owned())
        .collect();
    let mut isolate_indices = HashMap::with_capacity(isolates.len());
    for (index, isolate) in isolates.iter().enumerate() {
        if isolate.is_empty() || isolate_indices.insert(isolate.as_str(), index).is_some() {
            bail!(
                "Malformed pan-no-rec TSV! '{}' contains an empty or duplicate sample name",
                path.display()
            );
        }
    }

    let mut genes = Vec::new();
    let mut values = Vec::new();
    let mut seen_genes = HashMap::new();
    let mut frequency_values: Option<Vec<Option<(usize, usize)>>> = None;
    for (line_index, line) in lines.enumerate() {
        let line = line.trim_end_matches('\r');
        let fields: Vec<_> = line.split('\t').collect();
        if fields.first() == Some(&"#gene_frequency") {
            if frequency_values.is_some() {
                bail!(
                    "Malformed pan-no-rec TSV! '{}' contains multiple gene frequency rows",
                    path.display()
                );
            }
            let parsed = fields[1..]
                .iter()
                .map(|value| {
                    if *value == "NA" {
                        return Ok(None);
                    }
                    let (present, total) = value.split_once('/').with_context(|| {
                        format!(
                            "invalid gene frequency '{}' on TSV line {}",
                            value,
                            line_index + 2
                        )
                    })?;
                    let present = present.parse::<usize>().with_context(|| {
                        format!(
                            "invalid gene frequency '{}' on TSV line {}",
                            value,
                            line_index + 2
                        )
                    })?;
                    let total = total.parse::<usize>().with_context(|| {
                        format!(
                            "invalid gene frequency '{}' on TSV line {}",
                            value,
                            line_index + 2
                        )
                    })?;
                    if total == 0 || present > total {
                        bail!(
                            "invalid gene frequency '{}' on TSV line {}",
                            value,
                            line_index + 2
                        );
                    }
                    Ok(Some((present, total)))
                })
                .collect::<Result<Vec<_>>>()?;
            frequency_values = Some(parsed);
            continue;
        }
        if fields.len() != isolates.len() + 1 {
            bail!(
                "Malformed pan-no-rec TSV! '{}' line {} has {} columns; expected {}",
                path.display(),
                line_index + 2,
                fields.len(),
                isolates.len() + 1
            );
        }
        let gene = fields[0];
        if gene.is_empty() || seen_genes.insert(gene, ()).is_some() {
            bail!(
                "Malformed pan-no-rec TSV!'{}' line {} has an empty or duplicate gene name",
                path.display(),
                line_index + 2
            );
        }
        let row = fields[1..]
            .iter()
            .map(|value| match *value {
                "0" => Ok(false),
                "1" => Ok(true),
                _ => Err(anyhow::anyhow!(
                    "Malformed pan-no-rec TSV! '{}' line {} has non-binary value '{}'",
                    path.display(),
                    line_index + 2,
                    value
                )),
            })
            .collect::<Result<Vec<_>>>()?;
        genes.push(gene.to_owned());
        values.push(row);
    }

    let frequencies = match frequency_values {
        Some(frequencies) if frequencies.len() == genes.len() => frequencies,
        Some(frequencies) => bail!(
            "Malformed pan-no-rec TSV! '{}' gene frequency row has {} values; expected {}",
            path.display(),
            frequencies.len(),
            genes.len()
        ),
        None => vec![None; genes.len()],
    };

    Ok(Matrix {
        genes,
        isolates,
        values,
        frequencies,
    })
}

fn render_heatmap(matrix: &Matrix, tree: &TreeNode) -> Result<String> {
    let mut tip_names = Vec::new();
    collect_tip_names(tree, &mut tip_names)?;
    let mut seen_tips = HashMap::with_capacity(tip_names.len());
    for tip in &tip_names {
        if seen_tips.insert(tip.as_str(), ()).is_some() {
            bail!("Phylogeny contains duplicate tip name '{tip}'");
        }
    }

    let isolate_indices: HashMap<_, _> = matrix
        .isolates
        .iter()
        .enumerate()
        .map(|(index, isolate)| (isolate.as_str(), index))
        .collect();
    for tip in &tip_names {
        if !isolate_indices.contains_key(tip.as_str()) {
            bail!("Phylogeny tip label '{tip}' has no matching sample column in the TSV");
        }
    }

    let row_height = 18.0;
    let top = 18.0;
    let tree_left = 12.0;
    let tree_width = (heatmap_width(matrix.genes.len()) * 0.35).clamp(250.0, 1400.0);
    let label_width = tip_names
        .iter()
        .map(|tip| tip.chars().count())
        .max()
        .unwrap_or(0) as f64
        * 7.5
        + 14.0;
    let heatmap_x = tree_left + tree_width + label_width;
    let cell_width = cell_width(matrix.genes.len());
    let heatmap_width = matrix.genes.len() as f64 * cell_width;
    let height = top + tip_names.len() as f64 * row_height + 8.0;
    let width = heatmap_x + heatmap_width + 24.0;

    let mut max_depth = tree_max_depth(tree, 0.0);
    if max_depth == 0.0 {
        max_depth = tree_max_edges(tree, 0);
    }
    if max_depth == 0.0 {
        max_depth = 1.0;
    }
    let mut tip_index = 0usize;
    let positioned = position_tree(
        tree,
        0.0,
        max_depth,
        tree_left,
        tree_width,
        top,
        row_height,
        &mut tip_index,
    );

    let mut svg = String::new();
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{height:.0}" viewBox="0 0 {width:.0} {height:.0}">"#
    ));
    svg.push_str(r#"<rect width="100%" height="100%" fill="white"/>"#);
    draw_tree(&positioned, &mut svg);

    for (tip_row, tip) in tip_names.iter().enumerate() {
        let y = top + (tip_row as f64 + 0.5) * row_height;
        let x = tree_left + tree_width + 5.0;
        svg.push_str(&format!(
            r#"<text x="{x:.1}" y="{:.1}" text-anchor="start" font-family="sans-serif" font-size="11">{}</text>"#,
            y + 4.0,
            xml_escape(tip)
        ));
    }

    svg.push_str(&format!(
        r##"<rect x="{heatmap_x:.1}" y="{top:.1}" width="{heatmap_width:.1}" height="{:.1}" fill="#ffffff" stroke="#bbbbbb" stroke-width="0.5"/>"##,
        tip_names.len() as f64 * row_height
    ));
    for (gene_index, frequency) in matrix.frequencies.iter().enumerate() {
        if let Some((present, total)) = frequency {
            if (*present as u128) * 100 < (*total as u128) * 95 {
                let x = heatmap_x + gene_index as f64 * cell_width;
                svg.push_str(&format!(
                    r##"<rect x="{x:.1}" y="{top:.1}" width="{cell_width:.1}" height="{:.1}" fill="#eeeeee"/>"##,
                    tip_names.len() as f64 * row_height
                ));
            }
        }
    }
    let mut clade_ranges = Vec::new();
    let mut next_tip = 0;
    collect_clade_ranges(tree, &mut next_tip, &mut clade_ranges);
    for (gene_index, gene_values) in matrix.values.iter().enumerate() {
        let positive_tips: Vec<_> = tip_names
            .iter()
            .map(|tip| gene_values[isolate_indices[tip.as_str()]])
            .collect();
        let red_tips = monophyletic_positive_tips(&positive_tips, &clade_ranges);
        let x = heatmap_x + gene_index as f64 * cell_width;
        for (tip_row, &is_positive) in positive_tips.iter().enumerate() {
            if is_positive {
                let y = top + tip_row as f64 * row_height;
                let (color, opacity) = if red_tips[tip_row] {
                    ("#d73027", "0.65")
                } else {
                    ("#2166ac", "0.65")
                };
                svg.push_str(&format!(
                    r##"<rect x="{x:.1}" y="{y:.1}" width="{cell_width:.1}" height="{row_height:.1}" fill="{color}" fill-opacity="{opacity}"/>"##
                ));
            }
        }
    }

    svg.push_str("</svg>");
    Ok(svg)
}

fn collect_tip_names(node: &TreeNode, names: &mut Vec<String>) -> Result<()> {
    if node.children.is_empty() {
        let label = node
            .label
            .as_ref()
            .context("Phylogeny contains a tip without a name!")?;
        names.push(label.clone());
    } else {
        for child in &node.children {
            collect_tip_names(child, names)?;
        }
    }
    Ok(())
}

// Records each internal subtree as a half-open range in tree-tip order.
fn collect_clade_ranges(node: &TreeNode, next_tip: &mut usize, ranges: &mut Vec<(usize, usize)>) {
    let start = *next_tip;
    if node.children.is_empty() {
        *next_tip += 1;
    } else {
        for child in &node.children {
            collect_clade_ranges(child, next_tip, ranges);
        }
        let end = *next_tip;
        if end - start >= 2 {
            ranges.push((start, end));
        }
    }
}

// Marks positive tips that belong to a clade whose every tip is positive.
fn monophyletic_positive_tips(positive: &[bool], clades: &[(usize, usize)]) -> Vec<bool> {
    let mut positive_prefix = vec![0usize; positive.len() + 1];
    for (index, &is_positive) in positive.iter().enumerate() {
        positive_prefix[index + 1] = positive_prefix[index] + usize::from(is_positive);
    }

    let mut coverage = vec![0isize; positive.len() + 1];
    for &(start, end) in clades {
        if positive_prefix[end] - positive_prefix[start] == end - start {
            coverage[start] += 1;
            coverage[end] -= 1;
        }
    }

    let mut red_tips = vec![false; positive.len()];
    let mut covered = 0isize;
    for index in 0..positive.len() {
        covered += coverage[index];
        red_tips[index] = positive[index] && covered > 0;
    }
    red_tips
}

fn cell_width(gene_count: usize) -> f64 {
    if gene_count == 0 {
        4.0
    } else {
        (900.0 / gene_count as f64).clamp(2.0, 14.0)
    }
}

fn heatmap_width(gene_count: usize) -> f64 {
    gene_count as f64 * cell_width(gene_count)
}

fn tree_max_depth(node: &TreeNode, depth: f64) -> f64 {
    if node.children.is_empty() {
        depth
    } else {
        node.children
            .iter()
            .map(|child| tree_max_depth(child, depth + child.branch_length.unwrap_or(1.0)))
            .fold(depth, f64::max)
    }
}

fn tree_max_edges(node: &TreeNode, depth: usize) -> f64 {
    if node.children.is_empty() {
        depth as f64
    } else {
        node.children
            .iter()
            .map(|child| tree_max_edges(child, depth + 1))
            .fold(depth as f64, f64::max)
    }
}

#[allow(clippy::too_many_arguments)]
fn position_tree(
    node: &TreeNode,
    depth: f64,
    max_depth: f64,
    left: f64,
    width: f64,
    top: f64,
    row_height: f64,
    next_tip: &mut usize,
) -> PositionedNode {
    let x = left + width * depth / max_depth;
    if node.children.is_empty() {
        let y = top + (*next_tip as f64 + 0.5) * row_height;
        *next_tip += 1;
        PositionedNode {
            x,
            y,
            children: Vec::new(),
        }
    } else {
        let children: Vec<_> = node
            .children
            .iter()
            .map(|child| {
                position_tree(
                    child,
                    depth + child.branch_length.unwrap_or(1.0),
                    max_depth,
                    left,
                    width,
                    top,
                    row_height,
                    next_tip,
                )
            })
            .collect();
        let y = children.iter().map(|child| child.y).sum::<f64>() / children.len() as f64;
        PositionedNode { x, y, children }
    }
}

fn draw_tree(node: &PositionedNode, svg: &mut String) {
    if node.children.is_empty() {
        return;
    }
    let first_y = node.children.first().unwrap().y;
    let last_y = node.children.last().unwrap().y;
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="{first_y:.1}" x2="{:.1}" y2="{last_y:.1}" stroke="#444444" stroke-width="1"/>"##,
        node.x, node.x
    ));
    for child in &node.children {
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#444444" stroke-width="1"/>"##,
            node.x, child.y, child.x, child.y
        ));
        draw_tree(child, svg);
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

struct NewickParser {
    chars: Vec<char>,
    position: usize,
}

impl NewickParser {
    fn parse(input: &str) -> Result<TreeNode> {
        let mut parser = Self {
            chars: input.chars().collect(),
            position: 0,
        };
        parser.skip_ignored();
        let tree = parser.parse_subtree()?;
        parser.skip_ignored();
        if parser.consume(';') {
            parser.skip_ignored();
        }
        if parser.peek().is_some() {
            bail!("unexpected content after Newick tree");
        }
        Ok(tree)
    }

    fn parse_subtree(&mut self) -> Result<TreeNode> {
        self.skip_ignored();
        let mut children = Vec::new();
        if self.consume('(') {
            loop {
                children.push(self.parse_subtree()?);
                self.skip_ignored();
                if self.consume(',') {
                    continue;
                }
                if self.consume(')') {
                    break;
                }
                bail!("expected ',' or ')' in Newick tree");
            }
        }

        self.skip_ignored();
        let label = self.parse_label()?;
        self.skip_ignored();
        let branch_length = if self.consume(':') {
            self.skip_ignored();
            let value = self.parse_number()?;
            if !value.is_finite() || value < 0.0 {
                bail!("Newick branch lengths must be finite and non-negative");
            }
            Some(value)
        } else {
            None
        };

        Ok(TreeNode {
            label,
            branch_length,
            children,
        })
    }

    fn parse_label(&mut self) -> Result<Option<String>> {
        self.skip_ignored();
        if self.consume('\'') {
            let mut label = String::new();
            loop {
                match self.next() {
                    Some('\'') if self.peek() == Some('\'') => {
                        self.next();
                        label.push('\'');
                    }
                    Some('\'') => return Ok(Some(label)),
                    Some(character) => label.push(character),
                    None => bail!("unterminated quoted Newick label"),
                }
            }
        }

        let mut label = String::new();
        while let Some(character) = self.peek() {
            if character.is_whitespace()
                || matches!(character, '(' | ')' | ',' | ':' | ';' | '[' | ']')
            {
                break;
            }
            label.push(character);
            self.next();
        }
        Ok((!label.is_empty()).then_some(label))
    }

    fn parse_number(&mut self) -> Result<f64> {
        let mut number = String::new();
        while let Some(character) = self.peek() {
            if character.is_whitespace() || matches!(character, ',' | ')' | ';' | '[' | ']') {
                break;
            }
            number.push(character);
            self.next();
        }
        number
            .parse::<f64>()
            .with_context(|| format!("invalid Newick branch length '{number}'"))
    }

    fn skip_ignored(&mut self) {
        loop {
            while self.peek().is_some_and(char::is_whitespace) {
                self.next();
            }
            if !self.consume('[') {
                return;
            }
            let mut depth = 1usize;
            while let Some(character) = self.next() {
                match character {
                    '[' => depth += 1,
                    ']' => {
                        depth -= 1;
                        if depth == 0 {
                            break;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    fn consume(&mut self, expected: char) -> bool {
        if self.peek() == Some(expected) {
            self.next();
            true
        } else {
            false
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.position).copied()
    }

    fn next(&mut self) -> Option<char> {
        let character = self.peek()?;
        self.position += 1;
        Some(character)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn branch_lengths_below_one_use_the_full_tree_panel() {
        let matrix = Matrix {
            genes: vec!["gene_a".to_string()],
            isolates: vec!["alpha".to_string(), "beta".to_string()],
            values: vec![vec![false, true]],
            frequencies: vec![None],
        };
        let tree = NewickParser::parse("(alpha:0.0001,beta:0.0002);").unwrap();

        let svg = render_heatmap(&matrix, &tree).unwrap();

        assert!(
            svg.contains("x2=\"137.0\""),
            "SVG did not scale alpha branch: {svg}"
        );
        assert!(
            svg.contains("x2=\"262.0\""),
            "SVG did not scale beta branch: {svg}"
        );
        assert!(svg.contains("width=\"352\" height=\"62\""));
        assert!(!svg.contains("Recombinant genes by isolate"));
    }

    #[test]
    fn tree_panel_width_grows_with_heatmap_width() {
        let gene_count = 1_000;
        let matrix = Matrix {
            genes: (0..gene_count)
                .map(|index| format!("gene_{index}"))
                .collect(),
            isolates: vec!["alpha".to_string(), "beta".to_string()],
            values: vec![vec![false; 2]; gene_count],
            frequencies: vec![None; gene_count],
        };
        let tree = NewickParser::parse("(alpha:0.1,beta:0.2);").unwrap();

        let svg = render_heatmap(&matrix, &tree).unwrap();

        assert!(
            svg.contains("x2=\"712.0\""),
            "tree panel did not grow with heatmap width"
        );
    }

    #[test]
    fn png_is_rasterized_at_svg_dimensions_and_uses_sibling_name() {
        let matrix = Matrix {
            genes: vec!["gene_a".to_string()],
            isolates: vec!["alpha".to_string(), "beta".to_string()],
            values: vec![vec![false, true]],
            frequencies: vec![None],
        };
        let tree = NewickParser::parse("(alpha:0.1,beta:0.2);").unwrap();
        let svg = render_heatmap(&matrix, &tree).unwrap();

        let png = rasterize_svg(&svg).unwrap();

        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        assert_eq!(u32::from_be_bytes(png[16..20].try_into().unwrap()), 352);
        assert_eq!(u32::from_be_bytes(png[20..24].try_into().unwrap()), 62);
        assert_eq!(
            png_path_for(Path::new("plots/heatmap")),
            PathBuf::from("plots/heatmap.png")
        );
        assert_eq!(
            svg_path_for(Path::new("plots/heatmap")),
            PathBuf::from("plots/heatmap.svg")
        );
    }
}
