//! Exact-order segment metadata, independent of recombination inference.
use crate::blocks::{BlockGraph, LayoutResult};
use anyhow::{Context, Result, bail};
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Section {
    Backbone,
    Accessory,
    Unmapped,
}
impl Section {
    fn text(self) -> &'static str {
        match self {
            Self::Backbone => "backbone",
            Self::Accessory => "accessory",
            Self::Unmapped => "unmapped",
        }
    }
    fn parse(s: &str) -> Result<Self> {
        match s {
            "backbone" => Ok(Self::Backbone),
            "accessory" => Ok(Self::Accessory),
            "unmapped" => Ok(Self::Unmapped),
            _ => bail!("invalid segment section '{s}'"),
        }
    }
}
#[derive(Debug)]
struct Entry {
    gene: String,
    segment: String,
    section: Section,
}
pub(crate) struct Metadata {
    entries: Vec<Entry>,
}
impl Metadata {
    pub(crate) fn sequences(&self) -> Vec<crate::candidate::Sequence> {
        let mut sequences: Vec<crate::candidate::Sequence> = Vec::new();
        for (i, entry) in self.entries.iter().enumerate() {
            if sequences.last().is_none_or(|seq| seq.id != entry.segment) {
                sequences.push(crate::candidate::Sequence {
                    id: entry.segment.clone(),
                    start: i,
                    end: i + 1,
                });
            } else {
                sequences.last_mut().unwrap().end = i + 1;
            }
        }
        sequences
    }
    pub(crate) fn backbone_breaks(&self) -> Vec<usize> {
        self.entries
            .windows(2)
            .enumerate()
            .filter_map(|(i, p)| {
                (p[0].section == Section::Backbone
                    && p[1].section == Section::Backbone
                    && p[0].segment != p[1].segment)
                    .then_some(i + 1)
            })
            .collect()
    }
}
fn build(graph: &BlockGraph, layout: &LayoutResult, names: &[String]) -> Result<Metadata> {
    let mut ids: HashMap<&str, Vec<i64>> = HashMap::new();
    for node in &graph.nodes {
        if let Some(name) = graph.node_name(node.id) {
            ids.entry(name).or_default().push(node.id);
        }
    }
    let backbone: HashSet<_> = layout.node_order[..layout.backbone_len]
        .iter()
        .copied()
        .collect();
    let mut backbone_segment = HashMap::new();
    for (segment, &start) in layout.backbone_segment_starts.iter().enumerate() {
        let end = layout
            .backbone_segment_starts
            .get(segment + 1)
            .copied()
            .unwrap_or(layout.backbone_len);
        for &id in &layout.node_order[start..end] {
            backbone_segment.insert(id, segment);
        }
    }
    let edges: HashSet<_> = graph
        .edges
        .iter()
        .flat_map(|e| [(e.source, e.target), (e.target, e.source)])
        .collect();
    let mut seen = HashSet::new();
    let mut entries = Vec::new();
    let mut previous = None;
    let mut segment = 0;
    for name in names {
        if !seen.insert(name) {
            bail!("duplicate output feature '{name}' in layout segment export");
        }
        let id = match ids.get(name.as_str()) {
            None => None,
            Some(ids) if ids.len() == 1 => Some(ids[0]),
            Some(_) => bail!("feature '{name}' identifies multiple graph nodes"),
        };
        let connected = previous.zip(id).is_some_and(|(a, b)| {
            let same_backbone_segment = match (backbone_segment.get(&a), backbone_segment.get(&b)) {
                (Some(a), Some(b)) => a == b,
                _ => true,
            };
            a != b && same_backbone_segment && edges.contains(&(a, b))
        });
        if !connected {
            segment += 1;
        }
        let section = match id {
            None => Section::Unmapped,
            Some(id) if backbone.contains(&id) => Section::Backbone,
            Some(_) => Section::Accessory,
        };
        entries.push(Entry {
            gene: name.clone(),
            segment: format!("segment_{segment}"),
            section,
        });
        previous = id;
    }
    Ok(Metadata { entries })
}
pub(crate) fn write(
    path: &Path,
    graph: &BlockGraph,
    layout: &LayoutResult,
    names: &[String],
) -> Result<()> {
    let metadata = build(graph, layout, names)?;
    let mut writer =
        BufWriter::new(File::create(path).with_context(|| {
            format!("failed to create layout segment file '{}'", path.display())
        })?);
    writeln!(writer, "position\tgene\tsegment_id\tsection")?;
    for (i, entry) in metadata.entries.iter().enumerate() {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            i + 1,
            entry.gene,
            entry.segment,
            entry.section.text()
        )?;
    }
    writer.flush()?;
    Ok(())
}
pub(crate) fn read(path: &Path, names: &[String]) -> Result<Metadata> {
    let data = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read segment file '{}'", path.display()))?;
    parse(&data, names).with_context(|| format!("invalid segment file '{}'", path.display()))
}
fn parse(data: &str, names: &[String]) -> Result<Metadata> {
    let mut lines = data.lines();
    if lines.next().map(|s| s.trim_end_matches('\r')) != Some("position\tgene\tsegment_id\tsection")
    {
        bail!("expected segment header: position, gene, segment_id, section");
    }
    let mut entries = Vec::new();
    let mut completed = HashSet::new();
    let mut last = String::new();
    let mut seen_genes = HashSet::new();
    let mut accessory_started = false;
    for (i, line) in lines.enumerate() {
        let fields: Vec<_> = line.trim_end_matches('\r').split('\t').collect();
        if fields.len() != 4 {
            bail!("segment line {} must have four fields", i + 2);
        }
        if fields[0].parse::<usize>().ok() != Some(i + 1) {
            bail!(
                "segment positions must be consecutive and one-based at line {}",
                i + 2
            );
        }
        if names.get(i).map(String::as_str) != Some(fields[1]) {
            bail!(
                "segment feature '{}' at position {} does not match input TSV",
                fields[1],
                i + 1
            );
        }
        if !seen_genes.insert(fields[1]) || fields[2].is_empty() {
            bail!("duplicate feature or empty segment ID at line {}", i + 2);
        }
        if fields[2] != last {
            if completed.contains(fields[2]) {
                bail!(
                    "segment ID '{}' is repeated in non-contiguous runs",
                    fields[2]
                );
            }
            if !last.is_empty() {
                completed.insert(last.clone());
            }
            last = fields[2].to_owned();
        }
        let section = Section::parse(fields[3])?;
        if section == Section::Accessory {
            accessory_started = true;
        }
        if section == Section::Backbone && accessory_started {
            bail!("backbone feature appears after appended accessory section");
        }
        entries.push(Entry {
            gene: fields[1].to_owned(),
            segment: fields[2].to_owned(),
            section,
        });
    }
    if entries.len() != names.len() {
        bail!(
            "segment file has {} features but TSV has {}",
            entries.len(),
            names.len()
        );
    }
    Ok(Metadata { entries })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (BlockGraph, LayoutResult) {
        let graph=BlockGraph::from_gml_str(r#"graph [node [id 1 name "a" size 20] node [id 2 name "b" size 20] node [id 3 name "c" size 20] node [id 4 name "d" size 1] node [id 5 name "skip" size 20] edge [source 1 target 5] edge [source 5 target 2] edge [source 2 target 3] edge [source 3 target 4]]"#).unwrap();
        (
            graph,
            LayoutResult {
                node_order: vec![1, 5, 2, 3, 4],
                backbone_len: 4,
                backbone_segment_starts: vec![0],
            },
        )
    }
    #[test]
    fn filtered_order_sections_and_unmapped_singletons_roundtrip() {
        let (graph, layout) = fixture();
        let names = vec![
            "a".into(),
            "b".into(),
            "c".into(),
            "d".into(),
            "unknown1".into(),
            "unknown2".into(),
        ];
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("segments.tsv");
        write(&path, &graph, &layout, &names).unwrap();
        let data = read(&path, &names).unwrap();
        assert_eq!(data.backbone_breaks(), vec![1]);
        let sequences = data.sequences();
        assert_eq!(
            sequences
                .iter()
                .map(|s| (s.start, s.end))
                .collect::<Vec<_>>(),
            vec![(0, 1), (1, 4), (4, 5), (5, 6)]
        );
        assert_eq!(sequences[1].id, "segment_2"); // connected section transition remains one sequence

        assert_eq!(data.entries[0].section, Section::Backbone);
        assert_eq!(data.entries[3].section, Section::Accessory);
        assert_eq!(data.entries[4].section, Section::Unmapped);
        assert_ne!(data.entries[4].segment, data.entries[5].segment);
        assert_eq!(data.entries[2].segment, data.entries[3].segment);
        let all = vec![
            "a".into(),
            "skip".into(),
            "b".into(),
            "c".into(),
            "d".into(),
        ];
        assert!(
            build(&graph, &layout, &all)
                .unwrap()
                .backbone_breaks()
                .is_empty()
        );
    }
    #[test]
    fn ambiguous_graph_mapping_and_duplicate_export_features_fail() {
        let graph = BlockGraph::from_gml_str(
            r#"graph [node [id 1 name "a" size 1] node [id 2 name "a" size 1]]"#,
        )
        .unwrap();
        let layout = LayoutResult {
            node_order: vec![1, 2],
            backbone_len: 0,
            backbone_segment_starts: vec![],
        };
        assert!(build(&graph, &layout, &["a".into()]).is_err());
        let (graph, layout) = fixture();
        assert!(build(&graph, &layout, &["a".into(), "a".into()]).is_err());
    }
    #[test]
    fn explicit_backbone_gaps_survive_direct_edges_and_filtering() {
        let (graph, mut layout) = fixture();
        // skip--b and b--c are real edges, but these are separate layout blocks.
        layout.backbone_segment_starts = vec![0, 2, 3];
        let all: Vec<_> = ["a", "skip", "b", "c", "d"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(
            build(&graph, &layout, &all).unwrap().backbone_breaks(),
            vec![2, 3]
        );
        let filtered: Vec<_> = ["skip", "c"].iter().map(|s| s.to_string()).collect();
        assert_eq!(
            build(&graph, &layout, &filtered).unwrap().backbone_breaks(),
            vec![1]
        );
        let filtered: Vec<_> = ["b", "c"].iter().map(|s| s.to_string()).collect();
        assert_eq!(
            build(&graph, &layout, &filtered).unwrap().backbone_breaks(),
            vec![1]
        );
    }
    #[test]
    fn invalid_or_stale_files_fail() {
        let names = vec!["a".into(), "b".into(), "c".into()];
        let header = "position\tgene\tsegment_id\tsection\n";
        for rows in [
            "0\ta\ts1\tbackbone\n",
            "1\tb\ts1\tbackbone\n",
            "1\ta\ts1\tbackbone\n",
            "1\ta\ts1\tbackbone\n2\tb\ts2\tbackbone\n3\tc\ts1\tbackbone\n",
            "1\ta\ts1\twrong\n",
            "1\ta\ts1\taccessory\n2\tb\ts2\tbackbone\n3\tc\ts3\taccessory\n",
            "1\ta\t\tbackbone\n",
        ] {
            assert!(
                parse(&format!("{header}{rows}"), &names).is_err(),
                "accepted {rows}"
            );
        }
    }
    #[test]
    fn no_core_and_accessory_breaks_do_not_mark_backbone() {
        let (graph, mut layout) = fixture();
        layout.backbone_len = 0;
        let names = vec!["a".into(), "b".into(), "d".into()];
        let data = build(&graph, &layout, &names).unwrap();
        assert!(data.backbone_breaks().is_empty());
        assert!(data.entries.iter().all(|e| e.section == Section::Accessory));
    }
}
