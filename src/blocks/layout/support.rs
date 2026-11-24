//! Shared normalised GML edge membership.
use super::*;
use roaring::RoaringBitmap;
use std::collections::BTreeMap;

pub(super) type EdgeMembers = BTreeMap<(usize, usize), Option<RoaringBitmap>>;
pub(super) fn key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

// Invalid/missing records make the whole normalised edge unavailable. Never
// infer membership from edge sizes, or silently ignore a malformed duplicate.
pub(super) fn edge_members(graph: &BlockGraph) -> EdgeMembers {
    let mut result: EdgeMembers = BTreeMap::new();
    for edge in &graph.edges {
        let (Some(a), Some(b)) = (graph.node_index(edge.source), graph.node_index(edge.target))
        else {
            continue;
        };
        if a == b {
            continue;
        }
        let mut members = RoaringBitmap::new();
        let mut valid = true;
        for attribute in edge.attributes.iter().filter(|a| a.key == "members") {
            match &attribute.value {
                GmlValue::Integer(value) => {
                    if let Ok(value) = u32::try_from(*value) {
                        members.insert(value);
                    } else {
                        valid = false;
                    }
                }
                // NetworkX emits this marker before a one-element GML list.
                GmlValue::String(value) if value == "_networkx_list_start" => (),
                _ => valid = false,
            }
        }
        let parsed = (valid && !members.is_empty()).then_some(members);
        result
            .entry(key(a, b))
            .and_modify(|old| match (old.as_mut(), parsed.as_ref()) {
                (Some(old), Some(new)) => *old |= new,
                _ => *old = None,
            })
            .or_insert(parsed);
    }
    result
}
