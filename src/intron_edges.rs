// Adapted from SplAdder v3.1.1 editgraph.insert_intron_edges (BSD-3-Clause).
// See licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    graph::{Interval, SpliceGraph},
};
use anyhow::{Result, bail, ensure};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct IntronOptions {
    pub min_exon_len: i64,
    pub vicinity_region: i64,
    pub insert_intron_retention: bool,
    pub gene_merges: bool,
    pub append_new_terminal_exons: bool,
    pub append_new_terminal_exons_len: i64,
}

#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Inserted {
    pub intron_in_exon: usize,
    pub alt_53_prime: usize,
    pub exon_skip: usize,
    pub gene_merge: usize,
    pub new_terminal_exon: usize,
}

fn copy_exon(graph: &mut SpliceGraph, old: usize, exon: Interval, terminal: [bool; 2]) -> usize {
    let neighbors = graph.edges[old].clone();
    let diagonal = graph.connected(old, old);
    let new = graph.push_vertex(exon, terminal);
    for j in neighbors {
        graph.connect(new, j);
    }
    if diagonal {
        graph.connect(new, new);
    }
    new
}

/// Add missing splice sites using proximal exons, coverage, or a new terminal.
/// left=true means the donor has no matching exon end.
fn insert_missing(
    gene: &mut Gene,
    intron: Interval,
    all_introns: &[Interval],
    existing: &[usize],
    left: bool,
    options: IntronOptions,
    coverage: &mut impl FnMut(&Gene) -> Result<Vec<u64>>,
) -> Result<Option<bool>> {
    let site = intron[usize::from(!left)];
    let mut candidates: Vec<_> = gene
        .splicegraph
        .vertices
        .iter()
        .enumerate()
        .filter_map(|(i, v)| {
            let vicinity = if left {
                v[0] - options.vicinity_region <= site && v[1] + options.vicinity_region > site
            } else {
                v[0] - options.vicinity_region < site && v[1] + options.vicinity_region >= site
            };
            vicinity.then_some(i)
        })
        .collect();
    if candidates.is_empty() {
        let candidate = if left {
            gene.splicegraph.vertices.iter().rposition(|v| v[1] <= site)
        } else {
            gene.splicegraph.vertices.iter().position(|v| v[0] > site)
        };
        if let Some(i) = candidate {
            // gg = genes[i] aliases the gene in upstream, so these temporary
            // query bounds remain visible until gen_graphs resets them.
            if left {
                gene.start = gene.splicegraph.vertices[i][1];
                gene.stop = site;
            } else {
                gene.start = site;
                gene.stop = gene.splicegraph.vertices[i][1];
            }
            let track = coverage(gene)?;
            let fraction = track.iter().filter(|&&x| x > 10).count() as f64 / track.len() as f64;
            // Empty coverage gives NaN; the source rejects only fraction < .9.
            if fraction >= 0.9 || fraction.is_nan() {
                candidates.push(i);
            }
        }
    }
    let graph = &mut gene.splicegraph;
    if let Some(old) = candidates.into_iter().min_by_key(|&i| {
        let v = graph.vertices[i];
        if left {
            (v[0] - site).abs().min((v[1] - site - 1).abs())
        } else {
            (v[0] - site + 1).abs().min((v[1] - site).abs())
        }
    }) {
        let v = graph.vertices[old];
        let length = if left { site - v[0] } else { v[1] - site };
        if length >= options.min_exon_len {
            let mut terminals = graph.terminals[old];
            terminals[usize::from(left)] = false;
            let exon = if left { [v[0], site] } else { [site, v[1]] };
            let new = copy_exon(graph, old, exon, terminals);
            if left {
                graph.add_intron(&[new], false, existing, true);
            } else {
                graph.add_intron(existing, true, &[new], false);
            }
            return Ok(Some(false));
        }
    }
    if options.append_new_terminal_exons {
        let exon = if left {
            let start = (site - options.append_new_terminal_exons_len).max(0);
            [
                all_introns
                    .iter()
                    .map(|p| p[1])
                    .filter(|&p| p >= start && p < site)
                    .max()
                    .unwrap_or(start),
                site,
            ]
        } else {
            let stop = site + options.append_new_terminal_exons_len;
            [
                site,
                all_introns
                    .iter()
                    .map(|p| p[0])
                    .filter(|&p| p > site && p <= stop)
                    .min()
                    .unwrap_or(stop),
            ]
        };
        let new = graph.push_vertex(exon, [left, !left]);
        for &i in existing {
            if left && graph.terminals[i][0] && intron[1] < graph.vertices[i][1] {
                graph.vertices[i][0] = intron[1];
            } else if !left && graph.terminals[i][1] && intron[0] > graph.vertices[i][0] {
                graph.vertices[i][1] = intron[0];
            }
        }
        // Both source branches pass the existing nodes as the first argument.
        graph.add_intron(existing, true, &[new], false);
        if left {
            graph.uniquify();
        }
        return Ok(Some(true));
    }
    Ok(None)
}

/// Introns are supplied in source order for each gene's selected strand.
pub fn insert_edges(
    genes: &mut [Gene],
    introns: &[Vec<Interval>],
    options: IntronOptions,
    mut coverage: impl FnMut(&Gene) -> Result<Vec<u64>>,
) -> Result<Inserted> {
    ensure!(
        genes.len() == introns.len(),
        "intron lists must match genes"
    );
    let mut inserted = Inserted::default();
    for i in 0..genes.len() {
        for (j, &intron) in introns[i].iter().enumerate() {
            let graph = &mut genes[i].splicegraph;
            if j > 0 && graph.vertices.len() > 1 {
                graph.uniquify();
            }
            let left: Vec<_> = graph
                .vertices
                .iter()
                .enumerate()
                .filter_map(|(k, v)| (v[1] == intron[0]).then_some(k))
                .collect();
            let right: Vec<_> = graph
                .vertices
                .iter()
                .enumerate()
                .filter_map(|(k, v)| (v[0] == intron[1]).then_some(k))
                .collect();
            if left.is_empty() && right.is_empty() {
                if options.insert_intron_retention {
                    let contained: Vec<_> = graph
                        .vertices
                        .iter()
                        .enumerate()
                        .filter_map(|(k, v)| (intron[0] > v[0] && intron[1] < v[1]).then_some(k))
                        .collect();
                    for old in contained {
                        let v = graph.vertices[old];
                        let incoming: Vec<_> = graph.edges[old]
                            .iter()
                            .copied()
                            .filter(|&k| k <= old)
                            .collect();
                        let a =
                            graph.push_vertex([v[0], intron[0]], [graph.terminals[old][0], false]);
                        for k in incoming {
                            graph.connect(a, k);
                        }
                        let outgoing: Vec<_> = graph.edges[old]
                            .iter()
                            .copied()
                            .filter(|&k| k >= old)
                            .collect();
                        let b =
                            graph.push_vertex([intron[1], v[1]], [false, graph.terminals[old][1]]);
                        for k in outgoing {
                            graph.connect(b, k);
                        }
                        graph.connect(a, b);
                        inserted.intron_in_exon += 1;
                    }
                }
                continue;
            }
            let previous_match = left.is_empty()
                && i > 0
                && genes[i - 1].chr == genes[i].chr
                && genes[i - 1].strand == genes[i].strand
                && genes[i - 1]
                    .splicegraph
                    .vertices
                    .iter()
                    .any(|v| v[1] == intron[0]);
            let next_match = right.is_empty()
                && i + 1 < genes.len()
                && genes[i + 1].chr == genes[i].chr
                && genes[i + 1].strand == genes[i].strand
                && genes[i + 1]
                    .splicegraph
                    .vertices
                    .iter()
                    .any(|v| v[0] == intron[1]);
            if previous_match || next_match {
                if options.gene_merges {
                    // Upstream's initial np.c_[zeros((0,)), [a,b]] raises a
                    // shape error; a successful gene-merging result is undefined.
                    bail!(
                        "gene_merges encountered adjacent genes; SplAdder v3.1.1 fails constructing its merge index for this input"
                    );
                }
                continue;
            }
            if left.is_empty() {
                if let Some(terminal) = insert_missing(
                    &mut genes[i],
                    intron,
                    &introns[i],
                    &right,
                    true,
                    options,
                    &mut coverage,
                )? {
                    if terminal {
                        inserted.new_terminal_exon += 1;
                    } else {
                        inserted.alt_53_prime += 1;
                    }
                    continue;
                }
            } else if right.is_empty()
                && let Some(terminal) = insert_missing(
                    &mut genes[i],
                    intron,
                    &introns[i],
                    &left,
                    false,
                    options,
                    &mut coverage,
                )?
            {
                if terminal {
                    inserted.new_terminal_exon += 1;
                } else {
                    inserted.alt_53_prime += 1;
                }
                continue;
            }
            if left.is_empty() || right.is_empty() || left.len() > 20 || right.len() > 20 {
                continue;
            }
            let graph = &mut genes[i].splicegraph;
            for &a in &left {
                for &b in &right {
                    inserted.exon_skip += usize::from(!graph.connected(a, b));
                }
            }
            graph.add_intron(&left, true, &right, true);
        }
    }
    // The source's post-loop uniquify only applies to the last gene.
    if let Some(last) = genes.last_mut()
        && last.splicegraph.vertices.len() > 1
    {
        last.splicegraph.uniquify();
    }
    for gene in genes {
        ensure!(
            gene.splicegraph.vertices.iter().all(|v| v[0] <= v[1]),
            "intron insertion produced an invalid exon"
        );
        gene.splicegraph.sort();
    }
    Ok(inserted)
}
