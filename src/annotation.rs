// Adapted from SplAdder v3.1.1 init.py and classes/gene.py (BSD-3-Clause).
// See licenses/SplAdder-BSD.txt for the original copyright and license.
use crate::graph::{Interval, SegmentGraph, SpliceGraph};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Gene {
    pub name: String,
    pub start: i64,
    pub stop: i64,
    pub chr: String,
    pub strand: char,
    pub source: Option<String>,
    pub gene_type: Option<String>,
    pub symbol: Option<String>,
    pub transcripts: Vec<String>,
    pub exons: Vec<Vec<Interval>>,
    pub introns_anno: BTreeSet<Interval>,
    pub splicegraph: SpliceGraph,
    pub segmentgraph: SegmentGraph,
    pub is_alt: Option<bool>,
    pub is_alt_spliced: Option<bool>,
}

impl Gene {
    /// SplAdder's alternative-gene labels, using interval overlap instead of a
    /// gene-length bitmap. Preserve the source's positional-index semantics.
    pub fn label_alt(&mut self) {
        let graph = &self.splicegraph;
        let n = graph.vertices.len();
        if n < 2 {
            self.is_alt = Some(false);
            self.is_alt_spliced = Some(false);
            return;
        }
        let mut excluded = vec![false; n];
        for (i, exclude) in excluded.iter_mut().enumerate() {
            if graph.predecessors(i).next().is_none() {
                // np.where(~np.isin(np.where(equal)[0], i))[0] returns
                // POSITIONS in the equal-coordinate subset, not vertex IDs.
                let others: Vec<_> = graph
                    .vertices
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v[1] == graph.vertices[i][1])
                    .enumerate()
                    .filter(|(_, (j, _))| *j != i)
                    .map(|(position, _)| position)
                    .collect();
                if !others.is_empty() {
                    let simple = others.iter().all(|&j| {
                        graph.successors(j).eq(graph.successors(i))
                            && graph
                                .predecessors(j)
                                .all(|k| graph.vertices[k][1] <= graph.vertices[i][0])
                    });
                    if !simple {
                        continue;
                    }
                    *exclude = true;
                }
            }
            if graph.successors(i).next().is_none() {
                let others: Vec<_> = graph
                    .vertices
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| v[0] == graph.vertices[i][0])
                    .enumerate()
                    .filter(|(_, (j, _))| *j != i)
                    .map(|(position, _)| position)
                    .collect();
                if !others.is_empty() {
                    let simple = others.iter().all(|&j| {
                        graph.predecessors(j).eq(graph.predecessors(i))
                            // The original adds j (not j + 1) to slice indices.
                            && graph.successors(j).all(|k| graph.vertices[k - 1][0] >= graph.vertices[i][1])
                    });
                    if !simple {
                        continue;
                    }
                    *exclude = true;
                }
            }
        }
        let mut intervals = Vec::new();
        for i in 0..n {
            if excluded[i] {
                continue;
            }
            intervals.push(graph.vertices[i]);
            for j in graph.successors(i) {
                if !excluded[j] {
                    intervals.push([graph.vertices[i][1], graph.vertices[j][0]]);
                }
            }
        }
        intervals.retain(|v| v[0] < v[1]);
        intervals.sort_unstable();
        let mut stop = i64::MIN;
        let overlapping = intervals.iter().any(|v| {
            let overlaps = v[0] < stop;
            stop = stop.max(v[1]);
            overlaps
        });
        self.is_alt_spliced = Some(overlapping);
        self.is_alt = Some(overlapping || excluded.iter().any(|&x| x));
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct AnnotationFilters {
    pub overlap_genes: bool,
    pub overlap_exons: bool,
    pub overlap_transcripts: bool,
}

#[derive(Debug, Default)]
pub struct Annotation {
    pub genes: Vec<Gene>,
    /// Suffix and gene IDs, in the order used by upstream's exclusion reports.
    pub excluded: Vec<(String, Vec<String>)>,
}

struct Record {
    chr: String,
    source: String,
    feature: String,
    start: i64,
    stop: i64,
    strand: char,
    tags: HashMap<String, String>,
}

impl Record {
    fn required(&self, key: &str) -> Result<&str> {
        self.tags
            .get(key)
            .map(String::as_str)
            .with_context(|| format!("{} record lacks {key}", self.feature))
    }

    fn gene(&self, name: &str, gff: bool, inferred: bool) -> Result<Gene> {
        ensure!(
            matches!(self.strand, '+' | '-'),
            "gene {name} has unsupported strand {}",
            self.strand
        );
        Ok(Gene {
            name: name.to_owned(),
            start: self.start,
            stop: self.stop,
            chr: self.chr.clone(),
            strand: self.strand,
            source: Some(self.source.clone()),
            gene_type: if gff {
                Some(self.feature.clone())
            } else {
                self.tags
                    .get("gene_type")
                    .or_else(|| self.tags.get("gene_biotype"))
                    .cloned()
            },
            symbol: if inferred {
                None
            } else {
                self.tags.get("gene_name").cloned()
            },
            transcripts: Vec::new(),
            exons: Vec::new(),
            introns_anno: BTreeSet::new(),
            splicegraph: SpliceGraph::default(),
            segmentgraph: SegmentGraph::default(),
            is_alt: None,
            is_alt_spliced: None,
        })
    }
}

fn visit(path: &Path, gff: bool, mut callback: impl FnMut(Record) -> Result<()>) -> Result<()> {
    let file = File::open(path).with_context(|| format!("open annotation {}", path.display()))?;
    for (line_no, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        if gff && line.to_ascii_lowercase().starts_with("##fasta") {
            break;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<_> = line.trim_end().split('\t').collect();
        ensure!(
            fields.len() == 9,
            "{}:{}: expected 9 annotation columns",
            path.display(),
            line_no + 1
        );
        let mut tags = HashMap::new();
        for tag in fields[8].trim_matches(';').split(';') {
            let (key, value) = if gff {
                let (key, value) = tag.split_once('=').context("GFF3 attribute lacks '='")?;
                (key.trim(), value.split('=').next().unwrap().trim())
            } else {
                // Match get_tags_gtf: first space-delimited value, not an
                // invented interpretation of quoted multi-word attributes.
                let mut parts = tag.trim_matches(' ').split(' ');
                let key = parts.next().unwrap();
                let value = parts
                    .next()
                    .context("GTF attribute lacks a value")?
                    .trim_matches('"');
                (key, value)
            };
            tags.insert(key.to_owned(), value.to_owned());
        }
        let start = fields[3].parse::<i64>().map(|x| x - 1).unwrap_or(-1);
        let stop = fields[4].parse::<i64>().unwrap_or(-1);
        callback(Record {
            chr: fields[0].to_owned(),
            source: fields[1].to_owned(),
            feature: fields[2].to_owned(),
            start,
            stop,
            strand: fields[6].chars().next().unwrap_or('.'),
            tags,
        })
        .with_context(|| format!("{}:{}", path.display(), line_no + 1))?;
    }
    Ok(())
}

/// Two passes preserve Python dictionary insertion order, including genes
/// inferred from exon-only GTFs, without keeping the annotation text in memory.
pub fn read_annotation(path: &Path, filters: AnnotationFilters) -> Result<Annotation> {
    let extension = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let gff = match extension.as_str() {
        "gff" | "gff3" => true,
        "gtf" => false,
        _ => bail!(
            "annotation must end in .gtf, .gff or .gff3: {}",
            path.display()
        ),
    };
    let mut genes: Vec<Gene> = Vec::new();
    let mut gene_index: HashMap<String, usize> = HashMap::new();
    let mut transcript_parent: HashMap<String, String> = HashMap::new();
    visit(path, gff, |record| {
        let feature = record.feature.to_ascii_lowercase();
        let id = if gff {
            if matches!(feature.as_str(), "chromosome" | "contig" | "supercontig") {
                return Ok(());
            }
            match (record.tags.get("ID"), record.tags.get("Parent")) {
                (Some(id), Some(parent)) => {
                    transcript_parent.insert(id.clone(), parent.clone());
                    return Ok(());
                }
                (Some(id), None) => id.as_str(),
                _ => return Ok(()),
            }
        } else if matches!(
            feature.as_str(),
            "gene" | "pseudogene" | "unprocessed_pseudogene" | "transposable_element_gene"
        ) {
            record.required("gene_id")?
        } else {
            return Ok(());
        };
        let gene = record.gene(id, gff, false)?;
        if let Some(&i) = gene_index.get(id) {
            genes[i] = gene;
        } else {
            gene_index.insert(id.to_owned(), genes.len());
            genes.push(gene);
        }
        Ok(())
    })?;
    let mut inferred = false;
    let mut tx_indices: HashMap<(usize, String), usize> = HashMap::new();
    visit(path, gff, |record| {
        if (gff && !matches!(record.feature.as_str(), "exon" | "pseudogenic_exon"))
            || (!gff && !record.feature.eq_ignore_ascii_case("exon"))
        {
            return Ok(());
        }
        let (transcript, gene_id) = if gff {
            let tx = record.required("Parent")?;
            (
                tx,
                transcript_parent
                    .get(tx)
                    .context("exon parent has no transcript-to-gene mapping")?
                    .as_str(),
            )
        } else {
            (
                record.required("transcript_id")?,
                record.required("gene_id")?,
            )
        };
        let gene_i = if let Some(&i) = gene_index.get(gene_id) {
            i
        } else {
            ensure!(!gff, "missing gene parent {gene_id}");
            let i = genes.len();
            genes.push(record.gene(gene_id, false, true)?);
            gene_index.insert(gene_id.to_owned(), i);
            inferred = true;
            i
        };
        let gene = &mut genes[gene_i];
        let tx_i = *tx_indices
            .entry((gene_i, transcript.to_owned()))
            .or_insert_with(|| {
                let i = gene.transcripts.len();
                gene.transcripts.push(transcript.to_owned());
                gene.exons.push(Vec::new());
                i
            });
        gene.exons[tx_i].push([record.start, record.stop]);
        Ok(())
    })?;
    genes.par_iter_mut().for_each(|gene| {
        if inferred && !gene.exons.is_empty() {
            // Upstream applies this to ALL genes if any GTF gene was inferred.
            gene.start = *gene.exons.iter().flatten().flatten().min().unwrap();
            gene.stop = *gene.exons.iter().flatten().flatten().max().unwrap();
        }
        gene.splicegraph = SpliceGraph::from_transcripts(&gene.exons);
    });
    let mut result = Annotation {
        genes,
        excluded: Vec::new(),
    };
    let removed: BTreeSet<_> = result
        .genes
        .iter()
        .filter(|g| g.exons.is_empty())
        .map(|g| g.name.clone())
        .collect();
    result.exclude("no_exons", removed, false);
    if filters.overlap_genes {
        let mut order: Vec<_> = (0..result.genes.len()).collect();
        order.sort_by_key(|&i| {
            let g = &result.genes[i];
            (&g.chr, g.strand, g.start)
        });
        let mut previous: Option<usize> = None;
        let mut removed = BTreeSet::new();
        for i in order {
            let gene = &result.genes[i];
            if let Some(j) = previous {
                let prior = &result.genes[j];
                if (gene.chr.as_str(), gene.strand) == (prior.chr.as_str(), prior.strand)
                    && gene.start <= prior.stop
                {
                    removed.insert(gene.name.clone());
                    removed.insert(prior.name.clone());
                    if gene.stop > prior.stop {
                        previous = Some(i);
                    }
                    continue;
                }
            }
            previous = Some(i);
        }
        result.exclude("gene_overlap", removed, true);
    }
    if filters.overlap_exons {
        let mut owners = HashMap::new();
        let mut removed = BTreeSet::new();
        for gene in &result.genes {
            for exon in gene.exons.iter().flatten() {
                let first = owners
                    .entry((gene.chr.as_str(), *exon))
                    .or_insert(gene.name.as_str());
                if *first != gene.name {
                    removed.insert(first.to_string());
                    removed.insert(gene.name.clone());
                }
            }
        }
        result.exclude("exon_shared", removed, true);
    }
    if filters.overlap_transcripts {
        let removed = result
            .genes
            .iter()
            .filter(|gene| {
                gene.exons.iter().any(|exons| {
                    exons
                        .iter()
                        .enumerate()
                        .any(|(i, e)| exons[i + 1..].iter().any(|next| next[0] < e[1]))
                })
            })
            .map(|g| g.name.clone())
            .collect();
        result.exclude("exon_overlap", removed, true);
    }
    ensure!(
        !result.genes.is_empty(),
        "there are no valid genes left in the input annotation"
    );
    for gene in &mut result.genes {
        for exons in &gene.exons {
            let order = crate::sort::argsort_i64(&exons.iter().map(|e| e[0]).collect::<Vec<_>>());
            for pair in order.windows(2) {
                gene.introns_anno
                    .insert([exons[pair[0]][1], exons[pair[1]][0]]);
            }
        }
    }
    Ok(result)
}

impl Annotation {
    fn exclude(&mut self, suffix: &str, removed: BTreeSet<String>, sorted_report: bool) {
        if removed.is_empty() {
            return;
        }
        let mut ids = Vec::new();
        self.genes.retain(|g| {
            if removed.contains(&g.name) {
                ids.push(g.name.clone());
                false
            } else {
                true
            }
        });
        if sorted_report {
            ids.sort();
        }
        self.excluded.push((suffix.to_owned(), ids));
    }
}
