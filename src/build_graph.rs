// Adapted from SplAdder v3.1.1 core/gen_graphs.py and reads.get_intron_list.
// BSD-3-Clause; see licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    augment::{self, CassetteOptions, RetentionOptions},
    editgraph::{self, RemoveExons},
    intron_edges::{self, Inserted, IntronOptions},
    introns::{self, IntronLists},
    reads::{AlignmentReader, ReadFilter, ReadOptions},
    reference::Reference,
};
use anyhow::{Result, ensure};
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeSet,
    path::{Path, PathBuf},
};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GraphOptions {
    pub reads: ReadOptions,
    /// Upstream retains the original filter object after the first IR stage
    /// restores a copy into read_filter. Subsequent sample builds share this.
    #[serde(default)]
    pub retention_read_filter: Option<ReadFilter>,
    pub cassette: CassetteOptions,
    pub retention: RetentionOptions,
    pub intron_edges: IntronOptions,
    pub remove_exons: RemoveExons,
    pub insert_es: bool,
    pub insert_ir: bool,
    pub insert_ni: bool,
    pub remove_se: bool,
    pub introns_unstranded: bool,
    pub insert_intron_iterations: usize,
    /// None disables filtering; Some(false) is strict, Some(true) lenient.
    pub consensus: Option<bool>,
}

impl GraphOptions {
    pub fn confidence(level: u8, readlen: u64) -> Result<Self> {
        let filter = ReadFilter::confidence(level, readlen)?;
        Ok(Self {
            reads: ReadOptions {
                filter: Some(filter),
                ..Default::default()
            },
            retention_read_filter: None,
            cassette: CassetteOptions {
                min_cassette_cov: 5.0,
                min_cassette_region: 0.9,
                min_cassette_rel_diff: 0.5,
            },
            retention: RetentionOptions {
                min_retention_cov: [1.0, 2.0, 5.0, 10.0][level as usize],
                min_retention_region: if level < 2 { 0.75 } else { 0.9 },
                min_retention_rel_cov: if level < 2 { 0.1 } else { 0.2 },
                max_retention_rel_cov: if level == 0 { 2.0 } else { 1.2 },
                min_retention_max_exon_fold_diff: 4.0,
            },
            intron_edges: IntronOptions {
                min_exon_len: 50,
                vicinity_region: 40,
                insert_intron_retention: true,
                gene_merges: false,
                append_new_terminal_exons: true,
                append_new_terminal_exons_len: 200,
            },
            remove_exons: RemoveExons {
                terminal_short_len: 10,
                terminal_short_extend: 40,
                min_exon_len: 50,
                min_exon_len_remove: 10,
            },
            insert_es: true,
            insert_ir: true,
            insert_ni: true,
            remove_se: false,
            introns_unstranded: false,
            insert_intron_iterations: 5,
            consensus: None,
        })
    }
}

#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct GraphInserted {
    pub cassette_exon: usize,
    pub intron_retention: usize,
    #[serde(flatten)]
    pub edges: Inserted,
}

struct Evidence {
    readers: Vec<AlignmentReader>,
    contigs: BTreeSet<String>,
}

impl Evidence {
    fn open(bams: &[PathBuf], reference: Option<&Path>) -> Result<Self> {
        ensure!(!bams.is_empty(), "graph augmentation requires alignments");
        let readers: Vec<_> = bams
            .iter()
            .map(|path| AlignmentReader::open(path, reference))
            .collect::<Result<_>>()?;
        // init_regions stops after inspecting the first alignment file.
        let contigs = readers[0].contig_names().into_iter().collect();
        Ok(Self { readers, contigs })
    }

    fn coverage(&mut self, gene: &Gene, options: &ReadOptions) -> Result<Vec<u64>> {
        let mut options = options.clone();
        options.strand = Some(gene.strand);
        let mut track = vec![0; (gene.stop - gene.start) as usize];
        for reader in &mut self.readers {
            let evidence = reader.region(&gene.chr, gene.start, gene.stop, &options)?;
            for (total, count) in track.iter_mut().zip(evidence.coverage) {
                *total += count;
            }
        }
        Ok(track)
    }

    fn introns(
        &mut self,
        genes: &[Gene],
        options: &ReadOptions,
        unstranded: bool,
    ) -> Result<IntronLists> {
        let mut lists: IntronLists = vec![Default::default(); genes.len()];
        let mut options = options.clone();
        options.mapped = false;
        // Preserve upstream's compacted target index when an annotation
        // contig is absent from the first alignment's region list.
        for (index, gene) in genes
            .iter()
            .filter(|g| self.contigs.contains(&g.chr))
            .enumerate()
        {
            options.strand = Some(gene.strand);
            let start = (gene.start - 5000).max(1);
            let stop = gene.stop + 5000;
            let strand = usize::from(gene.strand == '-');
            let found = &mut lists[index][strand];
            for reader in &mut self.readers {
                let evidence = reader.junctions(&gene.chr, start, stop, &options)?;
                for (si, introns) in [evidence.introns_plus, evidence.introns_minus]
                    .into_iter()
                    .enumerate()
                {
                    if unstranded || si == strand {
                        found.extend(introns.into_iter().filter(|v| {
                            v[0] > start
                                && v[1] < stop
                                && options
                                    .filter
                                    .as_ref()
                                    .is_none_or(|f| v[2] >= f.mincount as i64)
                        }));
                    }
                }
            }
            // Direct BAM multi-sample mode retains separate support rows.
            found.sort();
        }
        Ok(lists)
    }
}

pub fn generate(
    mut genes: Vec<Gene>,
    bams: &[PathBuf],
    options: &mut GraphOptions,
    reference: Option<&Path>,
) -> Result<(Vec<Gene>, GraphInserted)> {
    ensure!(!genes.is_empty(), "no genes to augment");
    let mut evidence = Evidence::open(bams, reference)?;
    for gene in &mut genes {
        gene.splicegraph.sort();
        gene.label_alt();
        gene.splicegraph.update_terminals();
        gene.start = *gene
            .exons
            .iter()
            .flatten()
            .flatten()
            .min()
            .ok_or_else(|| anyhow::anyhow!("gene has no annotated exons"))?;
        gene.stop = *gene.exons.iter().flatten().flatten().max().unwrap();
    }
    let order = crate::sort::argsort_i64(&genes.iter().map(|g| g.stop).collect::<Vec<_>>());
    let mut slots: Vec<_> = genes.into_iter().map(Some).collect();
    genes = order
        .into_iter()
        .map(|i| slots[i].take().unwrap())
        .collect();
    genes.sort_by_key(|g| g.start);
    genes.sort_by(|a, b| a.chr.cmp(&b.chr));
    let mut introns: IntronLists = vec![Default::default(); genes.len()];
    if options.insert_es || options.insert_ir || options.insert_ni {
        introns = evidence.introns(&genes, &options.reads, options.introns_unstranded)?;
        if let Some(lenient) = options.consensus {
            introns::filter_consensus(
                &mut introns,
                &genes,
                &Reference::open(reference.ok_or_else(|| {
                    anyhow::anyhow!("consensus filtering requires a reference FASTA")
                })?)?,
                lenient,
            )?;
        }
        introns::filter_ambiguous(
            &mut introns,
            &genes,
            options.intron_edges.append_new_terminal_exons_len,
        )?;
        let mut filter = options
            .reads
            .filter
            .clone()
            .ok_or_else(|| anyhow::anyhow!("graph augmentation requires a read filter"))?;
        introns::make_feasible(&mut introns, &mut filter, |indices, filter| {
            let selected: Vec<_> = indices.iter().map(|&i| genes[i].clone()).collect();
            let mut reads = options.reads.clone();
            reads.filter = Some(filter.clone());
            evidence.introns(&selected, &reads, options.introns_unstranded)
        })?;
        options.reads.filter = Some(filter);
    }
    let mut inserted = GraphInserted::default();
    if options.insert_es {
        for (gene, lists) in genes.iter_mut().zip(&introns) {
            let track = evidence.coverage(gene, &options.reads)?;
            let pairs: Vec<_> = lists[usize::from(gene.strand == '-')]
                .iter()
                .map(|v| [v[0], v[1]])
                .collect();
            inserted.cassette_exon +=
                augment::insert_cassettes(gene, &pairs, &track, options.cassette)?;
        }
    }
    if options.insert_ir {
        if options.retention_read_filter.is_none() {
            options.retention_read_filter = options.reads.filter.clone();
        }
        let mut reads = options.reads.clone();
        reads.filter = options.retention_read_filter.clone();
        for gene in &mut genes {
            let track = evidence.coverage(gene, &reads)?;
            inserted.intron_retention +=
                augment::insert_retentions(gene, &track, options.retention)?;
        }
    }
    if options.remove_se {
        for gene in &mut genes {
            editgraph::remove_short_exons(gene, options.remove_exons);
        }
    }
    if options.insert_ni {
        let mut start = 0;
        while start < genes.len() {
            let stop = start + genes[start..].partition_point(|g| g.chr == genes[start].chr);
            let lists: Vec<_> = genes[start..stop]
                .iter()
                .zip(&introns[start..stop])
                .map(|(g, v)| {
                    v[usize::from(g.strand == '-')]
                        .iter()
                        .map(|p| [p[0], p[1]])
                        .collect()
                })
                .collect();
            for _ in 0..options.insert_intron_iterations {
                let before = genes[start..stop].to_vec();
                let added = intron_edges::insert_edges(
                    &mut genes[start..stop],
                    &lists,
                    options.intron_edges,
                    |gene| evidence.coverage(gene, &options.reads),
                )?;
                inserted.edges.intron_in_exon += added.intron_in_exon;
                inserted.edges.alt_53_prime += added.alt_53_prime;
                inserted.edges.exon_skip += added.exon_skip;
                inserted.edges.gene_merge += added.gene_merge;
                inserted.edges.new_terminal_exon += added.new_terminal_exon;
                for gene in &mut genes[start..stop] {
                    editgraph::merge_duplicate_exons(&mut gene.splicegraph);
                }
                if genes[start..stop] == before {
                    break;
                }
            }
            start = stop;
        }
    }
    for gene in &mut genes {
        gene.start = *gene.splicegraph.vertices.iter().flatten().min().unwrap();
        gene.stop = *gene.splicegraph.vertices.iter().flatten().max().unwrap();
        gene.label_alt();
    }
    Ok((genes, inserted))
}
