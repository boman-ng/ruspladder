// Adapted from SplAdder v3.1.1 core/gen_graphs.py and reads.get_intron_list.
// BSD-3-Clause; see licenses/SplAdder-BSD.txt.
use crate::{
    annotation::Gene,
    augment::{self, CassetteOptions, RetentionOptions},
    editgraph::{self, RemoveExons},
    intron_edges::{self, Inserted, IntronOptions},
    introns::{self, IntronLists},
    reads::{EvidenceReader, ReadFilter, ReadOptions},
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
    readers: Vec<EvidenceReader>,
    contigs: BTreeSet<String>,
    samples: usize,
}

impl Evidence {
    fn open(bams: &[PathBuf], reference: Option<&Path>, options: &ReadOptions) -> Result<Self> {
        ensure!(!bams.is_empty(), "graph augmentation requires alignments");
        // Source sparse stages reuse the first file's chromosome cache across
        // the complete multi-file loop.
        let paths = if options.sparse_confidence.is_some() {
            &bams[..1]
        } else {
            bams
        };
        let readers: Vec<_> = paths
            .iter()
            .map(|path| EvidenceReader::open(path, reference, options))
            .collect::<Result<_>>()?;
        // init_regions skips missing BAMs and stops at the first existing one.
        let contigs = if options.sparse_confidence.is_some() {
            match bams.iter().find(|p| p.exists()) {
                Some(path) => EvidenceReader::open(path, reference, options)?.contig_names()?,
                None => Vec::new(),
            }
        } else {
            readers[0].contig_names()?
        }
        .into_iter()
        .collect();
        Ok(Self {
            readers,
            contigs,
            samples: bams.len(),
        })
    }

    fn coverage(&mut self, gene: &Gene, options: &ReadOptions) -> Result<Vec<u64>> {
        let mut options = options.clone();
        options.strand = if options.sparse_confidence.is_some() {
            None
        } else {
            Some(gene.strand)
        };
        let mut track = vec![0; (gene.stop - gene.start) as usize];
        for reader in &mut self.readers {
            let evidence = reader.region(&gene.chr, gene.start, gene.stop, &options)?;
            ensure!(
                evidence.coverage.len() == track.len(),
                "graph exons extend beyond sparse coverage"
            );
            for (total, count) in track.iter_mut().zip(evidence.coverage) {
                *total += count;
            }
        }
        if options.sparse_confidence.is_some() && self.samples > 1 {
            // prep_sparse_bam writes collapsed, single-row uint32 coverage.
            for value in &mut track {
                *value = (*value as u32).wrapping_mul(self.samples as u32) as u64;
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
        if options.sparse_confidence.is_some() {
            return self.sparse_introns(genes, options, unstranded);
        }
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

    fn sparse_introns(
        &mut self,
        genes: &[Gene],
        options: &ReadOptions,
        unstranded: bool,
    ) -> Result<IntronLists> {
        let mut lists: IntronLists = vec![Default::default(); genes.len()];
        let selected: Vec<_> = genes
            .iter()
            .filter(|g| self.contigs.contains(&g.chr))
            .collect();
        let chromosomes: BTreeSet<_> = selected.iter().map(|g| g.chr.as_str()).collect();
        for chromosome in chromosomes {
            for strand in ['+', '-'] {
                for (gene_index, gene) in selected
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| g.chr == chromosome && g.strand == strand)
                {
                    let si = usize::from(strand == '-');
                    let mut reads = options.clone();
                    reads.strand = Some(strand);
                    let evidence = self.readers[0].junctions(
                        chromosome,
                        (gene.start - 5000).max(1),
                        gene.stop + 5000,
                        &reads,
                    )?;
                    let mut found = Vec::new();
                    for (s, rows) in [evidence.introns_plus, evidence.introns_minus]
                        .into_iter()
                        .enumerate()
                    {
                        if unstranded || s == si {
                            found.extend(rows);
                        }
                    }
                    let one_sample = found.clone();
                    for _ in 1..self.samples {
                        found.extend_from_slice(&one_sample);
                    }
                    found.sort();
                    let mut index = gene_index;
                    if self.samples > 1 {
                        let mut keep = vec![true; found.len()];
                        // Upstream accidentally reuses its gene index as the
                        // duplicate-intron loop variable. Preserve that target
                        // and report its observed out-of-bounds failure.
                        for i in 1..found.len() {
                            index = i;
                            if found[i][..2] == found[i - 1][..2] {
                                found[i][2] = (found[i][2] as u32)
                                    .wrapping_add(found[i - 1][2] as u32)
                                    as i64;
                                keep[i - 1] = false;
                            }
                        }
                        found = found
                            .into_iter()
                            .zip(keep)
                            .filter_map(|(row, keep)| keep.then_some(row))
                            .collect();
                    }
                    ensure!(
                        index < lists.len(),
                        "SplAdder sparse multi-BAM intron index {index} is out of bounds for {} genes",
                        lists.len()
                    );
                    lists[index][si] = found;
                }
            }
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
    let mut evidence = Evidence::open(bams, reference, &options.reads)?;
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
