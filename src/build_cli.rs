// Adapted from SplAdder v3.1.1 spladder_build.py (GPL-3.0-or-later).
use crate::prep::{alignments, prepare_annotation};
use crate::{
    analyze::{self, AnalysisOptions},
    annotation::{AnnotationFilters, Gene},
    build_graph::{self, GraphOptions},
    cache, count,
    count_io::{self, CountWriter},
    events::{self, EventType},
    expression,
    graph::SegmentGraph,
    output::{self, Format},
    reads::ReadOptions,
    verify::{AltPrimeOptions, MutexOptions, RetentionOptions, SkipOptions, VerifyOptions},
};
use anyhow::{Context, Result, ensure};
use clap::{ArgAction, Args};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

#[derive(Args, Debug)]
pub struct BuildArgs {
    #[arg(long = "bams", short = 'b')]
    pub bams: String,
    #[arg(long = "outdir", short = 'o')]
    pub outdir: PathBuf,
    #[arg(long = "annotation", short = 'a')]
    pub annotation: PathBuf,
    #[arg(long = "parallel", default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..=64))]
    pub parallel: u16,
    #[arg(long = "verbose", short = 'v', action = ArgAction::SetTrue)]
    pub verbose: bool,
    #[arg(long = "debug", short = 'd', action = ArgAction::SetTrue)]
    pub debug: bool,
    #[arg(long = "readlen", short = 'n', default_value_t = 50)]
    pub readlen: u64,
    #[arg(long = "primary-only", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_primary_only")]
    pub primary_only: bool,
    #[arg(long = "no-primary-only", action = ArgAction::SetTrue, overrides_with = "primary_only")]
    pub no_primary_only: bool,
    #[arg(long = "var-aware", action = ArgAction::SetTrue, overrides_with = "no_var_aware")]
    pub var_aware: bool,
    #[arg(long = "no-var-aware", action = ArgAction::SetTrue, overrides_with = "var_aware")]
    pub no_var_aware: bool,
    #[arg(long = "set-mm-tag", default_value = "NM")]
    pub mm_tag: String,
    #[arg(long = "labels", default_value = "-")]
    pub labels: String,
    #[arg(long = "filter-overlap-genes", action = ArgAction::SetTrue)]
    pub filter_overlap_genes: bool,
    #[arg(long = "filter-overlap-exons", action = ArgAction::SetTrue)]
    pub filter_overlap_exons: bool,
    #[arg(long = "filter-overlap-transcripts", action = ArgAction::SetTrue)]
    pub filter_overlap_transcripts: bool,
    #[arg(long = "filter-consensus", default_value = "")]
    pub filter_consensus: String,
    #[arg(long = "ignore-mismatches", action = ArgAction::SetTrue)]
    pub ignore_mismatches: bool,
    #[arg(long = "reference")]
    pub ref_genome: Option<PathBuf>,
    #[arg(long = "logfile", short = 'l', default_value = "-")]
    pub logfile: String,
    #[arg(long = "output-txt", action = ArgAction::SetTrue)]
    pub output_txt: bool,
    #[arg(long = "output-txt-conf", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_output_confirmed_txt")]
    pub output_confirmed_txt: bool,
    #[arg(long = "no-output-txt-conf", action = ArgAction::SetTrue, overrides_with = "output_confirmed_txt")]
    pub no_output_confirmed_txt: bool,
    #[arg(long = "output-gff3", action = ArgAction::SetTrue)]
    pub output_gff3: bool,
    #[arg(long = "output-gff3-conf", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_output_confirmed_gff3")]
    pub output_confirmed_gff3: bool,
    #[arg(long = "no-output-gff3-conf", action = ArgAction::SetTrue, overrides_with = "output_confirmed_gff3")]
    pub no_output_confirmed_gff3: bool,
    #[arg(long = "output-struc", action = ArgAction::SetTrue)]
    pub output_struc: bool,
    #[arg(long = "output-struc-conf", action = ArgAction::SetTrue)]
    pub output_confirmed_struc: bool,
    #[arg(long = "output-bed", action = ArgAction::SetTrue)]
    pub output_bed: bool,
    #[arg(long = "output-conf-bed", action = ArgAction::SetTrue)]
    pub output_confirmed_bed: bool,
    #[arg(long = "output-conf-tcga", action = ArgAction::SetTrue)]
    pub output_confirmed_tcga: bool,
    #[arg(long = "output-conf-icgc", action = ArgAction::SetTrue)]
    pub output_confirmed_icgc: bool,
    #[arg(long = "sparse-bam", action = ArgAction::SetTrue)]
    pub sparse_bam: bool,
    #[arg(long = "compress-text", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_compress_text")]
    pub compress_text: bool,
    #[arg(long = "no-compress-text", action = ArgAction::SetTrue, overrides_with = "compress_text")]
    pub no_compress_text: bool,
    #[arg(long = "tmp-dir")]
    pub tmpdir: Option<PathBuf>,
    #[arg(long = "confidence", short = 'c', default_value_t = 3, value_parser = clap::value_parser!(u8).range(0..=3))]
    pub confidence: u8,
    #[arg(long = "iterations", short = 'I', default_value_t = 5)]
    pub insert_intron_iterations: usize,
    #[arg(long = "merge-strat", short = 'M', default_value = "merge_graphs")]
    pub merge: String,
    #[arg(long = "chunked-merge", num_args = 4, action = ArgAction::Append)]
    pub chunked_merge: Vec<usize>,
    #[arg(long = "chunksize", default_value_t = 10)]
    pub chunksize: usize,
    #[arg(long = "insert-ir", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_insert_ir")]
    pub insert_ir: bool,
    #[arg(long = "no-insert-ir", action = ArgAction::SetTrue, overrides_with = "insert_ir")]
    pub no_insert_ir: bool,
    #[arg(long = "insert-es", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_insert_es")]
    pub insert_es: bool,
    #[arg(long = "no-insert-es", action = ArgAction::SetTrue, overrides_with = "insert_es")]
    pub no_insert_es: bool,
    #[arg(long = "insert-ni", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_insert_ni")]
    pub insert_ni: bool,
    #[arg(long = "no-insert-ni", action = ArgAction::SetTrue, overrides_with = "insert_ni")]
    pub no_insert_ni: bool,
    #[arg(long = "remove-se", action = ArgAction::SetTrue, overrides_with = "no_remove_se")]
    pub remove_se: bool,
    #[arg(long = "no-remove-se", action = ArgAction::SetTrue, overrides_with = "remove_se")]
    pub no_remove_se: bool,
    #[arg(long = "validate-sg", action = ArgAction::SetTrue, overrides_with = "no_validate_sg")]
    pub validate_sg: bool,
    #[arg(long = "no-validate-sg", action = ArgAction::SetTrue, overrides_with = "validate_sg")]
    pub no_validate_sg: bool,
    #[arg(long = "re-infer-sg", action = ArgAction::SetTrue, overrides_with = "no_infer_sg")]
    pub infer_sg: bool,
    #[arg(long = "no-re-infer-sg", action = ArgAction::SetTrue, overrides_with = "infer_sg")]
    pub no_infer_sg: bool,
    #[arg(long = "validate-sg-count", default_value_t = 10)]
    pub sg_min_edge_count: usize,
    #[arg(
        long = "event-types",
        default_value = "exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons"
    )]
    pub event_types: String,
    #[arg(long = "extract-ase", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_extract_as")]
    pub extract_as: bool,
    #[arg(long = "no-extract-ase", action = ArgAction::SetTrue, overrides_with = "extract_as")]
    pub no_extract_as: bool,
    #[arg(long = "ase-edge-limit", default_value_t = 500)]
    pub detect_edge_limit: usize,
    #[arg(long = "curate-alt-prime", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_curate_alt_prime")]
    pub curate_alt_prime: bool,
    #[arg(long = "no-curate-alt-prime", action = ArgAction::SetTrue, overrides_with = "curate_alt_prime")]
    pub no_curate_alt_prime: bool,
    #[arg(long = "quantify-graph", action = ArgAction::SetTrue, default_value_t = true, overrides_with = "no_quantify_graph")]
    pub quantify_graph: bool,
    #[arg(long = "no-quantify-graph", action = ArgAction::SetTrue, overrides_with = "quantify_graph")]
    pub no_quantify_graph: bool,
    #[arg(long = "use-anno-support", action = ArgAction::SetTrue, overrides_with = "no_use_anno_support")]
    pub use_anno_support: bool,
    #[arg(long = "no-use-anno-support", action = ArgAction::SetTrue, overrides_with = "use_anno_support")]
    pub no_use_anno_support: bool,
    #[arg(long = "psi-min-reads", default_value_t = 10)]
    pub psi_min_reads: usize,
    #[arg(long = "qmode", default_value = "all")]
    pub qmode: String,
}

impl BuildArgs {
    fn validated(&self) -> &'static str {
        if self.validate_sg && !self.no_validate_sg {
            ".validated"
        } else {
            ""
        }
    }
    fn graph(&self, tag: &str) -> PathBuf {
        self.outdir
            .join("spladder")
            .join(format!("genes_graph_conf{}.{}.hdf5", self.confidence, tag))
    }
    fn counts(&self, tag: &str) -> PathBuf {
        self.outdir.join("spladder").join(format!(
            "genes_graph_conf{}.{}.count.hdf5",
            self.confidence, tag
        ))
    }
    fn event_base(&self, tag: &str, kind: EventType) -> PathBuf {
        self.outdir
            .join(format!("{}_{}_C{}", tag, kind.as_str(), self.confidence))
    }
    fn verify(&self) -> Result<VerifyOptions> {
        let graph = GraphOptions::confidence(self.confidence, self.readlen)?;
        let skip = SkipOptions {
            min_skip_rel_cov: 0.05,
            min_non_skip_count: 3.0,
            min_skip_count: 3.0,
        };
        Ok(VerifyOptions {
            use_anno_support: self.use_anno_support && !self.no_use_anno_support,
            exon_skip: skip,
            mult_exon_skip: skip,
            intron_retention: RetentionOptions {
                min_retention_cov: graph.retention.min_retention_cov,
                min_retention_region: graph.retention.min_retention_region,
                min_retention_rel_cov: graph.retention.min_retention_rel_cov,
                min_non_retention_count: 3.0,
            },
            alt_prime: AltPrimeOptions {
                min_diff_rel_cov: 0.05,
                min_intron_count: 3.0,
            },
            mutex_exons: MutexOptions {
                min_skip_rel_cov: 0.05,
                min_conf_count: 2.0,
            },
        })
    }
    fn merge_graphs(&self, samples: &[String]) -> Result<()> {
        let mut output = self.graph(&self.merge);
        if output.exists() {
            return Ok(());
        }
        let mut paths: Vec<_> = samples.iter().map(|s| self.graph(s)).collect();
        if !self.chunked_merge.is_empty() {
            let [level, max_level, start, end] = self.chunked_merge[..4].try_into().unwrap();
            if level > 1 {
                let prefix = format!(
                    "genes_graph_conf{}.{}_level{}_chunk",
                    self.confidence,
                    self.merge,
                    level - 1
                );
                paths = fs::read_dir(self.outdir.join("spladder"))?
                    .map(|entry| entry.map(|e| e.path()))
                    .collect::<std::io::Result<Vec<_>>>()?
                    .into_iter()
                    .filter(|p| {
                        p.file_name()
                            .is_some_and(|n| n.to_string_lossy().starts_with(&prefix))
                            && p.extension().is_some_and(|n| n == "hdf5")
                    })
                    .collect();
            }
            paths.sort();
            let end = end.min(paths.len());
            if level == max_level {
                ensure!(
                    paths.len() <= self.chunksize,
                    "chunksize is {} but merge list has length {}",
                    self.chunksize,
                    paths.len()
                );
            } else {
                output = self.graph(&format!("{}_level{level}_chunk{start}_{end}", self.merge));
            }
            if output.exists() {
                return Ok(());
            }
            paths = paths
                .into_iter()
                .skip(start)
                .take(end.saturating_sub(start))
                .collect();
        }
        // do_merge_all is false in upstream default_settings, including merge_all.
        let genes = crate::merge::merge_samples(paths.iter().map(|p| cache::read_genes(p)))?;
        cache::write_genes(&output, &genes)
    }

    fn report_events(
        &self,
        tag: &str,
        kind: EventType,
        genes: &[Gene],
        samples: &[String],
        sample_idx: &[usize],
    ) -> Result<()> {
        let base = self.event_base(tag, kind);
        let events = cache::read_events(&with_suffix(&base, ".events.hdf5"))?;
        let counts = with_suffix(&base, ".counts.hdf5");
        let labels: Vec<_> = samples.iter().map(String::as_str).collect();
        let confirmed = if counts.exists() {
            let file = hdf5::File::open(&counts)?;
            if file.link_exists("conf_idx") {
                file.dataset("conf_idx")?
                    .read_raw::<u64>()?
                    .into_iter()
                    .map(|i| i as usize)
                    .collect()
            } else {
                Vec::new()
            }
        } else {
            cache::atomic_write(&counts, |path| {
                analyze::analyze(
                    genes,
                    &events,
                    kind,
                    &self.counts(&format!("{tag}{}", self.validated())),
                    path,
                    &AnalysisOptions {
                        samples: &labels,
                        sample_idx,
                        verify: self.verify()?,
                        psi_min_reads: self.psi_min_reads as f64,
                    },
                )
            })?
        };
        if events.is_empty() {
            return Ok(());
        }
        let gz = if self.compress_text && !self.no_compress_text {
            ".gz"
        } else {
            ""
        };
        let all: Vec<_> = (0..events.len()).collect();
        let selected_labels: Vec<_> = if self.merge == "single" {
            sample_idx.iter().map(|&i| labels[i]).collect()
        } else {
            labels.clone()
        };
        let outputs = [
            (self.output_txt, false, Format::Txt, format!(".txt{gz}")),
            (
                self.output_struc,
                false,
                Format::Structured,
                format!(".struc.txt{gz}"),
            ),
            (self.output_gff3, false, Format::Gff3, ".gff3".into()),
            (
                self.output_confirmed_gff3 && !self.no_output_confirmed_gff3,
                true,
                Format::Gff3,
                ".confirmed.gff3".into(),
            ),
            (
                self.output_confirmed_txt && !self.no_output_confirmed_txt,
                true,
                Format::Txt,
                format!(".confirmed.txt{gz}"),
            ),
            (
                self.output_confirmed_bed,
                true,
                Format::Bed,
                ".confirmed.bed".into(),
            ),
            (self.output_bed, false, Format::Bed, ".bed".into()),
            (
                self.output_confirmed_struc,
                true,
                Format::Structured,
                format!(".confirmed.struc.txt{gz}"),
            ),
            (
                self.output_confirmed_tcga,
                true,
                Format::Tcga,
                format!(".confirmed.tcga.txt{gz}"),
            ),
            (
                self.output_confirmed_icgc,
                true,
                Format::Icgc,
                format!(".confirmed.icgc.txt{gz}"),
            ),
        ];
        for (position, (enabled, conf, format, suffix)) in outputs.into_iter().enumerate() {
            // Source checks confirmation before all outputs except full TXT/struc.
            if !enabled || (position >= 2 && confirmed.is_empty()) {
                continue;
            }
            let path = with_suffix(&base, &suffix);
            if !path.exists() {
                output::write(
                    &path,
                    format,
                    &events,
                    &counts,
                    if conf { &selected_labels } else { &labels },
                    if conf { &confirmed } else { &all },
                )?;
            }
        }
        Ok(())
    }
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

pub fn run(options: &BuildArgs) -> Result<()> {
    ensure!(
        options.filter_consensus.is_empty() || options.ref_genome.is_some(),
        "consensus filtering requires --reference"
    );
    ensure!(
        ["single", "merge_graphs", "merge_bams", "merge_all"].contains(&options.merge.as_str()),
        "unknown merge strategy: {}",
        options.merge
    );
    ensure!(
        ["single", "all", "collect"].contains(&options.qmode.as_str()),
        "unknown quantification mode: {}",
        options.qmode
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.parallel as usize)
        .build()?;
    pool.install(|| run_build(options))
}

fn run_build(o: &BuildArgs) -> Result<()> {
    let (bams, samples) = alignments(&o.bams, o.sparse_bam)?;
    let kinds: Vec<EventType> = o
        .event_types
        .trim_matches(',')
        .split(',')
        .map(|name| {
            EventType::ALL
                .into_iter()
                .find(|k| k.as_str() == name)
                .with_context(|| format!("unknown event type: {name}"))
        })
        .collect::<Result<_>>()?;
    fs::create_dir_all(o.outdir.join("spladder"))?;
    fs::create_dir_all(o.tmpdir.clone().unwrap_or_else(|| o.outdir.join("tmp")))?;
    let annotation = prepare_annotation(
        &o.annotation,
        AnnotationFilters {
            overlap_genes: o.filter_overlap_genes,
            overlap_exons: o.filter_overlap_exons,
            overlap_transcripts: o.filter_overlap_transcripts,
        },
    )?;
    let chromosomes: BTreeMap<_, _> = annotation
        .iter()
        .map(|g| g.chr.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .enumerate()
        .map(|(i, chr)| (chr, i))
        .collect();
    let reads = ReadOptions {
        sparse_confidence: o.sparse_bam.then_some(o.confidence),
        primary_only: o.primary_only && !o.no_primary_only,
        var_aware: o.var_aware && !o.no_var_aware,
        no_mm: o.ignore_mismatches,
        mm_tag: o.mm_tag.clone(),
        ..ReadOptions::default()
    };
    let mut graph_options = GraphOptions::confidence(o.confidence, o.readlen)?;
    graph_options.reads = ReadOptions {
        filter: graph_options.reads.filter.take(),
        ..reads.clone()
    };
    graph_options.insert_ir = o.insert_ir && !o.no_insert_ir;
    graph_options.insert_es = o.insert_es && !o.no_insert_es;
    graph_options.insert_ni = o.insert_ni && !o.no_insert_ni;
    graph_options.remove_se = o.remove_se && !o.no_remove_se;
    graph_options.insert_intron_iterations = o.insert_intron_iterations;
    graph_options.consensus = match o.filter_consensus.as_str() {
        "" => None,
        "strict" => Some(false),
        "lenient" => Some(true),
        _ => anyhow::bail!("consensus must be strict or lenient"),
    };
    if o.merge != "merge_graphs" || !o.graph(&o.merge).exists() {
        if o.sparse_bam {
            crate::prep::prepare_summaries(
                &bams,
                &chromosomes.keys().cloned().collect::<Vec<_>>(),
                o.ref_genome.as_deref(),
                &graph_options.reads,
                o.confidence,
                o.parallel as usize,
            )?;
        }
        let mut generate = |tag: &str, paths: &[PathBuf]| -> Result<()> {
            let output = o.graph(tag);
            if !output.exists() {
                ensure!(
                    !o.infer_sg || o.no_infer_sg,
                    "SplAdder v3.1.1 --re-infer-sg calls undefined infer_splice_graph"
                );
                eprintln!("Generating graph: {tag}");
                let (genes, inserted) = build_graph::generate(
                    annotation.clone(),
                    paths,
                    &mut graph_options,
                    o.ref_genome.as_deref(),
                )?;
                cache::write_genes(&output, &genes)?;
                if o.logfile != "-" {
                    fs::write(
                        &o.logfile,
                        format!("{}\n", serde_json::to_string(&inserted)?),
                    )?;
                }
            }
            Ok(())
        };
        if ["single", "merge_graphs", "merge_all"].contains(&o.merge.as_str()) {
            for (sample, bam) in samples.iter().zip(&bams) {
                generate(sample, std::slice::from_ref(bam))?;
            }
        }
        if ["merge_bams", "merge_all"].contains(&o.merge.as_str()) {
            generate("merge_bams", &bams)?;
        }
        if ["merge_graphs", "merge_all"].contains(&o.merge.as_str()) {
            o.merge_graphs(&samples)?;
        }
    }
    drop(annotation);
    let merged_tag = format!("{}{}", o.merge, o.validated());
    if o.merge == "merge_graphs" && !o.validated().is_empty() && !o.graph(&merged_tag).exists() {
        let mut genes = cache::read_genes(&o.graph(&o.merge))?;
        crate::merge::filter_edge_support(
            &mut genes,
            o.sg_min_edge_count.min(samples.len()) as u64,
        )?;
        cache::write_genes(&o.graph(&merged_tag), &genes)?;
    }
    let indices: Vec<_> = if o.merge == "single" || o.qmode == "collect" {
        (0..samples.len()).collect()
    } else {
        vec![0]
    };
    if o.sparse_bam {
        crate::prep::prepare_summaries(
            &bams,
            &chromosomes.keys().cloned().collect::<Vec<_>>(),
            o.ref_genome.as_deref(),
            &reads,
            o.confidence,
            o.parallel as usize,
        )?;
    }
    let nonfinal_chunk = !o.chunked_merge.is_empty() && o.chunked_merge[0] < o.chunked_merge[1];
    if o.quantify_graph && !o.no_quantify_graph && !nonfinal_chunk {
        for &index in &indices {
            let tag = if o.merge == "single" {
                samples[index].clone()
            } else {
                merged_tag.clone()
            };
            let count_tag = if o.merge == "single" {
                tag.clone()
            } else if o.qmode == "single" {
                format!(
                    "{}.{}{}",
                    o.merge,
                    if o.merge == "merge_graphs" {
                        &samples[0]
                    } else {
                        "None"
                    },
                    o.validated()
                )
            } else {
                tag.clone()
            };
            let counts = o.counts(&count_tag);
            let graph_path = o.graph(&tag);
            let mut genes = cache::read_genes(&graph_path)?;
            if !counts.exists() {
                if o.merge == "merge_graphs" && o.qmode == "collect" {
                    let paths: Vec<_> = samples
                        .iter()
                        .map(|s| o.counts(&format!("{}.{s}{}", o.merge, o.validated())))
                        .collect();
                    cache::atomic_write(&counts, |path| count_io::collect(&paths, path))?;
                } else {
                    ensure!(!genes.is_empty(), "cannot quantify empty graph");
                    if genes[0].segmentgraph.segments.is_empty() {
                        for gene in &mut genes {
                            gene.segmentgraph = SegmentGraph::from_splice_graph(&gene.splicegraph);
                        }
                        cache::write_genes(&graph_path, &genes)?;
                    }
                    let selected: Vec<_> = if o.merge == "single" {
                        vec![index]
                    } else if o.merge == "merge_graphs" && o.qmode == "single" {
                        vec![0]
                    } else {
                        (0..bams.len()).collect()
                    };
                    eprintln!("Quantifying graph: {count_tag}");
                    cache::atomic_write(&counts, |path| {
                        let labels: Vec<_> = samples.iter().map(String::as_str).collect();
                        let mut writer =
                            CountWriter::create(path, &genes, &labels, selected.len())?;
                        for sample in selected {
                            writer.append(&count::count_sample(
                                &mut genes,
                                &bams[sample],
                                o.ref_genome.as_deref(),
                                &reads,
                            )?)?;
                        }
                        writer.finish()
                    })?;
                }
            }
            let expression_path = o.outdir.join("spladder").join(format!(
                "genes_graph_conf{}.{count_tag}.gene_exp.hdf5",
                o.confidence
            ));
            if !expression_path.exists() {
                cache::atomic_write(&expression_path, |path| {
                    expression::compute(
                        &genes,
                        &counts,
                        o.readlen as f64,
                        if o.merge == "single" {
                            Some(&[0])
                        } else {
                            None
                        },
                        Some(path),
                    )
                })?;
            }
        }
    }
    if o.extract_as && !o.no_extract_as {
        // Source collect_events stops after the first sample in single mode.
        let tag = if o.merge == "single" {
            &samples[0]
        } else {
            &o.merge
        };
        let gene_tag = if o.merge == "single" {
            tag.clone()
        } else {
            merged_tag.clone()
        };
        let graph_path = o.graph(&gene_tag);
        if kinds
            .iter()
            .any(|&k| !with_suffix(&o.event_base(tag, k), ".events.hdf5").exists())
        {
            let mut collected = if graph_path.exists() {
                events::collect(
                    &cache::read_genes(&graph_path)?,
                    &chromosomes,
                    o.detect_edge_limit,
                    o.curate_alt_prime && !o.no_curate_alt_prime,
                )
            } else {
                BTreeMap::new()
            };
            for &kind in &kinds {
                let path = with_suffix(&o.event_base(tag, kind), ".events.hdf5");
                if !path.exists() {
                    cache::write_events(&path, &collected.remove(&kind).unwrap_or_default())?;
                }
            }
        }
        if o.quantify_graph && !o.no_quantify_graph {
            for &index in &indices {
                let tag = if o.merge == "single" {
                    &samples[index]
                } else {
                    &o.merge
                };
                let genes = cache::read_genes(&o.graph(&format!("{tag}{}", o.validated())))?;
                let sample_idx: Vec<_> = if o.merge == "single" {
                    vec![index]
                } else {
                    (0..samples.len()).collect()
                };
                for &kind in &kinds {
                    o.report_events(tag, kind, &genes, &samples, &sample_idx)?;
                }
            }
        }
    }
    Ok(())
}
