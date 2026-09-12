// Adapted from SplAdder v3.1.1 spladder_prep.py and settings.py (BSD-3-Clause).
use crate::{
    annotation::{self, AnnotationFilters, Gene},
    annotation_loci::{self, AnnotationMode},
    cache,
    reads::{ReadFilter, ReadOptions},
    sparse,
};
use anyhow::{Context, Result, ensure};
use clap::{ArgAction, Args};
use regex::Regex;
use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

#[derive(Args, Debug)]
pub struct PrepArgs {
    #[arg(long = "bams", short = 'b', default_value = "-")]
    pub bams: String,
    #[arg(long = "sparse-bam", action = ArgAction::SetTrue)]
    pub sparse_bam: bool,
    #[arg(long = "readlen", short = 'n', default_value_t = 50)]
    pub readlen: u64,
    #[arg(long = "confidence", short = 'c', default_value_t = 3, value_parser = clap::value_parser!(u8).range(0..=3))]
    pub confidence: u8,
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
    #[arg(long = "ignore-mismatches", action = ArgAction::SetTrue)]
    pub ignore_mismatches: bool,
    #[arg(long = "reference")]
    pub ref_genome: Option<PathBuf>,
    #[arg(long = "annotation", short = 'a', default_value = "-")]
    pub annotation: PathBuf,
    /// Disambiguate GTF gene/transcript placements before graph construction.
    #[arg(long, value_enum, default_value_t = AnnotationMode::Spladder)]
    pub annotation_mode: AnnotationMode,
    #[arg(long = "filter-overlap-genes", action = ArgAction::SetTrue)]
    pub filter_overlap_genes: bool,
    #[arg(long = "filter-overlap-exons", action = ArgAction::SetTrue)]
    pub filter_overlap_exons: bool,
    #[arg(long = "filter-overlap-transcripts", action = ArgAction::SetTrue)]
    pub filter_overlap_transcripts: bool,
    #[arg(long = "verbose", short = 'v', action = ArgAction::SetTrue)]
    pub verbose: bool,
    #[arg(long = "parallel", default_value_t = 1, value_parser = clap::value_parser!(u16).range(1..=64))]
    pub parallel: u16,
    #[arg(long = "tmp-dir")]
    pub tmpdir: Option<PathBuf>,
}

pub fn prepare_annotation(path: &Path, filters: AnnotationFilters) -> Result<Vec<Gene>> {
    ensure!(
        path.is_file(),
        "annotation does not exist: {}",
        path.display()
    );
    if path.extension().is_some_and(|s| s == "hdf5") {
        return cache::read_genes(path);
    }
    let mut name = path.as_os_str().to_owned();
    name.push(".ruspladder.hdf5");
    let cached = PathBuf::from(name);
    if cached.exists() {
        return cache::read_genes(&cached);
    }
    let annotation = annotation::read_annotation(path, filters)?;
    cache::write_genes(&cached, &annotation.genes)?;
    for (suffix, names) in annotation.excluded {
        let mut report = path.as_os_str().to_owned();
        report.push(format!(".genes_excluded_{suffix}"));
        fs::write(report, names.join("\n") + "\n")?;
    }
    Ok(annotation.genes)
}

pub fn prepare_annotation_mode(
    path: &Path,
    filters: AnnotationFilters,
    mode: AnnotationMode,
) -> Result<Vec<Gene>> {
    if mode == AnnotationMode::Locus {
        prepare_annotation(&annotation_loci::normalize_gtf(path)?, filters)
    } else {
        prepare_annotation(path, filters)
    }
}

pub(crate) fn alignments(value: &str, sparse: bool) -> Result<(Vec<PathBuf>, Vec<String>)> {
    let mut names: Vec<String> = value
        .trim_matches(',')
        .split(',')
        .map(str::to_owned)
        .collect();
    if names.first().is_some_and(|p| p.ends_with(".txt")) {
        names = fs::read_to_string(&names[0])?
            .lines()
            .flat_map(|line| {
                line.split('#')
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .map(str::to_owned)
            })
            .collect();
    }
    ensure!(!names.is_empty(), "no alignment files supplied");
    let suffix = Regex::new(r"(?i)(\.bam|\.hdf5)|\.cram$")?;
    let mut samples = Vec::new();
    let mut paths = Vec::new();
    for name in names {
        let path = PathBuf::from(&name);
        let lower = name.to_lowercase();
        let cached_bam = sparse && lower.ends_with(".bam") && path.with_extension("hdf5").is_file();
        ensure!(
            path.is_file() || cached_bam,
            "alignment does not exist: {}",
            path.display()
        );
        if lower.ends_with(".bam") && path.is_file() {
            ensure!(
                Path::new(&format!("{name}.bai")).is_file(),
                "alignment is not indexed: {name}.bai"
            );
        } else if lower.ends_with(".cram") {
            ensure!(
                Path::new(&format!("{name}.crai")).is_file()
                    || path.with_extension("crai").is_file(),
                "alignment is not indexed: {name}"
            );
        }
        samples.push(
            suffix
                .replace_all(
                    path.file_name()
                        .context("alignment filename")?
                        .to_str()
                        .context("non-UTF8 alignment filename")?,
                    "",
                )
                .into_owned(),
        );
        paths.push(path);
    }
    Ok((paths, samples))
}

pub(crate) fn prepare_summaries(
    bams: &[PathBuf],
    chromosomes: &[String],
    reference: Option<&Path>,
    reads: &ReadOptions,
    confidence: u8,
    parallel: usize,
) -> Result<()> {
    for bam in bams {
        if !bam
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("bam") || e.eq_ignore_ascii_case("cram"))
        {
            continue;
        }
        let output = sparse::summary_path(bam, confidence, reads.filter.is_some())?;
        if output.exists() {
            continue;
        }
        let mut reads = reads.clone();
        // summarize_chr does not forward the CLI's ignore-mismatches switch.
        reads.no_mm = false;
        eprintln!("Preparing sparse alignments: {}", output.display());
        cache::atomic_write(&output, |path| {
            sparse::write_summary(
                bam,
                path,
                chromosomes,
                reference,
                &reads,
                sparse::SummaryOptions {
                    parallel,
                    window: 1048576,
                    unstranded: true,
                },
            )
        })?;
    }
    Ok(())
}

pub fn run(o: &PrepArgs) -> Result<()> {
    let bams = if o.bams == "-" {
        Vec::new()
    } else {
        alignments(&o.bams, o.sparse_bam)?.0
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(o.parallel as usize)
        .build()?;
    pool.install(|| {
        let genes = prepare_annotation_mode(
            &o.annotation,
            AnnotationFilters {
                overlap_genes: o.filter_overlap_genes,
                overlap_exons: o.filter_overlap_exons,
                overlap_transcripts: o.filter_overlap_transcripts,
            },
            o.annotation_mode,
        )?;
        let chromosomes: Vec<_> = genes
            .iter()
            .map(|g| g.chr.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        if o.verbose {
            eprintln!("Prepared {} genes", genes.len());
        }
        if o.sparse_bam && !bams.is_empty() {
            let mut reads = ReadOptions {
                filter: Some(ReadFilter::confidence(o.confidence, o.readlen)?),
                primary_only: o.primary_only && !o.no_primary_only,
                var_aware: o.var_aware && !o.no_var_aware,
                mm_tag: o.mm_tag.clone(),
                ..ReadOptions::default()
            };
            prepare_summaries(
                &bams,
                &chromosomes,
                o.ref_genome.as_deref(),
                &reads,
                o.confidence,
                o.parallel as usize,
            )?;
            reads.filter = None;
            prepare_summaries(
                &bams,
                &chromosomes,
                o.ref_genome.as_deref(),
                &reads,
                o.confidence,
                o.parallel as usize,
            )?;
        }
        Ok(())
    })
}
