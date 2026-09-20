use crate::{
    cache,
    events::EventType,
    hdf5io, quantify,
    statistics::{self, TestOptions},
    test_output::{self, Metadata},
    testing::{self, InputOptions},
};
use anyhow::{Context, Result, ensure};
use clap::Args;
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

#[derive(Args)]
pub struct TestArgs {
    #[arg(short = 'o', long)]
    pub outdir: PathBuf,
    #[arg(short = 'a', long = "conditionA")]
    pub condition_a: String,
    #[arg(short = 'b', long = "conditionB")]
    pub condition_b: String,
    #[arg(short = 'n', long, default_value_t = 50)]
    pub readlen: u64,
    #[arg(short='c',long,default_value_t=3,value_parser=clap::value_parser!(u8).range(0..=3))]
    pub confidence: u8,
    #[arg(short = 'M', long = "merge-strat", default_value = "merge_graphs")]
    pub merge: String,
    #[arg(
        long,
        default_value = "exon_skip,intron_retention,alt_3prime,alt_5prime,mult_exon_skip,mutex_exons"
    )]
    pub event_types: String,
    #[arg(long, overrides_with = "no_validate_sg")]
    pub validate_sg: bool,
    #[arg(long, overrides_with = "validate_sg")]
    pub no_validate_sg: bool,
    #[arg(short = 'C', long, default_value = "BH")]
    pub correction: String,
    #[arg(short = '0', long = "max-zero-frac", default_value_t = 0.5)]
    pub max_zero_fraction: f64,
    #[arg(long = "dpsi", default_value_t = 0.05)]
    pub min_dpsi: f64,
    #[arg(short = 'i', long, default_value_t = 10.0)]
    pub min_count: f64,
    #[arg(long, overrides_with = "no_cap_outliers")]
    pub cap_outliers: bool,
    #[arg(long, overrides_with = "cap_outliers")]
    pub no_cap_outliers: bool,
    #[arg(long, overrides_with = "no_cap_exp_outliers", default_value_t = true)]
    pub cap_exp_outliers: bool,
    #[arg(long, overrides_with = "cap_exp_outliers")]
    pub no_cap_exp_outliers: bool,
    #[arg(short = 'v', long)]
    pub verbose: bool,
    #[arg(short = 'd', long)]
    pub debug: bool,
    #[arg(short = 'D', long)]
    pub diagnose_plots: bool,
    #[arg(long)]
    pub timestamp: bool,
    #[arg(long = "labelA", default_value = "condA")]
    pub label_a: String,
    #[arg(long = "labelB", default_value = "condB")]
    pub label_b: String,
    #[arg(long, default_value = "-")]
    pub out_tag: String,
    #[arg(short = 'f', long, default_value = "png")]
    pub plot_format: String,
    #[arg(long,default_value_t=1,value_parser=clap::value_parser!(u16).range(1..=64))]
    pub parallel: u16,
    #[arg(long)]
    pub non_alt_norm: bool,
    #[arg(long)]
    pub high_memory: bool,
}

pub fn event_types(value: &str) -> Result<Vec<EventType>> {
    value
        .trim_matches(',')
        .split(',')
        .map(|name| {
            EventType::ALL
                .into_iter()
                .find(|k| k.as_str() == name)
                .ok_or_else(|| anyhow::anyhow!("unknown event type {name}"))
        })
        .collect()
}

fn conditions(value: &str) -> Result<Vec<String>> {
    let mut values: Vec<String> = value
        .trim_matches(',')
        .split(',')
        .map(str::to_owned)
        .collect();
    if values[0].to_lowercase().ends_with("txt") {
        let content = std::fs::read_to_string(&values[0])?;
        values = content
            .lines()
            .flat_map(|line| {
                line.split('#')
                    .next()
                    .unwrap()
                    .split_whitespace()
                    .map(str::to_owned)
            })
            .collect();
        ensure!(
            values.len() > 1,
            "SplAdder v3.1.1 cannot read a singleton condition list file"
        );
    }
    let suffix =
        regex::Regex::new(r"(\.[bB][aA][mM]|\.[hH][dD][fF]5)|\.[nN][pP][zZ]|\.[cC][rR][aA][mM]$")?;
    values
        .into_iter()
        .map(|name| {
            let name = suffix.replace_all(&name, "");
            Ok(Path::new(name.as_ref())
                .file_name()
                .context("empty condition sample")?
                .to_str()
                .context("sample label is not UTF-8")?
                .to_owned())
        })
        .collect()
}

fn columns(dataset: &hdf5::Dataset, indices: &[usize]) -> Result<Vec<Vec<f64>>> {
    let shape = dataset.shape();
    ensure!(
        shape.len() == 2 && indices.iter().all(|&i| i < shape[1]),
        "gene expression sample index mismatch"
    );
    let mut output = vec![vec![0.0; indices.len()]; shape[0]];
    let mut positions = BTreeMap::<usize, Vec<usize>>::new();
    for (position, &column) in indices.iter().enumerate() {
        positions.entry(column).or_default().push(position);
    }
    let selected: Vec<_> = positions.keys().copied().collect();
    let mut begin = 0;
    while begin < selected.len() {
        let mut end = begin + 1;
        while end < selected.len() && selected[end] == selected[end - 1] + 1 {
            end += 1;
        }
        let start = selected[begin];
        let values = dataset.read_slice_2d::<f64, _>((.., start..selected[end - 1] + 1))?;
        for &column in &selected[begin..end] {
            for &position in &positions[&column] {
                for (row, out) in output.iter_mut().enumerate() {
                    out[position] = values[[row, column - start]];
                }
            }
        }
        begin = end;
    }
    Ok(output)
}

pub fn run(options: &TestArgs) -> Result<()> {
    ensure!(
        !options.diagnose_plots,
        "visualization is outside ruspladder's current scope"
    );
    let kinds = event_types(&options.event_types)?;
    let correction = serde_json::from_value(serde_json::Value::String(options.correction.clone()))
        .context("unknown multiple-testing correction")?;
    let a = conditions(&options.condition_a)?;
    let b = conditions(&options.condition_b)?;
    let graph = options.outdir.join(format!(
        "spladder/genes_graph_conf{}.merge_graphs.hdf5",
        options.confidence
    ));
    ensure!(
        graph.exists(),
        "testing requires a merge_graphs native graph cache: {}",
        graph.display()
    );
    let tag = if options.non_alt_norm { ".non_alt" } else { "" };
    let mut name = format!("testing{tag}");
    // The source compares the boolean timestamp option with the string 'y',
    // so its public CLI never appends a timestamp.
    if options.label_a != "condA" && options.label_b != "condB" {
        name.push_str(&format!("_{}_vs_{}", options.label_a, options.label_b));
    }
    if options.out_tag != "-" {
        name.push_str(&format!("_{}", options.out_tag));
    }
    let output = options.outdir.join(name);
    std::fs::create_dir_all(&output)?;
    let validated = if options.validate_sg {
        ".validated"
    } else {
        ""
    };
    let expression_path = options.outdir.join(format!(
        "spladder/genes_graph_conf{}.{}{}.gene_exp.hdf5",
        options.confidence, options.merge, validated
    ));
    let file = hdf5::File::open(&expression_path)?;
    let samples = hdf5io::read_strings(&file.dataset("samples")?)?;
    let mut indices = Vec::new();
    for requested in a.iter().chain(&b) {
        let matched: Vec<_> = samples
            .iter()
            .enumerate()
            .filter_map(|(i, s)| (s == requested).then_some(i))
            .collect();
        ensure!(
            !matched.is_empty(),
            "condition sample {requested} is missing from gene expression"
        );
        indices.extend(matched);
    }
    let metadata = Metadata {
        samples: indices
            .iter()
            .map(|&i| {
                samples[i]
                    .split(':')
                    .nth(1)
                    .unwrap_or(&samples[i])
                    .to_owned()
            })
            .collect(),
        gene_ids: hdf5io::read_strings(&file.dataset("gene_ids")?)?,
        gene_symbols: hdf5io::read_strings(&file.dataset("gene_symbols")?)?,
    };
    let counts = columns(
        &file.dataset(if options.non_alt_norm {
            "raw_count_non_alt"
        } else {
            "raw_count"
        })?,
        &indices,
    )?;
    file.close()?;
    let expression = testing::prepare_expression(
        counts,
        indices.len(),
        options.cap_exp_outliers && !options.no_cap_exp_outliers,
    )?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(options.parallel as usize)
        .build()?;
    for kind in kinds {
        let event_tag = format!("merge_graphs_{}_C{}", kind.as_str(), options.confidence);
        let count_path = options.outdir.join(format!("{event_tag}.counts.hdf5"));
        let counted = hdf5::File::open(&count_path)?;
        if !counted.link_exists("conf_idx") {
            if options.verbose {
                eprintln!("Skipping {}: no events", kind.as_str());
            }
            continue;
        }
        counted.close()?;
        let quantified = quantify::from_counted_events(
            &count_path,
            &indices[..a.len()],
            &indices[a.len()..],
            kind,
        )?;
        ensure!(
            quantified.samples == metadata.samples,
            "gene/event sample labels differ"
        );
        if quantified.event_idx.is_empty() {
            continue;
        }
        let events = cache::read_events(&options.outdir.join(format!("{event_tag}.events.hdf5")))?;
        let Some(prepared) = testing::prepare_events(
            quantified,
            &expression,
            a.len(),
            InputOptions {
                cap_outliers: options.cap_outliers,
                max_zero_fraction: options.max_zero_fraction,
                min_dpsi: options.min_dpsi,
            },
        )?
        else {
            continue;
        };
        if options.verbose {
            eprintln!(
                "Testing {} {} events",
                prepared.event_idx.len(),
                kind.as_str()
            );
        }
        let result = pool.install(|| {
            statistics::run(
                &prepared.counts,
                &prepared.null,
                &prepared.alternative,
                &prepared.size_factors,
                &prepared.selected,
                TestOptions {
                    min_count: options.min_count,
                    max_zero_fraction: options.max_zero_fraction,
                },
            )
        })?;
        test_output::write_results(
            &output,
            &test_output::Options {
                confidence: options.confidence,
                kind,
                correction,
                group_a: a.len(),
                labels: [options.label_a.clone(), options.label_b.clone()],
            },
            &metadata,
            &events,
            &prepared,
            &result,
        )?;
    }
    Ok(())
}
