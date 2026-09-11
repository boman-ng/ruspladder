//! Rust-owned annotation and gene graph HDF5 cache. Dataset names are distinct
//! from SplAdder's public result files; Python pickle is not a cache interface.
use crate::annotation::Gene;
use crate::graph::{Interval, SegmentGraph, SpliceGraph};
use anyhow::{Context, Result, ensure};
use hdf5::{File, Group, H5Type, types::VarLenUnicode};
use std::collections::BTreeSet;
use std::path::Path;
use std::str::FromStr;

fn write_vector<T: H5Type>(group: &Group, name: &str, values: &[T]) -> Result<()> {
    group
        .new_dataset::<T>()
        .shape(values.len())
        .create(name)?
        .write_raw(values)?;
    Ok(())
}

fn write_pairs<T: H5Type + Copy>(group: &Group, name: &str, values: &[[T; 2]]) -> Result<()> {
    let flat: Vec<_> = values.iter().flatten().copied().collect();
    group
        .new_dataset::<T>()
        .shape((values.len(), 2))
        .create(name)?
        .write_raw(&flat)?;
    Ok(())
}

fn read_pairs<T: H5Type + Copy>(group: &Group, name: &str) -> Result<Vec<[T; 2]>> {
    let dataset = group.dataset(name)?;
    let shape = dataset.shape();
    ensure!(
        shape.len() == 2 && shape[1] == 2,
        "{name}: expected a two-column dataset"
    );
    Ok(dataset
        .read_raw::<T>()?
        .chunks_exact(2)
        .map(|p| [p[0], p[1]])
        .collect())
}

fn write_string(group: &Group, name: &str, value: &str) -> Result<()> {
    group
        .new_dataset::<VarLenUnicode>()
        .shape(())
        .create(name)?
        .write_scalar(&VarLenUnicode::from_str(value)?)?;
    Ok(())
}

fn read_string(group: &Group, name: &str) -> Result<String> {
    Ok(group
        .dataset(name)?
        .read_scalar::<VarLenUnicode>()?
        .as_str()
        .to_owned())
}

pub fn write_genes(path: &Path, genes: &[Gene]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("create graph cache {}", path.display()))?;
    write_string(&file, "format", "ruspladder-genes")?;
    file.new_dataset::<u64>()
        .shape(())
        .create("gene_count")?
        .write_scalar(&(genes.len() as u64))?;
    let root = file.create_group("genes")?;
    for (i, gene) in genes.iter().enumerate() {
        let group = root.create_group(&i.to_string())?;
        for (name, value) in [("name", gene.name.as_str()), ("chr", &gene.chr)] {
            write_string(&group, name, value)?;
        }
        write_string(&group, "strand", &gene.strand.to_string())?;
        for (name, value) in [
            ("gene_type", &gene.gene_type),
            ("symbol", &gene.symbol),
            ("source", &gene.source),
        ] {
            if let Some(value) = value {
                write_string(&group, name, value)?;
            }
        }
        write_vector(&group, "bounds", &[gene.start, gene.stop])?;
        write_vector(
            &group,
            "alt",
            &[
                gene.is_alt.map_or(-1, i8::from),
                gene.is_alt_spliced.map_or(-1, i8::from),
            ],
        )?;
        let transcripts: Vec<_> = gene
            .transcripts
            .iter()
            .map(|s| VarLenUnicode::from_str(s))
            .collect::<Result<_, _>>()?;
        write_vector(&group, "transcripts", &transcripts)?;
        let mut offsets = vec![0u64];
        for exons in &gene.exons {
            offsets.push(offsets.last().unwrap() + exons.len() as u64);
        }
        write_vector(&group, "transcript_offsets", &offsets)?;
        write_pairs(
            &group,
            "exons",
            &gene.exons.iter().flatten().copied().collect::<Vec<_>>(),
        )?;
        write_pairs(
            &group,
            "introns_anno",
            &gene.introns_anno.iter().copied().collect::<Vec<_>>(),
        )?;
        let splice = group.create_group("splice")?;
        write_pairs(&splice, "vertices", &gene.splicegraph.vertices)?;
        write_pairs(&splice, "terminals", &gene.splicegraph.terminals)?;
        let edges: Vec<_> = gene
            .splicegraph
            .edges
            .iter()
            .enumerate()
            .flat_map(|(i, neighbors)| {
                neighbors
                    .iter()
                    .filter(move |&&j| j >= i)
                    .map(move |&j| [i as u64, j as u64])
            })
            .collect();
        write_pairs(&splice, "edges", &edges)?;
        let segment = group.create_group("segment")?;
        write_pairs(&segment, "segments", &gene.segmentgraph.segments)?;
        write_pairs(
            &segment,
            "edges",
            &gene
                .segmentgraph
                .edges
                .iter()
                .map(|&[i, j]| [i as u64, j as u64])
                .collect::<Vec<_>>(),
        )?;
        let mut offsets = vec![0u64];
        for matches in &gene.segmentgraph.matches {
            offsets.push(offsets.last().unwrap() + matches.len() as u64);
        }
        write_vector(&segment, "match_offsets", &offsets)?;
        write_vector(
            &segment,
            "matches",
            &gene
                .segmentgraph
                .matches
                .iter()
                .flatten()
                .map(|&j| j as u64)
                .collect::<Vec<_>>(),
        )?;
    }
    file.close()?;
    Ok(())
}

pub fn read_genes(path: &Path) -> Result<Vec<Gene>> {
    let file = File::open(path).with_context(|| format!("open graph cache {}", path.display()))?;
    ensure!(
        read_string(&file, "format")? == "ruspladder-genes",
        "not a ruspladder graph cache"
    );
    let count = file.dataset("gene_count")?.read_scalar::<u64>()? as usize;
    let root = file.group("genes")?;
    let mut genes = Vec::with_capacity(count);
    for i in 0..count {
        let group = root.group(&i.to_string())?;
        let bounds = group.dataset("bounds")?.read_raw::<i64>()?;
        let alt = group.dataset("alt")?.read_raw::<i8>()?;
        ensure!(
            bounds.len() == 2 && alt.len() == 2,
            "gene {i}: malformed bounds/alt flags"
        );
        let exons: Vec<Interval> = read_pairs(&group, "exons")?;
        let offsets = group.dataset("transcript_offsets")?.read_raw::<u64>()?;
        let transcripts: Vec<_> = group
            .dataset("transcripts")?
            .read_raw::<VarLenUnicode>()?
            .into_iter()
            .map(|s| s.as_str().to_owned())
            .collect();
        ensure!(
            offsets.len() == transcripts.len() + 1,
            "gene {i}: transcript offsets disagree with names"
        );
        let exons = unpack(&exons, &offsets)?;
        let splice = group.group("splice")?;
        let vertices = read_pairs(&splice, "vertices")?;
        let terminals = read_pairs(&splice, "terminals")?;
        ensure!(
            vertices.len() == terminals.len(),
            "gene {i}: terminal/vertex lengths disagree"
        );
        let mut splicegraph = SpliceGraph {
            edges: vec![Vec::new(); vertices.len()],
            vertices,
            terminals,
        };
        for [a, b] in read_pairs::<u64>(&splice, "edges")? {
            ensure!(
                a < splicegraph.vertices.len() as u64 && b < splicegraph.vertices.len() as u64,
                "gene {i}: splice edge out of range"
            );
            splicegraph.connect(a as usize, b as usize);
        }
        let segment = group.group("segment")?;
        let segmentgraph = SegmentGraph {
            segments: read_pairs(&segment, "segments")?,
            edges: read_pairs::<u64>(&segment, "edges")?
                .into_iter()
                .map(|[i, j]| [i as usize, j as usize])
                .collect(),
            matches: unpack(
                &segment.dataset("matches")?.read_raw::<u64>()?,
                &segment.dataset("match_offsets")?.read_raw::<u64>()?,
            )?
            .into_iter()
            .map(|v| v.into_iter().map(|j| j as usize).collect())
            .collect(),
        };
        let optional = |name| -> Result<_> {
            if group.link_exists(name) {
                Ok(Some(read_string(&group, name)?))
            } else {
                Ok(None)
            }
        };
        let strand = read_string(&group, "strand")?;
        ensure!(
            strand == "+" || strand == "-",
            "gene {i}: invalid cached strand"
        );
        genes.push(Gene {
            name: read_string(&group, "name")?,
            chr: read_string(&group, "chr")?,
            source: optional("source")?,
            gene_type: optional("gene_type")?,
            symbol: optional("symbol")?,
            start: bounds[0],
            stop: bounds[1],
            strand: strand.chars().next().unwrap(),
            is_alt: (alt[0] >= 0).then_some(alt[0] != 0),
            is_alt_spliced: (alt[1] >= 0).then_some(alt[1] != 0),
            transcripts,
            exons,
            introns_anno: read_pairs::<i64>(&group, "introns_anno")?
                .into_iter()
                .collect::<BTreeSet<_>>(),
            splicegraph,
            segmentgraph,
        });
    }
    Ok(genes)
}

fn unpack<T: Clone>(values: &[T], offsets: &[u64]) -> Result<Vec<Vec<T>>> {
    ensure!(
        offsets.first() == Some(&0) && offsets.last() == Some(&(values.len() as u64)),
        "invalid ragged array offsets"
    );
    offsets
        .windows(2)
        .map(|pair| {
            ensure!(
                pair[0] <= pair[1] && pair[1] <= values.len() as u64,
                "invalid ragged array interval"
            );
            Ok(values[pair[0] as usize..pair[1] as usize].to_vec())
        })
        .collect()
}
