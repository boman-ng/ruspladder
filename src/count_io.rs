// Adapted from SplAdder v3.1.1 count.count_graph_coverage_wrapper (BSD-3-Clause).
use crate::{annotation::Gene, count::Counts, hdf5io};
use anyhow::{Result, ensure};
use hdf5::{Dataset, File};
use ndarray::ArrayView2;
use std::path::Path;

/// Match collect_single_quantification_results, including first-file dtypes,
/// extendable sample axes, S255 sample labels, and the absence of a strains link.
pub fn collect(paths: &[impl AsRef<Path>], output: &Path) -> Result<()> {
    ensure!(!paths.is_empty(), "no count files to collect");
    let out = File::create(output)?;
    let mut samples = Vec::new();
    let mut columns = 0;
    for (i, path) in paths.iter().enumerate() {
        let input = File::open(path)?;
        samples.extend(hdf5io::read_strings(&input.dataset("samples")?)?);
        if i == 0 {
            for key in [
                "gene_ids_segs",
                "edge_idx",
                "gene_ids_edges",
                "gene_names",
                "seg_len",
            ] {
                let source = input.dataset(key)?;
                let desc = source.dtype()?.to_descriptor()?;
                let dataset = out
                    .new_dataset_builder()
                    .empty_as(&desc)
                    .shape(source.shape())
                    .deflate(4)
                    .chunk_min_kb(64)
                    .create(key)?;
                match desc {
                    hdf5::types::TypeDescriptor::FixedAscii(_) => hdf5io::write_strings(
                        &dataset,
                        &hdf5io::read_strings(&source)?
                            .iter()
                            .map(String::as_str)
                            .collect::<Vec<_>>(),
                    )?,
                    hdf5::types::TypeDescriptor::Float(_) => {
                        dataset.write_raw(&source.read_raw::<f64>()?)?
                    }
                    hdf5::types::TypeDescriptor::Integer(_) => {
                        dataset.write_raw(&source.read_raw::<i64>()?)?
                    }
                    _ => anyhow::bail!("unsupported graph count metadata type"),
                }
            }
        }
        let added = input.dataset("segments")?.shape()[1];
        for key in ["edges", "segments", "seg_pos"] {
            let source = input.dataset(key)?;
            let shape = source.shape();
            ensure!(
                shape.len() == 2 && shape[1] == added,
                "count collection sample shape mismatch"
            );
            let target = if i == 0 {
                out.new_dataset_builder()
                    .empty_as(&source.dtype()?.to_descriptor()?)
                    .shape((shape[0], 0..))
                    .deflate(4)
                    .chunk((shape[0].clamp(1, 4096), 1))
                    .create(key)?
            } else {
                out.dataset(key)?
            };
            ensure!(
                shape[0] == target.shape()[0],
                "count collection row shape mismatch"
            );
            target.resize((shape[0], columns + added))?;
            for row in (0..shape[0]).step_by(65536) {
                let end = (row + 65536).min(shape[0]);
                let values = source.read_slice_2d::<f64, _>((row..end, ..))?;
                target.write_slice(values.view(), (row..end, columns..columns + added))?;
            }
        }
        columns += added;
    }
    let dataset = out
        .new_dataset_builder()
        .empty_as(&hdf5::types::TypeDescriptor::FixedAscii(255))
        .shape(0..)
        .deflate(4)
        .chunk((256,))
        .create("samples")?;
    dataset.resize(samples.len())?;
    hdf5io::write_strings(
        &dataset,
        &samples.iter().map(String::as_str).collect::<Vec<_>>(),
    )?;
    out.close()?;
    Ok(())
}

pub struct CountReader {
    file: File,
    segment_rows: Vec<Vec<usize>>,
    edge_rows: Vec<Vec<usize>>,
    edge_ids: Vec<u64>,
    samples: usize,
}

impl CountReader {
    pub fn open(path: &Path, genes: usize) -> Result<Self> {
        let file = File::open(path)?;
        let rows = |name| -> Result<Vec<Vec<usize>>> {
            let mut result = vec![Vec::new(); genes];
            for (row, id) in file
                .dataset(name)?
                .read_raw::<u64>()?
                .into_iter()
                .enumerate()
            {
                ensure!(id < genes as u64, "count gene index out of range");
                result[id as usize].push(row);
            }
            Ok(result)
        };
        Ok(Self {
            segment_rows: rows("gene_ids_segs")?,
            edge_rows: rows("gene_ids_edges")?,
            edge_ids: file.dataset("edge_idx")?.read_raw::<u64>()?,
            samples: file.dataset("segments")?.shape()[1],
            file,
        })
    }

    pub fn gene(&self, gene: usize, samples: &[usize]) -> Result<Vec<Counts>> {
        ensure!(
            gene < self.segment_rows.len() && samples.iter().all(|&s| s < self.samples),
            "count selection out of range"
        );
        let seg = &self.segment_rows[gene];
        let edge = &self.edge_rows[gene];
        let read_rows = |name: &str, rows: &[usize]| -> Result<Vec<Vec<f64>>> {
            let mut result = vec![Vec::with_capacity(rows.len()); samples.len()];
            let dataset = self.file.dataset(name)?;
            // Coalesce adjacent rows; reordered public count files remain readable.
            let mut offset = 0;
            while offset < rows.len() {
                let mut stop = offset + 1;
                while stop < rows.len() && rows[stop] == rows[stop - 1] + 1 {
                    stop += 1;
                }
                let block =
                    dataset.read_slice_2d::<f64, _>((rows[offset]..rows[stop - 1] + 1, ..))?;
                for (out, &sample) in result.iter_mut().zip(samples) {
                    out.extend(block.column(sample).iter().copied());
                }
                offset = stop;
            }
            Ok(result)
        };
        let segments = read_rows("segments", seg)?;
        let positions = read_rows("seg_pos", seg)?;
        let edges = read_rows("edges", edge)?;
        Ok(segments
            .into_iter()
            .zip(positions)
            .zip(edges)
            .map(|((segments, seg_pos), counts)| Counts {
                segments,
                seg_pos,
                edges: edge
                    .iter()
                    .zip(counts)
                    .map(|(&row, n)| [self.edge_ids[row], n as u64])
                    .collect(),
                edges_float: true,
            })
            .collect())
    }
}

/// Write one sample at a time. Only a single sample's count vectors are live.
pub struct CountWriter {
    file: File,
    segments: Dataset,
    seg_pos: Dataset,
    edges: Dataset,
    edge_ids: Vec<i64>,
    segment_sizes: Vec<usize>,
    edge_sizes: Vec<usize>,
    samples: usize,
    written: usize,
    edges_float: bool,
}

impl CountWriter {
    pub fn create(
        path: &Path,
        genes: &[Gene],
        labels: &[&str],
        sample_count: usize,
    ) -> Result<Self> {
        ensure!(
            !genes.is_empty() && sample_count > 0,
            "graph counting needs genes and samples"
        );
        let file = File::create(path)?;
        hdf5io::strings(&file, "samples", &[labels.len()], labels, false)?;
        file.link_soft("/samples", "strains")?;
        let names: Vec<_> = genes.iter().map(|g| g.name.as_str()).collect();
        hdf5io::strings(&file, "gene_names", &[genes.len(), 1], &names, false)?;
        let mut gene_segments = Vec::new();
        let mut gene_edges = Vec::new();
        let mut lengths = Vec::new();
        let mut edge_ids = Vec::new();
        let mut segment_sizes = Vec::new();
        let mut edge_sizes = Vec::new();
        for (i, gene) in genes.iter().enumerate() {
            let graph = &gene.segmentgraph;
            ensure!(
                !graph.segments.is_empty(),
                "count writer requires segment graphs"
            );
            segment_sizes.push(graph.segments.len());
            edge_sizes.push(graph.edges.len());
            gene_segments.extend(std::iter::repeat_n(i as i64, graph.segments.len()));
            gene_edges.extend(std::iter::repeat_n(i as i64, graph.edges.len()));
            lengths.extend(graph.segments.iter().map(|v| v[1] - v[0]));
            edge_ids.extend(
                graph
                    .edges
                    .iter()
                    .map(|&[a, b]| (a * graph.segments.len() + b) as i64),
            );
        }
        ensure!(
            !edge_ids.is_empty(),
            "SplAdder v3.1.1 cannot write count graphs with no edges"
        );
        hdf5io::array(
            &file,
            "gene_ids_segs",
            &[gene_segments.len(), 1],
            &gene_segments,
            false,
        )?;
        hdf5io::array(
            &file,
            "gene_ids_edges",
            &[gene_edges.len(), 1],
            &gene_edges,
            false,
        )?;
        hdf5io::array(&file, "seg_len", &[lengths.len(), 1], &lengths, false)?;
        let segments = file
            .new_dataset::<f64>()
            .shape((lengths.len(), sample_count))
            .create("segments")?;
        let seg_pos = file
            .new_dataset::<f64>()
            .shape((lengths.len(), sample_count))
            .create("seg_pos")?;
        let edges = file
            .new_dataset::<f64>()
            .shape((edge_ids.len(), sample_count))
            .create("edges")?;
        Ok(Self {
            file,
            segments,
            seg_pos,
            edges,
            edge_ids,
            segment_sizes,
            edge_sizes,
            samples: sample_count,
            written: 0,
            edges_float: false,
        })
    }

    pub fn append(&mut self, counts: &[Counts]) -> Result<()> {
        ensure!(
            counts.len() == self.segment_sizes.len() && self.written < self.samples,
            "count sample shape mismatch"
        );
        let mut segments = Vec::with_capacity(self.segments.shape()[0]);
        let mut positions = Vec::with_capacity(segments.capacity());
        let mut edges = Vec::with_capacity(self.edge_ids.len());
        for (i, count) in counts.iter().enumerate() {
            ensure!(
                count.segments.len() == self.segment_sizes[i]
                    && count.seg_pos.len() == self.segment_sizes[i]
                    && count.edges.len() == self.edge_sizes[i],
                "gene count shape mismatch"
            );
            self.edges_float |= !count.edges.is_empty() && count.edges_float;
            segments.extend_from_slice(&count.segments);
            positions.extend_from_slice(&count.seg_pos);
            for &[id, value] in &count.edges {
                ensure!(
                    id == self.edge_ids[edges.len()] as u64,
                    "count edge order mismatch"
                );
                edges.push(value as f64);
            }
        }
        for (dataset, values) in [
            (&self.segments, &segments),
            (&self.seg_pos, &positions),
            (&self.edges, &edges),
        ] {
            dataset.write_slice(
                ArrayView2::from_shape((values.len(), 1), values)?,
                (.., self.written..self.written + 1),
            )?;
        }
        self.written += 1;
        Ok(())
    }

    pub fn finish(self) -> Result<()> {
        ensure!(
            self.written == self.samples,
            "not all count samples have been written"
        );
        if self.edges_float {
            hdf5io::array(
                &self.file,
                "edge_idx",
                &[self.edge_ids.len()],
                &self.edge_ids.iter().map(|&i| i as f64).collect::<Vec<_>>(),
                false,
            )?;
        } else {
            // NumPy retains int64 only if every edge count is zero. HDF5's
            // default zero fill creates that matrix without allocating it.
            self.file.unlink("edges")?;
            self.file
                .new_dataset::<i64>()
                .shape((self.edge_ids.len(), self.samples))
                .create("edges")?;
            hdf5io::array(
                &self.file,
                "edge_idx",
                &[self.edge_ids.len()],
                &self.edge_ids,
                false,
            )?;
        }
        self.file.close()?;
        Ok(())
    }
}
