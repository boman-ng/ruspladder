//! Bounded-memory writer for SplAdder's public COO alignment summaries.
use crate::reads::{AlignmentReader, ReadOptions};
use anyhow::{Result, ensure};
use hdf5::{Dataset, File, H5Type};
use rayon::prelude::*;
use std::{collections::BTreeMap, path::Path};

fn append<T: H5Type>(dataset: &Dataset, values: &[T]) -> Result<()> {
    if !values.is_empty() {
        let offset = dataset.size();
        dataset.resize(offset + values.len())?;
        dataset.write_slice(values, offset..offset + values.len())?;
    }
    Ok(())
}

fn vector<T: H5Type>(file: &File, name: &str, compressed: bool) -> Result<Dataset> {
    let builder = file.new_dataset::<T>().shape(0..).chunk(65536);
    Ok(if compressed {
        builder.deflate(4).create(name)?
    } else {
        builder.create(name)?
    })
}

pub struct SummaryOptions {
    pub parallel: usize,
    pub window: usize,
    pub unstranded: bool,
}

pub fn write_summary(
    bam: &Path,
    output: &Path,
    chromosomes: &[String],
    reference: Option<&Path>,
    options: &ReadOptions,
    summary: SummaryOptions,
) -> Result<()> {
    let SummaryOptions {
        parallel,
        window,
        unstranded,
    } = summary;
    ensure!(
        parallel > 0 && window > 0,
        "threads and window must be positive"
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(parallel)
        .build()?;
    let mut readers: Vec<_> = (0..parallel)
        .map(|_| AlignmentReader::open(bam, reference))
        .collect::<Result<_>>()?;
    let lengths: BTreeMap<_, _> = readers[0].contigs().into_iter().collect();
    let file = File::create(output)?;
    for chromosome in chromosomes {
        let rows = vector::<u8>(&file, &format!("{chromosome}_reads_row"), true)?;
        let columns = vector::<i32>(&file, &format!("{chromosome}_reads_col"), true)?;
        let mut introns: [BTreeMap<[i64; 2], u32>; 2] = Default::default();
        let shape = if let Some(&length) = lengths.get(chromosome) {
            ensure!(
                length <= i32::MAX as u64,
                "contig exceeds the reference COO int32 coordinate range"
            );
            let data = vector::<u32>(&file, &format!("{chromosome}_reads_dat"), true)?;
            let nrows = if unstranded { 1 } else { 3 };
            for row in 0..nrows {
                for start in (0..length as usize).step_by(window * parallel) {
                    let results: Vec<_> = pool.install(|| {
                        readers
                            .par_iter_mut()
                            .enumerate()
                            .map(|(i, reader)| {
                                let a = (start + i * window).min(length as usize);
                                let b = (a + window).min(length as usize);
                                reader.coverage_window(
                                    chromosome,
                                    a as i64,
                                    b as i64,
                                    options,
                                    (!unstranded).then_some(row),
                                    row == 0,
                                )
                            })
                            .collect::<Result<_>>()
                    })?;
                    for chunk in results {
                        append(&rows, &vec![row; chunk.values.len()])?;
                        append(&columns, &chunk.columns)?;
                        append(&data, &chunk.values)?;
                        for (all, values) in introns.iter_mut().zip(chunk.introns) {
                            for (pair, n) in values {
                                let count = all.entry(pair).or_default();
                                *count = count.wrapping_add(n);
                            }
                        }
                    }
                }
            }
            [i64::from(nrows), length as i64]
        } else {
            // The original missing-contig COO starts from np.zeros((0, 1)),
            // so its empty data array is float64, unlike known contigs.
            vector::<f64>(&file, &format!("{chromosome}_reads_dat"), true)?;
            [0, 1]
        };
        let shape_dataset = file.new_dataset::<i64>().shape(2);
        // The multiprocessing prep collector compresses every dataset,
        // including the shape vector; serial prep leaves this one uncompressed.
        let shape_dataset = if parallel > 1 {
            shape_dataset.deflate(4)
        } else {
            shape_dataset
        };
        shape_dataset
            .create(format!("{chromosome}_reads_shp").as_str())?
            .write_raw(&shape)?;
        for (suffix, values) in ["p", "m"].into_iter().zip(introns) {
            let flat: Vec<u32> = values
                .into_iter()
                .flat_map(|([a, b], n)| [a as u32, b as u32, n])
                .collect();
            file.new_dataset::<u32>()
                .shape((flat.len() / 3, 3))
                .deflate(4)
                .create(format!("{chromosome}_introns_{suffix}").as_str())?
                .write_raw(&flat)?;
        }
    }
    file.close()?;
    Ok(())
}

/// Query upstream's sorted, row-major COO summaries without loading an entire
/// chromosome. HDF5 hyperslabs keep memory proportional to the requested gene.
pub struct SparseReader {
    file: File,
    contig: Option<(String, SparseContig)>,
}

struct SparseContig {
    length: usize,
    rows: Vec<std::ops::Range<usize>>,
    columns: Dataset,
    values: Dataset,
    introns: [Dataset; 2],
}

fn partition_point(
    mut start: usize,
    mut stop: usize,
    mut before: impl FnMut(usize) -> Result<bool>,
) -> Result<usize> {
    while start < stop {
        let middle = start + (stop - start) / 2;
        if before(middle)? {
            start = middle + 1;
        } else {
            stop = middle;
        }
    }
    Ok(start)
}

impl SparseReader {
    pub fn open(path: &Path) -> Result<Self> {
        Ok(Self {
            file: File::open(path)?,
            contig: None,
        })
    }

    pub fn contig_names(&self) -> Result<Vec<String>> {
        Ok(self
            .file
            .member_names()?
            .into_iter()
            .filter_map(|name| name.strip_suffix("_reads_shp").map(str::to_owned))
            .collect())
    }

    fn load_contig(&mut self, chromosome: &str) -> Result<&SparseContig> {
        if self
            .contig
            .as_ref()
            .is_none_or(|(name, _)| name != chromosome)
        {
            let shape = self
                .file
                .dataset(&format!("{chromosome}_reads_shp"))?
                .read_raw::<u64>()?;
            ensure!(shape.len() == 2, "invalid sparse coverage shape");
            let row_ids = self.file.dataset(&format!("{chromosome}_reads_row"))?;
            let mut rows = Vec::new();
            let mut start = 0;
            for row in 0..shape[0] {
                let stop = partition_point(start, row_ids.size(), |i| {
                    Ok(row_ids.read_slice_1d::<u64, _>(i..i + 1)?[0] <= row)
                })?;
                rows.push(start..stop);
                start = stop;
            }
            self.contig = Some((
                chromosome.into(),
                SparseContig {
                    length: shape[1] as usize,
                    rows,
                    columns: self.file.dataset(&format!("{chromosome}_reads_col"))?,
                    values: self.file.dataset(&format!("{chromosome}_reads_dat"))?,
                    introns: [
                        self.file.dataset(&format!("{chromosome}_introns_p"))?,
                        self.file.dataset(&format!("{chromosome}_introns_m"))?,
                    ],
                },
            ));
        }
        Ok(&self.contig.as_ref().unwrap().1)
    }

    pub fn query(
        &mut self,
        chromosome: &str,
        start: i64,
        stop: i64,
        options: &ReadOptions,
        coverage: bool,
        unstranded: bool,
    ) -> Result<crate::reads::ReadEvidence> {
        ensure!(start >= 0 && stop >= start, "invalid sparse query bounds");
        let contig = self.load_contig(chromosome)?;
        let mut result = crate::reads::ReadEvidence::default();
        if coverage {
            let len = if contig.rows.is_empty() {
                (stop - start) as usize
            } else {
                (stop as usize)
                    .min(contig.length)
                    .saturating_sub((start as usize).min(contig.length))
            };
            result.coverage = vec![0; len];
            for (row, range) in contig.rows.iter().enumerate() {
                if !unstranded
                    && contig.rows.len() > 1
                    && row != 0
                    && row != 1 + usize::from(options.strand == Some('-'))
                {
                    continue;
                }
                let first = partition_point(range.start, range.end, |i| {
                    Ok(contig.columns.read_slice_1d::<i64, _>(i..i + 1)?[0] < start)
                })?;
                let last = partition_point(first, range.end, |i| {
                    Ok(contig.columns.read_slice_1d::<i64, _>(i..i + 1)?[0] < stop)
                })?;
                if first == last {
                    continue;
                }
                let columns = contig.columns.read_slice_1d::<i64, _>(first..last)?;
                let values = contig.values.read_slice_1d::<u32, _>(first..last)?;
                for (&column, &value) in columns.iter().zip(&values) {
                    result.coverage[(column - start) as usize] += value as u64;
                }
            }
        }
        for (output, dataset) in [&mut result.introns_plus, &mut result.introns_minus]
            .into_iter()
            .zip(&contig.introns)
        {
            let n = dataset.shape()[0];
            let first = partition_point(0, n, |i| {
                Ok(dataset.read_slice_2d::<i64, _>((i..i + 1, 0..1))?[[0, 0]] <= start)
            })?;
            let last = partition_point(first, n, |i| {
                Ok(dataset.read_slice_2d::<i64, _>((i..i + 1, 0..1))?[[0, 0]] < stop)
            })?;
            if first == last {
                continue;
            }
            for row in dataset.read_slice_2d::<i64, _>((first..last, ..))?.rows() {
                if row[1] < stop
                    && options
                        .filter
                        .as_ref()
                        .is_none_or(|f| row[2] >= f.mincount as i64)
                {
                    output.push([row[0], row[1], row[2]]);
                }
            }
        }
        Ok(result)
    }
}

pub fn summary_path(path: &Path, confidence: u8, filtered: bool) -> Result<std::path::PathBuf> {
    let name = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("non-UTF8 sparse alignment filename"))?;
    if name.ends_with("hdf5") {
        return Ok(path.to_owned());
    }
    static SUFFIX: std::sync::LazyLock<regex::Regex> =
        std::sync::LazyLock::new(|| regex::Regex::new(r"(?i)\.bam|\.cram$").unwrap());
    let base = SUFFIX.replace_all(name, "");
    Ok(if filtered {
        format!("{base}.conf_{confidence}.filt.hdf5")
    } else {
        format!("{base}.hdf5")
    }
    .into())
}
