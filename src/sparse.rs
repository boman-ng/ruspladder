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
        file.new_dataset::<i64>()
            .shape(2)
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
    file.flush()?;
    Ok(())
}
