// Adapted from SplAdder v3.1.1 reads.py get_reads/filter_read and
// settings.py set_confidence_level (BSD-3-Clause). See licenses/SplAdder-BSD.txt.
use anyhow::{Context, Result, bail, ensure};
use rust_htslib::bam::{
    self, Read,
    record::{Aux, Cigar},
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReadFilter {
    pub intron: u64,
    pub exon_len: u64,
    pub mismatch: i64,
    pub mincount: u64,
}

impl ReadFilter {
    pub fn confidence(level: u8, readlen: u64) -> Result<Self> {
        ensure!(level <= 3, "confidence must be between 0 and 3");
        let (fraction, mismatch, mincount) = match level {
            0 => (0.10, 2.max((readlen as f64 * 0.03).floor() as i64), 1),
            1 => (0.15, 1.max((readlen as f64 * 0.02).floor() as i64), 2),
            2 => (0.20, 1.max((readlen as f64 * 0.01).floor() as i64), 2),
            _ => (0.25, 0, 2),
        };
        Ok(Self {
            intron: 350000,
            exon_len: (readlen as f64 * fraction).ceil() as u64,
            mismatch,
            mincount,
        })
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct ReadOptions {
    pub filter: Option<ReadFilter>,
    pub mapped: bool,
    pub spliced: bool,
    pub strand: Option<char>,
    pub primary_only: bool,
    pub var_aware: bool,
    pub no_mm: bool,
    pub mm_tag: String,
}

impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            filter: None,
            mapped: true,
            spliced: true,
            strand: None,
            primary_only: false,
            var_aware: false,
            no_mm: false,
            mm_tag: "NM".into(),
        }
    }
}

fn numeric_tag(record: &bam::Record, tag: &[u8]) -> Result<i64> {
    let value = record.aux(tag).with_context(|| {
        format!(
            "alignment {} lacks required {} tag; --ignore-mismatches disables mismatch filtering",
            String::from_utf8_lossy(record.qname()),
            String::from_utf8_lossy(tag)
        )
    })?;
    match value {
        Aux::I8(x) => Ok(x as i64),
        Aux::U8(x) => Ok(x as i64),
        Aux::I16(x) => Ok(x as i64),
        Aux::U16(x) => Ok(x as i64),
        Aux::I32(x) => Ok(x as i64),
        Aux::U32(x) => Ok(x as i64),
        _ => bail!("{} tag must be integral", String::from_utf8_lossy(tag)),
    }
}

fn strand_tag(record: &bam::Record) -> Option<char> {
    match record.aux(b"XS") {
        Ok(Aux::Char(value)) => Some(value as char),
        Ok(Aux::String(value)) if value.len() == 1 => value.chars().next(),
        // A present but non-strand XS tag is not treated as a missing tag.
        Ok(_) => Some('\0'),
        Err(_) => None,
    }
}

pub fn filter_read(record: &bam::Record, options: &ReadOptions) -> Result<bool> {
    if record.is_unmapped() || (options.primary_only && record.is_secondary()) {
        return Ok(true);
    }
    // Supplementary alignments are intentionally retained, as in upstream.
    let cigar = record.cigar();
    let is_spliced = cigar.iter().any(|op| matches!(op, Cigar::RefSkip(_)));
    if (is_spliced && !options.spliced) || (!is_spliced && !options.mapped) {
        return Ok(true);
    }
    if let Some(filter) = &options.filter {
        if options.var_aware {
            if numeric_tag(record, b"XM")? + numeric_tag(record, b"XG")? > filter.mismatch {
                return Ok(true);
            }
        } else if !options.no_mm
            && numeric_tag(record, options.mm_tag.as_bytes())? > filter.mismatch
        {
            return Ok(true);
        }
        if is_spliced {
            let mut segment = 0u64;
            for op in cigar.iter() {
                match op {
                    Cigar::Match(n) | Cigar::Del(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                        segment += *n as u64
                    }
                    Cigar::RefSkip(_) => {
                        if segment <= filter.exon_len {
                            return Ok(true);
                        }
                        segment = 0;
                    }
                    Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => {}
                    Cigar::Pad(_) => {
                        bail!("padded spliced CIGAR is unsupported by the reference filter")
                    }
                }
            }
            if segment <= filter.exon_len {
                return Ok(true);
            }
        }
    }
    if let (Some(required), Some(actual)) = (options.strand, strand_tag(record))
        && required != actual
    {
        return Ok(true);
    }
    Ok(false)
}

#[derive(Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct ReadEvidence {
    pub coverage: Vec<u64>,
    /// Genomic donor, acceptor, support. Not clipped to the query interval.
    pub introns_plus: Vec<[i64; 3]>,
    pub introns_minus: Vec<[i64; 3]>,
    pub read_count: usize,
}

pub struct AlignmentReader {
    reader: bam::IndexedReader,
}

/// One bounded coverage window for the public sparse BAM representation.
pub struct CoverageWindow {
    pub columns: Vec<i32>,
    pub values: Vec<u32>,
    pub introns: [BTreeMap<[i64; 2], u32>; 2],
}

impl AlignmentReader {
    pub fn contigs(&self) -> Vec<(String, u64)> {
        self.reader
            .header()
            .target_names()
            .iter()
            .enumerate()
            .map(|(i, name)| {
                (
                    String::from_utf8_lossy(name).into_owned(),
                    self.reader.header().target_len(i as u32).unwrap(),
                )
            })
            .collect()
    }

    /// Start/end difference accumulation, as used by mosdepth, restricted to a
    /// bounded indexed window. SplAdder includes deletions and overlapping mates.
    /// target_row=None collapses strands; Some(0/1/2) selects an XS row.
    pub fn coverage_window(
        &mut self,
        chromosome: &str,
        start: i64,
        stop: i64,
        options: &ReadOptions,
        target_row: Option<u8>,
        count_introns: bool,
    ) -> Result<CoverageWindow> {
        ensure!(
            start >= 0 && stop >= start && stop <= i32::MAX as i64,
            "unsupported sparse window bounds"
        );
        let mut result = CoverageWindow {
            columns: Vec::new(),
            values: Vec::new(),
            introns: Default::default(),
        };
        let Some(tid) = self.reader.header().tid(chromosome.as_bytes()) else {
            return Ok(result);
        };
        self.reader.fetch((tid, start, stop))?;
        let mut difference = Vec::<i64>::new();
        let mut record = bam::Record::new();
        while let Some(status) = self.reader.read(&mut record) {
            status?;
            if filter_read(&record, options)? {
                continue;
            }
            let strand = strand_tag(&record);
            let row = if strand == Some('-') {
                2
            } else if strand.is_some() {
                1
            } else {
                0
            };
            let owner = count_introns && record.pos() >= start && record.pos() < stop;
            let mut position = record.pos();
            for op in record.cigar().iter() {
                match op {
                    Cigar::RefSkip(n) => {
                        if owner {
                            let value = result.introns[usize::from(row == 2)]
                                .entry([position, position + *n as i64])
                                .or_default();
                            *value = value.wrapping_add(1);
                        }
                        position += *n as i64;
                    }
                    Cigar::Match(n) | Cigar::Del(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                        let end = position + *n as i64;
                        if target_row.is_none_or(|r| r == row) {
                            let a = position.max(start);
                            let b = end.min(stop);
                            if a < b {
                                if difference.is_empty() {
                                    difference.resize((stop - start) as usize + 1, 0);
                                }
                                difference[(a - start) as usize] += 1;
                                difference[(b - start) as usize] -= 1;
                            }
                        }
                        position = end;
                    }
                    _ => {}
                }
            }
        }
        if difference.is_empty() {
            return Ok(result);
        }
        let mut depth = 0i64;
        for (i, &delta) in difference[..difference.len() - 1].iter().enumerate() {
            depth += delta;
            let value = depth as u32;
            if value != 0 {
                result.columns.push((start + i as i64) as i32);
                result.values.push(value);
            }
        }
        Ok(result)
    }

    pub fn contig_names(&self) -> Vec<String> {
        self.reader
            .header()
            .target_names()
            .iter()
            .map(|name| String::from_utf8_lossy(name).into_owned())
            .collect()
    }
    pub fn open(path: &Path, reference: Option<&Path>) -> Result<Self> {
        let mut reader = bam::IndexedReader::from_path(path)
            .with_context(|| format!("open indexed alignment {}", path.display()))?;
        if let Some(reference) = reference {
            reader.set_reference(reference)?;
        }
        Ok(Self { reader })
    }

    /// Equivalent to get_reads(..., collapse=True). The handle is reusable
    /// across neighboring regions, avoiding repeated file/index opening.
    pub fn region(
        &mut self,
        chromosome: &str,
        start: i64,
        stop: i64,
        options: &ReadOptions,
    ) -> Result<ReadEvidence> {
        self.evidence(chromosome, start, stop, options, true)
    }

    pub fn junctions(
        &mut self,
        chromosome: &str,
        start: i64,
        stop: i64,
        options: &ReadOptions,
    ) -> Result<ReadEvidence> {
        self.evidence(chromosome, start, stop, options, false)
    }

    fn evidence(
        &mut self,
        chromosome: &str,
        start: i64,
        stop: i64,
        options: &ReadOptions,
        coverage: bool,
    ) -> Result<ReadEvidence> {
        ensure!(
            start >= 0 && stop >= start,
            "invalid region {chromosome}:{start}-{stop}"
        );
        ensure!(
            options.mm_tag.len() == 2,
            "mismatch tag must have two characters"
        );
        let mut result = ReadEvidence {
            coverage: if coverage {
                vec![0; (stop - start) as usize]
            } else {
                Vec::new()
            },
            ..Default::default()
        };
        // get_reads explicitly excludes MT (summarize_chr does not).
        if chromosome == "MT" || stop == start {
            return Ok(result);
        }
        let Some(tid) = self.reader.header().tid(chromosome.as_bytes()) else {
            return Ok(result);
        };
        self.reader.fetch((tid, start, stop))?;
        let mut plus = BTreeMap::<[i64; 2], i64>::new();
        let mut minus = BTreeMap::<[i64; 2], i64>::new();
        let mut record = bam::Record::new();
        while let Some(status) = self.reader.read(&mut record) {
            status.context("decode alignment record")?;
            if filter_read(&record, options)? {
                continue;
            }
            let introns = if strand_tag(&record) == Some('-') {
                &mut minus
            } else {
                &mut plus
            };
            let mut position = record.pos();
            for op in record.cigar().iter() {
                match op {
                    Cigar::RefSkip(n) => {
                        *introns.entry([position, position + *n as i64]).or_default() += 1;
                        position += *n as i64;
                    }
                    Cigar::Match(n) | Cigar::Del(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                        let end = position + *n as i64;
                        let a = position.max(start);
                        let b = end.min(stop);
                        if coverage && a < b {
                            for value in
                                &mut result.coverage[(a - start) as usize..(b - start) as usize]
                            {
                                *value += 1;
                            }
                        }
                        position = end;
                    }
                    _ => {}
                }
            }
            result.read_count += 1;
        }
        result.introns_plus = plus.into_iter().map(|([a, b], n)| [a, b, n]).collect();
        result.introns_minus = minus.into_iter().map(|([a, b], n)| [a, b, n]).collect();
        Ok(result)
    }
}
