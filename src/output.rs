// Adapted from SplAdder v3.1.1 alt_splice/write.py (BSD-3-Clause).
use crate::{
    analyze,
    events::{Event, EventType},
    graph::Interval,
    hdf5io,
};
use anyhow::{Result, ensure};
use flate2::{Compression, write::GzEncoder};
use serde::Deserialize;
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Clone, Copy, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Format {
    Txt,
    Structured,
    Bed,
    Gff3,
    Gtf,
    Tcga,
    Icgc,
}

fn float(x: f64, digits: usize) -> String {
    if x.is_nan() {
        "nan".into()
    } else {
        format!("{x:.digits$}")
    }
}

fn join(values: impl IntoIterator<Item = impl ToString>, separator: &str) -> String {
    values
        .into_iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(separator)
}

/// The three disjoint regions used by the flat alt-prime/retention output.
fn regions(event: &Event) -> Vec<Interval> {
    let a = &event.exons1;
    let b = &event.exons2;
    match event.event_type {
        EventType::ExonSkip | EventType::MultExonSkip => b.clone(),
        EventType::IntronRetention => vec![a[0], [a[0][1], a[1][0]], a[1]],
        EventType::MutexExons => vec![a[0], a[1], b[1], b[2]],
        _ if a[0] == b[0] => {
            let (long, short) = if a[1][0] < b[1][0] {
                (a[1], b[1])
            } else {
                (b[1], a[1])
            };
            vec![a[0], [long[0], short[0]], short]
        }
        _ => {
            let (short, long) = if a[0][1] < b[0][1] {
                (a[0], b[0])
            } else {
                (b[0], a[0])
            };
            vec![short, [short[1], long[1]], b[1]]
        }
    }
}

fn catalog_regions(event: &Event) -> Vec<Interval> {
    let a = &event.exons1;
    let b = &event.exons2;
    match event.event_type {
        EventType::IntronRetention => a.clone(),
        EventType::ExonSkip | EventType::MultExonSkip => b.clone(),
        EventType::MutexExons => vec![a[0], a[1], b[1], a[2]],
        _ if a[0] == b[0] => vec![a[0], a[1], b[1]],
        _ => vec![a[1], a[0], b[0]],
    }
}

fn assert_shorter(event: &Event) -> Result<()> {
    let len = |exons: &[Interval]| exons.iter().map(|v| v[1] - v[0]).sum::<i64>();
    ensure!(
        len(&event.exons1) <= len(&event.exons2),
        "output requires the shorter isoform first"
    );
    Ok(())
}

fn gff(out: &mut dyn Write, event: &Event, gtf: bool) -> Result<()> {
    assert_shorter(event)?;
    let kind = event.event_type.as_str();
    let name = format!("{kind}.{}", event.id);
    let gene = &event.gene_name[0];
    let start = event.exons1[0][0] + 1;
    let end = event.exons1.last().unwrap()[1];
    let novel = |mask| {
        if event.annotated & mask == mask {
            "N"
        } else {
            "Y"
        }
    };
    let prefix = |feature: &str, start, end| {
        format!(
            "{}\t{kind}\t{feature}\t{start}\t{end}\t.\t{}\t.\t",
            event.chr, event.strand
        )
    };
    if gtf {
        writeln!(
            out,
            "{}gene_id \"{name}\"; transcript_id \"{name}\"; gene_name \"{gene}\"; has_novel_junction \"{}\";",
            prefix("gene", start, end),
            novel(3)
        )?;
    } else {
        writeln!(
            out,
            "{}ID={name};GeneName=\"{gene}\";HasNovelJunction=\"{}\"",
            prefix("gene", start, end),
            novel(3)
        )?;
    }
    for (i, exons) in [&event.exons1, &event.exons2].into_iter().enumerate() {
        let iso = i + 1;
        if gtf {
            let feature = if i == 0 { "transcript" } else { "mRNA" };
            writeln!(
                out,
                "{}gene_id \"{name}\"; transcript_id \"{name}_iso{iso}\"; gene_name \"{gene}\"; has_novel_junction \"{}\"",
                prefix(feature, start, end),
                novel(1 << i)
            )?;
        } else {
            writeln!(
                out,
                "{}ID={name}_iso{iso};Parent={name};GeneName=\"{gene}\";HasNovelJunction=\"{}\"",
                prefix("mRNA", start, end),
                novel(1 << i)
            )?;
        }
        for (j, &[a, b]) in exons.iter().enumerate() {
            if gtf {
                writeln!(
                    out,
                    "{}gene_id \"{name}\"; transcript_id \"{name}_iso{iso}\"; exon_id \"{name}_iso{iso}_exon{}\"{}",
                    prefix("exon", a + 1, b),
                    j + 1,
                    if i == 1 { ";" } else { "" }
                )?;
            } else {
                writeln!(out, "{}Parent={name}_iso{iso}", prefix("exon", a + 1, b))?;
            }
        }
    }
    Ok(())
}

fn bed(out: &mut dyn Write, event: &Event) -> Result<bool> {
    assert_shorter(event)?;
    let a = &event.exons1;
    let b = &event.exons2;
    let start = a[0][0];
    let end = a.last().unwrap()[1];
    let mut exons = regions(event);
    let tag = match event.event_type {
        EventType::ExonSkip => "CA",
        EventType::IntronRetention => "IR",
        EventType::Alt3prime => "AA",
        EventType::Alt5prime => "AD",
        _ => {
            eprintln!(
                "SplAdder does not implement {} BED output",
                event.event_type.as_str()
            );
            return Ok(false);
        }
    };
    let inner = exons[1];
    if event.event_type == EventType::Alt5prime && a[0] == b[0] {
        // Preserve the source's overlapping third block on this path.
        exons[2] = if a[1][0] > b[1][0] { b[1] } else { a[1] };
    }
    writeln!(
        out,
        "{}\t{start}\t{end}\t{tag}-{tag}-{}-{}.0[L]\t0\t{}\t{start}\t{end}\t255,0,0\t3\t{}\t{}",
        event.chr,
        inner[0],
        inner[1],
        event.strand,
        join(exons.iter().map(|v| v[1] - v[0]), ","),
        join(exons.iter().map(|v| v[0] - start), ",")
    )?;
    Ok(true)
}

fn txt_header(out: &mut dyn Write, kind: EventType, samples: &[&str]) -> Result<()> {
    write!(
        out,
        "chrm\tstrand\tevent_id\tis_annotated\tgene_name\te1_start\te1_end"
    )?;
    if kind == EventType::MultExonSkip {
        write!(out, "\te2_starts\te2_ends\te3_start\te3_end")?;
    } else {
        write!(out, "\te2_start\te2_end\te3_start\te3_end")?;
        if kind == EventType::MutexExons {
            write!(out, "\te4_start\te4_end")?;
        }
    }
    for sample in samples {
        let fields = if kind == EventType::MultExonSkip {
            &[
                "e1_cov",
                "e2_cov",
                "e3_cov",
                "e1e2_conf",
                "sum_e2_conf",
                "num_e2",
                "e2e3_conf",
                "e1e3_conf",
            ][..]
        } else {
            &analyze::features(kind)[1..]
        };
        for field in fields {
            write!(out, "\t{sample}:{field}")?;
        }
        write!(out, "\t{sample}:psi")?;
    }
    writeln!(out)?;
    Ok(())
}

fn txt(
    out: &mut dyn Write,
    event: &Event,
    counts: &ndarray::Array2<f64>,
    psi: &[f64],
    samples: usize,
) -> Result<()> {
    let kind = event.event_type;
    write!(
        out,
        "{}\t{}\t{}.{}\t{}\t{}",
        event.chr,
        event.strand,
        kind.as_str(),
        event.id,
        event.annotated,
        event.gene_name[0]
    )?;
    let exons = regions(event);
    if kind == EventType::MultExonSkip {
        write!(
            out,
            "\t{}\t{}\t{}\t{}\t{}\t{}",
            exons[0][0] + 1,
            exons[0][1],
            join(exons[1..exons.len() - 1].iter().map(|v| v[0] + 1), ":"),
            join(exons[1..exons.len() - 1].iter().map(|v| v[1]), ":"),
            exons.last().unwrap()[0] + 1,
            exons.last().unwrap()[1]
        )?;
    } else {
        for &[a, b] in &exons {
            write!(out, "\t{}\t{b}", a + 1)?;
        }
    }
    for (s, &p) in psi.iter().enumerate().take(samples) {
        let c = counts.row(s);
        let valid = if matches!(
            kind,
            EventType::Alt3prime | EventType::Alt5prime | EventType::MultExonSkip
        ) {
            c[0] != 0.0
        } else {
            c[0] == 1.0
        };
        if !valid {
            let missing = match kind {
                EventType::ExonSkip | EventType::Alt3prime | EventType::Alt5prime => 5,
                EventType::IntronRetention => 4,
                EventType::MultExonSkip => 8,
                EventType::MutexExons => 9,
            };
            for _ in 0..missing {
                write!(out, "\t-1")?;
            }
            continue;
        }
        let coverage = if kind == EventType::MutexExons { 4 } else { 3 };
        for i in 1..=coverage {
            write!(out, "\t{}", float(c[i], 1))?;
        }
        if kind == EventType::IntronRetention {
            write!(out, "\t{}\t{}", c[4] as i64, float(c[5], 2))?;
        } else {
            let fields = if kind == EventType::MultExonSkip {
                vec![4, 7, 8, 5, 6]
            } else {
                (coverage + 1..c.len()).collect()
            };
            for i in fields {
                write!(out, "\t{}", c[i] as i64)?;
            }
        }
        write!(out, "\t{}", float(p, 2))?;
    }
    writeln!(out)?;
    Ok(())
}

fn icgc(out: &mut dyn Write, event: &Event, psi: &[f64], samples: usize) -> Result<()> {
    if psi.iter().all(|p| p.is_nan()) {
        return Ok(());
    }
    let kind = event.event_type;
    let short = match kind {
        EventType::ExonSkip => "ES",
        EventType::IntronRetention => "IR",
        EventType::Alt3prime => "A3",
        EventType::Alt5prime => "A5",
        EventType::MultExonSkip => "MES",
        EventType::MutexExons => "MEX",
    };
    let exons = catalog_regions(event);
    write!(
        out,
        "{}.{}\t{short}\t{}\t{}\t{}\t",
        kind.as_str(),
        event.id,
        event.chr,
        event.annotated,
        join(exons.iter().flatten(), ":")
    )?;
    let a = &event.exons1;
    let b = &event.exons2;
    let alt = match kind {
        EventType::IntronRetention => vec![a[0][1], a[1][0]],
        EventType::ExonSkip => b[1].to_vec(),
        EventType::MultExonSkip => b[1..b.len() - 1].iter().flatten().copied().collect(),
        EventType::MutexExons => vec![a[1][0], a[1][1], b[1][0], b[1][1]],
        _ => regions(event)[1].to_vec(),
    };
    if kind == EventType::MultExonSkip {
        write!(out, ":")?;
    }
    write!(out, "{}\t{}", join(alt, ":"), event.gene_name[0])?;
    for &p in &psi[..samples] {
        write!(out, "\t{}", float(p, 6))?;
    }
    writeln!(out)?;
    Ok(())
}

fn tcga(
    out: &mut dyn Write,
    event: &Event,
    counts: &ndarray::Array2<f64>,
    samples: usize,
) -> Result<()> {
    let kind = event.event_type;
    write!(
        out,
        "{}\t{}\t{}:",
        event.gene_name[0],
        kind.as_str(),
        event.chr
    )?;
    let exons = catalog_regions(event);
    for (i, v) in exons.iter().enumerate() {
        write!(out, ":{}-{}", v[0], v[1])?;
        if kind == EventType::MultExonSkip && i > 0 && i + 1 < exons.len() {
            write!(out, ":")?;
        }
    }
    for s in 0..samples {
        let c = counts.row(s);
        if c[0] != 1.0 {
            write!(out, "\tNA")?;
            continue;
        }
        let (num, denom, confirmation) = match kind {
            EventType::Alt3prime | EventType::Alt5prime => {
                let a = &event.exons1;
                let b = &event.exons2;
                (
                    if a[1][0] - a[0][1] < b[1][0] - b[0][1] {
                        c[4]
                    } else {
                        c[5]
                    },
                    c[4] + c[5],
                    c[4] + c[5],
                )
            }
            EventType::ExonSkip => (c[4] + c[5], c[4] + c[5] + 2.0 * c[6], c[4] + c[5] + c[6]),
            EventType::MultExonSkip => (
                c[4] + c[7] + c[5],
                c[4] + c[7] + c[5] + (2.0 + c[8]) * c[6],
                c[4] + c[7] + c[5] + c[6],
            ),
            EventType::IntronRetention => (c[4], 1.0, c[4]),
            EventType::MutexExons => (
                c[5] + c[6],
                c[5] + c[6] + c[7] + c[8],
                c[5] + c[6] + c[7] + c[8],
            ),
        };
        if confirmation < 10.0 {
            write!(out, "\tNA")?;
        } else {
            write!(out, "\t{}", float(num / denom, 1))?;
        }
    }
    writeln!(out)?;
    Ok(())
}

fn bytes_repr(value: &str) -> String {
    let quote = if value.contains('\'') && !value.contains('"') {
        '"'
    } else {
        '\''
    };
    let mut result = format!("b{quote}");
    for byte in value.bytes() {
        match byte {
            b'\\' => result.push_str("\\\\"),
            b'\n' => result.push_str("\\n"),
            b'\r' => result.push_str("\\r"),
            b'\t' => result.push_str("\\t"),
            b if b == quote as u8 => {
                result.push('\\');
                result.push(quote);
            }
            32..=126 => result.push(byte as char),
            _ => result.push_str(&format!("\\x{byte:02x}")),
        }
    }
    result.push(quote);
    result
}

fn structured(out: &mut dyn Write, event: &Event, psi: &[f64], samples: &[String]) -> Result<bool> {
    let a = &event.exons1;
    let b = &event.exons2;
    let kind = event.event_type;
    let name = format!("{}.{}", kind.as_str(), event.id);
    write!(
        out,
        "{}\tundefined\t{}\t{}\t{}\t.\t{}\t.\tgene_id \"{name}\"; transcript_id \"{name}\"; gene_name \"{}\"; is_annotated \"{}\";",
        event.chr,
        kind.as_str(),
        a[0][0],
        a.last().unwrap()[1],
        event.strand,
        event.gene_name[0],
        event.annotated
    )?;
    if kind == EventType::MultExonSkip {
        eprintln!("SplAdder does not implement mult_exon_skip structured output");
        return Ok(false);
    }
    let positive = event.strand == '+';
    let (structure, flanks, chain) = match kind {
        EventType::ExonSkip | EventType::MutexExons => {
            let flanks = if positive {
                format!("{}^,{}-", a[0][1], a.last().unwrap()[0] + 1)
            } else {
                format!("{}^,{}-", a.last().unwrap()[0] + 1, a[0][1])
            };
            let exon = |v: Interval| {
                if positive {
                    format!("{}-{}^", v[0] + 1, v[1])
                } else {
                    format!("{}-{}^", v[1], v[0] + 1)
                }
            };
            if kind == EventType::ExonSkip {
                ("0,1-2^", flanks, format!(",{}", exon(b[1])))
            } else {
                (
                    "1-2^,3-4^",
                    flanks,
                    format!("{},{}", exon(a[1]), exon(b[1])),
                )
            }
        }
        EventType::IntronRetention => {
            if positive {
                (
                    "0,1^2-",
                    format!("{}-,{}^", b[0][1] + 1, b[0][1]),
                    format!(",{}^{}-", a[0][1], a[1][0] + 1),
                )
            } else {
                (
                    "0,1^2-",
                    format!("{}-,{}^", b[0][1], b[0][0] + 1),
                    format!(",{}^{}-", a[1][0] + 1, a[0][1]),
                )
            }
        }
        _ if a[0] == b[0] => {
            if positive {
                (
                    "1-,2-",
                    format!(
                        "{}^,{}^",
                        a[0][1],
                        a.last().unwrap()[1].min(b.last().unwrap()[1])
                    ),
                    format!(
                        "{}-,{}-",
                        (a.last().unwrap()[0] + 1).min(b.last().unwrap()[0]) + 1,
                        a.last().unwrap()[0].max(b.last().unwrap()[0]) + 1
                    ),
                )
            } else {
                (
                    "1^,2^",
                    format!(
                        "{}-,{}-",
                        a.last().unwrap()[1].min(b.last().unwrap()[1]),
                        a[0][1]
                    ),
                    format!(
                        "{}^,{}^",
                        a.last().unwrap()[0].max(b.last().unwrap()[0]) + 1,
                        a.last().unwrap()[0].min(b.last().unwrap()[0]) + 1
                    ),
                )
            }
        }
        _ => {
            ensure!(a.last() == b.last(), "misconfigured alt-prime event");
            if positive {
                (
                    "1^,2^",
                    format!(
                        "{}-,{}-",
                        a[0][0].max(b[0][0]) + 1,
                        a.last().unwrap()[0] + 1
                    ),
                    format!("{}^,{}^", a[0][1].min(b[0][1]), a[0][1].max(b[0][1])),
                )
            } else {
                (
                    "1-,2-",
                    format!(
                        "{}^,{}^",
                        a.last().unwrap()[0] + 1,
                        a[0][0].max(b[0][0]) + 1
                    ),
                    format!("{}-,{}-", a[0][1].max(b[0][1]), a[0][1].min(b[0][1])),
                )
            }
        }
    };
    write!(
        out,
        " flanks \"{flanks}\"; structure \"{structure}\"; splice_chain \"{chain}\";"
    )?;
    for (i, &p) in psi.iter().enumerate() {
        if i > 0 {
            write!(out, ";")?;
        }
        write!(out, " psi_{} \"{}\"", bytes_repr(&samples[i]), float(p, 6))?;
    }
    writeln!(out, ";")?;
    Ok(true)
}

fn write_inner(
    out: &mut dyn Write,
    format: Format,
    events: &[Event],
    counts_path: &Path,
    samples: &[&str],
    indices: &[usize],
) -> Result<()> {
    let kind = events[0].event_type;
    if format == Format::Txt {
        txt_header(out, kind, samples)?;
    }
    if matches!(format, Format::Gff3 | Format::Gtf) {
        writeln!(out, "##gff-version 3")?;
    }
    if format == Format::Icgc {
        write!(
            out,
            "event_id\tevent_type\tevent_chr\tis_annotated\tevent_coordinates\talt_region_coordinates\tgene_name"
        )?;
    }
    if format == Format::Tcga {
        write!(out, "gene\teventtype\tcoordinates")?;
    }
    if matches!(format, Format::Icgc | Format::Tcga) {
        for sample in samples {
            write!(out, "\t{sample}")?;
        }
        writeln!(out)?;
    }
    let count_file = if matches!(format, Format::Gff3 | Format::Gtf | Format::Bed) {
        None
    } else {
        Some(hdf5::File::open(counts_path)?)
    };
    let stored_samples = if format == Format::Structured {
        hdf5io::read_strings(&count_file.as_ref().unwrap().dataset("samples")?)?
    } else {
        Vec::new()
    };
    for &index in indices {
        let event = events
            .get(index)
            .ok_or_else(|| anyhow::anyhow!("output event index out of range"))?;
        match format {
            Format::Gff3 | Format::Gtf => gff(out, event, format == Format::Gtf)?,
            Format::Bed => {
                if !bed(out, event)? {
                    break;
                }
            }
            _ => {
                let input = count_file.as_ref().unwrap();
                let psi = input
                    .dataset("psi")?
                    .read_slice_1d::<f64, _>((.., index))?
                    .to_vec();
                if format == Format::Structured {
                    if !structured(out, event, &psi, &stored_samples)? {
                        break;
                    }
                } else if format == Format::Icgc {
                    icgc(out, event, &psi, samples.len())?;
                } else {
                    let counts =
                        input
                            .dataset("event_counts")?
                            .read_slice_2d::<f64, _>((.., .., index))?;
                    if format == Format::Txt {
                        txt(out, event, &counts, &psi, samples.len())?;
                    } else {
                        tcga(out, event, &counts, samples.len())?;
                    }
                }
            }
        }
    }
    Ok(())
}

pub fn write(
    path: &Path,
    format: Format,
    events: &[Event],
    counts_path: &Path,
    samples: &[&str],
    indices: &[usize],
) -> Result<()> {
    if events.is_empty() {
        return Ok(());
    }
    let file = BufWriter::new(File::create(path)?);
    if path.extension().is_some_and(|e| e == "gz")
        && !matches!(format, Format::Gff3 | Format::Gtf | Format::Bed)
    {
        let mut out = GzEncoder::new(file, Compression::best());
        write_inner(&mut out, format, events, counts_path, samples, indices)?;
        out.finish()?.flush()?;
    } else {
        let mut out = file;
        write_inner(&mut out, format, events, counts_path, samples, indices)?;
        out.flush()?;
    }
    Ok(())
}
