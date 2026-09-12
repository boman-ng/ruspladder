//! GTF locus identity follows the GTF 2.2 gene_id contract. Biological IDs
//! can have several placements (NCBI GFF3); keep those coordinate frames apart.
//! See ANNOTATION_LOCUS_PROPOSAL.md for sources and compatibility boundaries.
use crate::cache;
use anyhow::{Context, Result, ensure};
use clap::ValueEnum;
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum AnnotationMode {
    #[default]
    Spladder,
    Locus,
}

impl AnnotationMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Spladder => "spladder",
            Self::Locus => "locus",
        }
    }
}

// A transcript belongs to a gene placement, so include its parent in its scope.
type Scope = (String, String, String);
type Placements = BTreeMap<String, BTreeSet<Scope>>;
type Names = BTreeMap<(String, Scope), String>;

fn tag_value<'a>(tag: &'a str, key: &str) -> Option<&'a str> {
    let rest = tag.trim_start().strip_prefix(key)?;
    if !rest.starts_with(char::is_whitespace) {
        return None;
    }
    let rest = rest.trim();
    rest.strip_prefix('"')?.strip_suffix('"')
}

fn attribute<'a>(attributes: &'a str, key: &str) -> Option<&'a str> {
    attributes.split(';').find_map(|tag| tag_value(tag, key))
}

fn visit(path: &Path, mut consume: impl FnMut(&str, &[&str]) -> Result<()>) -> Result<()> {
    let mut input = BufReader::new(File::open(path)?);
    let mut line = String::new();
    let mut number = 0;
    while input.read_line(&mut line)? != 0 {
        number += 1;
        if line.starts_with('#') || line.trim().is_empty() {
            consume(&line, &[])?;
        } else {
            let fields: Vec<_> = line.trim_end_matches(['\r', '\n']).split('\t').collect();
            ensure!(
                fields.len() == 9,
                "{}:{number}: expected 9 GTF fields",
                path.display()
            );
            consume(&line, &fields).with_context(|| format!("{}:{number}", path.display()))?;
        }
        line.clear();
    }
    Ok(())
}

fn names(placements: Placements) -> Names {
    let mut used: BTreeSet<String> = placements.keys().cloned().collect();
    let mut result = Names::new();
    for (id, scopes) in placements {
        if scopes.len() < 2 {
            continue;
        }
        for (index, scope) in scopes.into_iter().enumerate() {
            // Deterministic within the annotation, independent of record order.
            // Reserve original IDs too: a real identifier may already use this suffix.
            let mut name = format!("{id}__ruspladder_locus_{}", index + 1);
            while !used.insert(name.clone()) {
                name.push('_');
            }
            result.insert((id.clone(), scope), name);
        }
    }
    result
}

/// Persist a lossless (except disambiguated IDs) GTF for both native import and
/// an independent Python oracle. Like existing annotation caches, inputs are
/// immutable; use a new input path or remove its companion caches after editing.
pub fn normalize_gtf(path: &Path) -> Result<PathBuf> {
    ensure!(
        path.extension()
            .is_some_and(|s| s.eq_ignore_ascii_case("gtf")),
        "--annotation-mode locus requires a source GTF; use GFF3 feature ID/Parent input with the default mode"
    );
    let normalized = PathBuf::from(format!("{}.locus.gtf", path.display()));
    let mapping = PathBuf::from(format!("{}.loci.tsv", path.display()));
    if normalized.is_file() && mapping.is_file() {
        return Ok(normalized);
    }
    let mut genes = Placements::new();
    let mut transcripts = Placements::new();
    visit(path, |_, fields| {
        if fields.is_empty() {
            return Ok(());
        }
        let Some(gene) = attribute(fields[8], "gene_id") else {
            return Ok(());
        };
        let scope = (String::new(), fields[0].to_owned(), fields[6].to_owned());
        genes
            .entry(gene.to_owned())
            .or_default()
            .insert(scope.clone());
        if let Some(tx) = attribute(fields[8], "transcript_id") {
            transcripts.entry(tx.to_owned()).or_default().insert((
                gene.to_owned(),
                scope.1,
                scope.2,
            ));
        }
        Ok(())
    })?;
    let genes = names(genes);
    let transcripts = names(transcripts);
    cache::atomic_write(&mapping, |temporary| {
        let mut out = BufWriter::new(File::create(temporary)?);
        writeln!(
            out,
            "feature\toriginal_id\tlocus_id\tgene_id\tseqname\tstrand"
        )?;
        for (feature, names) in [("gene", &genes), ("transcript", &transcripts)] {
            for ((id, (parent, chr, strand)), name) in names {
                writeln!(out, "{feature}\t{id}\t{name}\t{parent}\t{chr}\t{strand}")?;
            }
        }
        out.flush()?;
        Ok(())
    })?;
    cache::atomic_write(&normalized, |temporary| {
        let mut out = BufWriter::new(File::create(temporary)?);
        visit(path, |line, fields| {
            if fields.is_empty() {
                out.write_all(line.as_bytes())?;
                return Ok(());
            }
            let Some(gene) = attribute(fields[8], "gene_id") else {
                out.write_all(line.as_bytes())?;
                return Ok(());
            };
            let chr = fields[0].to_owned();
            let strand = fields[6].to_owned();
            let new_gene = genes.get(&(
                gene.to_owned(),
                (String::new(), chr.clone(), strand.clone()),
            ));
            let new_tx = attribute(fields[8], "transcript_id")
                .and_then(|tx| transcripts.get(&(tx.to_owned(), (gene.to_owned(), chr, strand))));
            if new_gene.is_none() && new_tx.is_none() {
                out.write_all(line.as_bytes())?;
                return Ok(());
            }
            let offset = fields[..8].iter().map(|s| s.len() + 1).sum::<usize>();
            out.write_all(&line.as_bytes()[..offset])?;
            for tag in line[offset..].split_inclusive(';') {
                let body = tag.strip_suffix(';').unwrap_or(tag);
                let replacement = [("gene_id", new_gene), ("transcript_id", new_tx)]
                    .into_iter()
                    .find_map(|(key, replacement)| {
                        tag_value(body, key).and_then(|old| replacement.map(|new| (old, new)))
                    });
                if let Some((old, new)) = replacement {
                    let start = tag.find('"').unwrap() + 1;
                    write!(out, "{}{}{}", &tag[..start], new, &tag[start + old.len()..])?;
                } else {
                    out.write_all(tag.as_bytes())?;
                }
            }
            Ok(())
        })?;
        out.flush()?;
        Ok(())
    })?;
    eprintln!(
        "Normalized {} gene placements and {} transcript placements: {} (IDs: {})",
        genes.len(),
        transcripts.len(),
        normalized.display(),
        mapping.display()
    );
    Ok(normalized)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotation::{AnnotationFilters, read_annotation};
    use std::fs;

    #[test]
    fn keeps_coordinate_frames_and_transcript_placements_separate() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("input.gtf");
        let original = concat!(
            "# Synthetic exons at the audited ZNF84 placement bounds\n",
            "12\tRefSeq\texon\t133613872\t133613972\t.\t+\t.\tgene_id \"ZNF84\"; transcript_id \"tx\"; note \"kept\";\n",
            "12\tRefSeq\texon\t133639790\t133639890\t.\t+\t.\tgene_id \"ZNF84\"; transcript_id \"tx\";\n",
            "Un_gl000223\tRefSeq\texon\t42780\t42880\t.\t-\t.\tgene_id \"ZNF84\"; transcript_id \"tx\";\n",
            "Un_gl000223\tRefSeq\texon\t68470\t68570\t.\t-\t.\tgene_id \"ZNF84\"; transcript_id \"tx\";\n",
        );
        fs::write(&input, original)?;
        let normalized = normalize_gtf(&input)?;
        let annotation = read_annotation(&normalized, AnnotationFilters::default())?;
        assert_eq!(annotation.genes.len(), 2);
        let a = &annotation.genes[0];
        let b = &annotation.genes[1];
        assert_eq!(
            (&*a.chr, a.strand, a.start, a.stop),
            ("12", '+', 133613871, 133639890)
        );
        assert_eq!(
            (&*b.chr, b.strand, b.start, b.stop),
            ("Un_gl000223", '-', 42779, 68570)
        );
        assert_ne!(a.name, b.name);
        assert_ne!(a.transcripts, b.transcripts);
        let text = fs::read_to_string(&normalized)?;
        assert!(text.contains("note \"kept\";"));
        assert_eq!(text.lines().count(), original.lines().count());
        assert_eq!(fs::read_to_string(&input)?, original);
        // Record order cannot decide placement IDs (unlike first-record metadata).
        let reverse = dir.path().join("reverse.gtf");
        fs::write(
            &reverse,
            original.lines().rev().collect::<Vec<_>>().join("\n") + "\n",
        )?;
        normalize_gtf(&reverse)?;
        assert_eq!(
            fs::read(input.with_extension("gtf.loci.tsv"))?,
            fs::read(reverse.with_extension("gtf.loci.tsv"))?
        );
        Ok(())
    }

    #[test]
    fn retains_long_single_placement_and_all_unmodified_bytes() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("long.gtf");
        let original = concat!(
            "1\tX\tgene\t10\t9000000\t.\t+\t.\tgene_id \"long\";\r\n",
            "1\tX\texon\t10\t40\t.\t+\t.\tgene_id \"long\"; transcript_id \"a\";\r\n",
            "1\tX\texon\t8999900\t9000000\t.\t+\t.\tgene_id \"long\"; transcript_id \"b\";\r\n",
        );
        fs::write(&input, original)?;
        let normalized = normalize_gtf(&input)?;
        assert_eq!(fs::read_to_string(&normalized)?, original);
        let annotation = read_annotation(&normalized, AnnotationFilters::default())?;
        assert_eq!(annotation.genes.len(), 1);
        assert_eq!(annotation.genes[0].transcripts, ["a", "b"]);
        assert_eq!(annotation.genes[0].stop, 9000000);
        Ok(())
    }

    #[test]
    fn preserves_explicit_gene_rows_and_avoids_existing_id_collisions() -> Result<()> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("explicit.gtf");
        fs::write(
            &input,
            concat!(
                "1\tX\tgene\t10\t40\t.\t+\t.\tgene_id \"g\";\n",
                "2\tX\tgene\t100\t140\t.\t-\t.\tgene_id \"g\";\n",
                "1\tX\texon\t10\t40\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";\n",
                "2\tX\texon\t100\t140\t.\t-\t.\tgene_id \"g\"; transcript_id \"t\";\n",
                "3\tX\texon\t50\t60\t.\t+\t.\tgene_id \"g__ruspladder_locus_1\"; transcript_id \"third\";\n",
            ),
        )?;
        let annotation = read_annotation(&normalize_gtf(&input)?, AnnotationFilters::default())?;
        assert_eq!(annotation.genes.len(), 3);
        let names: BTreeSet<_> = annotation.genes.iter().map(|g| &*g.name).collect();
        assert_eq!(names.len(), 3);
        assert!(names.contains("g__ruspladder_locus_1"));
        assert_eq!(annotation.genes[1].chr, "2");
        assert_eq!(annotation.genes[1].start, 99);
        Ok(())
    }
}
