// Adapted from SplAdder v3.1.1 helpers.py (BSD-3-Clause).
use crate::{
    annotation::Gene,
    reads::ReadFilter,
    reference::{Reference, reverse_complement},
};
use anyhow::{Result, ensure};
use std::collections::BTreeMap;

pub type IntronLists = Vec<[Vec<[i64; 3]>; 2]>;

/// Count intervals using sorted starts/ends. This answers the same half-open
/// overlap query as intervaltree without materializing overlapping gene lists.
pub fn filter_ambiguous(introns: &mut IntronLists, genes: &[Gene], offset: i64) -> Result<()> {
    ensure!(
        introns.len() == genes.len(),
        "intron lists must match genes"
    );
    let mut bounds = BTreeMap::new();
    for gene in genes {
        ensure!(gene.start < gene.stop, "gene has an empty interval");
        let (starts, ends) = bounds
            .entry((gene.chr.as_str(), gene.strand))
            .or_insert_with(|| (Vec::<i64>::new(), Vec::<i64>::new()));
        starts.push(gene.start);
        ends.push(gene.stop);
    }
    for (starts, ends) in bounds.values_mut() {
        starts.sort_unstable();
        ends.sort_unstable();
    }
    for (gene, lists) in genes.iter().zip(introns) {
        for (strand, introns) in ['+', '-'].into_iter().zip(lists) {
            introns.retain(|intron| {
                let Some((starts, ends)) = bounds.get(&(gene.chr.as_str(), strand)) else {
                    return false;
                };
                let begin = intron[0] - offset;
                let end = intron[1] + offset;
                if begin >= end {
                    return false;
                }
                let opened = starts.partition_point(|&p| p < end);
                let closed = ends.partition_point(|&p| p <= begin);
                opened - closed == 1
            });
        }
    }
    Ok(())
}

pub fn filter_consensus(
    introns: &mut IntronLists,
    genes: &[Gene],
    reference: &Reference,
    lenient: bool,
) -> Result<()> {
    ensure!(
        introns.len() == genes.len(),
        "intron lists must match genes"
    );
    for si in 0..2 {
        for (gene, lists) in genes.iter().zip(introns.iter_mut()) {
            let mut kept = Vec::new();
            for &intron in &lists[si] {
                let left = reference.fetch(&gene.chr, intron[0], intron[0] + 2)?;
                let right = reference.fetch(&gene.chr, intron[1] - 2, intron[1])?;
                let (donor, acceptor) = if si == 0 {
                    (left, right)
                } else {
                    (reverse_complement(&right)?, reverse_complement(&left)?)
                };
                if acceptor == b"AG" && (donor == b"GT" || (lenient && donor == b"GC")) {
                    kept.push(intron);
                }
            }
            lists[si] = kept;
        }
    }
    Ok(())
}

/// Match the source's global filter mutation and re-query only infeasible genes.
pub fn make_feasible(
    introns: &mut IntronLists,
    filter: &mut ReadFilter,
    mut reload: impl FnMut(&[usize], &ReadFilter) -> Result<IntronLists>,
) -> Result<()> {
    loop {
        let indices: Vec<_> = introns
            .iter()
            .enumerate()
            .filter_map(|(i, lists)| lists.iter().any(|v| v.len() > 200).then_some(i))
            .collect();
        if indices.is_empty() {
            return Ok(());
        }
        filter.exon_len = (filter.exon_len + 4).min(36);
        filter.mincount *= 2;
        filter.mismatch = (filter.mismatch - 1).max(0);
        let new = reload(&indices, filter)?;
        ensure!(
            new.len() == indices.len(),
            "reloaded intron lists must match queried genes"
        );
        for (i, lists) in indices.into_iter().zip(new) {
            introns[i] = lists;
        }
    }
}
