#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""
Generates a bedMethyl file (modkit pileup format, 18 columns) from the MM/ML tags
of a BAM simulated by simulate_glioma_hg38.py. As modkit is not installed here, this script
directly reimplements the per-CpG-site aggregation from the SAM tags themselves -
no dependency on an external tool to test the rest of the chain.

Mimics `modkit pileup --cpg --combine-strands` (as run by Rapid-CNS2_nf): calls from
reads on both strands are aggregated on the C of the CpG on the + strand, and the
strand column is '.'.

The MM/ML tags are decoded with pysam (AlignedSegment.modified_bases), which handles the
orientation of reads aligned to the - strand: for those, the MM tag counts the Cs of the
ORIGINAL read (reverse complement of SEQ), so the modified base appears as a G in SEQ,
one base after the C of the CpG on the + strand.

Only columns 1,2,3 (chrom/start/end), 6 (strand) and 10,11 (coverage, % methylated)
are read by scr/methylation_classification.R; the others are filled in for consistency
with the standard modkit format.

Usage: python bedmethyl_from_bam.py --bam sim.bam --out sim.bedmethyl.bed
"""
import argparse
import pysam

ap = argparse.ArgumentParser(description=__doc__, formatter_class=type(
    'HelpFormatter', (argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter), {}))
ap.add_argument('--bam', required=True, help='BAM with MM/ML tags')
ap.add_argument('--out', required=True, help='bedMethyl file to write')
ap.add_argument('--min-cov', type=int, default=1, help='minimum coverage of a CpG site to be written')
a = ap.parse_args()

bam = pysam.AlignmentFile(a.bam, 'rb')
# site = (chrom, pos0 of the C of the CpG on the + strand) -> [n_total, n_methylated]
sites = {}

n_reads = n_tagged = 0
for r in bam.fetch(until_eof=True):
    n_reads += 1
    if r.is_unmapped or r.is_secondary or r.is_supplementary:
        continue
    if not r.has_tag('MM') or r.query_sequence is None:
        continue
    n_tagged += 1
    q2r = dict(r.get_aligned_pairs(matches_only=True))
    chrom = r.reference_name
    for (base, strand, mod), calls in r.modified_bases.items():
        if base != 'C' or mod != 'm':
            continue
        for qpos, ml in calls:
            # qpos is an index in SEQ: the C (forward read) or the G opposite the
            # modified C (reverse read). The CpG must be aligned without a gap.
            if r.is_reverse:
                rg, rc = q2r.get(qpos), q2r.get(qpos - 1)
            else:
                rc, rg = q2r.get(qpos), q2r.get(qpos + 1)
            if rc is None or rg is None or rg != rc + 1:
                continue
            d = sites.setdefault((chrom, rc), [0, 0])
            d[0] += 1
            if ml >= 128:
                d[1] += 1

print(f"{n_reads} reads read, {n_tagged} with MM tag, {len(sites)} CpG sites covered")

n_written = 0
with open(a.out, 'w') as f:
    for (chrom, pos0), (n, nmod) in sorted(sites.items()):
        if n < a.min_cov:
            continue
        pct = 100.0 * nmod / n
        # 18 modkit bedMethyl columns; only chrom,start,end,strand,Nvalid_cov,%mod matter
        # for methylation_classification.R (columns 1,2,3,6,10,11)
        row = [chrom, pos0, pos0 + 1, 'm', n, '.', pos0, pos0 + 1, '255,0,0',
               n, f'{pct:.2f}', nmod, n - nmod, 0, 0, 0, 0, 0]
        f.write('\t'.join(map(str, row)) + '\n')
        n_written += 1

print(f"{n_written} lines written to {a.out}")
