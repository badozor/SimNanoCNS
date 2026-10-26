#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""Checks the BAM produced by simulate_glioma_hg38.py (independently of the generator).

Usage: python verify_bam.py --bam sim/SIM_GBM_HG38.bam --data-dir Rapid-CNS2_nf/data [--ref hg38.fa | --standin]
       (with --standin, give the same --seed as the simulation)

Automatically picks, among the variants and wildtype sites of HOTSPOTS, those that are
actually covered in the BAM (to adapt to --profile gbm, astro or oligo).
"""
import argparse, os, sys, tempfile, statistics as st
import numpy as np
import pysam
from simulate_glioma_hg38 import RealRef, StandInRef, read_bed, HOTSPOTS

ap = argparse.ArgumentParser(description=__doc__, formatter_class=type(
    'HelpFormatter', (argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter), {}))
ap.add_argument('--bam', required=True, help='BAM produced by simulate_glioma_hg38.py')
ap.add_argument('--data-dir', required=True, help='data/ directory of Rapid-CNS2_nf')
ap.add_argument('--ref', help='hg38 FASTA used for the simulation')
ap.add_argument('--standin', action='store_true', help='the BAM was simulated with --standin')
ap.add_argument('--truth', help="truth_events.tsv (default: next to the BAM)")
ap.add_argument('--seed', type=int, default=2026, help='seed of the simulation (used with --standin)')
a = ap.parse_args()
if bool(a.ref) == bool(a.standin):
    sys.exit("Specify either --ref hg38.fa or --standin")

truth_path = a.truth or os.path.join(os.path.dirname(a.bam) or '.', 'truth_events.tsv')
truth_cnv = []  # (name, chrom, start, end, expected_ratio)
expected_bg = None
if os.path.exists(truth_path):
    with open(truth_path) as f:
        next(f)
        for line in f:
            name, typ, c, s0, e0, val = line.rstrip('\n').split('\t')
            if typ == 'copy_ratio':
                truth_cnv.append((name, c, int(s0), int(e0), float(val.split('=')[-1])))
            elif typ == 'background':
                expected_bg = float(val)
else:
    print(f"({truth_path} not found: CNV section skipped)")

bam = pysam.AlignmentFile(a.bam, 'rb')
panel = read_bed(f'{a.data_dir}/NPHD_panel_hg38.bed')
mgmt = read_bed(f'{a.data_dir}/mgmt_hg38.bed')[0]
genes = {r[3]: r for r in read_bed(f'{a.data_dir}/genes.bed')}

print('=== Header / format ===')
hdr = bam.header.to_dict()
ds = hdr['RG'][0]['DS']
print('basecall_model in @RG   :', [t for t in ds.split() if t.startswith('basecall_model=')][0])
print('contigs @SQ             :', len(hdr['SQ']), '| chr1 LN =', hdr['SQ'][0]['LN'])
n_mm = 0
for i, r in enumerate(bam.fetch(until_eof=True)):
    if i >= 2000:
        break
    n_mm += r.has_tag('MM')
print('reads with MM:Z (first 2000)    :', n_mm)

print('\n=== Mean coverage in the NPHD panel (samtools bedcov) ===')
with tempfile.NamedTemporaryFile('w', suffix='.bed', delete=False) as tmp:
    for c, s, e, _ in panel:
        tmp.write(f'{c}\t{s}\t{e}\n')
    bed_tmp = tmp.name
try:
    out = pysam.bedcov(bed_tmp, a.bam)        # samtools bedcov via pysam (no samtools binary needed)
finally:
    os.unlink(bed_tmp)
egfr_w = genes['EGFR']; cdk_w = genes['CDKN2A/B']
excl = [(egfr_w[0], egfr_w[1], egfr_w[2]), (cdk_w[0], cdk_w[1] - 50_000, cdk_w[2] + 50_000),
        (mgmt[0], mgmt[1] - 20_000, mgmt[2] + 20_000)]
for h in HOTSPOTS.values():
    excl.append((h['chrom'], h['pos0'] - 20_000, h['pos0'] + 20_000))
by = {}
for line in out.splitlines():
    c, s0, e0, tot = line.split('\t')[:4]; s0, e0 = int(s0), int(e0)
    if any(c == x[0] and s0 < x[2] and e0 > x[1] for x in excl):
        continue                      # focal event windows excluded from the per-chromosome calculation
    d = by.setdefault(c, [0, 0]); d[0] += int(tot); d[1] += e0 - s0
cov = {c: v[0] / v[1] for c, v in by.items() if v[1] > 0 and v[0] > 0}
print('Mean panel coverage per chromosome (outside event windows):')
print('  ' + '  '.join(f'{c}={cov[c]:.1f}x' for c in ('chr1', 'chr5', 'chr7', 'chr9', 'chr10') if c in cov))
print('  (not interpretable as a ratio: long reads overlap the neighbouring regions and the panel')
print('   density varies between chromosomes; the copy ratio is checked below, outside the panel)')

print('\n=== Copy number: OFF-panel reads per Mb (independent of the panel density) ===')
print('(works for a whole chromosome as well as for a partial arm: the window of each')
print(' event is taken directly from truth_events.tsv)')
import bisect
pan = {}
for c, s0, e0, _ in panel:
    pan.setdefault(c, []).append((s0, e0))
for c in pan:
    pan[c].sort()
pst = {c: [x[0] for x in v] for c, v in pan.items()}
pen = {c: [x[1] for x in v] for c, v in pan.items()}

off_starts = {}  # chrom -> sorted list of the start positions of the off-panel reads
bam.reset()      # ("excl" is NOT applied here: we want to be able to measure the density
                 #  INSIDE the focal event windows, not only around them)
for r in bam.fetch(until_eof=True):
    c, rs, re_ = r.reference_name, r.reference_start, r.reference_end
    i = bisect.bisect_right(pen.get(c, []), rs)
    if i < len(pst.get(c, [])) and pst[c][i] < re_:
        continue                       # overlaps the panel: skipped
    off_starts.setdefault(c, []).append(rs)
for c in off_starts:
    off_starts[c].sort()

def density(chrom, s0, e0):
    """Off-panel reads per Mb in [s0, e0)."""
    lst = off_starts.get(chrom, [])
    n = bisect.bisect_left(lst, e0) - bisect.bisect_left(lst, s0)
    return n / ((e0 - s0) / 1e6)

lens = {d['SN']: d['LN'] for d in hdr['SQ']}
event_chroms = {c for _, c, _, _, _ in truth_cnv}
neutral_names = [c for c in ('chr1', 'chr2', 'chr3', 'chr4', 'chr6', 'chr8', 'chr11', 'chr12')
                  if c not in event_chroms and c in lens and c in off_starts]
nb = st.median([density(c, 0, lens[c]) for c in neutral_names]) if neutral_names else None
LARGE_WINDOW = 5_000_000  # below this, the NPHD panel targets the window too densely for
                          # "off-panel reads" to remain a reliable estimator (see Focal events)
# chromosomes actually present in the BAM (simulate_glioma_hg38.py --chroms may restrict them):
# an event on a chromosome that was not simulated must be reported as such, not as ratio 0
simulated = {s.contig for s in bam.get_index_statistics() if s.mapped > 0}
if nb:
    print(f'(neutral baseline, median over {neutral_names}: {nb:.1f} reads/Mb)')
    for name, c, s0, e0, exp in truth_cnv:
        if e0 - s0 < LARGE_WINDOW:
            continue   # focal window: checked below (raw coverage, not "off-panel")
        if c not in simulated:
            print(f'{name:28s} {c}: chromosome not simulated (--chroms), not checked')
            continue
        d = density(c, s0, e0)
        print(f'{name:28s} {c}:{s0}-{e0}  {d:6.1f} reads/Mb   ratio = {d/nb:4.2f}   (expected {exp})')
else:
    print('(no neutral chromosome simulated: copy ratios not checked)')

def mean_cov(chrom, s, e, step=200):
    tot = n = 0
    for p in range(s, e, step):
        tot += sum(x[0] for x in bam.count_coverage(chrom, p, p + 1, quality_threshold=0)); n += 1
    return tot / max(n, 1)

print('\n=== Focal events (raw coverage: the "off-panel" density is not reliable here) ===')
# Reference chromosome for the baselines: a simulated chromosome WITHOUT any copy-number
# event (chr1 is not neutral for the oligo profile: 1p codeletion).
ref_chrom = neutral_names[0] if neutral_names else 'chr1'
bg = mean_cov(ref_chrom, 100_000_000, 100_200_000, 500)
panel_base = cov.get(ref_chrom, bg)   # panel depth on a neutral chromosome (computed above)
print(f'Off-target background {ref_chrom} (control) : {bg:.2f}x')
print(f'NPHD panel background {ref_chrom} (ref.)    : {panel_base:.1f}x')
focal = [t for t in truth_cnv if t[3] - t[2] < LARGE_WINDOW]
if not focal:
    print('(no focal event for this profile)')
for name, c, s0, e0, exp in focal:
    if c not in simulated:
        print(f'{name:35s} {c}: chromosome not simulated (--chroms), not checked')
        continue
    v = mean_cov(c, s0, e0, max(50, (e0 - s0) // 500))
    ratio = f'{v/panel_base:5.2f}' if panel_base else '  n/a'
    print(f'{name:35s} {c}:{s0}-{e0} (mean cov) : {v:6.2f}x   (ratio vs panel {ratio}, expected {exp})')
if focal:
    print('(the measured ratio can be lower than "expected": the mean covers the whole window,')
    print(' whereas the NPHD panel probes only cover part of the gene; the direction of the signal')
    print(' - clearly over- or under-represented relative to the reference panel - is the point to check)')

def counts(chrom, pos0):
    c = bam.count_coverage(chrom, pos0, pos0 + 1, quality_threshold=0)
    d = {b: c[i][0] for i, b in enumerate('ACGT')}
    return d, sum(d.values())

print('\n=== Hotspots (variants and wildtype sites) ===')
for name, h in HOTSPOTS.items():
    d, t = counts(h['chrom'], h['pos0'])
    if not t:
        continue   # hotspot not covered in this BAM (profile does not target it)
    alt_frac = d.get(h['alt_plus'], 0) / t
    non_ref = (t - d.get(h['ref_plus'], 0)) / t
    print(f"{name} {h['chrom']}:{h['pos0']+1} {d} depth {t}  "
          f"VAF({h['alt_plus']}) = {alt_frac:.2f}  non-ref = {non_ref:.2f}")

print('\n=== MGMT ===')
c, s, e = mgmt[0], mgmt[1], mgmt[2]
if c not in simulated:
    print(f'{c} not simulated (--chroms): MGMT not checked')
else:
    print(f'Mean depth {c}:{s}-{e} : {mean_cov(c, s, e, 20):.1f}x '
          '(--minimumMgmtCov of Rapid-CNS2: 5 by default, 10 in rapidcns2_gui_en.py and the USB run.sh)')
res = {'fwd': [], 'rev': []}; ctrl = []
for r in bam.fetch(c, s, e):
    if not r.has_tag('MM'):
        continue
    q2r = dict(r.get_aligned_pairs())
    for (_, _, _), lst in r.modified_bases.items():
        for qpos, ml in lst:
            rp = q2r.get(qpos)
            if rp is not None and s <= rp < e:
                res['rev' if r.is_reverse else 'fwd'].append(ml >= 128)
n_ctrl = 0
for pc, ps, pe, _ in panel:
    if pc != c or (ps < e + 20_000 and pe > s - 20_000):
        continue
    for r in bam.fetch(pc, ps, pe):
        if r.has_tag('MM'):
            for _, lst in r.modified_bases.items():
                ctrl += [ml >= 128 for _, ml in lst]
    n_ctrl += 1
    if n_ctrl >= 40:
        break
for k in ('fwd', 'rev'):
    if res[k]:
        print(f'Methylated sites (ML>=128), {k} reads: {np.mean(res[k]):.2f} (n={len(res[k])})')
if ctrl:
    exp_txt = f'~{expected_bg:.2f}' if expected_bg is not None else 'see truth_events.tsv'
    print(f'Control {c} (40 panel regions, outside MGMT) : {np.mean(ctrl):.2f} (n={len(ctrl)}) '
          f'(expected bimodal background {exp_txt}, slightly shifted by any transplanted HM450 probes)')

print('\n=== Alignment / sequence consistency (200 reads) ===')
ref = RealRef(a.ref) if a.ref else StandInRef(a.seed)
mis = tot = 0
bam.reset()
for i, r in enumerate(bam.fetch(until_eof=True)):
    if i >= 200:
        break
    if r.reference_name not in ref.lengths:
        continue
    seq = r.query_sequence
    rs = ref.fetch(r.reference_name, r.reference_start, r.reference_end)
    for qp, rp in r.get_aligned_pairs(matches_only=True):
        tot += 1
        mis += seq[qp] != rs[rp - r.reference_start]
if tot:
    print(f'Mismatch rate at aligned positions : {mis/tot*100:.1f}% (expected ~2%, sub 50% of 4%)')
    if mis / tot > 0.10 and a.standin:
        print('  -> far too high: with --standin, give verify_bam.py the --seed used for the simulation')
