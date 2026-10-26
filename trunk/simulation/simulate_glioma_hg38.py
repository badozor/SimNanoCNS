#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""
Simulate a Nanopore BAM (adaptive sampling, NPHD panel) of a diffuse glioma in hg38 coordinates,
to test Rapid-CNS2_nf. Three profiles:

  --profile astro   IDH-mutant astrocytoma (WHO CNS5, grade 2-3; grade 4 with --cdkn2a-del)
      IDH1 R132H (VAF = purity/2), TP53 R175H (VAF = purity/2), ATRX 1 bp frameshift
      (VAF = purity/2, XX donor), diploid genome (no +7/-10, no 1p/19q codeletion,
      no EGFR amplification, TERT wildtype), methylated MGMT promoter, G-CIMP background (hypermethylated).
  --profile oligo   IDH-mutant, 1p/19q-codeleted oligodendroglioma (WHO CNS5, grade 2-3)
      IDH1 R132H (VAF = purity/2), TERT promoter mutation C228T (VAF = purity/2, unlike the
      astrocytoma), 1p/19q codeletion (loss of one copy over the whole arm, ratio 0.5), CIC
      1 bp frameshift (VAF = purity/2), TP53 wildtype, ATRX untouched (no
      simulated site), methylated MGMT promoter, G-CIMP background (hypermethylated).
  --profile gbm     IDH-wildtype glioblastoma:
      +7/-10, EGFR amplification, CDKN2A/B deletion, TERT C228T (VAF 0.80), IDH1 wildtype,
      methylated MGMT promoter.

Usage (real reference, recommended):
  python simulate_glioma_hg38.py --profile oligo --ref /refs/hg38.fa --data-dir Rapid-CNS2_nf/data --out sim_oligo

Usage (plumbing test, SYNTHETIC sequence, unusable with the real hg38):
  python simulate_glioma_hg38.py --profile oligo --standin --data-dir Rapid-CNS2_nf/data --out sim_oligo_standin
"""

import argparse, os, re, sys, functools
from array import array
import numpy as np
import pysam

_NOT_ACGT = re.compile('[^ACGT]')  # detects N and any IUPAC ambiguity code (Y, R, W, S, K, M...)

HG38_LEN = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
    'chrM': 16569,
}
CHROMS = list(HG38_LEN)

# GRCh38 RefSeq accessions (stable since GRCh38.p1: patches p1-p14 only affect the
# alternate scaffolds, not these primary chromosome accessions). Allows the use of an NCBI
# FASTA (GCF_000001405.40_GRCh38.p14_genomic.fa) whose contigs are named NC_000001.11, etc.
REFSEQ_ACCESSION = {
    'chr1': 'NC_000001.11', 'chr2': 'NC_000002.12', 'chr3': 'NC_000003.12', 'chr4': 'NC_000004.12',
    'chr5': 'NC_000005.10', 'chr6': 'NC_000006.12', 'chr7': 'NC_000007.14', 'chr8': 'NC_000008.11',
    'chr9': 'NC_000009.12', 'chr10': 'NC_000010.11', 'chr11': 'NC_000011.10', 'chr12': 'NC_000012.12',
    'chr13': 'NC_000013.11', 'chr14': 'NC_000014.9', 'chr15': 'NC_000015.10', 'chr16': 'NC_000016.10',
    'chr17': 'NC_000017.11', 'chr18': 'NC_000018.10', 'chr19': 'NC_000019.10', 'chr20': 'NC_000020.11',
    'chr21': 'NC_000021.9', 'chr22': 'NC_000022.11', 'chrX': 'NC_000023.11', 'chrY': 'NC_000024.10',
    'chrM': 'NC_012920.1',
}

# Hotspots (0-based, + strand). ref_plus = reference base expected on the + strand.
HOTSPOTS = {
    'TERT_C228T': dict(chrom='chr5', pos0=1295112, ref_plus='G', alt_plus='A'),
    'IDH1_R132H': dict(chrom='chr2', pos0=208248387, ref_plus='C', alt_plus='T'),
    'TP53_R175H': dict(chrom='chr17', pos0=7675087, ref_plus='C', alt_plus='T'),
}
# Chromosome arms 1p and 19q (hg38): approximate centromere boundaries (middle of the
# UCSC cytoBand acrocentric region), sufficient to simulate a whole-arm loss.
ARM_1P = ('chr1', 0, 123_400_000)
ARM_19Q = ('chr19', 26_200_000, HG38_LEN['chr19'])

# vaf: 'het' = purity/2 (heterozygous variant), or a number
PROFILES = {
    'gbm': dict(label='IDH-wildtype glioblastoma', chrom_ratio={'chr7': 1.5, 'chr10': 0.5},
                events=['EGFR_amp', 'CDKN2A_del'], snvs={'TERT_C228T': 0.80}, wt=['IDH1_R132H'],
                del_gene=None, low_frac=0.5),
    'astro': dict(label='IDH-mutant astrocytoma', chrom_ratio={}, events=[],
                  snvs={'IDH1_R132H': 'het', 'TP53_R175H': 'het'}, wt=['TERT_C228T'],
                  del_gene='ATRX', low_frac=0.30),
    'oligo': dict(label='IDH-mutant, 1p/19q-codeleted oligodendroglioma', chrom_ratio={},
                  events=['1p_codel', '19q_codel'],
                  snvs={'IDH1_R132H': 'het', 'TERT_C228T': 'het'}, wt=['TP53_R175H'],
                  del_gene='CIC', low_frac=0.30),
}

MEAN_RL = 3100          # mean read length (lognormal mu=7.8, sigma=0.7)
ERR_RATE = 0.04         # basecalling errors (sup): 50% sub, 25% del, 25% ins
MGMT_BETA = 0.88
COMP = str.maketrans('ACGT', 'TGCA')


# ----------------------------------------------------------------- references
class RealRef:
    def __init__(self, path):
        if not os.path.exists(path + '.fai'):
            pysam.faidx(path)
        self.fa = pysam.FastaFile(path)
        avail = dict(zip(self.fa.references, self.fa.lengths))

        # 1) contigs already named chr1, chr2, ... (UCSC / Ensembl with prefix)
        if all(c in avail for c in CHROMS):
            self._actual = {c: c for c in CHROMS}
            self.lengths = {c: avail[c] for c in CHROMS}
            return

        # 2) contigs named by RefSeq accession (NCBI: NC_000001.11, ...)
        if all(REFSEQ_ACCESSION[c] in avail for c in CHROMS):
            self._actual = {c: REFSEQ_ACCESSION[c] for c in CHROMS}
            self.lengths = {c: avail[REFSEQ_ACCESSION[c]] for c in CHROMS}
            mismatched = [c for c in CHROMS if self.lengths[c] != HG38_LEN[c]]
            if mismatched:
                print(f"WARNING: length differs from the internal table for {mismatched[:5]} "
                      "(probably a sub-assembly or a different patch version); "
                      "the lengths read from the FASTA are used.", file=sys.stderr)
            print(f"Contigs detected by RefSeq accession (e.g. chr1 -> {REFSEQ_ACCESSION['chr1']}).")
            return

        missing = [c for c in CHROMS if c not in avail and REFSEQ_ACCESSION[c] not in avail]
        sys.exit(f"ERROR: hg38 contigs not found in {path} (neither 'chr*' names nor RefSeq accessions "
                 f"NC_0000xx): {missing[:5]}. Examples of contigs present: {list(avail)[:5]}")

    def fetch(self, chrom, s, e):
        return self.fa.fetch(self._actual[chrom], s, e).upper()


class StandInRef:
    """Deterministic synthetic sequence (GC 42%, CpG-depleted): NOT the real hg38."""
    B = 8192
    def __init__(self, seed):
        self.seed = seed
        self.lengths = dict(HG38_LEN)
        self.idx = {c: i for i, c in enumerate(CHROMS)}
        self.overrides = {}
        for h in HOTSPOTS.values():
            self.overrides.setdefault(h['chrom'], {})[h['pos0']] = h['ref_plus']
    @functools.lru_cache(maxsize=4096)
    def block(self, ci, b):
        rng = np.random.default_rng([self.seed, ci, b])
        a = rng.choice(np.frombuffer(b'ACGT', dtype=np.uint8), size=self.B, p=[.29, .21, .21, .29])
        cg = np.flatnonzero((a[:-1] == 67) & (a[1:] == 71))
        a[cg[rng.random(cg.size) < 0.75] + 1] = 65
        return a
    def finalize_overrides(self):
        """To be called once all overrides have been added (hotspots + HM450 probes): converts
        each dict {pos: base} into sorted arrays for a bisection search in
        fetch(), instead of a linear scan of ALL positions on every read - decisive
        with --mnp-targets, which adds about 200,000 forced positions in total."""
        self._sorted = {}
        for chrom, d in self.overrides.items():
            pos = np.array(sorted(d), dtype=np.int64)
            bases = np.array([ord(d[p]) for p in pos.tolist()], dtype=np.uint8)
            self._sorted[chrom] = (pos, bases)

    def fetch(self, chrom, s, e):
        ci = self.idx[chrom]; b0, b1 = s // self.B, (e - 1) // self.B
        arr = np.concatenate([self.block(ci, b) for b in range(b0, b1 + 1)])
        out = arr[s - b0 * self.B: e - b0 * self.B].copy()
        sorted_ov = getattr(self, '_sorted', None)
        if sorted_ov is not None and chrom in sorted_ov:
            pos, bases = sorted_ov[chrom]
            lo, hi = np.searchsorted(pos, [s, e])
            if hi > lo:
                out[pos[lo:hi] - s] = bases[lo:hi]
        else:
            for pos, base in self.overrides.get(chrom, {}).items():
                if s <= pos < e:
                    out[pos - s] = ord(base)
        return out.tobytes().decode()


# ------------------------------------------------------------------- BED / CN
def read_bed(path):
    rows = []
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 3 and not line.startswith(('#', 'track')):
            rows.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else ''))
    return rows


def build_model(data_dir, prof, cdkn2a_del):
    panel = read_bed(os.path.join(data_dir, 'NPHD_panel_hg38.bed'))
    mgmt = read_bed(os.path.join(data_dir, 'mgmt_hg38.bed'))[0]
    genes = {r[3]: r for r in read_bed(os.path.join(data_dir, 'genes.bed'))}
    defs = {
        'EGFR_amp': (genes['EGFR'][0], genes['EGFR'][1], genes['EGFR'][2], 9.0, 'EGFR_high_level_amplification'),
        'CDKN2A_del': (genes['CDKN2A/B'][0], genes['CDKN2A/B'][1] - 50_000, genes['CDKN2A/B'][2] + 50_000,
                       0.02, 'CDKN2A_B_homozygous_deletion'),
        '1p_codel': (ARM_1P[0], ARM_1P[1], ARM_1P[2], 0.5, '1p_codeletion'),
        '19q_codel': (ARM_19Q[0], ARM_19Q[1], ARM_19Q[2], 0.5, '19q_codeletion'),
    }
    names = list(prof['events']) + (['CDKN2A_del'] if cdkn2a_del and 'CDKN2A_del' not in prof['events'] else [])
    events = [defs[n] for n in names]
    return panel, mgmt, events, dict(prof['chrom_ratio'])


def frameshift_site(data_dir, gene):
    """Position (0-based) in the middle of the largest coding exon of the NM_ transcript with the longest CDS
    (refGene.txt.gz). A 1 bp deletion there causes a frameshift."""
    import gzip
    best = None
    for line in gzip.open(os.path.join(data_dir, 'refGene.txt.gz'), 'rt'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 13 or f[12] != gene or not f[1].startswith('NM_'):
            continue
        cds_s, cds_e = int(f[6]), int(f[7])
        if best is None or cds_e - cds_s > best[0]:
            best = (cds_e - cds_s, f[2], cds_s, cds_e, f[9], f[10], f[1])
    if best is None:
        sys.exit(f'ERROR: gene {gene} not found in refGene.txt.gz')
    _, chrom, cs, ce, starts, ends, nm = best
    ex = [(max(int(a), cs), min(int(b), ce)) for a, b in zip(starts.strip(',').split(','), ends.strip(',').split(','))]
    ex = [x for x in ex if x[1] > x[0]]
    s0, e0 = max(ex, key=lambda x: x[1] - x[0])
    return chrom, (s0 + e0) // 2, nm


def ratio_at(chrom, pos, events, chrom_ratio):
    for c, s, e, r, _ in events:
        if c == chrom and s <= pos < e:
            return r
    return chrom_ratio.get(chrom, 1.0)


# --------------------------------------------------------------------- reads
def revcomp(s):
    return s.translate(COMP)[::-1]


def add_errors(rng, frag):
    """Returns (seq, cigartuples, refpos per base of seq (-1 = insertion), NM)."""
    L = len(frag)
    k = rng.binomial(L - 2, ERR_RATE) if L > 4 else 0
    if k == 0:
        return frag, [(0, L)], np.arange(L), 0
    pos = np.sort(rng.choice(np.arange(1, L - 1), size=k, replace=False))
    typ = rng.choice(3, size=k, p=[.5, .25, .25])
    off = rng.integers(1, 4, size=k)
    ins = rng.integers(0, 4, size=k)
    idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    pieces, ops, rp = [], [], []
    prev = 0
    for p, t, o, i in zip(pos.tolist(), typ.tolist(), off.tolist(), ins.tolist()):
        if p > prev:
            pieces.append(frag[prev:p]); ops.append((0, p - prev)); rp.append(np.arange(prev, p))
        if t == 0:      # substitution
            pieces.append('ACGT'[(idx[frag[p]] + o) % 4]); ops.append((0, 1)); rp.append(np.array([p]))
        elif t == 1:    # deletion of reference base p
            ops.append((2, 1))
        else:           # insertion after base p
            pieces.append(frag[p] + 'ACGT'[i]); ops.append((0, 1)); ops.append((1, 1))
            rp.append(np.array([p, -1]))
        prev = p + 1
    if prev < L:
        pieces.append(frag[prev:]); ops.append((0, L - prev)); rp.append(np.arange(prev, L))
    merged = []
    for op, n in ops:
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + n)
        else:
            merged.append((op, n))
    return ''.join(pieces), merged, np.concatenate(rp), k


def cigar_from_rp(rp):
    """CIGAR (M/I/D) from the reference position of each read base (-1 = insertion)."""
    al = np.flatnonzero(rp >= 0)
    ev = [(int(q), 0, 1) for q in np.flatnonzero(rp < 0)]            # (index, order, length): I
    gap = np.diff(rp[al]) - 1
    for q0, g in zip(al[:-1][gap > 0].tolist(), gap[gap > 0].tolist()):
        ev.append((q0 + 0.5, 1, int(g)))                              # D after base q0
    ev.sort()
    ops, cur = [], 0
    def push(op, n):
        if n <= 0:
            return
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + n)
        else:
            ops.append((op, n))
    for q, kind, n in ev:
        if kind == 0:
            push(0, q - cur); push(1, 1); cur = q + 1
        else:
            push(0, int(q) + 1 - cur); push(2, n); cur = int(q) + 1
    push(0, len(rp) - cur)
    return ops


def methylation_tags(rng, seq, is_rev, refpos, chrom_i, chrom, mgmt, low_frac=0.5, mnp_targets=None):
    """MM/ML (5mC, CpG) compliant with the SAM spec: MM is expressed in the original
    orientation of the read (revcomp of SEQ for reads aligned to the - strand)."""
    L = len(seq)
    orig = revcomp(seq) if is_rev else seq
    arr = np.frombuffer(orig.encode(), dtype=np.uint8)
    c_idx = np.flatnonzero(arr == 67)
    if c_idx.size == 0:
        return None
    nxt = np.minimum(c_idx + 1, L - 1)
    is_cpg = (c_idx + 1 < L) & (arr[nxt] == 71)
    if not is_cpg.any():
        return None
    ranks = np.flatnonzero(is_cpg)
    skips = np.diff(np.concatenate([[-1], ranks])) - 1
    cpg_i = c_idx[is_cpg]
    seq_c = (L - 2 - cpg_i) if is_rev else cpg_i          # index of the C of the dinucleotide in SEQ
    seq_c = np.clip(seq_c, 0, L - 2)
    rc, rg = refpos[seq_c], refpos[seq_c + 1]
    valid = (rc >= 0) & (rg == rc + 1)
    # per-site methylation: deterministic bimodal background (hash of the position); MGMT methylated
    h = ((np.where(valid, rc, 0).astype(np.uint64) * np.uint64(2654435761)
          + np.uint64(chrom_i * 40503 + 12345)) & np.uint64(0xFFFFFFFF)).astype(np.float64) / 2**32
    beta = np.where(h < low_frac, 0.04, 0.90)
    beta = np.where(valid, beta, 0.5)
    if chrom == mgmt[0]:
        in_mgmt = valid & (rc >= mgmt[1]) & (rc < mgmt[2])
        beta = np.where(in_mgmt, MGMT_BETA, beta)
    if mnp_targets is not None:
        # HM450 probes (top_probes): real representative profile of the class (Capper et al.
        # betas, binarised) transplanted at the exact hg38 position of the probe -> the signal
        # carried by the BAM becomes that of a real sample of this tumour class.
        # Vectorised lookup (np.searchsorted): sorted positions prepared once per
        # chromosome in main(), no Python loop per CpG (~10x faster on a complete run).
        sub = mnp_targets.get(chrom)
        if sub is not None:
            pos_sorted, tgt_sorted = sub
            j = np.searchsorted(pos_sorted, rc)
            j = np.clip(j, 0, pos_sorted.size - 1)
            hit = valid & (pos_sorted[j] == rc)
            beta = np.where(hit, np.where(tgt_sorted[j] == 1, 0.92, 0.04), beta)
    meth = rng.random(beta.size) < beta
    ml = np.where(meth, rng.integers(200, 256, size=beta.size), rng.integers(0, 56, size=beta.size))
    # ML kept as a compact uint8 array (1 byte per CpG) until the read is written
    return "C+m," + ",".join(map(str, skips.tolist())) + ";", ml.astype(np.uint8)


# ---------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=type(
        'HelpFormatter', (argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter), {}))
    ap.add_argument('--profile', required=True, choices=list(PROFILES), help='tumour profile (required)')
    ap.add_argument('--ref', help='real hg38 FASTA (chr prefix or RefSeq accessions)')
    ap.add_argument('--standin', action='store_true', help='synthetic sequence (plumbing test)')
    ap.add_argument('--data-dir', required=True, help='data/ directory of Rapid-CNS2_nf')
    ap.add_argument('--out', default='sim_glioma', help='output folder')
    ap.add_argument('--sample-id', default=None, help='sample name, also the BAM file name; '
                    'default: SIM_<PROFILE>_HG38%(default).0s')
    ap.add_argument('--ontarget-depth', type=float, default=8.0,
                    help='depth per panel region before overlap; the effective panel coverage is ~2.5x higher (long reads)')
    ap.add_argument('--offtarget-depth', type=float, default=2.0, help='depth outside the panel (used for CNV calling)')
    ap.add_argument('--purity', type=float, default=0.80, help='tumour purity (heterozygous VAF = purity/2)')
    ap.add_argument('--cdkn2a-del', action='store_true', help='(astro) adds the CDKN2A/B homozygous deletion = grade 4')
    ap.add_argument('--chroms', default='', help='comma-separated list; default: chr1-22 and chrX (XX donor)%(default)s')
    ap.add_argument('--seed', type=int, default=2026, help='random seed (same seed and options = same BAM)')
    ap.add_argument('--skip-ref-check', action='store_true', help='skip the check of the hotspot reference bases (--ref)')
    ap.add_argument('--cnv-jitter-sd', type=float, default=0.0,
                    help="standard deviation of the multiplicative noise applied to each copy ratio "
                    "(trisomy 7, monosomy 10, EGFR amplification, deletions, 1p/19q codeletion) "
                    "to simulate inter-individual heterogeneity; 0 = disabled (previous "
                    "behaviour, fixed ratios per profile)")
    ap.add_argument('--mnp-targets', help='probe_targets_<profile>.tsv (01_export_probe_targets.R): '
                    'adds dedicated coverage and a real methylation signal at the ~100,000 HM450 probes '
                    "used by methylation_classification.R, for an interpretable classification")
    ap.add_argument('--mnp-depth', type=float, default=6.0, help='dedicated depth per HM450 probe')
    ap.add_argument('--mnp-window', type=int, default=150, help='half-window (bp) around each probe')
    ap.add_argument('--mgmt-depth-factor', type=float, default=4.0,
                    help='dedicated depth on the MGMT promoter = factor x --ontarget-depth, then '
                    'multiplied by the local copy ratio (0.5 with the monosomy 10 of gbm). Default 4 '
                    '-> ~18-22x on MGMT for gbm, above the --minimumMgmtCov threshold of Rapid-CNS2 '
                    '(5 by default, 10 in the launchers); the former fixed value 1.6 gave ~4-7x')
    a = ap.parse_args()
    if bool(a.ref) == bool(a.standin):
        sys.exit("Specify either --ref hg38.fa or --standin")
    if not 0 < a.purity <= 1:
        sys.exit(f"--purity must be in ]0, 1] (got {a.purity})")
    if a.chroms:
        unknown = [c for c in a.chroms.split(',') if c not in CHROMS]
        if unknown:
            sys.exit(f"--chroms: unknown chromosome(s) {unknown} (expected names such as chr1, chr10, chrX)")
    prof = PROFILES[a.profile]
    sample_id = a.sample_id or f'SIM_{a.profile.upper()}_HG38'

    ref = RealRef(a.ref) if a.ref else StandInRef(a.seed)
    panel, mgmt, events, chrom_ratio = build_model(a.data_dir, prof, a.cdkn2a_del)

    if a.cnv_jitter_sd > 0:
        # Reproducible multiplicative noise (derived from --seed): each individual of a
        # cohort gets a slightly different CNV amplitude, instead of a fixed ratio
        # identical for all individuals of the same profile.
        jrng = np.random.default_rng([a.seed, 777])

        def jitter(r):
            return float(r * np.clip(1 + jrng.normal(0, a.cnv_jitter_sd), 0.5, 1.5))

        chrom_ratio = {c: jitter(r) for c, r in chrom_ratio.items()}
        events = [(c, s, e, jitter(r), n) for c, s, e, r, n in events]
        print(f"CNV jitter (sd={a.cnv_jitter_sd}) applied: "
              + ", ".join(f"{c}={r:.2f}" for c, r in chrom_ratio.items())
              + (", " if chrom_ratio and events else "")
              + ", ".join(f"{n}={r:.2f}" for _, _, _, r, n in events))

    mnp_targets = None   # {chrom: {pos0: 0|1}}
    mnp_windows = []     # (chrom, start, end) dedicated coverage windows
    if a.mnp_targets:
        mnp_targets = {}
        n_probes = 0
        with open(a.mnp_targets) as f:
            next(f)  # header
            for line in f:
                probe, chrom, p0, p1, strand, target = line.rstrip('\n').split('\t')
                if chrom not in CHROMS or chrom in ('chrX', 'chrY', 'chrM'):
                    continue   # methylation_classification.R excludes chrX/chrY
                p0, p1, target = int(p0), int(p1), int(target)
                mnp_targets.setdefault(chrom, {})[p0] = target
                mnp_windows.append((chrom, max(0, p0 - a.mnp_window), p1 + a.mnp_window))
                n_probes += 1
                if isinstance(ref, StandInRef):
                    # The stand-in sequence deliberately depletes CpGs (see StandInRef);
                    # without this, the real hg38 coordinate of a probe has almost no
                    # chance of falling on a real "CG" in the synthetic sequence, and the
                    # methylation transplant can never apply. We therefore force the
                    # CG dinucleotide exactly at each targeted probe position.
                    ref.overrides.setdefault(chrom, {})[p0] = 'C'
                    ref.overrides.setdefault(chrom, {})[p0 + 1] = 'G'
        print(f"{n_probes} HM450 probes loaded from {a.mnp_targets} "
              f"(dedicated coverage {a.mnp_depth}x, window ±{a.mnp_window} bp)")
        for c in list(mnp_targets):
            d = mnp_targets[c]
            pos = np.array(sorted(d), dtype=np.int64)
            tgt = np.array([d[p] for p in pos.tolist()], dtype=np.int8)
            mnp_targets[c] = (pos, tgt)   # (sorted positions, 0/1 targets) for vectorised lookup
    snvs = {n: (a.purity / 2 if v == 'het' else v) for n, v in prof['snvs'].items()}
    by_chrom_snv = {}
    for n, v in snvs.items():
        by_chrom_snv.setdefault(HOTSPOTS[n]['chrom'], []).append((n, HOTSPOTS[n], v))
    del_site = None
    if prof['del_gene']:
        dc, dp, dnm = frameshift_site(a.data_dir, prof['del_gene'])
        del_site = dict(name=f"{prof['del_gene']}_frameshift_1bp_del", chrom=dc, pos0=dp, vaf=a.purity / 2, nm=dnm)
        print(f"Site {del_site['name']}: {dc}:{dp+1} ({dnm}, middle of the largest coding exon)")

    if a.ref and not a.skip_ref_check:
        for n in list(snvs) + prof['wt']:
            h = HOTSPOTS[n]
            b = ref.fetch(h['chrom'], h['pos0'], h['pos0'] + 1)
            if b != h['ref_plus']:
                sys.exit(f"ERROR: {n} {h['chrom']}:{h['pos0']+1} ref base={b}, expected {h['ref_plus']}. "
                         "Check the genome version (hg38) or use --skip-ref-check.")

    chroms = a.chroms.split(',') if a.chroms else [c for c in CHROMS if c not in ('chrM', 'chrY')]
    os.makedirs(a.out, exist_ok=True)
    bam_path = os.path.join(a.out, f'{sample_id}.bam')
    dorado = 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0'
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': c, 'LN': int(ref.lengths[c])} for c in CHROMS],
        'RG': [{'ID': 'sim_run', 'PL': 'ONT', 'SM': sample_id,
                'DS': f'basecall_model={dorado} modbase_models={dorado}_5mCG_5hmCG@v2'}],
        'CO': [f"SIMULATED DATA ({prof['label']}; simulate_glioma_hg38.py). Reference sequence: " +
               ('real hg38 from --ref' if a.ref else 'SYNTHETIC STAND-IN, NOT real hg38')],
    }
    out = pysam.AlignmentFile(bam_path, 'wb', header=header)
    tid = {c: i for i, c in enumerate(CHROMS)}

    # extra targets (dedicated coverage)
    extra = [(mgmt[0], mgmt[1], mgmt[2], a.ontarget_depth * a.mgmt_depth_factor, 'mgmt')]
    for n in list(snvs) + prof['wt']:
        h = HOTSPOTS[n]
        extra.append((h['chrom'], h['pos0'] - 300, h['pos0'] + 300, a.ontarget_depth, n))
    if del_site:
        extra.append((del_site['chrom'], del_site['pos0'] - 300, del_site['pos0'] + 300, a.ontarget_depth, 'del'))
    for c, s, e, r, n in events:
        if r > 1:  # background coverage in the amplified window
            extra.append((c, s, e, 1.0, n))

    if isinstance(ref, StandInRef):
        ref.finalize_overrides()

    mnp_by_chrom = {}
    for c, s, e in mnp_windows:
        mnp_by_chrom.setdefault(c, []).append((s, e))

    n_reads = 0
    for chrom in chroms:
        ci = CHROMS.index(chrom); L = int(ref.lengths[chrom])
        rng = np.random.default_rng([a.seed, ci])
        rng_m = np.random.default_rng([a.seed, ci, 99])
        jobs = []
        cuts = sorted({0, L, *[x for c, s, e, _, _ in events if c == chrom for x in (max(0, s), min(L, e))]})
        for s, e in zip(cuts[:-1], cuts[1:]):
            jobs.append((s, e, a.offtarget_depth * ratio_at(chrom, (s + e) // 2, events, chrom_ratio), False))
        for c, s, e, _ in panel:
            if c == chrom:
                jobs.append((s, e, a.ontarget_depth * ratio_at(chrom, (s + e) // 2, events, chrom_ratio), False))
        for c, s, e, d, _ in extra:
            if c == chrom:
                jobs.append((max(0, s), min(L, e), d * ratio_at(chrom, (s + e) // 2, events, chrom_ratio), False))
        for s, e in mnp_by_chrom.get(chrom, []):
            s, e = max(0, s), min(L, e)
            # narrow window (±mnp_window) -> dedicated short reads: same useful coverage
            # on the probe, far fewer bp "wasted" outside it (hence faster).
            jobs.append((s, e, a.mnp_depth * ratio_at(chrom, (s + e) // 2, events, chrom_ratio), True))

        recs = []
        for s, e, depth, short in jobs:
            if depth <= 0 or e <= s:
                continue
            mean_rl = 500 if short else MEAN_RL
            n = rng.poisson(depth * (e - s + mean_rl) / mean_rl)
            if n == 0:
                continue
            if short:
                lens = np.clip(rng.lognormal(6.1, 0.4, n), 250, 2000).astype(int)
            else:
                lens = np.clip(rng.lognormal(7.8, 0.7, n), 300, 40000).astype(int)
            starts = rng.integers(s - lens + 1, e)
            for st, rl in zip(starts.tolist(), lens.tolist()):
                st = max(0, st); en = min(L, st + rl)
                if en - st < 300:
                    continue
                frag = ref.fetch(chrom, st, en)
                if _NOT_ACGT.search(frag):
                    continue  # N or IUPAC ambiguity code (Y, R, W, S, K, M...): read skipped
                for _, h, vaf in by_chrom_snv.get(chrom, []):
                    if st <= h['pos0'] < en and rng.random() < vaf:
                        i = h['pos0'] - st
                        frag = frag[:i] + h['alt_plus'] + frag[i + 1:]
                didx = None
                if del_site and chrom == del_site['chrom'] and st + 20 < del_site['pos0'] < en - 20 \
                        and rng.random() < del_site['vaf']:
                    didx = del_site['pos0'] - st
                    frag = frag[:didx] + frag[didx + 1:]
                is_rev = bool(rng.random() < 0.5)
                seq, _, rp, nm = add_errors(rng, frag)
                if didx is not None:
                    rp = np.where(rp >= 0, rp + (rp >= didx), -1); nm += 1
                rp = np.where(rp >= 0, rp + st, -1)
                # MM/ML computed right away: the per-base reference position array (rp,
                # 8 bytes per base) is NOT kept for every read of the chromosome until the
                # sort, which divided the peak memory by several (several GB per job on
                # chr1, a problem with build_cohort.py --jobs 16).
                mm = methylation_tags(rng_m, seq, is_rev, rp, ci, chrom, mgmt, prof['low_frac'], mnp_targets)
                recs.append((st, seq, cigar_from_rp(rp), nm, is_rev, mm))
        recs.sort(key=lambda r: r[0])
        for st, seq, cigar, nm, is_rev, mm in recs:
            n_reads += 1
            r = pysam.AlignedSegment(out.header)
            r.query_name = f'sim-{ci:02d}-{n_reads:07d}'
            r.flag = 16 if is_rev else 0
            r.reference_id = tid[chrom]
            r.reference_start = st
            r.mapping_quality = 60
            r.cigartuples = cigar
            r.query_sequence = seq
            r.query_qualities = pysam.qualitystring_to_array('5' * len(seq))
            r.set_tag('NM', int(nm))
            r.set_tag('RG', 'sim_run')
            if mm:
                r.set_tag('MM', mm[0], value_type='Z')
                r.set_tag('ML', array('B', mm[1].tobytes()))   # B:C array, as written by dorado
            out.write(r)
        print(f'{chrom}: {len(recs)} reads', flush=True)
    out.close()
    pysam.index(bam_path)

    bg = prof['low_frac'] * 0.04 + (1 - prof['low_frac']) * 0.90
    with open(os.path.join(a.out, 'truth_events.tsv'), 'w') as f:
        f.write('event\ttype\tchrom\tstart\tend\tvalue\n')
        for c, r in chrom_ratio.items():
            f.write(f'{c}_copy_ratio\tcopy_ratio\t{c}\t0\t{HG38_LEN[c]}\t{r}\n')
        for c, s, e, r, n in events:
            f.write(f'{n}\tcopy_ratio\t{c}\t{s}\t{e}\t{r}\n')
        f.write(f'MGMT_promoter_methylated\tmethylation\t{mgmt[0]}\t{mgmt[1]}\t{mgmt[2]}\t{MGMT_BETA}\n')
        f.write(f'background_methylated_fraction\tbackground\t.\t0\t0\t{bg:.3f}\n')
        for n, v in snvs.items():
            h = HOTSPOTS[n]
            f.write(f"{n}\tsnv\t{h['chrom']}\t{h['pos0']}\t{h['pos0']+1}\t{h['ref_plus']}>{h['alt_plus']};vaf={v:.3f}\n")
        for n in prof['wt']:
            h = HOTSPOTS[n]
            f.write(f"{n}_wildtype\twildtype\t{h['chrom']}\t{h['pos0']}\t{h['pos0']+1}\t{h['ref_plus']}\n")
        if del_site:
            f.write(f"{del_site['name']}\tdel\t{del_site['chrom']}\t{del_site['pos0']}\t{del_site['pos0']+1}\tvaf={del_site['vaf']:.3f}\n")
    print(f'Done: {bam_path} ({n_reads} reads)')


if __name__ == '__main__':
    main()
