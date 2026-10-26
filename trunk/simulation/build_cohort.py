#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""
Builds a simulated cohort of N individuals (one BAM per individual), split as
evenly as possible between the gbm / astro / oligo profiles (33/33/33 for N=99).

Each individual gets:
  - a unique seed (reproducible: seed_base + index)
  - a tumour purity drawn at random within --purity-range (inter-individual variability)
  - if --cnv-jitter-sd > 0: a CNV amplitude (trisomy 7, monosomy 10, EGFR amplification,
    deletions, 1p/19q codeletion) slightly different per individual, instead of the fixed
    ratio of the profile
  - if --astro-grade4-fraction > 0: a fraction of the "astro" individuals additionally gets
    --cdkn2a-del (CDKN2A/B homozygous deletion, grade 4 criterion per WHO 2021) AND a
    methylation profile transplanted from a real A_IDH_HG sample (rather than A_IDH) -
    methylation/CNV consistency, not just an isolated deletion added
  - if --mnp-targets-dir is given (outputs of 02_export_cohort_probe_targets.R): the methylation
    profile of a real, DISTINCT sample of the Capper cohort for its class, hence real
    inter-individual diversity rather than copies of the same representative

Parallelised over --jobs processes for the generation (each individual is independent): useful
since generating a single individual with dedicated MNP probes takes several minutes.

With --classify: then classifies each BAM (bedMethyl + methylation_classification_local.R,
real Capper et al. data) and computes the correct classification rate of the cohort.
Separate step, parallelised over --classify-jobs (kept low by default: each classification
can use several GB of RAM).

All the scripts called here (simulate_glioma_hg38.py, bedmethyl_from_bam.py,
methylation_classification_local.R) must be in the same directory as this one.

Usage (real hg38, recommended, with realistic classification, 16 cores):
  Rscript 02_export_cohort_probe_targets.R Rapid-CNS2_nf/data cohort_probes 33
  python3 build_cohort.py --n 99 --ref /path/hg38.fa --data-dir Rapid-CNS2_nf/data \\
      --out-dir cohort_bams --mnp-targets-dir cohort_probes --jobs 16 \\
      --cnv-jitter-sd 0.15 --astro-grade4-fraction 0.3 \\
      --classify --classify-jobs 2 --num-trees 2000

Usage (plumbing test, stand-in sequence, fast, without MNP probes):
  python3 build_cohort.py --n 9 --standin --data-dir Rapid-CNS2_nf/data \\
      --out-dir cohort_bams_test --jobs 3
"""
import argparse
import csv
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROFILES = ['gbm', 'astro', 'oligo']
EXPECTED_DX = {'gbm': 'GBM_RTK_II', 'astro': 'A_IDH', 'oligo': 'O_IDH'}
GRADE4_DX = 'A_IDH_HG'


def parse_profile_counts(spec, n):
    """'gbm=50,astro=30,oligo=19' -> {'gbm': 50, 'astro': 30, 'oligo': 19}; None if spec is empty."""
    if not spec:
        return None
    counts = {p: 0 for p in PROFILES}
    for item in spec.split(','):
        name, _, val = item.partition('=')
        name = name.strip()
        if name not in counts or not val.strip().isdigit():
            sys.exit(f"--profile-counts: invalid item '{item}' (expected e.g. gbm=50,astro=30,oligo=19)")
        counts[name] = int(val)
    if sum(counts.values()) != n:
        sys.exit(f"--profile-counts: the sum of the group sizes ({sum(counts.values())}) must equal --n ({n})")
    return counts


def build_assignment(n, seed, fixed_counts=None):
    """Splits n individuals between the 3 profiles as evenly as possible, shuffled order."""
    base = n // 3
    counts = {p: base for p in PROFILES}
    for p in PROFILES[:n - base * 3]:
        counts[p] += 1  # distributes the remainder (non-zero only if n % 3 != 0)
    if fixed_counts:
        counts = dict(fixed_counts)
    assignment = []
    per_profile_idx = {p: 0 for p in PROFILES}
    for p in PROFILES:
        for _ in range(counts[p]):
            per_profile_idx[p] += 1
            assignment.append((p, per_profile_idx[p]))
    rng = random.Random(seed)
    rng.shuffle(assignment)
    return assignment, counts


# ---------------------------------------------------------------------- step 1: generation
def run_one(args_tuple):
    (idx, profile, profile_idx, a) = args_tuple
    ind_id = f"IND{idx:03d}_{profile.upper()}"
    out_dir = os.path.join(a.out_dir, ind_id)
    seed = a.seed_base + idx
    purity = round(random.Random(seed).uniform(*a.purity_range), 3)
    grade4 = (profile == 'astro' and random.Random(seed + 500_000).random() < a.astro_grade4_fraction)
    expected_dx = GRADE4_DX if grade4 else EXPECTED_DX[profile]

    cmd = [sys.executable, os.path.join(SCRIPT_DIR, 'simulate_glioma_hg38.py'),
           '--profile', profile, '--data-dir', a.data_dir, '--out', out_dir,
           '--sample-id', ind_id, '--seed', str(seed), '--purity', str(purity)]
    cmd += ['--ref', a.ref] if a.ref else ['--standin']
    if a.chroms:
        cmd += ['--chroms', a.chroms]
    if a.cnv_jitter_sd > 0:
        cmd += ['--cnv-jitter-sd', str(a.cnv_jitter_sd)]
    if grade4:
        cmd += ['--cdkn2a-del']

    mnp_file = None
    if a.mnp_targets_dir:
        mnp_profile = 'astro_hg' if grade4 else profile
        mnp_file = os.path.join(a.mnp_targets_dir, f'probe_targets_{mnp_profile}_{profile_idx:02d}.tsv')
        if os.path.exists(mnp_file):
            cmd += ['--mnp-targets', mnp_file, '--mnp-depth', str(a.mnp_depth)]
        else:
            mnp_file = None  # file missing (n_per_class too small in the R export): continue without it

    t0 = time.time()
    log_path = os.path.join(a.out_dir, f'{ind_id}.log')
    os.makedirs(a.out_dir, exist_ok=True)
    with open(log_path, 'w') as logf:
        rc = subprocess.run(cmd, stdout=logf, stderr=subprocess.STDOUT).returncode
    elapsed = time.time() - t0
    bam_path = os.path.join(out_dir, f'{ind_id}.bam')
    ok = rc == 0 and os.path.exists(bam_path)
    return dict(individual=ind_id, profile=profile, profile_index=profile_idx, seed=seed,
                purity=purity, grade4=grade4, expected_dx=expected_dx, mnp_file=mnp_file or '',
                bam=bam_path if ok else '', status='OK' if ok else f'FAILED(rc={rc}, see {log_path})',
                elapsed_s=round(elapsed, 1))


# ------------------------------------------------------------------ step 2: classification
def classify_one(args_tuple):
    (row, a) = args_tuple
    ind = row['individual']
    out_dir = os.path.dirname(row['bam'])
    bedmethyl = os.path.join(out_dir, f'{ind}.bedmethyl.bed')
    classif_dir = os.path.join(a.classify_out_dir, ind)
    os.makedirs(classif_dir, exist_ok=True)
    log_path = os.path.join(classif_dir, f'{ind}_classify.log')

    with open(log_path, 'w') as logf:
        rc1 = subprocess.run([sys.executable, os.path.join(SCRIPT_DIR, 'bedmethyl_from_bam.py'),
                               '--bam', row['bam'], '--out', bedmethyl],
                              stdout=logf, stderr=subprocess.STDOUT).returncode
        if rc1 != 0:
            return dict(individual=ind, **_classify_fail(row, f'bedmethyl_from_bam.py rc={rc1}'))
        rc2 = subprocess.run(
            ['Rscript', os.path.join(SCRIPT_DIR, 'methylation_classification_local.R'),
             '-s', ind, '-o', classif_dir, '-i', bedmethyl,
             '-p', os.path.join(a.data_dir, 'top_probes_hm450.Rdata'),
             '-a', os.path.join(a.data_dir, 'HM450.hg38.manifest.gencode.v22.Rdata'),
             '-d', os.path.join(a.data_dir, 'capper_top_100k_betas_binarised.Rdata'),
             '-t', str(a.classify_threads), '-n', str(a.num_trees)],
            stdout=logf, stderr=subprocess.STDOUT).returncode

    result_path = os.path.join(classif_dir, f'{ind}_calibrated_classification.tsv')
    if rc2 != 0 or not os.path.exists(result_path):
        return dict(individual=ind, **_classify_fail(row, f'methylation_classification_local.R rc={rc2}, see {log_path}'))

    fields = {}
    for line in open(result_path):
        if ':' not in line:
            continue          # blank or unexpected line
        k, v = line.rstrip('\n').split(':', 1)
        fields[k.strip()] = v.strip()
    predicted = fields.get('Predicted Class', '')
    return dict(individual=ind, profile=row['profile'], expected_dx=row['expected_dx'],
                predicted_dx=predicted, correct=(predicted == row['expected_dx']),
                calibrated_score=fields.get('Calibrated Score', ''),
                n_features=fields.get('Number of features', ''), status='OK')


def _classify_fail(row, reason):
    return dict(profile=row['profile'], expected_dx=row['expected_dx'], predicted_dx='',
                correct=False, calibrated_score='', n_features='', status=f'FAILED({reason})')


# --------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=type(
        'HelpFormatter', (argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter), {}))
    ap.add_argument('--n', type=int, default=99, help='total number of individuals (default 99)')
    ap.add_argument('--profile-counts', default='',
                    help="group sizes, e.g. gbm=50,astro=30,oligo=19 (sum = --n); "
                         "omitted = equal split%(default)s")
    ap.add_argument('--ref', help='real hg38 FASTA (chr prefix or RefSeq accessions)')
    ap.add_argument('--standin', action='store_true', help='synthetic sequence (plumbing test)')
    ap.add_argument('--data-dir', required=True, help='data/ directory of Rapid-CNS2_nf')
    ap.add_argument('--out-dir', default='cohort_bams', help='output folder (one subfolder per individual)')
    ap.add_argument('--mnp-targets-dir', help='outputs of 02_export_cohort_probe_targets.R '
                    '(real methylation diversity); omitted = random background as before')
    ap.add_argument('--mnp-depth', type=float, default=6.0, help='dedicated depth per HM450 probe')
    ap.add_argument('--chroms', default='', help='restrict the chromosomes, comma-separated (quick tests); default: chr1-22 and chrX%(default)s')
    ap.add_argument('--purity-range', type=float, nargs=2, default=[0.60, 0.95],
                    metavar=('MIN', 'MAX'), help='tumour purity drawn at random per individual')
    ap.add_argument('--cnv-jitter-sd', type=float, default=0.0,
                    help='inter-individual variability of the copy ratios (0 = fixed profile ratios)')
    ap.add_argument('--astro-grade4-fraction', type=float, default=0.0,
                    help='fraction of astro individuals getting --cdkn2a-del + A_IDH_HG methylation '
                    '(WHO 2021 grade 4 subgroup; 0 = disabled, all astro remain A_IDH)')
    ap.add_argument('--seed-base', type=int, default=2026, help='seed of individual i = seed-base + i')
    ap.add_argument('--jobs', type=int, default=4, help='individuals generated in parallel')
    ap.add_argument('--classify', action='store_true',
                    help='classifies each BAM after generation (bedMethyl + real random forest)')
    ap.add_argument('--classify-out-dir', default='cohort_classification', help='classification output folder')
    ap.add_argument('--classify-jobs', type=int, default=2,
                    help='parallel classifications (each can use several GB of RAM: keep low)')
    ap.add_argument('--num-trees', type=int, default=2000,
                    help='random forest trees (20000 in production = much more accurate but '
                    'much slower; reduced value here for cohort screening)')
    ap.add_argument('--classify-threads', type=int, default=2, help='R threads per individual')
    a = ap.parse_args()
    if bool(a.ref) == bool(a.standin):
        sys.exit("Specify either --ref hg38.fa or --standin")
    lo, hi = a.purity_range
    if not 0 < lo <= hi <= 1:
        sys.exit(f"--purity-range: expected 0 < MIN <= MAX <= 1 (got {lo} {hi})")
    if a.n < 1 or a.jobs < 1 or a.classify_jobs < 1:
        sys.exit("--n, --jobs and --classify-jobs must be >= 1")
    if not 0 <= a.astro_grade4_fraction <= 1:
        sys.exit("--astro-grade4-fraction must be between 0 and 1")
    if not a.mnp_targets_dir and (a.classify or a.astro_grade4_fraction > 0):
        sys.exit("--classify and --astro-grade4-fraction require --mnp-targets-dir "
                 "(real methylation signal; without it the classification or the "
                 "A_IDH/A_IDH_HG distinction is meaningless)")

    assignment, counts = build_assignment(a.n, a.seed_base, parse_profile_counts(a.profile_counts, a.n))
    print(f"Cohort of {a.n} individuals: " + ", ".join(f"{p}={c}" for p, c in counts.items()))
    if a.mnp_targets_dir:
        needed = [p for p in PROFILES if counts.get(p)]
        if a.astro_grade4_fraction > 0 and counts.get('astro'):
            needed.append('astro_hg')
        missing = [p for p in needed if not os.path.exists(
            os.path.join(a.mnp_targets_dir, f'probe_targets_{p}_01.tsv'))]
        if missing:
            print(f"WARNING: no probe file found for {missing} in {a.mnp_targets_dir} "
                  "- run 02_export_cohort_probe_targets.R first (with astro_hg if "
                  "--astro-grade4-fraction > 0). These individuals will be generated without a real signal.")

    os.makedirs(a.out_dir, exist_ok=True)
    jobs = [(idx + 1, p, pidx, a) for idx, (p, pidx) in enumerate(assignment)]

    print(f"\n== Generating the {len(jobs)} BAMs ({a.jobs} in parallel) ==")
    results = []
    t_start = time.time()
    with ProcessPoolExecutor(max_workers=a.jobs) as pool:
        futures = {pool.submit(run_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futures), 1):
            r = fut.result()
            results.append(r)
            g4 = ' [grade4]' if r['grade4'] else ''
            print(f"[{i}/{len(jobs)}] {r['individual']}{g4} -> {r['status']} ({r['elapsed_s']}s)", flush=True)

    manifest_path = os.path.join(a.out_dir, 'cohort_manifest.tsv')
    with open(manifest_path, 'w', newline='') as f:
        fields = ['individual', 'profile', 'profile_index', 'seed', 'purity', 'grade4',
                  'expected_dx', 'mnp_file', 'bam', 'status', 'elapsed_s']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in sorted(results, key=lambda r: r['individual']):
            w.writerow(r)
    n_ok = sum(1 for r in results if r['status'] == 'OK')
    print(f"\n{n_ok}/{len(jobs)} individuals generated successfully in "
          f"{(time.time() - t_start) / 60:.1f} min. Manifest: {manifest_path}")
    if n_ok < len(jobs):
        print("See the individual .log files in --out-dir for the failures.")

    if not a.classify:
        return

    print(f"\n== Classifying the {n_ok} generated BAMs ({a.classify_jobs} in parallel, "
          f"{a.num_trees} trees) ==")
    os.makedirs(a.classify_out_dir, exist_ok=True)
    ok_rows = [r for r in results if r['status'] == 'OK']
    t_start = time.time()
    class_results = []
    with ProcessPoolExecutor(max_workers=a.classify_jobs) as pool:
        futures = {pool.submit(classify_one, (r, a)): r for r in ok_rows}
        for i, fut in enumerate(as_completed(futures), 1):
            cr = fut.result()
            class_results.append(cr)
            mark = '✓' if cr.get('correct') else ('✗' if cr['status'] == 'OK' else '!')
            print(f"[{i}/{len(ok_rows)}] {cr['individual']} expected={cr['expected_dx']} "
                  f"predicted={cr.get('predicted_dx','')} {mark} ({cr['status']})", flush=True)

    summary_path = os.path.join(a.classify_out_dir, 'cohort_classification_summary.tsv')
    with open(summary_path, 'w', newline='') as f:
        fields = ['individual', 'profile', 'expected_dx', 'predicted_dx', 'correct',
                  'calibrated_score', 'n_features', 'status']
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in sorted(class_results, key=lambda r: r['individual']):
            w.writerow(r)

    n_done = sum(1 for r in class_results if r['status'] == 'OK')
    n_correct = sum(1 for r in class_results if r.get('correct'))
    print(f"\n{n_done}/{len(ok_rows)} classifications completed in "
          f"{(time.time() - t_start) / 60:.1f} min.")
    if n_done:
        print(f"Correct classification rate: {n_correct}/{n_done} = {100*n_correct/n_done:.1f}%")
    print(f"Summary: {summary_path}")


if __name__ == '__main__':
    main()
