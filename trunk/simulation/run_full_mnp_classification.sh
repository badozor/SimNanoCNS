#!/bin/bash
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# Complete pipeline: simulation -> bedMethyl -> local methylation classification (random forest,
# real Capper et al. data), for a given profile of simulate_glioma_hg38.py.
#
# Usage:
#   ./run_full_mnp_classification.sh --profile gbm|astro|oligo --data-dir Rapid-CNS2_nf/data \
#       (--ref hg38.fa | --standin) [--num-trees 20000] [--threads 16]
#
# Examples:
#   ./run_full_mnp_classification.sh --profile astro --data-dir Rapid-CNS2_nf/data \
#       --ref /path/to/hg38.fa --num-trees 20000 --threads 16
#   ./run_full_mnp_classification.sh --profile gbm --data-dir Rapid-CNS2_nf/data --standin
#
# Outputs (in the current directory): sim_<profile>_mnp/ (BAM, bedMethyl, truth_events.tsv)
# and classif_<profile>/ (classification). The probe files go to <script dir>/probe_targets/.
#
# MEMORY WARNING: with the complete Capper dataset (2801 samples x up to
# 10,000 probes), the classification step can use more than 4 GB of RAM even with a
# reduced number of trees. On a machine with limited memory, provide swap:
#   fallocate -l 4G /swapfile && chmod 600 /swapfile && mkswap /swapfile && swapon /swapfile
set -euo pipefail

PROFILE="" DATA_DIR="" REF="" STANDIN=0 NUM_TREES=20000 THREADS=16
while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)   PROFILE="$2"; shift 2 ;;
        --data-dir)  DATA_DIR="$2"; shift 2 ;;
        --ref)       REF="$2"; shift 2 ;;
        --standin)   STANDIN=1; shift ;;
        --num-trees) NUM_TREES="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        -h|--help)   sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "Unknown argument: $1 (see --help)"; exit 1 ;;
    esac
done

if [[ -z "$PROFILE" || -z "$DATA_DIR" ]]; then
    echo "Error: --profile and --data-dir are required (see --help)"; exit 1
fi
case "$PROFILE" in
    gbm|astro|oligo) ;;
    *) echo "Error: --profile must be gbm, astro or oligo (got '$PROFILE')"; exit 1 ;;
esac
if [[ "$STANDIN" -eq 1 && -n "$REF" ]]; then
    echo "Error: --ref and --standin are mutually exclusive"; exit 1
fi
if [[ "$STANDIN" -eq 0 && -z "$REF" ]]; then
    echo "Error: specify --ref hg38.fa or --standin"; exit 1
fi
if [[ "$STANDIN" -eq 1 ]]; then REF_ARGS=(--standin); else REF_ARGS=(--ref "$REF"); fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="sim_${PROFILE}_mnp"
SAMPLE="SIM_${PROFILE^^}_HG38"

echo "== 1/4: export of the real representative methylation profile (R) =="
mkdir -p "$SCRIPT_DIR/probe_targets"
Rscript "$SCRIPT_DIR/01_export_probe_targets.R" "$DATA_DIR" "$SCRIPT_DIR/probe_targets"

echo "== 2/4: BAM simulation with methylation transplant at the HM450 probes =="
python3 "$SCRIPT_DIR/simulate_glioma_hg38.py" --profile "$PROFILE" "${REF_ARGS[@]}" --data-dir "$DATA_DIR" \
    --out "$OUT" --mnp-targets "$SCRIPT_DIR/probe_targets/probe_targets_${PROFILE}.tsv"

echo "== 3/4: bedMethyl from the BAM =="
python3 "$SCRIPT_DIR/bedmethyl_from_bam.py" --bam "$OUT/${SAMPLE}.bam" --out "$OUT/${SAMPLE}.bedmethyl.bed"

echo "== 4/4: classification (random forest, $NUM_TREES trees, $THREADS threads) =="
mkdir -p "classif_${PROFILE}"
Rscript "$SCRIPT_DIR/methylation_classification_local.R" \
    -s "$SAMPLE" -o "classif_${PROFILE}" \
    -i "$OUT/${SAMPLE}.bedmethyl.bed" \
    -p "$DATA_DIR/top_probes_hm450.Rdata" \
    -a "$DATA_DIR/HM450.hg38.manifest.gencode.v22.Rdata" \
    -d "$DATA_DIR/capper_top_100k_betas_binarised.Rdata" \
    -t "$THREADS" -n "$NUM_TREES"

echo "== Result =="
cat "classif_${PROFILE}/${SAMPLE}_calibrated_classification.tsv"
