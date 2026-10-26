#!/usr/bin/env bash
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# =============================================================================
# align_fastq_rapidcns2_en.sh
#
# Preparation of Nanopore data (FASTQ) for the Rapid-CNS2 pipeline:
#   1. merge the FASTQ files (a single file or a folder of .fastq / .fastq.gz)
#   2. check for methylation tags (MM/ML) in the FASTQ headers
#   3. align to hg38 with dorado aligner (minimap2, CPU)
#   4. sort + index (samtools)
#   5. checks: integrity, statistics, contig names, MM/ML tags
#
# Output: <outdir>/<sample>.sorted.bam (+ .bai), ready to be given to the pipeline.
#
# Usage:
#   ./align_fastq_rapidcns2_en.sh -i <fastq|folder> -r <hg38.fa> -s <sample> \
#                                 [-o <outdir>] [-t <threads>] [-m "<minimap2 options>"]
#
# Example:
#   ./align_fastq_rapidcns2_en.sh -i fastq_pass/ -r hg38.fa -s GBM01 -o results -t 24
#
# Dependencies: dorado (aligner subcommand), samtools (>= 1.16), gzip; bgzip if available.
# =============================================================================

set -euo pipefail

# ----------------------------- parameters ------------------------------------
INPUT=""
REF=""
SAMPLE=""
OUTDIR="rapidcns2_input"
THREADS=16
MM2_OPTS=""          # e.g. "-x lr:hq"; empty = dorado's default preset

usage() { sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while getopts "i:r:s:o:t:m:h" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MM2_OPTS="$OPTARG" ;;
    *) usage ;;
  esac
done

log()  { echo "[$(date '+%F %T')] $*" >&2; }
die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

# ----------------------------- sanity checks ---------------------------------
[[ -n "$INPUT" && -n "$REF" && -n "$SAMPLE" ]] || usage
[[ -e "$INPUT" ]] || die "input not found: $INPUT"
[[ -f "$REF"   ]] || die "reference not found: $REF"
[[ "$THREADS" =~ ^[0-9]+$ && "$THREADS" -ge 1 ]] || die "-t must be a positive integer"

for tool in dorado samtools gzip; do
  command -v "$tool" >/dev/null 2>&1 || die "$tool is not in the PATH"
done

mkdir -p "$OUTDIR/tmp"
MERGED="$OUTDIR/${SAMPLE}.merged.fastq.gz"
BAM="$OUTDIR/${SAMPLE}.sorted.bam"
QC="$OUTDIR/${SAMPLE}.qc.txt"

log "dorado   : $(dorado --version 2>&1 | head -n 1)"
log "samtools : $(samtools --version | head -n 1)"

# ----------------------------- 1. merge the FASTQ files ----------------------
# Single file: used as is. Folder: all FASTQ files are merged
# (gzip -cdf reads compressed and uncompressed files alike).
if [[ -d "$INPUT" ]]; then
  mapfile -t FASTQS < <(find "$INPUT" -type f \
      \( -name '*.fastq' -o -name '*.fq' -o -name '*.fastq.gz' -o -name '*.fq.gz' \) | sort)
  (( ${#FASTQS[@]} > 0 )) || die "no FASTQ file found in $INPUT"
  log "Merging ${#FASTQS[@]} FASTQ files -> $MERGED"
  if command -v bgzip >/dev/null 2>&1; then
    gzip -cdf "${FASTQS[@]}" | bgzip -@ "$THREADS" > "$MERGED"
  else
    gzip -cdf "${FASTQS[@]}" | gzip -c > "$MERGED"
  fi
  READS="$MERGED"
else
  READS="$INPUT"
fi

# ----------------------------- 2. methylation tags ---------------------------
# Rapid-CNS2 needs methylation (classification, MGMT). In a FASTQ it only
# exists if the MM/ML tags are in the headers (FASTQ produced by
# "samtools fastq -T '*'" or equivalent). A standard MinKNOW FASTQ has none:
# in that case, start again from the POD5 files or the unaligned BAM.
HAS_MM=$( (set +o pipefail; gzip -cdf "$READS" | head -n 4000 | awk 'NR%4==1' | grep -c 'MM:Z:') || true )
if [[ "${HAS_MM:-0}" -eq 0 ]]; then
  warn "no MM:Z: tag in the FASTQ headers."
  warn "The BAM will NOT contain methylation: the methylation modules of"
  warn "Rapid-CNS2 (classification, MGMT) will not be able to run."
  warn "Only CNV / variants / SV will be usable."
else
  log "MM/ML tags detected in the FASTQ headers (methylation kept)."
fi

# ----------------------------- 3-4. alignment + sort -------------------------
SORT_THREADS=$(( THREADS / 4 > 1 ? THREADS / 4 : 1 ))
ALIGN_ARGS=( "$REF" "$READS" -t "$THREADS" )
if [[ -n "$MM2_OPTS" ]]; then ALIGN_ARGS+=( --mm2-opts "$MM2_OPTS" ); fi

log "Alignment (dorado aligner) then sort (samtools sort) -> $BAM"
dorado aligner "${ALIGN_ARGS[@]}" \
  | samtools sort -@ "$SORT_THREADS" -m 2G -T "$OUTDIR/tmp/${SAMPLE}" -o "$BAM" -

log "Indexing"
samtools index -@ "$SORT_THREADS" "$BAM"

# ----------------------------- 5. checks -------------------------------------
log "Checks -> $QC"
samtools quickcheck "$BAM" || die "corrupted or truncated BAM: $BAM"

{
  echo "== Sample    : $SAMPLE"
  echo "== Reference : $REF"
  echo
  echo "== samtools flagstat"
  samtools flagstat -@ "$SORT_THREADS" "$BAM"
  echo
  echo "== Aligned reads per chromosome (chr1-22, X, Y)"
  samtools idxstats "$BAM" | awk '$1 ~ /^(chr)?([0-9]+|X|Y)$/ {print $1"\t"$3}'
  echo
  echo "== Reads carrying MM tags (among the first 2000 alignments)"
  (set +o pipefail; samtools view "$BAM" | head -n 2000 | grep -c 'MM:Z:') || true
} > "$QC"

# The panel BED files and hg38 annotations use the "chr" prefix.
# (contig list captured first: "grep -q" closing the pipe early would trip pipefail)
CONTIGS=$(samtools idxstats "$BAM" | cut -f1)
if ! grep -qx 'chr7' <<< "$CONTIGS"; then
  warn "contig 'chr7' missing: the reference does not use the 'chr' prefix."
  warn "The hg38 BED/annotation files of the pipeline will not match."
fi

rmdir "$OUTDIR/tmp" 2>/dev/null || true
# The merged FASTQ is no longer needed; uncomment to delete it:
# [[ "$READS" == "$MERGED" ]] && rm -f "$MERGED"

log "Done."
log "BAM: $BAM"
log "QC : $QC"
