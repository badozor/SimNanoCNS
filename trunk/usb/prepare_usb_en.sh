#!/usr/bin/env bash
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# =============================================================================
# prepare_usb_en.sh: builds a self-contained USB stick for Rapid-CNS2_nf
#
# The stick holds: the pipeline, Nextflow (+ Java), the containers, the hg38
# reference, ANNOVAR (+ the 5 databases of the default protocol), the AnnotSV
# annotations, the graphical interface and launchers (run.sh, run_gui.sh).
# Nextflow's "work" directory stays on the host's disk.
#
# The graphical interface (rapidcns2_gui_en.py) is looked for next to this script,
# then in ../rapidcns2/ (layout of the NanoPoreSeq project); or set GUI_PY.
#
# Usage: ./prepare_usb_en.sh [options] /media/$USER/MY_STICK
# Run with --dry-run first to check the sizes.
# =============================================================================
set -euo pipefail

# ----------------------------------------------------------------------------
# Sources (P3Ultra machine). Edit here or override with environment variables.
# ----------------------------------------------------------------------------
PIPELINE_DIR=${PIPELINE_DIR:-/export/scratch3/gitProject/Rapid-CNS2_nf}
REF_FA=${REF_FA:-/export/scratch3/data/Annotations/UCSC_HG38/hg38.fa}
ANNOVAR_DIR=${ANNOVAR_DIR:-/export/scratch3/data/Annotations/humandb/annovar}
ANNOVAR_DB=${ANNOVAR_DB:-$ANNOVAR_DIR/humandb}
ANNOTSV_DIR=${ANNOTSV_DIR:-/export/scratch3/gitProject/AnnotSV/share/AnnotSV}
SCRIPT_DIR=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
if [ -z "${GUI_PY:-}" ]; then
    GUI_PY=$SCRIPT_DIR/rapidcns2_gui_en.py
    [ -f "$GUI_PY" ] || GUI_PY=$SCRIPT_DIR/../rapidcns2/rapidcns2_gui_en.py
fi

# ANNOVAR databases of the default protocol (refGene,cytoBand,clinvar_20240917,avsnp151,1000g2015aug_eur)
ANNOVAR_REQUIRED=(hg38_refGene.txt hg38_refGeneMrna.fa hg38_cytoBand.txt
                  hg38_clinvar_20240917.txt hg38_avsnp151.txt hg38_EUR.sites.2015_08.txt)
ANNOVAR_OPTIONAL=(hg38_refGeneVersion.txt hg38_clinvar_20240917.txt.idx
                  hg38_avsnp151.txt.idx hg38_EUR.sites.2015_08.txt.idx)

# Fallback images if they cannot be read from the pipeline (README v3.0.x)
FALLBACK_IMAGES=(areebapatel/rapid_cns:3.0.2
                 quay.io/biocontainers/ont-modkit:0.6.4--h7f49ad2_0
                 hkubal/clair3:v2.0.2
                 quay.io/biocontainers/severus:1.7--pyhdfd78af_0
                 quay.io/biocontainers/savana:1.3.8--pyhdfd78af_0)

# ----------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $0 [options] USB_STICK_MOUNT_POINT

Options:
  --engine apptainer|docker|both  container format on the stick
                                  (default: apptainer if installed, otherwise docker)
  --no-annovar                    do not copy ANNOVAR (stick handed to a third party:
                                  ANNOVAR may not be redistributed)
  --with-samples DIR              also copy the BAMs (+ .bai) from DIR into samples/
  --no-java                       do not bundle Java (the host must have Java >= 17)
  --dry-run                       show what would be copied and the sizes, without writing
  -h, --help                      this help

The stick must be formatted as ext4 (or xfs/btrfs), not exFAT/FAT32.
EOF
}

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33mWARNING: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }
human() { numfmt --to=iec --suffix=B "$1" 2>/dev/null || echo "$1 B"; }
sz()    { du -sbL "$@" 2>/dev/null | awk '{s+=$1} END {print s+0}'; }
# File name Nextflow uses in its Singularity/Apptainer cache
sif_name() { local n=${1#docker://}; n=${n//\//-}; n=${n//:/-}; echo "$n.img"; }

ENGINE=""; WITH_ANNOVAR=1; SAMPLES_DIR=""; WITH_JAVA=1; DRY=0; USB=""
while [ $# -gt 0 ]; do
    case "$1" in
        --engine)       [ $# -ge 2 ] || die "--engine needs a value"; ENGINE=$2; shift 2 ;;
        --no-annovar)   WITH_ANNOVAR=0; shift ;;
        --with-samples) [ $# -ge 2 ] || die "--with-samples needs a directory"
                        SAMPLES_DIR=$(readlink -f "$2"); shift 2 ;;
        --no-java)      WITH_JAVA=0; shift ;;
        --dry-run)      DRY=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        -*)             die "unknown option: $1 (see --help)" ;;
        *)              USB=$(readlink -f "$1"); shift ;;
    esac
done
[ -n "$USB" ] || { usage; exit 1; }
R="$USB/rapidcns2"

# ----------------------------------------------------------------------------
say "Checks"
# ----------------------------------------------------------------------------
[ -d "$USB" ] || die "$USB does not exist (stick not mounted?)"
FSTYPE=$(findmnt -no FSTYPE --target "$USB" 2>/dev/null || echo unknown)
case "$FSTYPE" in
    ext4|ext3|xfs|btrfs|f2fs) echo "File system: $FSTYPE ✔" ;;
    *) if [ "${RCNS_SKIP_FS_CHECK:-0}" = 1 ]; then warn "file system $FSTYPE accepted (RCNS_SKIP_FS_CHECK=1)"
       else die "the stick is formatted as '$FSTYPE'. Nextflow needs symbolic links and Unix permissions,
       and FAT32 limits files to 4 GB. Reformat it as ext4 (ERASES ALL ITS CONTENT):
         lsblk                                # identify the stick, e.g. /dev/sdX1
         sudo umount /dev/sdX1
         sudo mkfs.ext4 -m 0 -L RAPIDCNS2 /dev/sdX1
         # replug the stick, then: sudo chown \$USER: /media/\$USER/RAPIDCNS2"
       fi ;;
esac
[ "$DRY" = 1 ] || [ -w "$USB" ] || die "no write permission on $USB: sudo chown $USER: $USB"
[ -z "$SAMPLES_DIR" ] || [ -d "$SAMPLES_DIR" ] || die "--with-samples: directory not found: $SAMPLES_DIR"

[ -f "$PIPELINE_DIR/main.nf" ] || die "pipeline not found: $PIPELINE_DIR/main.nf"
[ -f "$REF_FA" ]               || die "reference not found: $REF_FA"
[ -f "$REF_FA.fai" ]           || die "missing index: samtools faidx $REF_FA"
[ -d "$ANNOTSV_DIR/Annotations_Human" ] || die "Annotations_Human missing from $ANNOTSV_DIR"
[ -f "$GUI_PY" ]               || warn "graphical interface not found ($GUI_PY): it will not be copied"
command -v nextflow >/dev/null || die "nextflow not in the PATH"
command -v rsync    >/dev/null || die "rsync required (sudo apt install rsync)"

if [ "$WITH_ANNOVAR" = 1 ]; then
    for f in "${ANNOVAR_REQUIRED[@]}"; do
        [ -f "$ANNOVAR_DB/$f" ] || die "ANNOVAR database missing or still compressed: $ANNOVAR_DB/$f"
    done
    [ -f "$ANNOVAR_DIR/table_annovar.pl" ] || die "table_annovar.pl not found in $ANNOVAR_DIR"
fi

if [ -z "$ENGINE" ]; then
    if command -v apptainer >/dev/null || command -v singularity >/dev/null; then ENGINE=apptainer
    elif command -v docker >/dev/null; then ENGINE=docker
    else die "neither Apptainer/Singularity nor Docker on this machine"; fi
fi
case "$ENGINE" in apptainer|docker|both) ;; *) die "--engine: apptainer, docker or both" ;; esac
APPT=$(command -v apptainer || command -v singularity || true)
if [ "$ENGINE" = apptainer ] || [ "$ENGINE" = both ]; then
    [ -n "$APPT" ] || die "apptainer/singularity not found to build the .sif images"
fi
if [ "$ENGINE" = docker ] || [ "$ENGINE" = both ]; then
    command -v docker >/dev/null || die "docker not found"
fi
echo "Container format: $ENGINE"

# ----------------------------------------------------------------------------
say "Container images used by the pipeline"
# ----------------------------------------------------------------------------
mapfile -t IMAGES < <(
    grep -rhoE "container[[:space:]]*=?[[:space:]]*['\"][^'\"]+['\"]" \
        "$PIPELINE_DIR"/nextflow.config "$PIPELINE_DIR"/*.nf \
        "$PIPELINE_DIR"/conf "$PIPELINE_DIR"/nextflow 2>/dev/null \
    | sed -E "s/.*['\"]([^'\"]+)['\"]$/\1/" | sed 's#^docker://##' \
    | grep -v '[$]' | grep -E '^[a-zA-Z0-9]' | sort -u )
if [ "${#IMAGES[@]}" -eq 0 ]; then
    warn "no image read from the pipeline, using the README list"
    IMAGES=("${FALLBACK_IMAGES[@]}")
fi
if grep -rqE "container[[:space:]]*=?[[:space:]]*['\"][^'\"]*[$]" \
        "$PIPELINE_DIR"/nextflow.config "$PIPELINE_DIR"/conf "$PIPELINE_DIR"/nextflow 2>/dev/null; then
    warn "some images are defined through a variable and could not be listed: check with run.sh --check"
fi
printf '  %s\n' "${IMAGES[@]}"

# ----------------------------------------------------------------------------
say "Estimating the required space"
# ----------------------------------------------------------------------------
S_PIPE=$(sz "$PIPELINE_DIR")
S_REF=$(sz "$REF_FA"*)
S_ASV=$(sz "$ANNOTSV_DIR")
S_ANV=0
if [ "$WITH_ANNOVAR" = 1 ]; then
    files=(); for f in "${ANNOVAR_REQUIRED[@]}" "${ANNOVAR_OPTIONAL[@]}"; do
        [ -f "$ANNOVAR_DB/$f" ] && files+=("$ANNOVAR_DB/$f"); done
    S_ANV=$(sz "${files[@]}" "$ANNOVAR_DIR"/*.pl)
fi
S_SMP=0; [ -n "$SAMPLES_DIR" ] && S_SMP=$(sz "$SAMPLES_DIR")
S_IMG=0
if command -v docker >/dev/null; then
    for i in "${IMAGES[@]}"; do
        s=$(docker image inspect -f '{{.Size}}' "$i" 2>/dev/null || echo 2000000000)
        S_IMG=$((S_IMG + s))
    done
else
    S_IMG=$(( ${#IMAGES[@]} * 2000000000 ))   # ~2 GB per image, estimate
fi
[ "$ENGINE" = both ] && S_IMG=$((S_IMG * 2))
S_JAVA=$(( WITH_JAVA * 400000000 ))
TOTAL=$((S_PIPE + S_REF + S_ASV + S_ANV + S_SMP + S_IMG + S_JAVA + 500000000))
AVAIL=$(df -B1 --output=avail "$USB" | tail -1 | tr -d ' ')
printf '  %-28s %10s\n' "Pipeline" "$(human "$S_PIPE")" "hg38 reference" "$(human "$S_REF")" \
       "ANNOVAR" "$(human "$S_ANV")" "AnnotSV" "$(human "$S_ASV")" \
       "Containers (approx.)" "$(human "$S_IMG")" "Java + Nextflow" "$(human "$S_JAVA")" \
       "Samples" "$(human "$S_SMP")"
printf '  %-28s %10s\n' "TOTAL (approx.)" "$(human "$TOTAL")" "Free on the stick" "$(human "$AVAIL")"
[ "$TOTAL" -lt "$AVAIL" ] || warn "free space looks insufficient"
if [ "$DRY" = 1 ]; then echo; echo "--dry-run: nothing was written."; exit 0; fi
[ "$TOTAL" -lt "$AVAIL" ] || die "not enough space on the stick"

mkdir -p "$R"/{bin,pipeline,refs,annotsv,containers/sif,containers/docker,nxf_home,conf,gui,samples,results}

# ----------------------------------------------------------------------------
say "Nextflow and Java"
# ----------------------------------------------------------------------------
NF_VER=$(nextflow -version 2>/dev/null | grep -oE 'version [0-9]+\.[0-9]+\.[0-9]+' | awk '{print $2}' | head -1 || true)
echo "Nextflow version on this machine: ${NF_VER:-unknown}"
# The self-contained binary (all dependencies included) is ideal offline. It is called
# "-dist" in recent releases and "-all" in older ones: both are tried.
NF_OK=0
if [ -n "$NF_VER" ]; then
    for suffix in dist all; do
        if curl -fsSL -o "$R/bin/nextflow.tmp" \
                "https://github.com/nextflow-io/nextflow/releases/download/v${NF_VER}/nextflow-${NF_VER}-${suffix}"; then
            mv "$R/bin/nextflow.tmp" "$R/bin/nextflow"
            echo "nextflow-${NF_VER}-${suffix} downloaded ✔"
            NF_OK=1; break
        fi
    done
fi
if [ "$NF_OK" = 0 ]; then
    rm -f "$R/bin/nextflow.tmp"
    warn "could not download the self-contained binary: copying the local launcher"
    cp "$(readlink -f "$(command -v nextflow)")" "$R/bin/nextflow"
fi
chmod +x "$R/bin/nextflow"

if [ "$WITH_JAVA" = 1 ]; then
    JAVA_BIN=$(readlink -f "$(command -v java)" 2>/dev/null || true)
    if [ -n "$JAVA_BIN" ]; then
        JHOME=$(dirname "$(dirname "$JAVA_BIN")")
        rsync -a --delete "$JHOME/" "$R/jre/"
        echo "Java copied from $JHOME ✔"
    else
        warn "java not found: the host will have to provide Java >= 17"
    fi
fi

# ----------------------------------------------------------------------------
say "Pipeline"
# ----------------------------------------------------------------------------
rsync -a --delete --exclude 'work/' --exclude '.nextflow*' --exclude 'results*/' \
      "$PIPELINE_DIR/" "$R/pipeline/"
PIPE_COMMIT=$(git -C "$PIPELINE_DIR" rev-parse --short HEAD 2>/dev/null || echo "n/a")
echo "Rapid-CNS2_nf copied (commit $PIPE_COMMIT) ✔"

# Nextflow plugins the pipeline may declare, cached on the stick
mapfile -t PLUGINS < <(grep -hoE "id[[:space:]]+['\"]nf-[^'\"]+['\"]" "$R/pipeline/nextflow.config" 2>/dev/null \
                       | sed -E "s/.*['\"](nf-[^'\"]+)['\"]/\1/" | sort -u)
(
    export NXF_HOME="$R/nxf_home"
    [ -x "$R/jre/bin/java" ] && export JAVA_HOME="$R/jre"
    "$R/bin/nextflow" -version >/dev/null 2>&1 || warn "first Nextflow launch failed"
    for p in "${PLUGINS[@]}"; do
        "$R/bin/nextflow" plugin install "$p" && echo "plugin $p ✔" || warn "plugin $p not installed"
    done
)

# ----------------------------------------------------------------------------
say "Containers"
# ----------------------------------------------------------------------------
: > "$R/containers/images.tsv"
for img in "${IMAGES[@]}"; do
    if [ "$ENGINE" = apptainer ] || [ "$ENGINE" = both ]; then
        out="$R/containers/sif/$(sif_name "$img")"
        if [ -s "$out" ]; then echo "already present: $(basename "$out")"
        else
            echo "building $(basename "$out")"
            if command -v docker >/dev/null && docker image inspect "$img" >/dev/null 2>&1; then
                "$APPT" build "$out.tmp" "docker-daemon://$img" || "$APPT" pull "$out.tmp" "docker://$img"
            else
                "$APPT" pull "$out.tmp" "docker://$img"
            fi
            mv "$out.tmp" "$out"
        fi
    fi
    if [ "$ENGINE" = docker ] || [ "$ENGINE" = both ]; then
        tgz="$(sif_name "$img")"; tgz="${tgz%.img}.tar.gz"
        if [ -s "$R/containers/docker/$tgz" ]; then echo "already present: $tgz"
        else
            docker image inspect "$img" >/dev/null 2>&1 || docker pull "$img"
            echo "exporting $img"
            docker save "$img" | gzip -1 > "$R/containers/docker/$tgz.tmp"
            mv "$R/containers/docker/$tgz.tmp" "$R/containers/docker/$tgz"
        fi
        printf '%s\t%s\n' "$img" "$tgz" >> "$R/containers/images.tsv"
    fi
done

# ----------------------------------------------------------------------------
say "Reference, AnnotSV, ANNOVAR"
# ----------------------------------------------------------------------------
REF_BASE=$(basename "$REF_FA")
rsync -aL --info=progress2 "$REF_FA" "$REF_FA".* "$R/refs/" 2>/dev/null || rsync -aL "$REF_FA"* "$R/refs/"
rsync -a --info=progress2 "$ANNOTSV_DIR/" "$R/annotsv/"
if [ "$WITH_ANNOVAR" = 1 ]; then
    mkdir -p "$R/annovar/humandb"
    rsync -a "$ANNOVAR_DIR"/*.pl "$R/annovar/"
    for f in "${ANNOVAR_REQUIRED[@]}" "${ANNOVAR_OPTIONAL[@]}"; do
        [ -f "$ANNOVAR_DB/$f" ] && rsync -aL --info=progress2 "$ANNOVAR_DB/$f" "$R/annovar/humandb/"
    done
else
    echo "ANNOVAR not copied (--no-annovar): users will have to point to their own installation"
fi

# ----------------------------------------------------------------------------
say "Graphical interface and samples"
# ----------------------------------------------------------------------------
[ -f "$GUI_PY" ] && cp "$GUI_PY" "$R/gui/rapidcns2_gui.py"
if [ -n "$SAMPLES_DIR" ]; then
    rsync -a --info=progress2 --include='*/' --include='*.bam' --include='*.bai' --exclude='*' \
          --prune-empty-dirs "$SAMPLES_DIR/" "$R/samples/"
fi

# ----------------------------------------------------------------------------
say "Launchers and configuration"
# ----------------------------------------------------------------------------
cat > "$R/env.sh" <<EOF
# shellcheck shell=bash
# env.sh: shared environment, sourced by run.sh and run_gui.sh
R=\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)
export RCNS_USB="\$R"
export PATH="\$R/bin:\$PATH"
if [ -x "\$R/jre/bin/java" ]; then export JAVA_HOME="\$R/jre"; export PATH="\$JAVA_HOME/bin:\$PATH"; fi
export NXF_HOME="\$R/nxf_home"
export NXF_OFFLINE=true
export NXF_DISABLE_CHECK_LATEST=true
# Nextflow >= 26.04 uses the strict syntax by default, not yet followed by Rapid-CNS2_nf 3.0:
# legacy parser unless already set (NXF_SYNTAX_PARSER=v2 once the pipeline is migrated)
export NXF_SYNTAX_PARSER="\${NXF_SYNTAX_PARSER:-v1}"
export NXF_SINGULARITY_CACHEDIR="\$R/containers/sif"
export NXF_APPTAINER_CACHEDIR="\$R/containers/sif"

# Container engine: can be forced with RCNS_ENGINE=singularity|docker
if [ -z "\${RCNS_ENGINE:-}" ]; then
    if { command -v apptainer || command -v singularity; } >/dev/null 2>&1 \\
       && compgen -G "\$R/containers/sif/*.img" >/dev/null; then RCNS_ENGINE=singularity
    elif command -v docker >/dev/null 2>&1 && compgen -G "\$R/containers/docker/*.tar.gz" >/dev/null; then
        RCNS_ENGINE=docker
    fi
fi
export RCNS_ENGINE="\${RCNS_ENGINE:-none}"

# Values read by the graphical interface (RCNS_ prefix + field name in upper case)
export RCNS_PIPELINE="\$R/pipeline/main.nf"
export RCNS_COHORT="\${RCNS_COHORT:-\$R/samples}"
export RCNS_RUNDIR="\${RCNS_WORKDIR:-\$HOME/rapidcns2_work}"   # work directory: host disk
export RCNS_OUTDIR="\${RCNS_OUTDIR:-\$R/results/{id}}"
export RCNS_CONFIG="\$R/conf/usb.config"
export RCNS_REF="\$R/refs/$REF_BASE"
export RCNS_ANNOVARPATH="\${RCNS_ANNOVARPATH:-\$R/annovar}"
export RCNS_ANNOVARDB="\${RCNS_ANNOVARDB:-\$R/annovar/humandb}"
export RCNS_ANNOTSVANNOT="\$R/annotsv"
export RCNS_PROFILE="local,\$RCNS_ENGINE"
export RCNS_BINDPATHS="\$R"
# wf-human-variation starts a nested Nextflow run that downloads its workflow and
# containers: impossible offline, so it is disabled by default on the stick.
export RCNS_HUMANVAR="\${RCNS_HUMANVAR:-false}"
export RCNS_MAXCPUS="\${RCNS_MAXCPUS:-\$(nproc)}"
export RCNS_MAXMEMORY="\${RCNS_MAXMEMORY:-\$(awk '/MemTotal/ {printf "%d", \$2/1048576*0.8}' /proc/meminfo)}"

rcns_check_engine() {
    if [ "\$RCNS_ENGINE" = none ]; then
        echo "ERROR: no usable container engine (Apptainer/Singularity or Docker)" >&2
        echo "       for the containers stored on the stick." >&2
        return 1
    fi
}

# Loads into Docker the stick's images that are missing on the host
rcns_load_docker() {
    [ "\$RCNS_ENGINE" = docker ] || return 0
    while IFS=\$'\\t' read -r img tgz; do
        [ -n "\$img" ] || continue
        if ! docker image inspect "\$img" >/dev/null 2>&1; then
            echo "Loading \$img into Docker…"
            gunzip -c "\$R/containers/docker/\$tgz" | docker load
        fi
    done < "\$R/containers/images.tsv"
}
EOF

cat > "$R/conf/usb.config" <<'EOF'
// usb.config: stick-specific configuration (loaded with -c by run.sh and the interface)
// The stick is mounted into the containers through --containerBindPaths (added by
// run.sh and the interface), which nextflow.config turns into -B or -v.
def usb = System.getenv('RCNS_USB')

singularity.cacheDir = "${usb}/containers/sif"
apptainer.cacheDir   = "${usb}/containers/sif"

// Steps that fail on simulated BAMs (no SV, no heterozygous SNP).
// COMMENT OUT for real data, so that errors do not go unnoticed.
process {
    withName: 'annotSv.*' { errorStrategy = 'ignore' }
    withName: 'savanaCna' { errorStrategy = 'ignore' }
}
EOF

cat > "$R/run.sh" <<'EOF'
#!/usr/bin/env bash
# run.sh: runs Rapid-CNS2 from the stick
#   ./run.sh --check                      check the installation
#   ./run.sh BAM [ID] [nextflow options]  analyse one sample
set -euo pipefail
source "$(cd "$(dirname "$(readlink -f "$0")")" && pwd)/env.sh"

if [ "${1:-}" = --check ]; then
    ok=1
    chk() { if eval "$2" >/dev/null 2>&1; then echo "  ✔ $1"; else echo "  ✘ $1"; ok=0; fi; }
    echo "Stick: $R"
    echo "Container engine: $RCNS_ENGINE"
    chk "usable container engine"          "[ '$RCNS_ENGINE' != none ]"
    chk "Java"                             "java -version"
    chk "Nextflow (offline)"               "nextflow -version"
    chk "pipeline"                         "[ -f '$RCNS_PIPELINE' ]"
    chk "reference + index"                "[ -f '$RCNS_REF' ] && [ -f '$RCNS_REF.fai' ]"
    chk "ANNOVAR (table_annovar.pl)"       "[ -f '$RCNS_ANNOVARPATH/table_annovar.pl' ]"
    chk "ANNOVAR databases"                "[ -f '$RCNS_ANNOVARDB/hg38_refGene.txt' ]"
    chk "AnnotSV"                          "[ -d '$RCNS_ANNOTSVANNOT/Annotations_Human' ]"
    chk "work directory writable"          "mkdir -p '$RCNS_RUNDIR'"
    if [ -d "$RCNS_RUNDIR" ]; then
        echo "Work directory: $RCNS_RUNDIR ($(df -h --output=avail "$RCNS_RUNDIR" | tail -1 | tr -d ' ') free; allow ~150 GB/sample)"
    fi
    echo "Resources: $RCNS_MAXCPUS CPUs, $RCNS_MAXMEMORY GB"
    echo "Samples found: $(find "$RCNS_COHORT" -maxdepth 2 -name '*.bam' 2>/dev/null | wc -l)"
    [ "$ok" = 1 ] && echo "Everything is in place." || { echo "Some items are missing."; exit 1; }
    exit 0
fi

[ $# -ge 1 ] || { sed -n '2,5p' "$0"; exit 1; }
BAM=$(readlink -f "$1"); shift
ID=$(basename "$BAM" .bam)
if [ $# -ge 1 ] && [ "${1#-}" = "$1" ]; then ID=$1; shift; fi
[ -f "$BAM" ] || { echo "BAM not found: $BAM" >&2; exit 1; }

rcns_check_engine
rcns_load_docker
mkdir -p "$RCNS_RUNDIR"; cd "$RCNS_RUNDIR"
OUT=${RCNS_OUTDIR//\{id\}/$ID}

exec nextflow run "$RCNS_PIPELINE" -ansi-log false \
    --input "$BAM" --id "$ID" --outDir "$OUT" \
    --ref "$RCNS_REF" \
    --annovarPath "$RCNS_ANNOVARPATH" --annovarDB "$RCNS_ANNOVARDB" \
    --annotsvAnnot "$RCNS_ANNOTSVANNOT" \
    --containerBindPaths "$RCNS_BINDPATHS" \
    --minimumMgmtCov 10 --mnpFlex true --runHumanVariation "$RCNS_HUMANVAR" \
    --savanaCnaArgs '--overrule_cellularity 0.8' \
    --maxCpus "$RCNS_MAXCPUS" --maxMemory "$RCNS_MAXMEMORY" \
    -profile "$RCNS_PROFILE" -c "$RCNS_CONFIG" -resume "$@"
EOF

cat > "$R/run_gui.sh" <<'EOF'
#!/usr/bin/env bash
# run_gui.sh: opens the graphical interface, preconfigured for the stick
set -euo pipefail
source "$(cd "$(dirname "$(readlink -f "$0")")" && pwd)/env.sh"
[ -f "$R/gui/rapidcns2_gui.py" ] || { echo "graphical interface missing from the stick: use ./run.sh" >&2; exit 1; }
python3 -c 'import tkinter' 2>/dev/null || { echo "tkinter missing: sudo apt install python3-tk" >&2; exit 1; }
rcns_check_engine
rcns_load_docker
mkdir -p "$RCNS_RUNDIR"; cd "$RCNS_RUNDIR"
exec python3 "$R/gui/rapidcns2_gui.py"
EOF
chmod +x "$R/run.sh" "$R/run_gui.sh"

cat > "$R/README.txt" <<EOF
Rapid-CNS2 on a USB stick
=========================
Built on $(date '+%Y-%m-%d') on $(hostname)
Pipeline: Rapid-CNS2_nf, commit $PIPE_COMMIT
Nextflow: ${NF_VER:-?}    Containers: $ENGINE
ANNOVAR included: $([ "$WITH_ANNOVAR" = 1 ] && echo yes || echo no)

Requirements on the host computer
  - Linux x86_64
  - Apptainer/Singularity (recommended) or Docker
  - python3-tk for the graphical interface
  - ~150 GB free for the work directory (default \$HOME/rapidcns2_work,
    can be changed: export RCNS_WORKDIR=/path/to/disk)

Usage
  ./run.sh --check                     check the installation on this host
  ./run_gui.sh                         graphical interface
  ./run.sh samples/X/X.bam X           command line
Results: results/<id>/   Report: results/<id>/report/<id>_Rapid-CNS2_report_full.html

Notes
  - wf-human-variation is disabled: it downloads its own workflow and
    containers, which is impossible offline (with network: export RCNS_HUMANVAR=true).
  - Offline, automatic upload to MNP-Flex is not possible;
    the file mnpflex/<id>.MNPFlex.input.bed can be uploaded by hand.
  - conf/usb.config ignores AnnotSV and SAVANA failures (simulated BAMs):
    comment these lines out for real data.
  - ANNOVAR is free for academic and non-profit use and may not be
    redistributed: do not hand a stick containing it to a third party.
EOF

sync
say "Done"
echo "Stick ready in $R ($(du -sh "$R" | cut -f1))."
echo "Check: $R/run.sh --check"
