#!/usr/bin/env bash
# =============================================================================
# make_usb_stick_en.sh: creates a Rapid-CNS2 USB stick from start to finish
#
#   1. checks that the device is a removable USB disk, never a system disk
#   2. asks for confirmation, then ERASES it: one GPT partition formatted as ext4
#      (or only reformats the partition given, or keeps it with --keep-format)
#   3. mounts it and makes the current user its owner
#   4. runs prepare_usb_en.sh (same folder) to copy everything onto it
#   5. runs the stick's run.sh --check, then optionally ejects the stick
#
# Usage:
#   ./make_usb_stick_en.sh --list
#   ./make_usb_stick_en.sh [options] DEVICE [-- prepare_usb_en.sh options]
#
# DEVICE: whole disk (/dev/sdX: new partition table, ONE partition)
#         or partition (/dev/sdX1: only this partition is reformatted)
#         or mount folder of a stick already formatted (/media/$USER/RAPIDCNS2:
#         nothing is erased, as with --keep-format)
#
# Examples:
#   ./make_usb_stick_en.sh --list
#   ./make_usb_stick_en.sh --dry-run /dev/sdb
#   ./make_usb_stick_en.sh /dev/sdb -- --engine apptainer --with-samples ../cohort_bams
#   ./make_usb_stick_en.sh --keep-format --eject /dev/sdb1 -- --no-annovar
#
# Run it as your normal user: sudo is called only for the privileged steps.
# Needs: lsblk, findmnt, wipefs, mkfs.ext4, parted (or sgdisk), sudo;
#        udisksctl is used to mount and eject when available.
#
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
# =============================================================================
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
PREPARE="$SCRIPT_DIR/prepare_usb_en.sh"

LABEL="RAPIDCNS2"; KEEP=0; YES=0; ALLOW_NONREM=0; DRY=0; EJECT=0; LIST=0; MP_OPT=""; DEV_ARG=""
PREP_ARGS=()

usage() { sed -n '2,31p' "$0" | sed 's/^# \{0,1\}//'; cat <<EOF
Options:
  --list                  list the disks and their removable/USB status, then exit
  --label NAME            file system label (default: RAPIDCNS2, at most 16 characters)
  --keep-format           keep the existing ext4/xfs/btrfs file system of the partition
  --yes                   no interactive confirmation (the safety checks still apply)
  --allow-non-removable   accept a disk that the system does not report as removable/USB
  --mountpoint DIR        mount point if udisksctl is not available
                          (default: /media/\$USER/LABEL)
  --eject                 unmount and power off the stick at the end
  --dry-run               run the checks and show the actions, without changing anything
  -h, --help              this help
Everything after "--" is passed to prepare_usb_en.sh (e.g. --engine, --no-annovar,
--with-samples DIR, --no-java).
EOF
}

say()  { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33mWARNING: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }
# runs a command, or only shows it with --dry-run
run()  { if [ "$DRY" = 1 ]; then printf '  [dry-run] %s\n' "$*"; else "$@"; fi; }
# name of the first partition: /dev/sdb -> /dev/sdb1, /dev/nvme0n1 or /dev/loop0 -> ...p1
part1() { case "$1" in *[0-9]) echo "${1}p1" ;; *) echo "${1}1" ;; esac; }

while [ $# -gt 0 ]; do
    case "$1" in
        --list)                LIST=1; shift ;;
        --label)               [ $# -ge 2 ] || die "--label needs a value"; LABEL=$2; shift 2 ;;
        --keep-format)         KEEP=1; shift ;;
        --yes)                 YES=1; shift ;;
        --allow-non-removable) ALLOW_NONREM=1; shift ;;
        --mountpoint)          [ $# -ge 2 ] || die "--mountpoint needs a value"; MP_OPT=$2; shift 2 ;;
        --eject)               EJECT=1; shift ;;
        --dry-run)             DRY=1; shift ;;
        -h|--help)             usage; exit 0 ;;
        --)                    shift; PREP_ARGS=("$@"); break ;;
        -*)                    die "unknown option: $1 (see --help)" ;;
        *)                     [ -z "$DEV_ARG" ] || die "only one device can be given"; DEV_ARG=$1; shift ;;
    esac
done

if [ "$LIST" = 1 ]; then
    echo "Disks (TRAN=usb or RM=1: removable; MOUNTPOINT: currently mounted):"
    lsblk -o NAME,SIZE,TYPE,FSTYPE,LABEL,MODEL,TRAN,RM,MOUNTPOINT
    exit 0
fi
[ -n "$DEV_ARG" ] || { usage; exit 1; }
[[ "$LABEL" =~ ^[A-Za-z0-9_-]{1,16}$ ]] || die "--label: 1 to 16 characters among letters, digits, _ and -"
[ -x "$PREPARE" ] || [ -f "$PREPARE" ] || die "prepare_usb_en.sh not found next to this script ($PREPARE)"

# ---------------------------------------------------------------------------- user and sudo
if [ "$(id -u)" = 0 ]; then
    [ -z "${SUDO_USER:-}" ] || die "run this script as your normal user, not with sudo: \
it calls sudo itself, and prepare_usb_en.sh must run as you (Docker/Apptainer, Nextflow)"
    SUDO=(); OWNER=root
else
    command -v sudo >/dev/null || die "sudo is required for formatting and mounting"
    SUDO=(sudo); OWNER=$(id -un)
fi
LOG="make_usb_stick_$(date +%Y%m%d_%H%M%S).log"
[ "$DRY" = 1 ] || exec > >(tee -a "$LOG") 2>&1

# ------------------------------------------------------------------------------- checks
say "Device checks"
for tool in lsblk findmnt wipefs mkfs.ext4; do
    command -v "$tool" >/dev/null || die "$tool not found (package util-linux / e2fsprogs)"
done
DEV=$(readlink -f "$DEV_ARG")
if [ -d "$DEV" ]; then
    # a folder: accepted only if it is the mount point of a partition, which is then kept
    SRC=$(findmnt -nro SOURCE --mountpoint "$DEV" 2>/dev/null | head -1 || true)
    [ -n "$SRC" ] || die "$DEV_ARG is a folder, not a mount point: give the stick device \
(e.g. /dev/sdb, see --list) or the folder where the stick is mounted"
    echo "$DEV_ARG is the mount point of $SRC: this partition is kept, nothing is erased"
    [ -n "$MP_OPT" ] || MP_OPT=$DEV
    DEV=$(readlink -f "$SRC"); KEEP=1
fi
[ -b "$DEV" ] || die "$DEV_ARG is neither a block device nor a mount point (see --list)"
TYPE=$(lsblk -dno TYPE "$DEV")
case "$TYPE" in
    disk|loop) DISK=$DEV; PART="" ;;
    part)      DISK=/dev/$(lsblk -no PKNAME "$DEV" | head -1); PART=$DEV ;;
    *)         die "$DEV is of type '$TYPE': give a whole disk (/dev/sdX) or a partition (/dev/sdX1)" ;;
esac
[ -b "$DISK" ] || die "cannot find the disk holding $DEV"
if [ "$KEEP" = 1 ] && [ -z "$PART" ]; then
    die "--keep-format needs a partition (e.g. $(part1 "$DISK")), not a whole disk"
fi

# 1) never a disk used by the system: any mount outside /media, /run/media, /mnt and the
#    --mountpoint folder, or any swap, means the disk belongs to the system
MP_ABS=""; [ -z "$MP_OPT" ] || MP_ABS=$(readlink -m "$MP_OPT")
while read -r name mp; do
    [ -n "$mp" ] || continue
    mp=$(printf '%b' "$mp")                     # lsblk -r escapes spaces as \x20
    case "$mp" in
        /media/*|/run/media/*|/mnt|/mnt/*) continue ;;
    esac
    [ -n "$MP_ABS" ] && { [ "$mp" = "$MP_ABS" ] || [ "${mp#"$MP_ABS"/}" != "$mp" ]; } && continue
    die "$name is mounted on '$mp': $DISK is used by the system, it will not be touched"
done < <(lsblk -nrpo NAME,MOUNTPOINT "$DISK")

# 2) a removable / USB disk, unless explicitly overridden
RM=$(lsblk -dno RM "$DISK" | tr -d ' ')
TRAN=$(lsblk -dno TRAN "$DISK" | tr -d ' ')
if [ "$TRAN" != usb ] && [ "$RM" != 1 ]; then
    [ "$ALLOW_NONREM" = 1 ] || die "$DISK is neither USB nor removable (TRAN='${TRAN:-?}', RM=$RM). \
If you are sure, add --allow-non-removable"
    warn "$DISK is neither USB nor removable: accepted (--allow-non-removable)"
fi

SIZE_B=$(lsblk -dbno SIZE "$DISK" | tr -d ' ')
SIZE_H=$(lsblk -dno SIZE "$DISK" | tr -d ' ')
[ "$SIZE_B" -ge $((64 * 1024 ** 3)) ] || warn "$DISK is only $SIZE_H: a complete stick needs ~100 GB or more \
(run prepare_usb_en.sh --dry-run to know the exact size)"
lsblk -o NAME,SIZE,TYPE,FSTYPE,LABEL,MODEL,TRAN,MOUNTPOINT "$DISK"

if [ "$KEEP" = 1 ]; then
    # file system type: lsblk (udev), else findmnt if mounted, else blkid
    FST=$(lsblk -dno FSTYPE "$PART" | tr -d ' ')
    [ -n "$FST" ] || FST=$(findmnt -nro FSTYPE "$PART" 2>/dev/null | head -1 || true)
    [ -n "$FST" ] || FST=$("${SUDO[@]}" blkid -o value -s TYPE "$PART" 2>/dev/null || true)
    case "$FST" in ext4|ext3|xfs|btrfs|f2fs) echo "Existing file system kept: $FST" ;;
        *) die "$PART is formatted as '${FST:-none}': Nextflow needs ext4/xfs/btrfs (drop --keep-format)" ;;
    esac
    ACTION="keep $PART ($FST) and add Rapid-CNS2 to it"
elif [ -n "$PART" ]; then
    ACTION="ERASE the partition $PART and format it as ext4 (label $LABEL)"
else
    ACTION="ERASE the whole disk $DISK ($SIZE_H): one GPT partition, ext4, label $LABEL"
fi

# --------------------------------------------------------------------------- confirmation
say "Planned action: $ACTION"
if [ "$DRY" = 0 ] && [ "$KEEP" = 0 ] && [ "$YES" = 0 ]; then
    [ -t 0 ] || die "no terminal to confirm: add --yes to run without confirmation"
    printf 'ALL DATA on %s will be lost. Type the device name (%s) to confirm: ' "$DISK" "$(basename "$DISK")"
    read -r answer
    [ "$answer" = "$(basename "$DISK")" ] || die "confirmation not given: nothing was changed"
fi
[ "$DRY" = 1 ] || "${SUDO[@]}" true    # asks for the sudo password once, before the long steps

# ------------------------------------------------------------------- unmount the stick
say "Unmounting"
while read -r name mp; do
    [ -n "$mp" ] || continue
    if command -v udisksctl >/dev/null && [ "$DRY" = 0 ] && udisksctl unmount -b "$name" >/dev/null 2>&1; then
        echo "  $name unmounted"
    else
        run "${SUDO[@]}" umount "$name"
    fi
done < <(lsblk -nrpo NAME,MOUNTPOINT "$DISK" | if [ "$KEEP" = 1 ]; then grep -v "^$PART "; else cat; fi)

# ----------------------------------------------------------------- partition and format
if [ "$KEEP" = 0 ]; then
    if [ -z "$PART" ]; then
        say "New partition table on $DISK"
        while read -r child; do run "${SUDO[@]}" wipefs -aq "$child"; done < <(lsblk -nrpo NAME,TYPE "$DISK" | awk '$2=="part"{print $1}')
        run "${SUDO[@]}" wipefs -aq "$DISK"
        if command -v parted >/dev/null; then
            run "${SUDO[@]}" parted -s "$DISK" mklabel gpt mkpart "$LABEL" ext4 1MiB 100%
        elif command -v sgdisk >/dev/null; then
            run "${SUDO[@]}" sgdisk -Z -n 1:1MiB:0 -t 1:8300 -c "1:$LABEL" "$DISK"
        else
            die "parted or sgdisk is required to partition the disk (sudo apt install parted)"
        fi
        if [ "$DRY" = 1 ]; then
            PART="<first partition of $DISK>"
        else
            command -v partprobe >/dev/null && "${SUDO[@]}" partprobe "$DISK" 2>/dev/null || true
            command -v udevadm >/dev/null && "${SUDO[@]}" udevadm settle 2>/dev/null || true
            for _ in $(seq 1 20); do
                PART=$(lsblk -nrpo NAME,TYPE "$DISK" | awk '$2=="part"{print $1; exit}')
                [ -n "$PART" ] && [ -b "$PART" ] && break
                sleep 0.5
            done
            [ -n "$PART" ] && [ -b "$PART" ] || die "the new partition of $DISK did not appear: unplug and replug the stick, then rerun with its partition (e.g. $(part1 "$DISK"))"
        fi
    fi
    say "Formatting $PART as ext4 (label $LABEL)"
    run "${SUDO[@]}" wipefs -aq "$PART"
    run "${SUDO[@]}" mkfs.ext4 -F -q -m 0 -L "$LABEL" "$PART"
fi

# -------------------------------------------------------------------------------- mount
say "Mounting"
MP=""
if [ "$DRY" = 1 ]; then
    MP=${MP_OPT:-/media/$OWNER/$LABEL}
    run mount "$PART" "$MP"
    echo "  (udisksctl, if present, may choose another mount point)"
else
    MP=$(findmnt -nro TARGET "$PART" 2>/dev/null | head -1 || true)
    if [ -z "$MP" ] && [ -z "$MP_OPT" ] && command -v udisksctl >/dev/null; then
        out=$(udisksctl mount -b "$PART" 2>/dev/null || true)       # "Mounted /dev/sdb1 at /media/u/LABEL"
        MP=$(printf '%s' "$out" | sed -n 's/^Mounted .* at \(.*\)$/\1/p' | sed 's/\.$//')
    fi
    if [ -z "$MP" ]; then
        MP=${MP_OPT:-/media/$OWNER/$LABEL}
        "${SUDO[@]}" mkdir -p "$MP"
        "${SUDO[@]}" mount "$PART" "$MP"
    fi
    [ -d "$MP" ] || die "mounting $PART failed"
fi
[ "$DRY" = 1 ] || echo "  $PART mounted on $MP"
run "${SUDO[@]}" chown "$OWNER": "$MP"

# ---------------------------------------------------------------- content of the stick
say "Copying Rapid-CNS2 onto the stick (prepare_usb_en.sh)"
if [ "$DRY" = 1 ]; then
    run bash "$PREPARE" "${PREP_ARGS[@]}" "$MP"
    echo "  (to estimate the space needed: $PREPARE --dry-run <any ext4 mount point>)"
else
    bash "$PREPARE" "${PREP_ARGS[@]}" "$MP"
fi

say "Check of the stick"
if [ "$DRY" = 1 ]; then
    run "$MP/rapidcns2/run.sh" --check
else
    "$MP/rapidcns2/run.sh" --check || warn "some items are missing on the stick (see above)"
    sync
fi

# --------------------------------------------------------------------------------- end
if [ "$EJECT" = 1 ]; then
    say "Ejecting"
    if command -v udisksctl >/dev/null && [ "$DRY" = 0 ]; then
        udisksctl unmount -b "$PART" && udisksctl power-off -b "$DISK" || warn "eject failed: unmount it manually"
    else
        run "${SUDO[@]}" umount "$MP"
    fi
    echo "The stick can be removed."
elif [ "$DRY" = 0 ]; then
    echo
    echo "Stick ready: $MP/rapidcns2   (unmount before removing it: umount $MP)"
fi
[ "$DRY" = 1 ] && echo "--dry-run: nothing was changed." || echo "Log: $LOG"
