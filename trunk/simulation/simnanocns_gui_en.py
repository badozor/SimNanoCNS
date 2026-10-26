#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""
SimNanoCNS - graphical interface (Streamlit) to prepare and launch simulations of Nanopore BAM
files (NPHD panel, hg38) used to test Rapid-CNS2_nf. English version of simnanocns_gui.py.

Two modes:
  * Single sample  -> simulate_glioma_hg38.py (+ verify_bam.py, optional)
  * Cohort         -> build_cohort.py (+ optional Capper probe export and classification)

Run:
    pip install "streamlit>=1.37" pandas
    streamlit run simnanocns_gui_en.py

Written for the English scripts (src_en/): simulate_glioma_hg38.py, verify_bam.py,
build_cohort.py, bedmethyl_from_bam.py, methylation_classification_local.R,
02_export_cohort_probe_targets.R.
From another machine:  ssh -L 8501:localhost:8501 P3Ultra   then http://localhost:8501

The interface does no computation itself: it builds the command line, shows it, then runs it
in the background (log in <working directory>/gui_logs/). The simulation keeps running if the
browser tab is closed; it stops if the Streamlit server is stopped.
"""
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import time
from datetime import datetime
from pathlib import Path

import streamlit as st

CONFIG_PATH = Path.home() / ".simnanocns_gui_en.json"
DEV = "/export/scratch3/gitProject/SimNanoCNS/trunk/dev"

DEFAULTS = {
    # --- paths (sidebar) ---
    "workdir": DEV,
    "script_dir": f"{DEV}/src_en",
    "python": "python3",
    "ref_mode": "Real hg38",
    "ref": "/export/scratch3/data/Annotations/UCSC_HG38/hg38.fa",
    "data_dir": f"{DEV}/Rapid-CNS2_nf/data",
    # --- single sample ---
    "s_profile": "gbm",
    "s_sample_id": "",
    "s_out": "",
    "s_purity": 0.80,
    "s_cdkn2a_del": False,
    "s_ontarget": 8.0,
    "s_offtarget": 2.0,
    "s_jitter": 0.0,
    "s_seed": 2026,
    "s_chroms": [],
    "s_skip_ref_check": False,
    "s_mnp_targets": "",
    "s_mnp_depth": 6.0,
    "s_mnp_window": 150,
    "s_verify": True,
    "s_extra": "",
    # --- cohort ---
    "c_n": 99,
    "c_w_gbm": 1.0,
    "c_w_astro": 1.0,
    "c_w_oligo": 1.0,
    "c_out_dir": "cohort_bams",
    "c_jobs": 16,
    "c_jitter": 0.15,
    "c_set_purity": False,
    "c_purity_range": [0.60, 0.95],
    "c_grade4": 0.30,
    "c_mnp_dir": "cohort_probes",
    "c_mnp_depth": 6.0,
    "c_chroms": [],
    "c_classify_out": "cohort_classification",
    "c_export": False,
    "c_export_script": "02_export_cohort_probe_targets.R",
    "c_classify": False,
    "c_classify_jobs": 2,
    "c_classify_threads": 8,
    "c_num_trees": 2000,
    "c_set_seed": False,
    "c_seed": 2026,
    "c_extra": "",
}

PROFILES = {
    "gbm": {
        "label": "IDH-wildtype glioblastoma",
        "methylation_class": "GBM_RTK_II",
        "alterations": "+7 / -10, EGFR amplification, CDKN2A/B homozygous deletion, "
                       "TERT C228T (fixed VAF 0.80), IDH1 wildtype, methylated MGMT promoter",
    },
    "astro": {
        "label": "IDH-mutant astrocytoma",
        "methylation_class": "A_IDH (A_IDH_HG if grade 4)",
        "alterations": "IDH1 R132H, TP53 R175H, ATRX frameshift (1 bp), diploid genome, "
                       "TERT wildtype, methylated MGMT, G-CIMP background; grade 4 with CDKN2A/B deletion",
    },
    "oligo": {
        "label": "IDH-mutant, 1p/19q-codeleted oligodendroglioma",
        "methylation_class": "O_IDH",
        "alterations": "IDH1 R132H, TERT C228T, 1p/19q codeletion (ratio 0.5), CIC frameshift "
                       "(1 bp), TP53 wildtype, methylated MGMT, G-CIMP background",
    },
}
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


# ----------------------------------------------------------------------------------------------
# Pure functions (testable without Streamlit)
# ----------------------------------------------------------------------------------------------
def resolve(path, base):
    """Absolute path; relative paths are resolved against `base`."""
    p = Path(os.path.expanduser(str(path)))
    return p if p.is_absolute() else Path(base) / p


def sample_id_of(p):
    return p["s_sample_id"].strip() or f"SIM_{p['s_profile'].upper()}_HG38"


def single_out_of(p):
    return p["s_out"].strip() or f"sim_{p['s_profile']}"


def classify_out_of(p):
    return p["c_classify_out"].strip() or "cohort_classification"


def ref_args(p):
    return ["--standin"] if p["ref_mode"] == "standin" else ["--ref", p["ref"]]


def build_single_cmds(p, verify_flags=None):
    """Commands for one sample: simulation, then verification if requested."""
    sd = Path(p["script_dir"])
    cmd = [p["python"], str(sd / "simulate_glioma_hg38.py"), "--profile", p["s_profile"],
           *ref_args(p), "--data-dir", p["data_dir"], "--out", single_out_of(p),
           "--sample-id", sample_id_of(p),
           "--purity", f"{p['s_purity']:.2f}",
           "--ontarget-depth", f"{p['s_ontarget']:g}",
           "--offtarget-depth", f"{p['s_offtarget']:g}",
           "--seed", str(int(p["s_seed"]))]
    if p["s_profile"] == "astro" and p["s_cdkn2a_del"]:
        cmd.append("--cdkn2a-del")
    if p["s_jitter"] > 0:
        cmd += ["--cnv-jitter-sd", f"{p['s_jitter']:g}"]
    if p["s_chroms"]:
        cmd += ["--chroms", ",".join(c for c in CHROMS if c in p["s_chroms"])]
    if p["s_skip_ref_check"] and p["ref_mode"] != "standin":
        cmd.append("--skip-ref-check")
    if p["s_mnp_targets"].strip():
        cmd += ["--mnp-targets", p["s_mnp_targets"].strip()]
        if p["s_mnp_depth"] != 6.0:
            cmd += ["--mnp-depth", f"{p['s_mnp_depth']:g}"]
        if int(p["s_mnp_window"]) != 150:
            cmd += ["--mnp-window", str(int(p["s_mnp_window"]))]
    cmd += shlex.split(p["s_extra"])
    cmds = [cmd]
    if p["s_verify"]:
        out = single_out_of(p)
        v = [p["python"], str(sd / "verify_bam.py"), "--bam", f"{out}/{sample_id_of(p)}.bam",
             "--data-dir", p["data_dir"], *ref_args(p)]
        if p["ref_mode"] == "standin" and (verify_flags is None or "--seed" in verify_flags):
            # the stand-in sequence is derived from the seed: verify_bam.py must rebuild the
            # same one, otherwise the mismatch check reports ~75% instead of ~2%
            v += ["--seed", str(int(p["s_seed"]))]
        if verify_flags is None or "--truth" in verify_flags:
            v += ["--truth", f"{out}/truth_events.tsv"]
        cmds.append(v)
    return cmds


GROUPS = ["gbm", "astro", "oligo"]


def equal_counts(n):
    """Default split of build_cohort.py: n // 3 each, remainder given to gbm then astro."""
    n = int(n)
    base, rest = divmod(n, 3)
    return {g: base + (i < rest) for i, g in enumerate(GROUPS)}


def group_counts(n, weights):
    """Group sizes for n individuals from relative weights (largest remainder)."""
    n = int(n)
    w = [max(0.0, float(weights[g])) for g in GROUPS]
    tot = sum(w)
    if tot <= 0:
        raise ValueError("at least one group must have a positive weight")
    if len(set(w)) == 1:
        return equal_counts(n)
    exact = [n * x / tot for x in w]
    counts = [int(e) for e in exact]
    order = sorted(range(3), key=lambda i: (-(exact[i] - counts[i]), i))
    for i in order[:n - sum(counts)]:
        counts[i] += 1
    return dict(zip(GROUPS, counts))


def counts_of(p):
    return group_counts(p["c_n"], {g: p[f"c_w_{g}"] for g in GROUPS})


def per_class(p):
    """Number of methylation profiles to export per class: the size of the largest group."""
    return max(counts_of(p).values())


def build_cohort_cmds(p, cohort_flags=None):
    """Commands for a cohort: optional probe export, then build_cohort.py."""
    sd = Path(p["script_dir"])
    cmds = []
    mnp = p["c_mnp_dir"].strip()
    if p["c_export"] and mnp:
        cmds.append(["Rscript", str(resolve(p["c_export_script"], sd)), p["data_dir"], mnp,
                     str(per_class(p))])
    cmd = [p["python"], str(sd / "build_cohort.py"), "--n", str(int(p["c_n"])), *ref_args(p),
           "--data-dir", p["data_dir"], "--out-dir", p["c_out_dir"].strip() or "cohort_bams",
           "--jobs", str(int(p["c_jobs"]))]
    counts = counts_of(p)
    if counts != equal_counts(p["c_n"]):
        cmd += ["--profile-counts", ",".join(f"{g}={counts[g]}" for g in GROUPS)]
    if mnp:
        cmd += ["--mnp-targets-dir", mnp]
        if p["c_mnp_depth"] != 6.0:
            cmd += ["--mnp-depth", f"{p['c_mnp_depth']:g}"]
    if p["c_chroms"]:
        cmd += ["--chroms", ",".join(c for c in CHROMS if c in p["c_chroms"])]
    if p["c_jitter"] > 0:
        cmd += ["--cnv-jitter-sd", f"{p['c_jitter']:g}"]
    if p["c_set_purity"]:
        lo, hi = p["c_purity_range"]
        cmd += ["--purity-range", f"{lo:.2f}", f"{hi:.2f}"]
    if p["c_grade4"] > 0:
        cmd += ["--astro-grade4-fraction", f"{p['c_grade4']:g}"]
    if p["c_set_seed"]:
        cmd += ["--seed-base", str(int(p["c_seed"]))]
    if p["c_classify"]:
        cmd += ["--classify", "--classify-jobs", str(int(p["c_classify_jobs"])),
                "--classify-threads", str(int(p["c_classify_threads"])),
                "--num-trees", str(int(p["c_num_trees"]))]
        if classify_out_of(p) != "cohort_classification":
            cmd += ["--classify-out-dir", classify_out_of(p)]
    cmd += shlex.split(p["c_extra"])
    cmds.append(cmd)
    return cmds


def cohort_errors(p):
    errs = []
    try:
        counts = counts_of(p)
    except ValueError as e:
        return [f"Group split: {e}."]
    if counts["astro"] == 0 and p["c_grade4"] > 0:
        errs.append("The grade 4 fraction is non-zero but the cohort contains no "
                    "astrocytoma: set the fraction to 0.")
    if not p["c_mnp_dir"].strip():
        if p["c_grade4"] > 0:
            errs.append("The grade 4 astrocytoma fraction requires an MNP probe directory "
                        "(without real methylation, A_IDH and A_IDH_HG cannot be told apart). "
                        "Set the directory or set the fraction to 0.")
        if p["c_classify"]:
            errs.append("Classification requires an MNP probe directory.")
    return errs


def used_flags(cmd):
    return [a for a in cmd if re.fullmatch(r"--[a-z0-9][a-z0-9-]*", a)]


def parse_help_flags(text):
    return set(re.findall(r"(?<![\w-])--[a-z0-9][a-z0-9-]*", text))


def unknown_flags(cmd, known):
    """Command options missing from the script help (None if the help is unavailable)."""
    if known is None:
        return None
    return [f for f in used_flags(cmd) if f not in known]


def shell_line(cmds):
    return " && \\\n".join(shlex.join(c) for c in cmds)


def preflight(p, mode):
    """List of (level, message); level in {'ok', 'warn', 'err'}."""
    out = []
    wd = Path(p["workdir"])
    out.append(("ok" if wd.is_dir() else "err", f"Working directory: {wd}"))
    sd = Path(p["script_dir"])
    needed = ["simulate_glioma_hg38.py"]
    if mode == "single" and p["s_verify"]:
        needed.append("verify_bam.py")
    if mode == "cohort":
        needed.append("build_cohort.py")
        if p["c_classify"]:
            needed += ["bedmethyl_from_bam.py", "methylation_classification_local.R"]
    for s in needed:
        out.append(("ok" if (sd / s).is_file() else "err", f"Script {s} in {sd}"))
    dd = Path(p["data_dir"])
    out.append(("ok" if (dd / "NPHD_panel_hg38.bed").is_file() else "err",
                f"Panel NPHD_panel_hg38.bed in {dd}"))
    if p["ref_mode"] == "standin":
        out.append(("warn", "Stand-in sequence (--standin): BAM unusable with Rapid-CNS2_nf "
                            "and a real hg38, plumbing test only"))
    else:
        ref = Path(p["ref"])
        out.append(("ok" if ref.is_file() else "err", f"Reference: {ref}"))
        if ref.is_file() and not Path(str(ref) + ".fai").is_file():
            out.append(("warn", "Reference .fai index missing: the simulator creates it on first "
                                "use (write access to the reference folder needed)"))
    py = p["python"]
    out.append(("ok" if shutil.which(py) or Path(py).is_file() else "err", f"Interpreter {py}"))
    if mode == "cohort" and (p["c_export"] or p["c_classify"]):
        out.append(("ok" if shutil.which("Rscript") else "err", "Rscript in PATH"))
    if mode == "cohort" and p["c_export"]:
        ex = resolve(p["c_export_script"], sd)
        out.append(("ok" if ex.is_file() else "err", f"Probe export script: {ex}"))
    if mode == "cohort" and (p["c_export"] or p["c_classify"]):
        for f in ("top_probes_hm450.Rdata", "HM450.hg38.manifest.gencode.v22.Rdata",
                  "capper_top_100k_betas_binarised.Rdata"):
            out.append(("ok" if (dd / f).is_file() else "err", f"Reference data {f} in {dd}"))
    if mode == "cohort" and p["c_mnp_dir"].strip() and not p["c_export"]:
        m = resolve(p["c_mnp_dir"].strip(), wd)
        out.append(("ok" if m.is_dir() else "err", f"MNP probe directory: {m}"))
        if m.is_dir():
            # build_cohort.py looks for probe_targets_<group>_<NN>.tsv, NN = 01..group size
            try:
                need = dict(counts_of(p))
            except ValueError:
                need = {}
            if p["c_grade4"] > 0 and need.get("astro"):
                need["astro_hg"] = need["astro"]
            for g, k in need.items():
                if k and not (m / f"probe_targets_{g}_{k:02d}.tsv").is_file():
                    have = len(list(m.glob(f"probe_targets_{g}_[0-9]*.tsv")))
                    out.append(("warn", f"Only {have} probe file(s) for {g} in {m}, {k} needed: the "
                                        "others will be simulated without a real methylation signal "
                                        "(tick the probe export, or rerun 02_export_cohort_probe_targets.R)"))
    if mode == "single" and p["s_mnp_targets"].strip():
        m = resolve(p["s_mnp_targets"].strip(), wd)
        out.append(("ok" if m.is_file() else "err", f"MNP probe file: {m}"))
    target = resolve(single_out_of(p) if mode == "single" else p["c_out_dir"], wd)
    if target.exists() and any(target.iterdir()):
        out.append(("warn", f"Output directory {target} already exists and is not empty"))
    if wd.is_dir():
        free = shutil.disk_usage(wd).free / 1e9
        out.append(("ok" if free > 50 else "warn", f"Free disk space: {free:.0f} GB"))
    return out


def tail(path, n=200):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - 200_000))
            return "\n".join(fh.read().decode("utf-8", "replace").splitlines()[-n:])
    except OSError:
        return ""


# ----------------------------------------------------------------------------------------------
# Background job
# ----------------------------------------------------------------------------------------------
@st.cache_resource
def registry():
    """One job at a time, shared across page reloads."""
    return {"job": None}


@st.cache_data(ttl=300, show_spinner=False)
def help_text(python, script, mtime):
    try:
        r = subprocess.run([python, script, "--help"], capture_output=True, text=True, timeout=60)
        return r.stdout + r.stderr if r.returncode == 0 else None
    except (OSError, subprocess.TimeoutExpired):
        return None


def script_help(p, name):
    path = Path(p["script_dir"]) / name
    if not path.is_file():
        return None
    return help_text(p["python"], str(path), path.stat().st_mtime)


def launch(p, cmds, label, mode, out_dir, expected, classify_out=None):
    wd = Path(p["workdir"])
    logs = wd / "gui_logs"
    logs.mkdir(exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log = logs / f"{stamp}_{label}.log"
    params = logs / f"{stamp}_{label}.params.json"
    params.write_text(json.dumps({"commands": [shlex.join(c) for c in cmds],
                                  "parameters": p}, indent=2, ensure_ascii=False))
    line = "set -eo pipefail; " + " && ".join(shlex.join(c) for c in cmds)
    fh = open(log, "wb")
    fh.write(("# " + shell_line(cmds) + "\n\n").encode())
    fh.flush()
    env = dict(os.environ, PYTHONUNBUFFERED="1")
    proc = subprocess.Popen(["bash", "-c", line], cwd=wd, stdout=fh, stderr=subprocess.STDOUT,
                            start_new_session=True, env=env)
    fh.close()   # the child process keeps its own copy of the descriptor
    registry()["job"] = dict(proc=proc, log=str(log), params=str(params), label=label, mode=mode,
                             start=time.time(), end=None, out_dir=str(resolve(out_dir, wd)),
                             expected=expected, workdir=str(wd), stopped=False,
                             classify_out=str(resolve(classify_out, wd)) if classify_out else None)


def stop_job():
    job = registry()["job"]
    if job and job["proc"].poll() is None:
        job["stopped"] = True
        try:
            os.killpg(os.getpgid(job["proc"].pid), signal.SIGTERM)
        except ProcessLookupError:
            pass


def job_running():
    job = registry()["job"]
    return bool(job) and job["proc"].poll() is None


# ----------------------------------------------------------------------------------------------
# Interface
# ----------------------------------------------------------------------------------------------
def load_saved():
    try:
        return {k: v for k, v in json.loads(CONFIG_PATH.read_text()).items() if k in DEFAULTS}
    except (OSError, ValueError):
        return {}


LEGACY = {"hg38 réel": "Real hg38", "Séquence de secours (--standin)": "Stand-in sequence (--standin)",
          "espace": "space"}


def norm(k, v):
    """Accept values saved or exported by the French version of the interface."""
    if isinstance(v, str):
        v = LEGACY.get(v, v)
    return tuple(v) if k == "c_purity_range" else v


def init_state():
    if st.session_state.get("_init"):
        return
    for k, v in {**DEFAULTS, **load_saved()}.items():
        st.session_state[k] = norm(k, v)
    st.session_state["_init"] = True


def params():
    p = {k: st.session_state[k] for k in DEFAULTS}
    p["c_purity_range"] = list(p["c_purity_range"])
    p["ref_mode"] = "standin" if p["ref_mode"].startswith("Stand-in") else "hg38"
    return p


def apply_config(cfg):
    for k, v in cfg.items():
        if k in DEFAULTS:
            st.session_state[k] = norm(k, v)


def cb_import():
    up = st.session_state.get("_upload")
    if up is None:
        return
    try:
        cfg = json.loads(up.getvalue().decode())
        apply_config(cfg.get("gui_parameters", cfg))
        st.session_state["_msg"] = "Configuration loaded."
    except (ValueError, AttributeError):
        st.session_state["_msg"] = "Unreadable configuration file."


def cb_reset():
    apply_config(DEFAULTS)


def exportable():
    return {k: (list(v) if isinstance(v, tuple) else v)
            for k, v in ((k, st.session_state[k]) for k in DEFAULTS)}


def sidebar():
    sb = st.sidebar
    sb.header("Paths")
    sb.text_input("Working directory", key="workdir",
                  help="Relative paths (outputs, probes) are resolved from here.")
    sb.text_input("Script directory", key="script_dir",
                  help="Contains simulate_glioma_hg38.py, build_cohort.py, verify_bam.py.")
    sb.text_input("Python interpreter", key="python",
                  help="The one from the environment where pysam and numpy are installed.")
    sb.text_input("Rapid-CNS2_nf data/ directory", key="data_dir")
    sb.radio("Reference sequence", ["Real hg38", "Stand-in sequence (--standin)"],
             key="ref_mode")
    sb.text_input("hg38 FASTA", key="ref", disabled=st.session_state["ref_mode"] != "Real hg38")
    sb.divider()
    sb.header("Configuration")
    if sb.button("Save as defaults"):
        CONFIG_PATH.write_text(json.dumps(exportable(), indent=2, ensure_ascii=False))
        sb.success(f"Saved to {CONFIG_PATH}")
    sb.download_button("Export (JSON)", json.dumps(exportable(), indent=2, ensure_ascii=False),
                       file_name="simnanocns_config.json", mime="application/json")
    sb.file_uploader("Import (JSON)", type="json", key="_upload", on_change=cb_import)
    sb.button("Reset", on_click=cb_reset)
    if st.session_state.get("_msg"):
        sb.info(st.session_state.pop("_msg"))


def show_checks(checks):
    icon = {"ok": ":material/check_circle:", "warn": ":material/warning:", "err": ":material/error:"}
    for lvl, msg in checks:
        if lvl == "ok":
            st.caption(f"{icon[lvl]} {msg}")
        elif lvl == "warn":
            st.warning(msg, icon=icon[lvl])
        else:
            st.error(f"Missing or invalid: {msg}", icon=icon[lvl])
    return not any(lvl == "err" for lvl, _ in checks)


def show_command(p, cmds, scripts):
    """Command preview + check of the options against each script's help."""
    st.subheader("Command")
    st.code(f"cd {shlex.quote(p['workdir'])}\n" + shell_line(cmds), language="bash")
    ok = True
    for cmd in cmds:
        name = Path(cmd[1]).name
        if name not in scripts:
            continue
        bad = unknown_flags(cmd, scripts[name])
        if bad is None:
            st.caption(f"Options of {name} not checked (script help unavailable).")
        elif bad:
            ok = False
            st.error(f"{name} does not know: {', '.join(bad)}. The script directory probably "
                     "holds an older version: point it to the English scripts (src_en), "
                     "or check the script help below.")
    return ok


def flags_for(p, name):
    txt = script_help(p, name)
    return (parse_help_flags(txt) if txt else None), txt


def tab_single(p):
    left, right = st.columns([3, 2], gap="large")
    with left:
        st.subheader("Tumour profile")
        st.radio("Profile", list(PROFILES), key="s_profile", horizontal=True,
                 format_func=lambda k: f"{k}: {PROFILES[k]['label']}", label_visibility="collapsed")
        prof = st.session_state["s_profile"]
        st.info(f"**Simulated alterations**: {PROFILES[prof]['alterations']}.  \n"
                f"**Expected methylation class** (with MNP probes): {PROFILES[prof]['methylation_class']}")
        c1, c2 = st.columns(2)
        c1.text_input("Sample ID", key="s_sample_id",
                      placeholder=f"SIM_{prof.upper()}_HG38")
        c2.text_input("Output directory", key="s_out", placeholder=f"sim_{prof}")
        st.subheader("Tumour")
        c1, c2 = st.columns(2)
        c1.slider("Tumour purity", 0.10, 1.00, step=0.05, key="s_purity")
        vaf = st.session_state["s_purity"] / 2
        note = f"VAF of heterozygous variants: {vaf:.2f}"
        if prof == "gbm":
            note += " (TERT C228T stays at 0.80)"
        c1.caption(note)
        c2.number_input("Copy-ratio variability (SD)", 0.0, 0.5, step=0.05,
                        key="s_jitter", help="--cnv-jitter-sd; 0 = fixed profile ratios.")
        c2.checkbox("CDKN2A/B homozygous deletion (grade 4)", key="s_cdkn2a_del",
                    disabled=prof != "astro", help="astro profile only.")
        st.subheader("Sequencing")
        c1, c2, c3 = st.columns(3)
        c1.number_input("On-target depth (x)", 1.0, 100.0, step=1.0, key="s_ontarget",
                        help="Per panel region before overlap; the effective panel coverage "
                             "is about 2.5 times higher.")
        c2.number_input("Off-target depth (x)", 0.0, 10.0, step=0.5, key="s_offtarget",
                        help="CNVpytor fails below about 1 to 2x (script default: 2x).")
        c3.number_input("Random seed", 0, 2**31 - 1, step=1, key="s_seed")
        with st.expander("Advanced options"):
            st.multiselect("Simulated chromosomes (empty = chr1-22 and chrX)", CHROMS, key="s_chroms")
            st.text_input("MNP probe file (--mnp-targets)", key="s_mnp_targets",
                          placeholder="probe_targets/probe_targets_gbm.tsv",
                          help="Transplanted real methylation profile; required for an "
                               "interpretable classification.")
            c1, c2 = st.columns(2)
            c1.number_input("Depth per HM450 probe (x)", 1.0, 50.0, step=1.0, key="s_mnp_depth",
                            help="--mnp-depth; only used with an MNP probe file.")
            c2.number_input("Half-window around each probe (bp)", 50, 1000, step=50,
                            key="s_mnp_window", help="--mnp-window; only used with an MNP probe file.")
            st.checkbox("Skip the reference base check (--skip-ref-check)",
                        key="s_skip_ref_check")
            st.checkbox("Verify the BAM after the simulation (verify_bam.py)", key="s_verify")
            st.text_input("Extra options", key="s_extra")
    p = params()
    with right:
        sim_flags, sim_help = flags_for(p, "simulate_glioma_hg38.py")
        ver_flags, _ = flags_for(p, "verify_bam.py")
        try:
            cmds = build_single_cmds(p, ver_flags)
        except ValueError as e:
            st.error(f"Unreadable extra options: {e}")
            return
        ok_cmd = show_command(p, cmds, {"simulate_glioma_hg38.py": sim_flags,
                                        "verify_bam.py": ver_flags})
        st.subheader("Checks")
        ok = show_checks(preflight(p, "single")) and ok_cmd
        if st.button("Run the simulation", type="primary", disabled=not ok or job_running(),
                     key="go_single"):
            launch(p, cmds, sample_id_of(p), "single", single_out_of(p), 1)
            st.rerun()
        if job_running():
            st.caption("A job is already running (see the Run tab).")
        if sim_help:
            with st.expander("Help of simulate_glioma_hg38.py"):
                st.code(sim_help, language=None)


def tab_cohort(p):
    left, right = st.columns([3, 2], gap="large")
    with left:
        st.subheader("Cohort")
        c1, c2, c3 = st.columns(3)
        c1.number_input("Number of individuals", 1, 2000, step=3, key="c_n")
        c2.text_input("Output directory", key="c_out_dir")
        c3.number_input("Parallel simulations", 1, 256, step=1, key="c_jobs")
        n = int(st.session_state["c_n"])
        st.subheader("Split between tumour groups")
        cols = st.columns(3)
        for col, g in zip(cols, GROUPS):
            col.number_input(f"Weight {g}", 0.0, 1000.0, step=1.0, key=f"c_w_{g}",
                             help=f"{PROFILES[g]['label']}. Relative weight: 1/1/1 = equal shares, "
                                  "2/1/1 = half gbm, 50/30/20 = percentages; 0 excludes the group.")
        try:
            counts = group_counts(n, {g: st.session_state[f"c_w_{g}"] for g in GROUPS})
            for col, g in zip(cols, GROUPS):
                col.metric(f"Individuals {g}", counts[g], f"{100 * counts[g] / n:.0f}%",
                           delta_color="off")
            n_g4 = round(counts["astro"] * st.session_state["c_grade4"])
            st.caption("The script shuffles the order of individuals."
                       + (f" About {n_g4} of {counts['astro']} astrocytoma(s) are grade 4."
                          if n_g4 else ""))
        except ValueError:
            pass  # reported in the checks
        st.subheader("Variability between individuals")
        c1, c2 = st.columns(2)
        c1.number_input("Copy-ratio variability (SD)", 0.0, 0.5, step=0.05,
                        key="c_jitter", help="--cnv-jitter-sd; 0.15 gives about +/-15%.")
        c2.slider("Fraction of grade 4 astrocytomas", 0.0, 1.0, step=0.05, key="c_grade4",
                  help="CDKN2A/B deletion and A_IDH_HG methylation; requires MNP probes.")
        st.checkbox("Set the tumour purity range (otherwise script default: 0.60 to 0.95)",
                    key="c_set_purity")
        if st.session_state["c_set_purity"]:
            st.slider("Purity range", 0.10, 1.00, step=0.05, key="c_purity_range")
        st.subheader("Methylation and classification")
        st.text_input("MNP probe directory (--mnp-targets-dir)", key="c_mnp_dir",
                      help="Real profiles from the Capper cohort, one per individual. Empty = random "
                           "methylation background.")
        try:
            n_exp = per_class(params())
        except ValueError:
            n_exp = n
        st.checkbox(f"Export the probes first ({n_exp} per class) with the R script",
                    key="c_export")
        if st.session_state["c_export"]:
            st.text_input("Export script", key="c_export_script",
                          help="Relative to the script directory, or an absolute path.")
        st.checkbox("Classify each BAM after generation (--classify)", key="c_classify")
        if st.session_state["c_classify"]:
            c1, c2, c3 = st.columns(3)
            c1.number_input("Parallel classifications", 1, 32, step=1, key="c_classify_jobs",
                            help="Several GB of RAM per classification.")
            c2.number_input("Threads per classification", 1, 256, step=1,
                            key="c_classify_threads")
            c3.number_input("Number of trees", 100, 50000, step=500, key="c_num_trees",
                            help="2000 for screening, 20000 for a final validation "
                                 "(about 10 times longer).")
        with st.expander("Advanced options"):
            c1, c2 = st.columns(2)
            c1.checkbox("Set the base seed", key="c_set_seed")
            c2.number_input("Seed", 0, 2**31 - 1, step=1, key="c_seed",
                            disabled=not st.session_state["c_set_seed"])
            st.multiselect("Simulated chromosomes (empty = chr1-22 and chrX; restrict for quick tests)",
                           CHROMS, key="c_chroms")
            c1, c2 = st.columns(2)
            c1.number_input("Depth per HM450 probe (x)", 1.0, 50.0, step=1.0, key="c_mnp_depth",
                            help="--mnp-depth; only used with an MNP probe directory.")
            c2.text_input("Classification output directory", key="c_classify_out",
                          help="--classify-out-dir; only used with classification.")
            st.caption("Off-target depth: build_cohort.py has no option for it; the default of "
                       "simulate_glioma_hg38.py applies (2x).")
            st.text_input("Extra options", key="c_extra")
    p = params()
    with right:
        coh_flags, coh_help = flags_for(p, "build_cohort.py")
        try:
            cmds = build_cohort_cmds(p, coh_flags)
        except ValueError as e:
            st.error(f"Invalid cohort parameters: {e}")
            return
        ok_cmd = show_command(p, cmds, {"build_cohort.py": coh_flags})
        st.subheader("Checks")
        errs = cohort_errors(p)
        for e in errs:
            st.error(e)
        ok = show_checks(preflight(p, "cohort")) and ok_cmd and not errs
        if st.button("Run the cohort", type="primary", disabled=not ok or job_running(),
                     key="go_cohort"):
            launch(p, cmds, f"cohort_n{int(p['c_n'])}", "cohort", p["c_out_dir"], int(p["c_n"]),
                   classify_out_of(p) if p["c_classify"] else None)
            st.rerun()
        if job_running():
            st.caption("A job is already running (see the Run tab).")
        if coh_help:
            with st.expander("Help of build_cohort.py"):
                st.code(coh_help, language=None)


def show_tsv(path, title):
    import pandas as pd
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception as e:  # partial file while being written, unexpected format
        st.caption(f"{path}: unreadable ({e})")
        return None
    st.markdown(f"**{title}**  \n`{path}`")
    st.dataframe(df, hide_index=True)
    return df


def show_results(job):
    out = Path(job["out_dir"])
    if not out.is_dir():
        return
    if job["mode"] == "single":
        truth = out / "truth_events.tsv"
        if truth.is_file():
            show_tsv(truth, "Simulated events (ground truth)")
    else:
        man = out / "cohort_manifest.tsv"
        if man.is_file():
            dfm = show_tsv(man, "Cohort manifest")
            if dfm is not None and "status" in dfm.columns:
                failed = dfm[dfm["status"].astype(str) != "OK"]
                if len(failed):
                    st.error(f"{len(failed)} individual(s) failed: " + ", ".join(failed["individual"].astype(str))
                             + f". See the matching .log files in {out}.")
        if job.get("classify_out"):
            summ = Path(job["classify_out"]) / "cohort_classification_summary.tsv"
            if summ.is_file() and summ.stat().st_mtime >= job["start"]:
                df = show_tsv(summ, "Classification")
                if df is not None and {"correct", "status"} <= set(df.columns):
                    done = df[df["status"].astype(str) == "OK"]
                    good = done["correct"].astype(str).str.lower().eq("true").sum()
                    c1, c2 = st.columns(2)
                    c1.metric("Correct classifications", f"{good} / {len(done)}")
                    c2.metric("Failed classifications", f"{len(df) - len(done)}")
            else:
                st.warning(f"No classification summary found in {job['classify_out']}.")
    bams = sorted(out.rglob("*.bam"))
    if bams:
        with st.expander(f"{len(bams)} BAM file(s)"):
            st.code("\n".join(str(b) for b in bams[:500]), language=None)


@st.fragment(run_every=2)
def tab_run():
    job = registry()["job"]
    if not job:
        st.info("No job launched since the interface started.")
        return
    rc = job["proc"].poll()
    if rc is not None and job["end"] is None:
        job["end"] = time.time()
    elapsed = (job["end"] or time.time()) - job["start"]
    c1, c2, c3 = st.columns([2, 1, 1])
    c1.markdown(f"**{job['label']}**  \nLog: `{job['log']}`")
    c2.metric("Duration", f"{int(elapsed // 3600):d} h {int(elapsed % 3600 // 60):02d} min {int(elapsed % 60):02d} s")
    out = Path(job["out_dir"])
    done = sum(1 for b in out.rglob("*.bam") if b.stat().st_mtime >= job["start"]) if out.is_dir() else 0
    c3.metric("BAMs produced", f"{done} / {job['expected']}")
    if rc is None:
        st.progress(min(1.0, done / max(1, job["expected"])), text="Simulation running")
        st.button("Stop", on_click=stop_job, key="stop")
    elif job["stopped"]:
        st.warning("Job stopped by the user.")
    elif rc == 0:
        st.success("Job finished without error.")
    else:
        st.error(f"Job finished with exit code {rc} (see the end of the log).")
    st.code(tail(job["log"]) or "(empty log)", language=None)
    if rc is not None:
        show_results(job)


def main():
    st.set_page_config(page_title="SimNanoCNS", layout="wide")
    init_state()
    st.title("SimNanoCNS: simulation setup")
    st.caption("Simulated Nanopore BAMs (adaptive sampling, NPHD panel, hg38) to test Rapid-CNS2_nf")
    sidebar()
    p = params()
    t1, t2, t3 = st.tabs(["Single sample", "Cohort", "Run"])
    with t1:
        tab_single(p)
    with t2:
        tab_cohort(p)
    with t3:
        tab_run()


if __name__ == "__main__":
    main()
