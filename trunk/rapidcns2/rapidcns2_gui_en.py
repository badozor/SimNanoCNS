#!/usr/bin/env python3
# Developed with the help of AI (Claude, Anthropic), under active review by Pierre Bady.
"""Graphical interface to run Rapid-CNS2_nf on a sample chosen from a menu.

Usage: python3 rapidcns2_gui_en.py
Requirements: Python 3 and tkinter (sudo apt install python3-tk if missing).
"""
import os
import queue
import re
import time
import shlex
import shutil
import subprocess
import threading
import webbrowser
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext

# ----------------------------------------------------------------------------
# Default values: adapt them here once and for all
# ----------------------------------------------------------------------------
DEFAULTS = {
    "pipeline": "/export/scratch3/gitProject/Rapid-CNS2_nf/main.nf",
    "cohort": "/export/scratch3/gitProject/SimNanoCNS/trunk/dev/cohort_test",
    "rundir": "/export/scratch3/gitProject/SimNanoCNS/trunk/dev/test_pipeline",
    "outdir": "results_{id}",          # relative to the launch directory; {id} = sample
    "config": "custom.config",         # relative to the launch directory; ignored if missing
    "ref": "/export/scratch3/data/Annotations/UCSC_HG38/hg38.fa",
    "annovarPath": "/export/scratch3/data/Annotations/humandb/annovar",
    "annovarDB": "/export/scratch3/data/Annotations/humandb/annovar/humandb",
    "annotsvAnnot": "/export/scratch3/gitProject/AnnotSV/share/AnnotSV",
    "profile": "local,docker",
    "minimumMgmtCov": "10",
    "maxCpus": "16",
    "maxMemory": "64",
    "purity": "0.8",                   # --overrule_cellularity for SAVANA; empty = not passed
}

# Environment overrides (used by run_gui.sh on the USB stick):
# RCNS_<UPPERCASE KEY>, e.g. RCNS_REF, RCNS_ANNOVARDB, RCNS_MAXCPUS.
for _k in DEFAULTS:
    DEFAULTS[_k] = os.environ.get("RCNS_" + _k.upper(), DEFAULTS[_k])
HUMANVAR_DEFAULT = os.environ.get("RCNS_HUMANVAR", "true").lower() == "true"
BIND_PATHS = os.environ.get("RCNS_BINDPATHS", "")   # --containerBindPaths if set

# Nextflow >= 26.04 parses with the strict syntax by default, which Rapid-CNS2_nf 3.0 does not
# follow yet (import in nextflow.config, statements outside the workflow in main.nf): the legacy
# parser is used unless NXF_SYNTAX_PARSER is already set (e.g. to v2 once the pipeline is migrated).
NXF_ENV = {"NXF_SYNTAX_PARSER": os.environ.get("NXF_SYNTAX_PARSER", "v1")}


# Rapid-CNS2 v3.0 processes for one sample (taken from a complete run).
# Used to estimate the total number of steps for the progress bar.
PIPELINE_TASKS = [
    "prepareBam", "subsetBam", "mosdepth", "methylationCalls", "methylationClassification",
    "checkMgmtCoverage", "mgmtPromoterMethyartist", "mgmtPred", "copyNumberVariants",
    "cnvAnnotated", "detectClair3Model", "clair3", "recodeVCF", "convert2annovar",
    "tableAnnovar", "filterReport", "igv_reports", "structuralVariants", "annotSvSniffles",
    "severus", "annotSvSeverus", "savanaTo", "savanaCna", "savanaCnvPlot",
    "fusionsSniffles", "fusionsSeverus", "fusionsSavana", "reportRendering",
]
HUMANVAR_TASKS = ["human_variation_snp", "human_variation_sv"]
MNPFLEX_TASKS = ["mnpFlex"]
DONE_STATUS = {"COMPLETED", "CACHED", "FAILED", "ABORTED"}
# "[ab/123456] Submitted process > WORKFLOW:SUB:taskName (1)": the process name may carry
# workflow prefixes separated by ":"; only the last segment is kept (see task_name()).
LOG_TASK_RE = re.compile(r"\]\s+(Submitted|Cached) process > ([A-Za-z0-9_:]+)")


def task_name(raw):
    """'WF:SUB:taskName (1)' -> 'taskName' (same convention as count_trace())."""
    return raw.split(" (")[0].split(":")[-1]


def count_trace(trace_file):
    """Finished processes (any outcome) according to a Nextflow trace file."""
    done = set()
    try:
        with open(trace_file) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_name, i_status = header.index("name"), header.index("status")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(i_name, i_status) and f[i_status] in DONE_STATUS:
                    done.add(task_name(f[i_name]))
    except (OSError, ValueError):
        pass
    return done


def find_bams(cohort_dir):
    """Return {id: path} for the BAMs in <cohort>/*.bam and <cohort>/*/*.bam."""
    d = Path(cohort_dir).expanduser()
    if not d.is_dir():
        return {}
    bams = sorted(list(d.glob("*.bam")) + list(d.glob("*/*.bam")))
    return {p.stem: p for p in bams}


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Rapid-CNS2: pipeline launcher")
        self.geometry("980x780")
        self.proc = None
        self.log_queue = queue.Queue()
        self.samples = {}
        self.vars = {k: tk.StringVar(value=v) for k, v in DEFAULTS.items()}
        self.sample = tk.StringVar()
        self.mnpflex = tk.BooleanVar(value=True)
        self.humanvar = tk.BooleanVar(value=HUMANVAR_DEFAULT)
        self.resume = tk.BooleanVar(value=True)
        self.status = tk.StringVar(value="Ready")
        self.progress_text = tk.StringVar(value="")
        self.start_time = None
        self.expected = set()
        self.log_done = set()
        self.current_task = ""
        self._build()
        self.refresh_samples()
        self.after(100, self._poll_log)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ------------------------------------------------------------------ UI
    def _build(self):
        pad = {"padx": 6, "pady": 3}

        # --- Sample
        f1 = ttk.LabelFrame(self, text="Sample")
        f1.pack(fill="x", padx=8, pady=6)
        f1.columnconfigure(1, weight=1)
        ttk.Label(f1, text="Cohort directory").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(f1, textvariable=self.vars["cohort"]).grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(f1, text="Browse…", command=self._browse_cohort).grid(row=0, column=2, **pad)
        ttk.Button(f1, text="Refresh", command=self.refresh_samples).grid(row=0, column=3, **pad)
        ttk.Label(f1, text="Sample").grid(row=1, column=0, sticky="w", **pad)
        self.combo = ttk.Combobox(f1, textvariable=self.sample, state="readonly")
        self.combo.grid(row=1, column=1, sticky="ew", **pad)
        self.combo.bind("<<ComboboxSelected>>", lambda e: self._update_preview())
        self.count_label = ttk.Label(f1, text="")
        self.count_label.grid(row=1, column=2, columnspan=2, sticky="w", **pad)

        # --- Paths
        f2 = ttk.LabelFrame(self, text="Paths")
        f2.pack(fill="x", padx=8, pady=6)
        f2.columnconfigure(1, weight=1)
        rows = [
            ("Pipeline (main.nf)", "pipeline", "file"),
            ("Launch directory", "rundir", "dir"),
            ("Output directory", "outdir", None),
            ("Config file (-c)", "config", None),
            ("hg38 reference", "ref", "file"),
            ("ANNOVAR", "annovarPath", "dir"),
            ("ANNOVAR humandb", "annovarDB", "dir"),
            ("AnnotSV annotations", "annotsvAnnot", "dir"),
        ]
        for i, (label, key, kind) in enumerate(rows):
            ttk.Label(f2, text=label).grid(row=i, column=0, sticky="w", **pad)
            ttk.Entry(f2, textvariable=self.vars[key]).grid(row=i, column=1, sticky="ew", **pad)
            if kind:
                ttk.Button(f2, text="…", width=3,
                           command=lambda k=key, t=kind: self._browse(k, t)).grid(row=i, column=2, **pad)

        # --- Options
        f3 = ttk.LabelFrame(self, text="Options")
        f3.pack(fill="x", padx=8, pady=6)
        opts = [("Profile", "profile", 16), ("Min. MGMT cov.", "minimumMgmtCov", 5),
                ("CPU max", "maxCpus", 5), ("Max memory (GB)", "maxMemory", 5),
                ("SAVANA purity", "purity", 5)]
        for i, (label, key, width) in enumerate(opts):
            ttk.Label(f3, text=label).grid(row=0, column=2 * i, sticky="w", **pad)
            ttk.Entry(f3, textvariable=self.vars[key], width=width).grid(row=0, column=2 * i + 1, **pad)
        ttk.Checkbutton(f3, text="MNP-Flex", variable=self.mnpflex,
                        command=self._update_preview).grid(row=1, column=0, columnspan=2, sticky="w", **pad)
        ttk.Checkbutton(f3, text="Human variation", variable=self.humanvar,
                        command=self._update_preview).grid(row=1, column=2, columnspan=2, sticky="w", **pad)
        ttk.Checkbutton(f3, text="Resume (-resume)", variable=self.resume,
                        command=self._update_preview).grid(row=1, column=4, columnspan=2, sticky="w", **pad)

        # --- Buttons
        f4 = ttk.Frame(self)
        f4.pack(fill="x", padx=8, pady=4)
        self.btn_run = ttk.Button(f4, text="▶ Run", command=self.run)
        self.btn_run.pack(side="left", padx=4)
        self.btn_stop = ttk.Button(f4, text="■ Stop", command=self.stop, state="disabled")
        self.btn_stop.pack(side="left", padx=4)
        ttk.Button(f4, text="Open report", command=self.open_report).pack(side="left", padx=4)
        ttk.Button(f4, text="Copy command", command=self.copy_command).pack(side="left", padx=4)
        ttk.Label(f4, textvariable=self.status).pack(side="right", padx=8)

        # --- Progress
        f5 = ttk.Frame(self)
        f5.pack(fill="x", padx=8, pady=(2, 4))
        style = ttk.Style(self)
        style.configure("ok.Horizontal.TProgressbar", background="#2e8b57")
        style.configure("err.Horizontal.TProgressbar", background="#c0392b")
        self.pbar = ttk.Progressbar(f5, mode="determinate", maximum=1, value=0)
        self.pbar.pack(fill="x", padx=4)
        ttk.Label(f5, textvariable=self.progress_text).pack(anchor="w", padx=4, pady=(2, 0))

        # --- Command + log
        ttk.Label(self, text="Command").pack(anchor="w", padx=10)
        self.preview = tk.Text(self, height=5, wrap="word", font=("monospace", 9))
        self.preview.pack(fill="x", padx=8, pady=2)
        ttk.Label(self, text="Log").pack(anchor="w", padx=10)
        self.log = scrolledtext.ScrolledText(self, wrap="none", font=("monospace", 9), state="disabled")
        self.log.pack(fill="both", expand=True, padx=8, pady=(2, 8))

        for v in self.vars.values():
            v.trace_add("write", lambda *a: self._update_preview())

    # ------------------------------------------------------------ helpers
    def _browse(self, key, kind):
        start = self.vars[key].get() or "/"
        if kind == "dir":
            path = filedialog.askdirectory(initialdir=start)
        else:
            path = filedialog.askopenfilename(initialdir=str(Path(start).parent))
        if path:
            self.vars[key].set(path)

    def _browse_cohort(self):
        path = filedialog.askdirectory(initialdir=self.vars["cohort"].get() or "/")
        if path:
            self.vars["cohort"].set(path)
            self.refresh_samples()

    def refresh_samples(self):
        self.samples = find_bams(self.vars["cohort"].get())
        ids = list(self.samples)
        self.combo["values"] = ids
        self.count_label.config(text=f"{len(ids)} BAM file(s) found")
        if self.sample.get() not in ids:
            self.sample.set(ids[0] if ids else "")
        self._update_preview()

    def _outdir(self):
        return self.vars["outdir"].get().replace("{id}", self.sample.get())

    def build_command(self):
        v = {k: var.get().strip() for k, var in self.vars.items()}
        sid = self.sample.get()
        cmd = ["nextflow", "run", v["pipeline"], "-ansi-log", "false",
               "--input", str(self.samples.get(sid, "")),
               "--id", sid,
               "--outDir", self._outdir(),
               "--ref", v["ref"],
               "--annovarPath", v["annovarPath"],
               "--annovarDB", v["annovarDB"],
               "--annotsvAnnot", v["annotsvAnnot"],
               "--minimumMgmtCov", v["minimumMgmtCov"],
               "--mnpFlex", "true" if self.mnpflex.get() else "false",
               "--runHumanVariation", "true" if self.humanvar.get() else "false"]
        if BIND_PATHS:
            cmd += ["--containerBindPaths", BIND_PATHS]
        if v["purity"]:
            # value with a space: Nextflow accepts it even though it starts with "--"
            cmd += ["--savanaCnaArgs", f"--overrule_cellularity {v['purity']}"]
        cmd += ["-profile", v["profile"]]
        if v["config"] and (Path(v["rundir"]) / v["config"]).is_file():
            cmd += ["-c", v["config"]]
        cmd += ["--maxCpus", v["maxCpus"], "--maxMemory", v["maxMemory"]]
        if self.resume.get():
            cmd.append("-resume")
        return cmd

    @staticmethod
    def _command_line(cmd):
        """Command as it can be pasted in a terminal, with the Nextflow environment."""
        env = " ".join(f"{k}={shlex.quote(v)}" for k, v in NXF_ENV.items())
        return env + " " + " ".join(shlex.quote(c) for c in cmd)

    def _update_preview(self):
        self.preview.delete("1.0", "end")
        if self.sample.get():
            self.preview.insert("1.0", self._command_line(self.build_command()))

    def copy_command(self):
        self.clipboard_clear()
        self.clipboard_append(self.preview.get("1.0", "end").strip())
        self.status.set("Command copied")

    def _append_log(self, text):
        self.log.config(state="normal")
        self.log.insert("end", text)
        self.log.see("end")
        self.log.config(state="disabled")

    # ------------------------------------------------------------- actions
    def run(self):
        sid = self.sample.get()
        rundir = Path(self.vars["rundir"].get())
        errors = []
        if not sid or sid not in self.samples:
            errors.append("No sample selected.")
        elif not self.samples[sid].is_file():
            errors.append(f"BAM not found: {self.samples[sid]}")
        if not Path(self.vars["pipeline"].get()).is_file():
            errors.append("main.nf not found.")
        if not rundir.is_dir():
            errors.append("The launch directory does not exist.")
        if shutil.which("nextflow") is None:
            errors.append("nextflow is not in the PATH.")
        if errors:
            messagebox.showerror("Cannot start", "\n".join(errors))
            return

        cmd = self.build_command()
        self._append_log("\n$ " + self._command_line(cmd) + "\n\n")
        self.proc = subprocess.Popen(cmd, cwd=rundir, env={**os.environ, **NXF_ENV}, stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT, text=True, bufsize=1)
        threading.Thread(target=self._reader, args=(self.proc,), daemon=True).start()
        self.btn_run.config(state="disabled")
        self.btn_stop.config(state="normal")
        self.status.set(f"Running: {sid}")
        self._start_progress()

    def _reader(self, proc):
        for line in proc.stdout:
            self.log_queue.put(line)
        proc.wait()
        self.log_queue.put(("__done__", proc.returncode))

    def _poll_log(self):
        try:
            while True:
                item = self.log_queue.get_nowait()
                if isinstance(item, tuple):
                    code = item[1]
                    self.proc = None
                    self.btn_run.config(state="normal")
                    self.btn_stop.config(state="disabled")
                    self.status.set("Completed successfully" if code == 0 else f"Failed (code {code})")
                    self._append_log(f"\n--- End (exit code {code}) ---\n")
                    self._finish_progress(code)
                else:
                    self._append_log(item)
                    m = LOG_TASK_RE.search(item)
                    if m:
                        self.current_task = task_name(m.group(2))
                        if m.group(1) == "Cached":
                            self.log_done.add(self.current_task)
        except queue.Empty:
            pass
        self.after(100, self._poll_log)

    def stop(self):
        if self.proc and messagebox.askyesno("Stop", "Stop the running pipeline?"):
            self.proc.terminate()   # SIGTERM: Nextflow shuts its tasks down cleanly
            self.status.set("Stop requested…")

    # ------------------------------------------------------------ progress
    def _outdir_path(self):
        outdir = Path(self._outdir())
        return outdir if outdir.is_absolute() else Path(self.vars["rundir"].get()) / outdir

    def _start_progress(self):
        self.expected = set(PIPELINE_TASKS)
        if self.humanvar.get():
            self.expected |= set(HUMANVAR_TASKS)
        if self.mnpflex.get():
            self.expected |= set(MNPFLEX_TASKS)
        self.log_done = set()
        self.current_task = ""
        self.start_time = time.time()
        self.pbar.config(maximum=len(self.expected), value=0, style="Horizontal.TProgressbar")
        self.progress_text.set("Starting Nextflow…")
        self.after(2000, self._poll_progress)

    def _latest_trace(self):
        """Trace file of this run: the newest one in pipeline_info/, created after launch."""
        info = self._outdir_path() / "pipeline_info"
        traces = [p for p in info.glob("*_trace.txt")
                  if p.stat().st_mtime >= self.start_time - 5] if info.is_dir() else []
        return max(traces, key=lambda p: p.stat().st_mtime) if traces else None

    def _done_tasks(self):
        trace = self._latest_trace()
        done = count_trace(trace) if trace else set()
        return (done | self.log_done) & self.expected

    def _elapsed(self):
        s = int(time.time() - self.start_time)
        return f"{s // 3600:d}h{s % 3600 // 60:02d}m{s % 60:02d}s"

    def _poll_progress(self):
        if self.proc is None:          # run finished: _finish_progress took over
            return
        n, total = len(self._done_tasks()), len(self.expected)
        self.pbar.config(value=n)
        cur = f"   ·   running: {self.current_task}" if self.current_task else ""
        self.progress_text.set(f"{n} / {total} steps ({100 * n // total} %){cur}   ·   {self._elapsed()}")
        self.after(2000, self._poll_progress)

    def _finish_progress(self, code):
        if self.start_time is None:
            return
        done, total = len(self._done_tasks()), len(self.expected)
        if code == 0:
            # some steps do not always run (e.g. savanaCnvPlot when savanaCna fails)
            self.pbar.config(value=total, style="ok.Horizontal.TProgressbar")
            self.progress_text.set(f"Completed: {done} steps run   ·   {self._elapsed()}")
        else:
            self.pbar.config(value=done, style="err.Horizontal.TProgressbar")
            last = f" (last step: {self.current_task})" if self.current_task else ""
            self.progress_text.set(f"Stopped at {done} / {total} steps{last}   ·   {self._elapsed()}")

    def open_report(self):
        sid = self.sample.get()
        rep_dir = self._outdir_path() / "report"
        for kind in ("full", "lite"):
            report = rep_dir / f"{sid}_Rapid-CNS2_report_{kind}.html"
            if report.is_file():
                webbrowser.open(report.as_uri())
                return
        messagebox.showinfo("Report", f"Report not found in:\n{rep_dir}")

    def _on_close(self):
        if self.proc and not messagebox.askyesno(
                "Quit", "A pipeline is running. Stop it and quit?"):
            return
        if self.proc:
            self.proc.terminate()
        self.destroy()


if __name__ == "__main__":
    App().mainloop()
