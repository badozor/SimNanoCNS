#!/usr/bin/env python3
"""Generate the SimNanoCNS hex sticker logo (SVG + optional PNG).

Usage:
    python3 make_logo.py [seed]          # default seed 7
    rsvg-convert -w 1040 simnanocns_hex.svg -o simnanocns_hex.png

Inputs : brain_contour.npy (brain silhouette, x/y points) in the same folder.
Outputs: simnanocns_hex.svg in the same folder.
Requires: numpy, scipy, matplotlib.
"""
import os
import math, sys
import numpy as np
from scipy.interpolate import splprep, splev
from scipy.spatial import Delaunay
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 7
rng = np.random.default_rng(seed)
W, H = 520, 600
cx, cy, R = 260, 300, 288
hw = R * math.sqrt(3) / 2
hexpts = [(cx, cy - R), (cx + hw, cy - R / 2), (cx + hw, cy + R / 2),
          (cx, cy + R), (cx - hw, cy + R / 2), (cx - hw, cy - R / 2)]
hexd = "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in hexpts) + " Z"

# ---------- brain outline (lateral view, frontal lobe on the left) ----------
raw = np.load(os.path.join(HERE, "brain_contour.npy"))
x0i, y0i = raw[:, 0].min(), raw[:, 1].min()
sc = 350 / (raw[:, 0].max() - x0i)
raw = np.c_[(raw[:, 0] - x0i) * sc + 85, (raw[:, 1] - y0i) * sc + 78]
raw = raw[::4]
tck, _ = splprep([np.r_[raw[:, 0], raw[0, 0]], np.r_[raw[:, 1], raw[0, 1]]], s=len(raw) * 1.5, per=True)
ox, oy = splev(np.linspace(0, 1, 1500), tck)
outline = np.c_[ox, oy]
print("bbox", outline.min(0), outline.max(0))
poly = Path(outline)

def along(curve, n):
    seg = np.sqrt((np.diff(curve, axis=0) ** 2).sum(1))
    cum = np.r_[0, np.cumsum(seg)]
    return np.array([curve[min(np.searchsorted(cum, t), len(curve) - 1)]
                     for t in np.linspace(0, cum[-1], n, endpoint=False)])

# contour nodes + a jittered inner ring -> dense rim like the reference
cont = along(outline, 32)
cont += rng.normal(0, 0.8, cont.shape)
pts = list(cont)

def ok(p, dmin):
    return min(np.hypot(*(p - q)) for q in pts) >= dmin

minx, miny = outline.min(0); maxx, maxy = outline.max(0)
# rim band (denser), then interior (sparser)
for radius, dmin, n_try in [(-4, 34, 2000), (-30, 44, 14000)]:
    for _ in range(n_try):
        p = np.array([rng.uniform(minx, maxx), rng.uniform(miny, maxy)])
        if radius == -4:
            if not poly.contains_point(p, radius=-4) or poly.contains_point(p, radius=-20):
                continue
        elif not poly.contains_point(p, radius=radius):
            continue
        if ok(p, dmin):
            pts.append(p)
pts = np.array(pts)

tri = Delaunay(pts)
edges = set()
for s in tri.simplices:
    for i in range(3):
        edges.add(tuple(sorted((s[i], s[(i + 1) % 3]))))
good = []
for a_, b_ in edges:
    p, q = pts[a_], pts[b_]
    d = np.hypot(*(p - q))
    if d < 100 and poly.contains_point((p + q) / 2, radius=3):
        good.append((p, q, d))
# always close the outline: link consecutive contour nodes
have = {(round(p[0],1),round(p[1],1),round(q[0],1),round(q[1],1)) for p,q,_ in good}
for i in range(len(cont)):
    p, q = pts[i], pts[(i + 1) % len(cont)]
    k1 = (round(p[0],1),round(p[1],1),round(q[0],1),round(q[1],1))
    k2 = (round(q[0],1),round(q[1],1),round(p[0],1),round(p[1],1))
    if k1 not in have and k2 not in have:
        good.append((p, q, min(np.hypot(*(p - q)), 99)))
# a few long-range "axon" links inside
inner = [p for p in pts if poly.contains_point(p, radius=-30)]
for _ in range(0):
    i, j = rng.choice(len(inner), 2, replace=False)
    p, q = inner[i], inner[j]
    d = np.hypot(*(p - q))
    if 60 < d < 140:
        good.append((p, q, d))

# ---------- SVG ----------
out = []; a = out.append
a(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">')
a('<defs>')
a(f'<clipPath id="hex"><path d="{hexd}"/></clipPath>')
a('<linearGradient id="bg" x1="0" y1="0" x2="0" y2="1"><stop offset="0" stop-color="#06163A"/><stop offset="0.55" stop-color="#0B3566"/><stop offset="1" stop-color="#0A2A4A"/></linearGradient>')
a('<filter id="blur" x="-50%" y="-50%" width="200%" height="200%"><feGaussianBlur stdDeviation="22"/></filter>')
a('<filter id="soft" x="-20%" y="-20%" width="140%" height="140%"><feGaussianBlur stdDeviation="1.2" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter>')
a('<filter id="glow" x="-50%" y="-50%" width="200%" height="200%"><feGaussianBlur stdDeviation="3.5" result="b"/><feMerge><feMergeNode in="b"/><feMergeNode in="b"/><feMergeNode in="SourceGraphic"/></feMerge></filter>')
a('<radialGradient id="core" cx="50%" cy="50%" r="50%"><stop offset="0" stop-color="#C6F28A" stop-opacity="0.85"/><stop offset="0.45" stop-color="#6DBE33" stop-opacity="0.45"/><stop offset="1" stop-color="#0B3566" stop-opacity="0"/></radialGradient>')
a('</defs>')
a(f'<path d="{hexd}" fill="url(#bg)"/>')
a('<g clip-path="url(#hex)">')

# blurred bokeh background
a('<g filter="url(#blur)">')
for _ in range(16):
    x, y = rng.uniform(0, W), rng.uniform(0, H)
    r = rng.uniform(30, 80)
    col = rng.choice(["#1B6FA8", "#2A9FB8", "#0E4C86", "#3BB3C9"])
    a(f'<circle cx="{x:.0f}" cy="{y:.0f}" r="{r:.0f}" fill="{col}" fill-opacity="{rng.uniform(0.25, 0.55):.2f}"/>')
a('</g>')

# DNA double helix in background, bottom-left -> top-right
A = np.array([62.0, 425.0]); B = np.array([446.0, 108.0])
_c = (A + B) / 2; _th = math.radians(10)  # -10 deg (math convention) = clockwise on screen
_rot = np.array([[math.cos(_th), -math.sin(_th)], [math.sin(_th), math.cos(_th)]])
A = _c + _rot @ (A - _c); B = _c + _rot @ (B - _c)
L = np.hypot(*(B - A)); uax = (B - A) / L; nax = np.array([-uax[1], uax[0]])
amp, P, Sa = 24.0, 100.0, 55.0
ts = np.linspace(0, L, 1600)
# S-shaped axis: one full sine period of lateral offset along the straight diagonal
axis = A + np.outer(ts, uax) + np.outer(Sa * np.sin(2 * np.pi * ts / L), nax)
tan = np.gradient(axis, axis=0); tan /= np.linalg.norm(tan, axis=1)[:, None]
nrm = np.c_[-tan[:, 1], tan[:, 0]]
# arc length along curved axis to keep a regular helix pitch
arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(axis, axis=0), axis=1))]
def strand(phase):
    return axis + nrm * (amp * np.sin(2 * np.pi * arc / P + phase))[:, None]
s1, s2 = strand(0), strand(np.pi)
a('<g stroke-linecap="round">')
for k, ar in enumerate(np.arange(P / 20, arc[-1], P / 10)):
    idx = np.searchsorted(arc, ar)
    idx = min(idx, len(arc) - 1)
    p1, p2 = s1[idx], s2[idx]
    col = ["#5FE1F5", "#FF6FA0", "#FFD166", "#8FB3DE"][k % 4]
    a(f'<line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="{col}" stroke-width="3" stroke-opacity="0.38"/>')
for st in (s1, s2):
    a(f'<path d="M{" L".join(f"{x:.1f},{y:.1f}" for x, y in st)}" fill="none" stroke="#8DD94A" stroke-width="5" stroke-opacity="0.75"/>')
a('</g>')

# glowing core inside the brain
a(f'<path d="M{" L".join(f"{x:.1f},{y:.1f}" for x, y in outline[::6])} Z" fill="#6DBE33" fill-opacity="0.10"/>')
a('<ellipse cx="262" cy="190" rx="165" ry="105" fill="url(#core)"/>')

# electric "lightning" filaments
a('<g filter="url(#glow)" stroke="#B5F06A" fill="none" stroke-linecap="round" stroke-linejoin="round">')
for _ in range(14):
    x, y = rng.uniform(170, 360), rng.uniform(158, 258)
    ang = rng.uniform(0, 2 * math.pi)
    d = f"M{x:.1f},{y:.1f}"
    for _ in range(rng.integers(5, 10)):
        ang += rng.normal(0, 0.7)
        step = rng.uniform(6, 14)
        nx, ny = x + step * math.cos(ang), y + step * math.sin(ang)
        if not poly.contains_point((nx, ny), radius=-25):
            break
        x, y = nx, ny
        d += f" L{x:.1f},{y:.1f}"
    a(f'<path d="{d}" stroke-width="{rng.uniform(1.0, 2.0):.1f}" stroke-opacity="0.9"/>')
a('</g>')

# mesh edges: thin, slightly curved white threads
a('<g filter="url(#soft)" stroke="#FFFFFF" fill="none" stroke-linecap="round">')
for p, q, d in good:
    m = (p + q) / 2
    nvec = np.array([-(q - p)[1], (q - p)[0]]) / max(d, 1e-6)
    m = m
    op = 0.85 - 0.35 * min(d / 100, 1)
    if d > 75:
        op = 0.22
    a(f'<path d="M{p[0]:.1f},{p[1]:.1f} Q{m[0]:.1f},{m[1]:.1f} {q[0]:.1f},{q[1]:.1f}" stroke-width="1.2" stroke-opacity="{op:.2f}"/>')
a('</g>')

# nodes
a('<g fill="#FFFFFF" filter="url(#soft)">')
for p in pts:
    a(f'<circle cx="{p[0]:.1f}" cy="{p[1]:.1f}" r="{rng.uniform(2.0, 3.2):.1f}"/>')
a('</g>')

# text
a('<text x="260" y="464" text-anchor="middle" font-family="DejaVu Sans, Helvetica, Arial, sans-serif" font-size="55" font-weight="bold" fill="#F1FAEE" letter-spacing="0">Sim<tspan fill="#8DD94A">Nano</tspan>CNS</text>')
a('<text x="260" y="500" text-anchor="middle" font-family="DejaVu Sans, Helvetica, Arial, sans-serif" font-size="19" fill="#A9CDEB" letter-spacing="1.5">nanopore BAM simulator</text>')
a('</g>')
a(f'<path d="{hexd}" fill="none" stroke="#8DD94A" stroke-width="12" stroke-linejoin="round"/>')
a('</svg>')
open(os.path.join(HERE, "simnanocns_hex.svg"), "w").write("\n".join(out))
print(len(pts), "nodes", len(good), "edges")
