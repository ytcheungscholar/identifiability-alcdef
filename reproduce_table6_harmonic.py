# -*- coding: utf-8 -*-
"""In-sample rate for the fourth condition, to fill the missing row of Table 6.

For every campaign of the LONO sample, fit

    H0  = {f0, 2f0}             (4 params + one offset per night)
    H34 = {f0, 2f0, 3f0, 4f0}   (8 params + one offset per night)

and count the fraction whose in-sample F exceeds the exact F_crit(4, N-K-8).
Hold-out improvement (46.7%) comes from the existing LONO run.

Output: out_insample_c34_20261121.csv
"""
import csv
import os
import time

import numpy as np
from scipy.stats import f as fdist

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
OUT = os.path.join(H, "out_insample_c34_20261121.csv")
KEYS = os.path.join(H, "out_lono_v5_c34.csv")


def parse(path):
    jd, mag = [], []
    for line in open(path, "r", errors="ignore"):
        if line.startswith("DATA="):
            p = line[5:].strip().split("|")
            if len(p) >= 2:
                try:
                    j = float(p[0]); m = float(p[1])
                except ValueError:
                    continue
                if j > 2.0e6 and -5 < m < 30:
                    jd.append(j); mag.append(m)
    return np.array(jd), np.array(mag)


def nights_of(t, gap=7.2):
    o = np.argsort(t); s = t[o]
    cut = np.where(np.diff(s) > gap)[0] + 1
    return [g for g in np.split(o, cut) if len(g) >= 3]


def design(t, f, kind):
    w = 2 * np.pi * f * t
    c1, s1 = np.cos(w), np.sin(w)
    c2, s2 = np.cos(2 * w), np.sin(2 * w)
    if kind == "H0":
        return np.column_stack([c1, s1, c2, s2])
    return np.column_stack([c1, s1, c2, s2,
                            np.cos(3 * w), np.sin(3 * w),
                            np.cos(4 * w), np.sin(4 * w)])


def rss(t, y, g, f, kind):
    X = design(t, f, kind).copy()
    yv = y.copy()
    for gg in g:
        X[gg] = X[gg] - X[gg].mean(axis=0)
        yv[gg] = yv[gg] - yv[gg].mean()
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    r = yv - X @ b
    return float(r @ r)


rows = [r for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"),
                                      encoding="utf-8"))
        if int(float(r["nights"])) >= 5]
feat = {(r["file"], r["app"]): r for r in rows}
keys = [(r["file"], r["app"]) for r in
        csv.DictReader(open(KEYS, encoding="utf-8"))]
print("campaigns to process:", len(keys), flush=True)

t0 = time.time()
n = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "K", "npt", "F_in", "F_crit", "reject"])
    for base, appk in keys:
        p = os.path.join(ROOT, base)
        if not os.path.exists(p) or (base, appk) not in feat:
            continue
        jd, mag = parse(p)
        if jd.size < 30:
            continue
        t = (jd - jd.min()) * 24.0
        o = np.argsort(t); s = t[o]
        cut = np.where(np.diff(s) / 24.0 > 90.0)[0] + 1
        k = int(appk)
        parts = np.split(o, cut)
        if k >= len(parts):
            continue
        ts = t[parts[k]]; ys = mag[parts[k]]
        g = nights_of(ts)
        if len(g) < 5:
            continue
        f = 1.0 / float(feat[(base, appk)]["P"])
        y = ys - np.median(ys)
        r0 = rss(ts, y, g, f, "H0")
        r1 = rss(ts, y, g, f, "H34")
        if r0 is None or r1 is None or r1 <= 0:
            continue
        dof = max(ts.size - 8 - len(g), 1)
        F = ((r0 - r1) / 4.0) / (r1 / dof)
        Fc = float(fdist.ppf(0.95, 4, dof))
        w.writerow([base, k, len(g), int(ts.size), round(F, 4), round(Fc, 4),
                    int(F > Fc)])
        n += 1
        if n % 1000 == 0:
            print("  %d  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

L = list(csv.DictReader(open(OUT, encoding="utf-8")))
F = np.array([float(r["F_in"]) for r in L])
R = np.array([int(r["reject"]) for r in L])
print("\nfourth condition, in sample:  %d/%d = %.1f%%  (%.1f min)"
      % (R.sum(), len(R), 100 * R.mean(), (time.time() - t0) / 60))
print("median F_in %.2f  median F_crit %.3f"
      % (np.median(F), np.median([float(r["F_crit"]) for r in L])))
