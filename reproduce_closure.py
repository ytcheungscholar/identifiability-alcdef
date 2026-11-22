# -*- coding: utf-8 -*-
"""Closure runs: the fourth condition {3f0,4f0} evaluated at the SAME frequency
and with the SAME folds as the existing three-condition numbers.

mode 'freqopt'  : archive campaigns, per-fold frequency refined on training
                  nights only (identical golden-section path to
                  freqopt_v5_20261121.py, so f_hat matches the existing CSV)
                  -> out_freqopt_c34.csv  (file, app, d_c34)
mode 'u3'       : positive control, f_trial = 2 f_cat, fixed frequency
                  -> out_u3_c34.csv      (file, app, d_c34)

Only H0 and C34 are fitted (d_c34 = d_H0 - d_C34), because F_in / d_half /
d_quad / d_gam already exist in the reference CSVs.  The campaign universe and
the fold guards are copied verbatim from the reference scripts, so the merge on
(file, app) is exact.

Usage: python closure_c34_20261121.py freqopt|u3
"""
import csv
import os
import sys
import time

import numpy as np

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
KMIN, MIN_DOF, COND_MAX = 5, 5, 1e6
MODE = sys.argv[1] if len(sys.argv) > 1 else "freqopt"
if MODE == "freqopt":
    OUT = os.path.join(H, "out_freqopt_c34.csv")
elif MODE == "u3":
    OUT = os.path.join(H, "out_u3_c34.csv")
else:
    raise SystemExit("mode must be freqopt or u3")


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


def design(t, f, kind, u):
    w = 2 * np.pi * f * t
    c1, s1 = np.cos(w), np.sin(w)
    c2, s2 = np.cos(2 * w), np.sin(2 * w)
    if kind == "H0":
        return np.column_stack([c1, s1, c2, s2])
    if kind == "C34":
        return np.column_stack([c1, s1, c2, s2,
                                np.cos(3 * w), np.sin(3 * w),
                                np.cos(4 * w), np.sin(4 * w)])
    raise ValueError(kind)


def train_rss(tr_t, tr_y, tr_g, f, kind, u, pmax=8):
    if tr_t.size - len(tr_g) - pmax < MIN_DOF:
        return None
    X = design(tr_t, f, kind, u).copy()
    yv = tr_y.copy()
    for gg in tr_g:
        X[gg] = X[gg] - X[gg].mean(axis=0)
        yv[gg] = yv[gg] - yv[gg].mean()
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    if np.linalg.cond(X) > COND_MAX:
        return None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    r = yv - X @ b
    return float(r @ r), b


def holdout_rss(ho_t, ho_y, b, f, kind, u):
    Xh = design(ho_t, f, kind, u)
    r = ho_y - Xh @ b
    r = r - r.mean()
    return float(r @ r)


def golden(fun, lo, hi, it=14):
    gr = (np.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = fun(c), fun(d)
    for _ in range(it):
        if fc is None or fd is None:
            return None
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fun(d)
    return 0.5 * (a + b)


rows = [r for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"),
                                      encoding="utf-8"))
        if int(float(r["nights"])) >= KMIN]
feat = {(r["file"], r["app"]): r for r in rows}
keys = sorted(feat)
print("mode %s : campaigns in the K>=5 universe: %d" % (MODE, len(keys)), flush=True)

t0 = time.time()
n = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "K", "npt", "d_c34", "d_h0_ho"])
    for base, appk in keys:
        p = os.path.join(ROOT, base)
        if not os.path.exists(p):
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
        idx = parts[k]
        ts = t[idx]; ys = mag[idx]
        g = nights_of(ts)
        if len(g) < KMIN:
            continue
        P = float(feat[(base, appk)]["P"])
        f = (2.0 / P) if MODE == "u3" else (1.0 / P)
        y = ys - np.median(ys)
        u = (ts - ts.mean()) / max(ts.max() - ts.min(), 1e-9)
        span = float(ts.max() - ts.min())
        half_win = 0.5 / span
        acc = {"H0": 0.0, "C34": 0.0}
        ok = True
        for h in range(len(g)):
            tr = [gg for j, gg in enumerate(g) if j != h]
            ho = g[h]
            tridx = np.concatenate(tr)
            pos = np.full(ts.size, -1, dtype=int)
            pos[tridx] = np.arange(tridx.size)
            trg = [pos[gg] for gg in tr]
            trT = ts[tridx]; trY = y[tridx]; trU = u[tridx]
            if trT.size - len(trg) - 8 < MIN_DOF:
                ok = False
                break
            if MODE == "u3":
                fh = f
            else:
                def obj(ff):
                    v = train_rss(trT, trY, trg, ff, "H0", trU)
                    return None if v is None else v[0]
                fh = golden(obj, f - half_win, f + half_win)
                if fh is None or not np.isfinite(fh):
                    fh = f
            got = {}
            for kind in ("H0", "C34"):
                v = train_rss(trT, trY, trg, fh, kind, trU)
                if v is None:
                    got = None
                    break
                got[kind] = holdout_rss(ts[ho], y[ho], v[1], fh, kind, u[ho])
            if got is None:
                ok = False
                break
            for kk in acc:
                acc[kk] += got[kk]
        if not ok:
            continue
        w.writerow([base, k, len(g), int(ts.size),
                    round(acc["H0"] - acc["C34"], 5), round(acc["H0"], 5)])
        n += 1
        if n % 100 == 0:
            print("  %d done  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

print("evaluated %d  (%.1f min) -> %s" % (n, (time.time() - t0) / 60, OUT), flush=True)
