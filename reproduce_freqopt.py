# -*- coding: utf-8 -*-
"""Step 2: does gamma12 survive when the period is re-optimised inside each fold?

Threat being tested: cos(2*pi*(f+df)*t) ~ cos - 2*pi*df*t*sin, so a fixed
catalogue frequency that is slightly wrong looks exactly like the gamma12
basis.  Fix: in every LONO fold, refine f on the TRAINING nights only
(f_hat in f_cat +/- 1/(2 T_train)), then evaluate every model at f_hat on the
held-out night.  Never optimise on the full campaign -- that would leak the
held-out phase.

Usage:  python freqopt_v5_20261121.py [nsample] [seed]
Output: out_lono_v5_freqopt.csv
"""
import csv
import glob
import os
import sys
import time

import numpy as np

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
OUT = os.path.join(H, "out_lono_v5_freqopt.csv")
KMIN, MIN_DOF, COND_MAX = 5, 5, 1e6
NS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 12345
MODE = sys.argv[3] if len(sys.argv) > 3 else "opt"      # opt | fix
if MODE == "fix":
    OUT = os.path.join(H, "out_lono_v5_fixref_20261121.csv")


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
    if kind == "H1":
        return np.column_stack([np.cos(w / 2), np.sin(w / 2),
                                np.cos(1.5 * w), np.sin(1.5 * w), c1, s1, c2, s2])
    if kind == "Cq":
        return np.column_stack([c1, s1, c2, s2, u, u ** 2, u ** 3, u ** 4])
    if kind == "Cg":
        return np.column_stack([c1, s1, c2, s2, u * c1, u * s1, u * c2, u * s2])
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
    r = r - r.mean()          # the held-out night keeps its own free zero point
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
rng = np.random.default_rng(SEED)
if NS < len(keys):
    pick = rng.choice(len(keys), size=NS, replace=False)
    keys = [keys[i] for i in pick]
keys = sorted((str(a), str(b)) for a, b in keys)
print("campaigns to process: %d" % len(keys), flush=True)

t0 = time.time()
n = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "K", "npt", "span_d", "f_shift_rel",
                "F_in", "d_half", "d_quad", "d_gam", "rss_h0_ho", "d_half_old",
                "d_gam_old"])
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
        f = 1.0 / float(feat[(base, appk)]["P"])
        y = ys - np.median(ys)
        u = (ts - ts.mean()) / max(ts.max() - ts.min(), 1e-9)
        span = float(ts.max() - ts.min())
        half_win = 0.5 / span
        acc = {kk: 0.0 for kk in ("H0", "H1", "Cq", "Cg")}
        shifts, ok = [], True
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

            def obj(ff):
                v = train_rss(trT, trY, trg, ff, "H0", trU)
                return None if v is None else v[0]

            if MODE == "fix":
                fh = f
                shifts.append(0.0)
            else:
                fh = golden(obj, f - half_win, f + half_win)
                if fh is None or not np.isfinite(fh):
                    fh = f
                shifts.append((fh - f) / f)
            got = {}
            for kind in ("H0", "H1", "Cq", "Cg"):
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
        fbar = float(np.median([f * (1 + s) for s in shifts]))
        gall = list(range(len(g)))
        v0 = train_rss(ts, y, gall, fbar, "H0", u)
        v1 = train_rss(ts, y, gall, fbar, "H1", u)
        Fin = 0.0
        if v0 and v1 and v1[0] > 0:
            dof1 = max(ts.size - 8 - len(g), 1)
            Fin = ((v0[0] - v1[0]) / 4.0) / (v1[0] / dof1)
        w.writerow([base, k, len(g), int(ts.size), round(span, 3),
                    round(float(np.median(shifts)), 6),
                    round(Fin, 4),
                    round(acc["H0"] - acc["H1"], 5),
                    round(acc["H0"] - acc["Cq"], 5),
                    round(acc["H0"] - acc["Cg"], 5),
                    round(acc["H0"], 5), "", ""])
        n += 1
        if n % 50 == 0:
            print("  %d done  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

print("evaluated %d  (%.1f min) -> %s" % (n, (time.time() - t0) / 60, OUT), flush=True)

L = list(csv.DictReader(open(OUT, encoding="utf-8")))
dh = np.array([float(r["d_half"]) for r in L])
dq = np.array([float(r["d_quad"]) for r in L])
dg = np.array([float(r["d_gam"]) for r in L])
print("with per-fold frequency refinement (n=%d):" % len(L))
for nm, v in (("{f/2,3f/2}", dh), ("quartic", dq), ("gamma12", dg)):
    print("  %-12s improves %.1f%%  median %+.4f" % (nm, 100 * (v > 0).mean(), np.median(v)))
print("  paired gam>half %.1f%%  half>gam %.1f%%"
      % (100 * (dg > dh).mean(), 100 * (dh > dg).mean()))
