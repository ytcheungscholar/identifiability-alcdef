# -*- coding: utf-8 -*-
"""Positive control on real data: pretend the catalogue reports half the period.

For LCDB U=3 objects (high-confidence, so P_cat is very likely the true period)
we feed the screen f_trial = 2 f_cat, i.e. exactly the situation in which the
catalogue locked onto P/2.  Then

    H0 = {2f, 4f}                 (the half-period model)
    H1 = {f, 2f, 3f, 4f}          (the truth, since harmonics sit at f..4f)

so a pass = the screen recovers the true period.  This is a recall measurement
that uses real noise and real phase coverage, with no synthetic drift at all.

Output: out_u3_control_20261121.csv
"""
import csv
import glob
import os
import time

import numpy as np

H = r"C:\Users\jacky\.openclaw\workspace\ac_repo_plan"
ROOT = r"C:\Users\jacky\.openclaw\workspace\alcdef_data"
OUT = os.path.join(H, "out_u3_control_20261121.csv")
KMIN, MIN_DOF, COND_MAX = 5, 5, 1e6
FACTOR = 2.0          # trial frequency = FACTOR * f_catalogue


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


def fit_insample(t, y, g, f, kind, u):
    X = design(t, f, kind, u)
    yv = y - y.mean()
    for gg in g:
        X[gg] = X[gg] - X[gg].mean(axis=0)
        yv[gg] = yv[gg] - yv[gg].mean()
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None, None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    r = yv - X @ b
    return float(r @ r), b


def feval_score(tr_t, tr_y, tr_g, ho_t, ho_y, f, kind, u_tr, u_ho, pmax):
    nt = tr_t.size
    if nt - len(tr_g) - pmax < MIN_DOF:
        return None
    X = design(tr_t, f, kind, u_tr)
    yv = tr_y.copy()
    for gg in tr_g:
        X[gg] = X[gg] - X[gg].mean(axis=0)
        yv[gg] = yv[gg] - yv[gg].mean()
    if np.linalg.matrix_rank(X) < X.shape[1] or np.linalg.cond(X) > COND_MAX:
        return None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    Xh = design(ho_t, f, kind, u_ho)
    r = ho_y - Xh @ b
    r = r - r.mean()
    return float(r @ r)


rows = [r for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"),
                                      encoding="utf-8"))
        if int(float(r["nights"])) >= KMIN and str(r.get("U", "")).strip() == "3"]
feat = {(r["file"], r["app"]): r for r in rows}
keep = {k[0] for k in feat}
print("U=3 campaigns with K>=%d: %d  (files %d)" % (KMIN, len(feat), len(keep)), flush=True)

t0 = time.time()
n = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "num", "K", "npt", "span_d", "f0L", "F_in",
                "d_half", "d_quad", "d_gam"])
    for p in sorted(glob.glob(os.path.join(ROOT, "*.txt"))):
        base = os.path.basename(p)
        if base not in keep:
            continue
        try:
            jd, mag = parse(p)
        except Exception:
            continue
        if jd.size < 30:
            continue
        t = (jd - jd.min()) * 24.0
        o = np.argsort(t); s = t[o]
        cut = np.where(np.diff(s) / 24.0 > 90.0)[0] + 1
        for k, idx in enumerate(np.split(o, cut)):
            key = (base, str(k))
            if key not in feat:
                continue
            ts = t[idx]; ys = mag[idx]
            g = nights_of(ts)
            if len(g) < KMIN:
                continue
            f = FACTOR / float(feat[key]["P"])
            y = ys - np.median(ys)
            u = (ts - ts.mean()) / max(ts.max() - ts.min(), 1e-9)
            r0, _ = fit_insample(ts, y, g, f, "H0", u)
            r1, _ = fit_insample(ts, y, g, f, "H1", u)
            if r0 is None or r1 is None or r1 <= 0:
                continue
            dof1 = max(ts.size - 8 - len(g), 1)
            F_in = ((r0 - r1) / 4.0) / (r1 / dof1)
            d, ok, sk = {}, True, 0
            for kind in ("H0", "H1", "Cq", "Cg"):
                tot = 0.0
                for h in range(len(g)):
                    tr = [gg for j, gg in enumerate(g) if j != h]
                    mat = np.zeros(ts.size, bool); mat[g[h]] = True
                    tridx = np.concatenate(tr)
                    pos = np.full(ts.size, -1, dtype=int)
                    pos[tridx] = np.arange(tridx.size)
                    trg = [pos[gg] for gg in tr]
                    v = feval_score(ts[tridx], y[tridx], trg, ts[mat], y[mat],
                                    f, kind, u[tridx], u[mat], 8)
                    if v is None:
                        sk += 1; tot = np.nan
                    else:
                        tot += v
                if not np.isfinite(tot):
                    ok = False
                    break
                d[kind] = tot
            if not ok:
                continue
            w.writerow([base, k, feat[key]["num"], len(g), int(ts.size),
                        feat[key]["span_d"], feat[key]["f0L"], round(F_in, 4),
                        round(d["H0"] - d["H1"], 5),
                        round(d["H0"] - d["Cq"], 5),
                        round(d["H0"] - d["Cg"], 5)])
            n += 1
        if n and n % 500 == 0:
            print("  %d  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

L = list(csv.DictReader(open(OUT, encoding="utf-8")))
F = np.array([float(r["F_in"]) for r in L])
dh = np.array([float(r["d_half"]) for r in L])
dq = np.array([float(r["d_quad"]) for r in L])
dg = np.array([float(r["d_gam"]) for r in L])
print("\nU=3 positive control, f_trial = 2 f_cat  (n=%d, %.1f min)" % (len(L), (time.time() - t0) / 60))
for nm, v in (("{f/2,3f/2}", dh), ("quartic", dq), ("gamma12", dg)):
    print("  %-12s improves %.1f%%" % (nm, 100 * (v > 0).mean()))
PS = (F > 2.432) & (dh > 0) & (dh > dq) & (dh > dg)
print("  full screen pass (recall) %.1f%%  (%d/%d)" % (100 * PS.mean(), PS.sum(), len(L)))
print("  in-sample F>2.432 alone    %.1f%%" % (100 * (F > 2.432).mean()))
print("->", OUT)
