# -*- coding: utf-8 -*-
"""v5 §4.2: unified 4 -> 8 nested test on the physical sample with K >= 5.

Definitions (identical everywhere in v5)
  H0  = {f, 2f}                      4 parameters, + K nightly offsets
  H1  = {f/2, f, 3f/2, 2f}           8 parameters, + K nightly offsets  (dof = 4)
  Cq  = H0 + {u, u^2, u^3, u^4}      8 parameters, 4 extra  (flexibility control)
  Cg  = H0 + {u cos w, u sin w, u cos 2w, u sin 2w}  8 parameters, 4 extra
        (drift of the fundamental and of the second harmonic)

Sample: physical selection with K >= 5 nights (8,896 campaigns). No permutation
cut, no peak-to-peak cut -- this is the unbiased domain.

Guards (per held-out fold): the training block must satisfy
    N_train - (K-1) - p_max >= 5,  full column rank,  cond < 1e6
otherwise the fold is skipped and counted; lstsq is never allowed to return a
minimum-norm pseudo-solution silently.

Output: out_lono_v5.csv  (file, app, K, npt, F_in, d_half, d_quad, d_gam, folds_skipped)
"""
import os
import csv, os, glob, time
import numpy as np

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
OUT = os.path.join(H, "out_lono_v5_c34.csv")
KMIN = 5
MIN_DOF = 5
COND_MAX = 1e6


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
        X = np.column_stack([c1, s1, c2, s2])
    elif kind == "H1":
        X = np.column_stack([np.cos(w / 2), np.sin(w / 2),
                             np.cos(1.5 * w), np.sin(1.5 * w), c1, s1, c2, s2])
    elif kind == "Cq":
        X = np.column_stack([c1, s1, c2, s2, u, u ** 2, u ** 3, u ** 4])
    elif kind == "C34":
        X = np.column_stack([c1, s1, c2, s2,
                             np.cos(3 * w), np.sin(3 * w),
                             np.cos(4 * w), np.sin(4 * w)])
    elif kind == "Cg":
        X = np.column_stack([c1, s1, c2, s2,
                             u * c1, u * s1, u * c2, u * s2])
    else:
        raise ValueError(kind)
    return X


def fit_insample(t, y, g, f, kind, u):
    X = design(t, f, kind, u)
    yv = y - y.mean()
    for gg in g:
        X[gg] = X[gg] - X[gg].mean(axis=0)
        yv[gg] = yv[gg] - yv[gg].mean()
    p = X.shape[1]
    if np.linalg.matrix_rank(X) < p:
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
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    if np.linalg.cond(X) > COND_MAX:
        return None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    Xh = design(ho_t, f, kind, u_ho)
    r = ho_y - Xh @ b
    r = r - r.mean()
    return float(r @ r)


rows = [r for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"), encoding="utf-8"))
        if int(float(r["nights"])) >= KMIN]
print("physics sample with K>=%d: %d" % (KMIN, len(rows)), flush=True)
keep = set()
feat = {(r["file"], r["app"]): r for r in rows}
for k in feat:
    keep.add(k[0])

t0 = time.time()
n = 0
skip_fold_tot = 0
skip_camp = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "K", "npt", "F_in", "d_half", "d_quad", "d_gam", "folds_skipped", "d_c34"])
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
            f = 1.0 / float(feat[key]["P"])
            y = ys - np.median(ys)
            u = (ts - ts.mean()) / max(ts.max() - ts.min(), 1e-9)
            r0, _ = fit_insample(ts, y, g, f, "H0", u)
            r1, _ = fit_insample(ts, y, g, f, "H1", u)
            if r0 is None or r1 is None or r1 <= 0:
                skip_camp += 1
                continue
            dof1 = max(ts.size - 8 - len(g), 1)
            F_in = ((r0 - r1) / 4.0) / (r1 / dof1)
            d = {}
            ok = True
            skipped = 0
            for kind in ("H0", "H1", "Cq", "Cg", "C34"):
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
                        skipped += 1
                        tot += np.nan
                    else:
                        tot += v
                if not np.isfinite(tot):
                    ok = False
                    break
                d[kind] = tot
            skip_fold_tot += skipped
            if not ok:
                skip_camp += 1
                continue
            w.writerow([base, k, len(g), int(ts.size), round(F_in, 4),
                        round(d["H0"] - d["H1"], 5),
                        round(d["H0"] - d["Cq"], 5),
                        round(d["H0"] - d["Cg"], 5), skipped,
                        round(d["H0"] - d["C34"], 5)])
            n += 1
        if n and n % 1000 == 0:
            print("  %d campaigns  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

print("evaluated %d campaigns (%d skipped), folds skipped %d, time %.1f min"
      % (n, skip_camp, skip_fold_tot, (time.time() - t0) / 60), flush=True)
print("->", OUT, flush=True)

L = list(csv.DictReader(open(OUT, encoding="utf-8")))
F = np.array([float(r["F_in"]) for r in L])
dh = np.array([float(r["d_half"]) for r in L])
dq = np.array([float(r["d_quad"]) for r in L])
dg = np.array([float(r["d_gam"]) for r in L])
d34 = np.array([float(r["d_c34"]) for r in L])
for name, v in (("{f/2,3f/2}", dh), ("quartic", dq), ("gamma12", dg), ("{3f,4f}", d34)):
    print("  %-12s improves %5d (%.1f%%)  median gain %+.4f"
          % (name, int((v > 0).sum()), 100 * np.mean(v > 0), np.median(v)))
PS = (F > 2.432) & (dh > 0) & (dh > dq) & (dh > dg)
print("  pass (F_in>2.432 & gain>0 & beats both): %d = %.1f%%" % (PS.sum(), 100 * PS.mean()))
PS2 = (F > 2.432) & (dh > 0) & (dh > dq) & (dh > dg) & (dh > d34)
print("  pass if the test must ALSO beat {3f,4f}: %d = %.2f%%" % (PS2.sum(), 100 * PS2.mean()))
print("  of the %d passing campaigns, %d still beat {3f,4f}"
      % (PS.sum(), (PS & (dh > d34)).sum()))
