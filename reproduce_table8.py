# -*- coding: utf-8 -*-
"""Table 8, single source of truth: the exact profiled kernel on ALL cadences.

M = I - sum_k kappa_k P_Zk - v v^T/(1^T v),   v = (I - sum_k kappa_k P_Zk) 1
with kappa_k = N_k / (N_k + lambda), lambda = (sigma/sigma_zp)^2.

O(N) closed form (verified against the explicit N x N matrix to 1e-13):
  X^T M X = (X^T X - sum_k kappa_k N_k xbar_k xbar_k^T)
            - (sum_k w_k xbar_k)(sum_k w_k xbar_k)^T / sum_k w_k
  w_k = N_k lambda/(N_k + lambda)
sigma_zp -> inf gives w=0, kappa=1 (hard demeaning); sigma_zp -> 0 gives
kappa=0, w=N_k (single global offset).
"""
import csv
import os
import time

import numpy as np
import pandas as pd

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
OUT = os.path.join(H, "out_table8_full_20261121.csv")
SZP = [None, 0.020, 0.001, 0.0]
BINS = [(0, .15), (.15, .30), (.30, .50), (.50, 1.0), (1.0, 1e9)]
LAB = ["<0.15", "0.15-0.30", "0.30-0.50", "0.50-1", ">1"]


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


def r0_kernel(t, groups, f0, sig, szp):
    N = t.size
    cols = []
    for m in (1, 2, 4):
        w = 2 * np.pi * m * f0 * t
        cols += [np.sqrt(2) * np.cos(w), np.sqrt(2) * np.sin(w)]
    X = np.column_stack(cols)
    XTX = X.T @ X
    if szp is None:
        lam = 0.0
    elif szp == 0.0:
        lam = np.inf
    else:
        lam = (sig / szp) ** 2
    if szp is None:
        # hard demeaning
        S = XTX.copy()
        for gg in groups:
            xb = X[gg].sum(axis=0)
            S -= np.outer(xb, xb) / len(gg)
        G = S / N
    elif szp == 0.0:
        s = X.sum(axis=0)
        G = (XTX - np.outer(s, s) / N) / N
    else:
        kaps = np.array([len(g) for g in groups], float)
        kaps = kaps / (kaps + lam)
        S = XTX.copy()
        acc = np.zeros(6)
        wtot = 0.0
        for kk, gg in enumerate(groups):
            nk = len(gg)
            xbar = X[gg].mean(axis=0)
            S -= kaps[kk] * nk * np.outer(xbar, xbar)
            wk = nk * lam / (nk + lam)
            acc += wk * xbar
            wtot += wk
        if wtot > 1e-12:
            S -= np.outer(acc, acc) / wtot
        G = S / N

    def blk(i, j):
        return G[2 * i:2 * i + 2, 2 * j:2 * j + 2]

    G24 = np.block([[blk(1, 1), blk(1, 2)], [blk(2, 1), blk(2, 2)]])
    G12 = np.block([[blk(0, 0), blk(0, 1)], [blk(1, 0), blk(1, 1)]])
    G1_24 = np.hstack([blk(0, 1), blk(0, 2)])
    G24_1 = np.vstack([blk(1, 0), blk(2, 0)])
    G4_12 = np.hstack([blk(2, 0), blk(2, 1)])
    G12_4 = np.vstack([blk(0, 2), blk(1, 2)])
    S1 = blk(0, 0) - G1_24 @ np.linalg.pinv(G24) @ G24_1
    S4 = blk(2, 2) - G4_12 @ np.linalg.pinv(G12) @ G12_4
    a1, a4 = 0.5 * float(np.trace(S1)), 0.5 * float(np.trace(S4))
    if a1 <= 0 or a4 <= 0:
        return np.nan, np.nan
    km = np.nan
    if szp not in (None, 0.0):
        lam = (sig / szp) ** 2
        nks = np.array([len(g) for g in groups], float)
        km = float(np.median(nks / (nks + lam)))
    return float(np.sqrt(a1 / a4)), km


rows = list(csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"),
                                encoding="utf-8")))
by_file = {}
for r in rows:
    by_file.setdefault(r["file"], []).append(r)
print("campaigns", len(rows), "files", len(by_file), flush=True)

t0 = time.time()
out = []
for fi, (base, grp) in enumerate(sorted(by_file.items())):
    p = os.path.join(ROOT, base)
    if not os.path.exists(p):
        continue
    jd, mag = parse(p)
    if jd.size < 30:
        continue
    t_all = (jd - jd.min()) * 24.0
    o = np.argsort(t_all); s = t_all[o]
    cut = np.where(np.diff(s) / 24.0 > 90.0)[0] + 1
    parts = np.split(o, cut)
    for r in grp:
        k = int(r["app"])
        if k >= len(parts):
            continue
        ts = t_all[parts[k]]
        if (ts.max() - ts.min()) > 400 * 24:
            continue
        g = nights_of(ts)
        if len(g) < 3:
            continue
        f0 = 1.0 / float(r["P"])
        sig = max(float(r["sig"]), 1e-3)
        rec = dict(file=base, app=k, num=r["num"], U=r.get("U", ""),
                   f0L=float(r["f0L"]), span=float(r["span_d"]), sig=sig,
                   npt=int(ts.size), K=len(g))
        for sz in SZP:
            tag = "hard" if sz is None else ("none" if sz == 0.0 else "%.3f" % sz)
            val, km = r0_kernel(ts, g, f0, sig, sz)
            rec["r0_" + tag] = val
            if tag == "0.020":
                rec["kappa"] = km
        out.append(rec)
    if (fi + 1) % 1000 == 0:
        print("  %d/%d files  %d campaigns  %.1f min"
              % (fi + 1, len(by_file), len(out), (time.time() - t0) / 60), flush=True)

d = pd.DataFrame(out)
d.to_csv(OUT, index=False)
print("done %d campaigns in %.1f min -> %s" % (len(d), (time.time() - t0) / 60, OUT))

arch = pd.read_csv(os.path.join(H, "out_zp_baseline_20261121.csv"))
arch["app"] = arch["app"].astype(str)
d2 = d.copy(); d2["app"] = d2["app"].astype(str)
m = d2.merge(arch[["file", "app", "r0n"]], on=["file", "app"], how="left")
diff = (m.r0_hard - m.r0n).abs()
print("hard vs archive: median %.4f  p99 %.4f  max %.4f"
      % (diff.median(), diff.quantile(.99), diff.max()))
bad = m.loc[diff > 0.05, ["file", "app", "f0L", "npt", "K", "r0_hard", "r0n"]]
print("campaigns differing by >0.05: %d" % len(bad))
print(bad.head(8).to_string(index=False))

d["b"] = pd.cut(d.f0L, [x[0] for x in BINS] + [1e9], labels=LAB, right=False)
d["clean"] = d.U.fillna("").astype(str).str.strip().isin(["3", "3-"]) & (d.span <= 30)
gb = d.groupby("b", observed=True).agg(
    n=("sig", "size"), sig=("sig", "median"), Nk=("npt", "median"),
    hard=("r0_hard", "median"), p020=("r0_0.020", "median"),
    p001=("r0_0.001", "median"), kap=("kappa", "median"))
gb["gain%"] = 100 * (gb.p020 / gb.hard - 1)
print(gb.round(3).to_string())
print("\nclean subset (U>=3, span<=30 d):")
gc = d[d.clean].groupby("b", observed=True).agg(
    n=("sig", "size"), hard=("r0_hard", "median"), p020=("r0_0.020", "median"))
gc["gain%"] = 100 * (gc.p020 / gc.hard - 1)
print(gc.round(3).to_string())
