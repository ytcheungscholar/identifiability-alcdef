# -*- coding: utf-8 -*-
"""Table 2 v2: dense/refined frequency search, two noise regimes.

The 180-point grid of the first version was too coarse for M=4 at 6 nights:
the four-harmonic peak is only 1/(4T) wide. Here the scan is 600 points plus a
golden-section refinement, which is enough to remove the grid artefact.

Regime A (high S/N):  6 nights x 30 points, sigma = 0.05 mag
Regime B (low S/N):   5 nights x 12 points, sigma = 0.15 mag
Both: P ~ U(6.5,11) h, zero points N(0,0.25^2), 2% tolerance, 40 trials/cell.
"""
import csv
import os
import time

import numpy as np

H = r"C:\Users\jacky\.openclaw\workspace\ac_repo_plan"
OUT = os.path.join(H, "out_table2_n200.csv")
SIGMA_ZP = 0.25
A1, A2, A3 = 0.10, 0.25, 0.03
R4S = [0.2, 0.5, 0.8]
MS = [1, 2, 3, 4]
NTRIAL = 200
TOL = 0.02
FGRID = np.linspace(1 / 16.0, 1 / 4.0, 600)


def gen_trial(rng, K, NPTS, night_h, r4):
    P = rng.uniform(6.5, 11.0)
    ph = rng.uniform(0, 2 * np.pi, 4)
    t, g, idx = [], [], 0
    for k in range(K):
        start = k * 24.0 + rng.uniform(0, 1.0)
        t.append(start + np.linspace(0, night_h, NPTS))
        g.append(np.arange(idx, idx + NPTS))
        idx += NPTS
    t = np.concatenate(t)
    y = np.zeros_like(t)
    for m, a in enumerate([A1, A2, A3, r4 * A1]):
        y += a * np.cos(2 * np.pi * (m + 1) * t / P + ph[m])
    z = rng.normal(0, SIGMA_ZP, K)
    for k, gg in enumerate(g):
        y[gg] += z[k]
    return t, g, y, P


def fitted_period(t, g, y, f, M, offsets):
    cols = []
    for m in range(1, M + 1):
        w = 2 * np.pi * m * f * t
        cols += [np.cos(w), np.sin(w)]
    X = np.column_stack(cols)
    if offsets:
        Z = np.zeros((t.size, len(g)))
        for j, gg in enumerate(g):
            Z[gg, j] = 1.0
        X = np.hstack([X, Z])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    return float(r @ r)


def best_period(t, g, y, M, offsets):
    rss = [fitted_period(t, g, y, f, M, offsets) for f in FGRID]
    i = int(np.argmin(rss))
    lo = FGRID[max(i - 1, 0)]
    hi = FGRID[min(i + 1, len(FGRID) - 1)]
    gr = (np.sqrt(5) - 1) / 2
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = fitted_period(t, g, y, c, M, offsets), fitted_period(t, g, y, d, M, offsets)
    for _ in range(14):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = fitted_period(t, g, y, c, M, offsets)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = fitted_period(t, g, y, d, M, offsets)
    return 0.5 * (a + b)


rows = []
t0 = time.time()
for regime, (K, NPTS, night_h, sigma) in (("A", (6, 30, 4.0, 0.05)),
                                          ("B", (5, 12, 4.0, 0.15))):
    rng = np.random.default_rng(2026 if regime == "A" else 2027)
    for r4 in R4S:
        hit = {(M, off): 0 for M in MS for off in (0, 1)}
        for _ in range(NTRIAL):
            t, g, y, P = gen_trial(rng, K, NPTS, night_h, r4)
            y = y + rng.normal(0, sigma, t.size)
            for off in (0, 1):
                for M in MS:
                    Pf = 1.0 / best_period(t, g, y, M, off)
                    if abs(Pf / P - 1.0) < TOL:
                        hit[(M, off)] += 1
        for M in MS:
            p0 = 100.0 * hit[(M, 0)] / NTRIAL
            p1 = 100.0 * hit[(M, 1)] / NTRIAL
            rows.append(dict(regime=regime, r4=r4, M=M, no_off=p0, nightly=p1))
            print("%s r4=%.1f M=%d  none %5.1f  nightly %5.1f  d=%+5.1f  (%.1f min)"
                  % (regime, r4, M, p0, p1, p1 - p0, (time.time() - t0) / 60), flush=True)

with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("->", OUT)
