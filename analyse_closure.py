# -*- coding: utf-8 -*-
"""Closure analysis for v8.

1. four-condition screen at the refined frequency (archive)
2. four-condition screen on the U=3 positive control
3. cluster-bootstrap 95% CIs for the Table-7 rates (three- and four-condition)
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import f as fdist

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
DEN = 6657.0          # denominator used in the paper (guarded campaigns)


def screen(df, fcol="F_in", extra=None, denom=DEN):
    d = df.copy()
    d["dof"] = np.maximum(d.npt - 8 - d.K, 1)
    Fc = fdist.ppf(0.95, 4, d.dof)
    p3 = (d[fcol] > Fc) & (d.d_half > 0) & (d.d_half > d.d_quad) & (d.d_half > d.d_gam)
    out = {"n": len(d), "pass3": int(p3.sum()), "rate3": 100 * p3.sum() / denom,
           "rate3_sample": 100 * p3.mean()}
    if extra is not None:
        p4 = p3 & (d.d_half > d[extra])
        out.update({"pass4": int(p4.sum()), "rate4": 100 * p4.sum() / denom,
                    "rate4_sample": 100 * p4.mean(),
                    "beaten": int((p3 & (d.d_half <= d[extra])).sum())})
    return out, p3


# ---------- 1. archive, catalogue frequency vs refined frequency ----------
cat = pd.read_csv(os.path.join(H, "out_lono_v5_c34.csv"))
print("=== archive, catalogue frequency (sanity: should be 332 / 4.99% / 3.35%)")
r, _ = screen(cat, extra="d_c34")
print(r)

fo = pd.read_csv(os.path.join(H, "out_lono_v5_freqopt.csv"))
fo34 = pd.read_csv(os.path.join(H, "out_freqopt_c34.csv"))
fom = fo.merge(fo34[["file", "app", "d_c34"]], on=["file", "app"], how="inner")
print("\n=== archive, per-fold refined frequency   merged n=%d (freqopt n=%d)"
      % (len(fom), len(fo)))
r, _ = screen(fom, extra="d_c34")
print(r)

# ---------- 2. U=3 positive control ----------
u = pd.read_csv(os.path.join(H, "out_u3_control_20261121.csv"))
u34 = pd.read_csv(os.path.join(H, "out_u3_c34.csv"))
um = u.merge(u34[["file", "app", "d_c34"]], on=["file", "app"], how="inner")
print("\n=== U=3 positive control, f_trial = 2 f_cat   n=%d (sanity: 1,264 / 55.5%%)"
      % len(um))
NU3 = float(len(um))
r, _ = screen(um, extra="d_c34", denom=NU3)
print(r)

# ---------- 3. cluster bootstrap on the Table-7 rates ----------
LAB = dict(A="A true 2P, stationary      ", B="B true P, amplitude modul. ",
           C="C true P, four-harmonic drift", D="D true 2P + drift         ",
           E="E true P, stationary (null) ", F="F strict H0, A3=A4=0     ")
rng = np.random.default_rng(20261121)
for suf, a1 in (("_a05", 0.05), ("_a10", 0.10)):
    z = np.load(os.path.join(H, "calib_v5b%s.npz" % suf))
    n = len(z["acc1_A"])
    print("\n=== Table 7 cluster bootstrap, A1 = %.2f mag, %d campaigns x 8 reps" % (a1, n))
    for arm in "ABCDEF":
        for tag, ac, rp in (("3-cond", z["acc1_" + arm], z["rep1_" + arm]),
                            ("4-cond", z["acc4_" + arm], z["rep1_" + arm])):
            a = ac.sum(); b = rp.sum()
            pt = 100 * a / max(b, 1)
            bs = np.empty(2000)
            for i in range(2000):
                idx = rng.integers(0, n, n)
                bs[i] = 100 * ac[idx].sum() / max(rp[idx].sum(), 1)
            lo, hi = np.percentile(bs, [2.5, 97.5])
            print("  %s %s  %5.1f%%  [%4.1f, %5.1f]   %d/%d"
                  % (LAB[arm], tag, pt, lo, hi, int(a), int(b)))

# ---------- 4. Table 8 medians per six bins, for the figure caption ----------
t8 = pd.read_csv(os.path.join(H, "out_table8_full_20261121.csv"))
B6 = [(0, .05), (.05, .15), (.15, .30), (.30, .50), (.50, 1.0), (1.0, 1e9)]
print("\n=== Table 8, six bins (medians), for the Fig. 4 caption")
for a, b in B6:
    s = t8[(t8.f0L >= a) & (t8.f0L < b)]
    print("  %5.2f-%-5.2f n=%4d  sigma=%6.3f  kappa=%6.3f  hard=%5.3f  p02=%5.3f  p001=%5.3f"
          % (a, b, len(s), s.sig.median(), s.kappa.median(), s.r0_hard.median(),
             s["r0_0.020"].median(), s["r0_0.001"].median()))
