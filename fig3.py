# -*- coding: utf-8 -*-
"""v9 Fig. 3 -- panel (a) now carries the fourth extension.

(a) four extensions of the same 4 -> 8 model on the same campaigns: in-sample
    flag rate (exact per-campaign F_crit) against held-out improvement.
    The higher-harmonic control {3f0,4f0} is the one that also predicts the
    withheld nights, which is why it becomes the fourth condition.
(b) calibration at the two reference amplitudes under the three conditions,
    with the archive rate marked (unchanged).
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import f as fdist

H = r"C:\Users\jacky\.openclaw\workspace\ac_repo_plan"
OUT = r"C:\Users\jacky\.openclaw\workspace\v2_manuscript\figures\fig3_lono.png"

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 8.5,
                     "xtick.labelsize": 7.4, "ytick.labelsize": 7.2,
                     "legend.fontsize": 6.8})

ins = pd.read_csv(H + r"\out_insample_v5.csv")
lon = pd.read_csv(H + r"\out_lono_v5_c34.csv")
insc = pd.read_csv(H + r"\out_insample_c34_20261121.csv")
for d in (ins, lon, insc):
    d["app"] = d["app"].astype(str)
lon["Fcrit"] = fdist.ppf(0.95, 4, (lon.npt - lon.K - 8).clip(lower=1))
m = ins.merge(lon[["file", "app", "Fcrit", "d_half", "d_quad", "d_gam", "d_c34"]],
              on=["file", "app"], how="inner")
m = m.merge(insc[["file", "app", "reject"]], on=["file", "app"], how="inner")

names = [r"$\{f_0/2,3f_0/2\}$", "quartic", r"$\gamma_{12}$", r"$\{3f_0,4f_0\}$"]
insr = [100 * float((m[c] > m.Fcrit).mean()) for c in ("F_half", "F_quad", "F_gam")]
insr.append(100 * float(m["reject"].mean()))
held = [100 * float((m[c] > 0).mean()) for c in ("d_half", "d_quad", "d_gam", "d_c34")]
print("n =", len(m))
print("in-sample", np.round(insr, 1), " held-out", np.round(held, 1))

fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.95))
n = len(names)
xx = np.arange(n); w = 0.36
ax[0].bar(xx - w / 2, insr, w, color="#C44E52", label=r"in-sample $F>F_{\rm crit}$")
ax[0].bar(xx + w / 2, held, w, color="#4C72B0", label="improves held-out nights")
for i in range(n):
    ax[0].text(xx[i] - w / 2, insr[i] + 1.4, "%.1f" % insr[i], ha="center", fontsize=6.8)
    ax[0].text(xx[i] + w / 2, held[i] + 1.4, "%.1f" % held[i], ha="center", fontsize=6.8)
ax[0].set_xticks(xx); ax[0].set_xticklabels(names)
ax[0].set_ylabel("per cent of campaigns")
ax[0].set_ylim(0, 115)
ax[0].grid(alpha=0.22, axis="y")
ax[0].legend(loc="upper left", framealpha=0.9)
ax[0].set_title("(a) same 4$\\to$8 model, %s campaigns" % format(len(m), ","))

arms = ["A", "E", "F", "B", "C", "D"]
labels = ["$2P$\nstat.", "$P$\nstat.", "strict\n$H_0$", "$P$\namp-mod", "$P$\ndrift", "$2P$\n+drift"]
res = {}
for suf, a1v in (("_a05", 0.05), ("_a10", 0.10)):
    d = np.load(H + r"\calib_v5%s.npz" % suf, allow_pickle=True)
    r1 = [100 * d["acc1_" + a].sum() / d["rep1_" + a].sum() for a in arms]
    res[a1v] = r1
    print("A1=%.2f" % a1v, np.round(r1, 1))

x2 = np.arange(len(arms)); w2 = 0.38
ax[1].bar(x2 - w2 / 2, res[0.05], w2, color="#8C8C8C", label=r"$A_1=0.05$ mag")
ax[1].bar(x2 + w2 / 2, res[0.10], w2, color="#4C72B0", label=r"$A_1=0.10$ mag")
for i in range(len(arms)):
    ax[1].text(x2[i] - w2 / 2, res[0.05][i] + 0.9, "%.1f" % res[0.05][i], ha="center", fontsize=6.4)
    ax[1].text(x2[i] + w2 / 2, res[0.10][i] + 0.9, "%.1f" % res[0.10][i], ha="center", fontsize=6.4)
ax[1].axhline(4.99, color="#C44E52", ls="--", lw=0.9)
ax[1].set_xlim(-0.62, 5.45)
ax[1].set_xticks(x2); ax[1].set_xticklabels(labels)
ax[1].set_ylabel("accepted [%]")
ax[1].set_ylim(0, 78)
ax[1].grid(alpha=0.22, axis="y")
ax[1].legend(loc="upper center", framealpha=0.9, ncol=2, bbox_to_anchor=(0.5, 0.96))
ax[1].set_title(r"(b) calibration, six arms")

fig.tight_layout()
fig.savefig(OUT, dpi=220)
print("wrote", OUT)
