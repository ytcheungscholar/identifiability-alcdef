# -*- coding: utf-8 -*-
"""v8 Figure 4, six bins -- identical to the bins of Table 8, with the 0/0
regime of the lowest bin shaded so it is not read as a recovery measurement."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = r"C:\Users\jacky\.openclaw\workspace\ac_repo_plan"
OUT = r"C:\Users\jacky\.openclaw\workspace\v2_manuscript\figures\fig4.png"

d = pd.read_csv(H + r"\out_table8_full_20261121.csv")
BINS = [(0, .05), (.05, .15), (.15, .30), (.30, .50), (.50, 1.0), (1.0, 1e9)]
lab = ["$<0.05$", "0.05\u20130.15", "0.15\u20130.30", "0.30\u20130.50",
       "0.50\u20131", "$>1$"]

plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9.5,
                     "legend.fontsize": 7.6})
fig, ax = plt.subplots(figsize=(5.2, 3.5))
x = np.arange(len(BINS))
series = [("r0_hard", "hard demeaning ($\\kappa=1$)", "#4C72B0", "o-"),
          ("r0_0.020", "prior $\\sigma_{\\rm zp}=0.02$ mag", "#55A868", "s-"),
          ("r0_0.001", "prior $\\sigma_{\\rm zp}=0.001$ mag", "#C44E52", "^-")]
for col, name, c, mk in series:
    med, lo, hi, ns = [], [], [], []
    for a, b in BINS:
        s = d[(d.f0L >= a) & (d.f0L < b)][col].dropna()
        ns.append(len(s))
        med.append(s.median()); lo.append(s.quantile(.1)); hi.append(s.quantile(.9))
    ax.fill_between(x, lo, hi, color=c, alpha=0.09, lw=0)
    ax.plot(x, med, mk, ms=5, color=c, label=name)
    print(name, ["%.3f" % m for m in med], "n", ns)

# the 0/0 regime: shade and mark it, do not let it look like a recovery point
ax.axvspan(-0.5, 0.5, color="0.7", alpha=0.28, lw=0)
ax.text(0.0, 0.415, "$0/0$ regime", ha="center", va="bottom", fontsize=7.0,
        color="0.25", rotation=90)

ax.set_xticks(x)
ax.set_xticklabels(lab, fontsize=8)
ax.set_xlim(-0.5, len(BINS) - 0.5)
ax.set_xlabel(r"$f_0L$")
ax.set_ylabel(r"$r_0^*$")
ax.set_ylim(0.4, 1.06)
ax.grid(alpha=0.25)
ax.legend(loc="lower right", framealpha=0.92, bbox_to_anchor=(1.0, 0.02))
ax.set_title("Absolute calibration helps only where coverage is short", fontsize=9)
fig.tight_layout()
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
