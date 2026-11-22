# -*- coding: utf-8 -*-
"""v6 Fig. 2: a priori audit with the nightly-zero-point cost separated out.

(a) fraction of cadences whose minimum identifiable amplitude exceeds a fixed
    reference amplitude, against coverage.  Solid: the paper's design (K free
    nightly offsets).  Dashed: the same cadences with a single global offset,
    i.e. the pure window/photometric baseline.  The gap between the two green
    curves is what the nightly zero points cost.
(b) r0* median with its 10-90 band, both designs.
lambda_crit is the exact per-campaign value (noncentral F, 80% power).
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(H, "fig2.png")

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "legend.fontsize": 6.9, "xtick.labelsize": 8, "ytick.labelsize": 8})

d = pd.read_csv(H + r"\out_zp_baseline_20261121.csv")
d["app"] = d["app"].astype(str)
meta = pd.read_csv(H + r"\out_full_lcdb_20261119.csv")
meta["app"] = meta["app"].astype(str)
meta["U"] = meta["U"].fillna("").astype(str).str.strip()
d = d.merge(meta[["file", "app", "U", "span_d"]], on=["file", "app"], how="left",
            suffixes=("", "_m"))
sp = d["span_d_m"] if "span_d_m" in d.columns else d["span_d"]
d["clean"] = d.U.isin(["3", "3-"]) & (sp.astype(float) <= 30.0)
print("rows", len(d), " clean subset", int(d.clean.sum()))

BINS = [(0, .05), (.05, .10), (.10, .15), (.15, .25), (.25, .35), (.35, .50), (.50, 1.0), (1.0, 1e9)]
LAB = [r"$<0.05$", r"$0.05$--$0.10$", r"$0.10$--$0.15$", r"$0.15$--$0.25$",
       r"$0.25$--$0.35$", r"$0.35$--$0.50$", r"$0.50$--$1.0$", r"$>1.0$"]
x = np.arange(len(BINS))

fracn = {0.02: [], 0.05: [], 0.10: []}
fracg = {0.05: []}
fracc = []
medn, lon, hin, medg, log, hig, ns = [], [], [], [], [], [], []
for a, b in BINS:
    m = (d.f0L >= a) & (d.f0L < b)
    ns.append(int(m.sum()))
    for ref in fracn:
        v = d.loc[m, "An_%.2f" % ref]
        fracn[ref].append(100 * np.mean(v.dropna() > ref))
    v = d.loc[m, "Ag_0.05"]
    fracg[0.05].append(100 * np.mean(v.dropna() > 0.05))
    vc = d.loc[m & d.clean, "An_0.05"]
    fracc.append(100 * np.mean(vc.dropna() > 0.05) if len(vc.dropna()) else np.nan)
    print("    clean n in bin: %d" % len(vc.dropna()))
    rn = d.loc[m, "r0n"].dropna()
    rg = d.loc[m, "r0g"].dropna()
    medn.append(rn.median()); lon.append(rn.quantile(.1)); hin.append(rn.quantile(.9))
    medg.append(rg.median()); log.append(rg.quantile(.1)); hig.append(rg.quantile(.9))
    print("  %.2f-%.2f n=%5d night>0.05 %5.1f global %5.1f | r0n %.2f r0g %.2f"
          % (a, b, ns[-1], fracn[0.05][-1], fracg[0.05][-1], medn[-1], medg[-1]))

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.5))
ax[0].fill_between(x, fracg[0.05], fracn[0.05], color="#55A868", alpha=0.13, lw=0,
                   label="cost of nightly zero points")
for ref, c, mk in ((0.02, "#4C72B0", "o"), (0.05, "#55A868", "s"), (0.10, "#C44E52", "^")):
    ax[0].plot(x, fracn[ref], mk + "-", ms=4.5, color=c,
               label=r"$A_{\rm ref}=%.2f$, nightly $c_k$" % ref)
ax[0].plot(x, fracg[0.05], "s--", ms=4.5, color="#55A868", mfc="white",
           label=r"$A_{\rm ref}=0.05$, one global offset")
ax[0].plot(x, fracc, "D-.", ms=4.0, color="0.15", lw=1.4,
           label=r"$A_{\rm ref}=0.05$, clean subset ($U\geq3$, $\leq30$ d)")
ax[0].set_xticks(x); ax[0].set_xticklabels(LAB, rotation=45, ha="right")
ax[0].set_xlabel(r"$f_0L$")
ax[0].set_ylabel("cadences below the floor [%]")
ax[0].set_ylim(0, 112)
ax[0].grid(alpha=0.25)
h_, l_ = ax[0].get_legend_handles_labels()
fig.legend(h_, l_, loc="lower center", ncol=3, fontsize=6.2, frameon=False,
           bbox_to_anchor=(0.5, 0.0), handlelength=1.8)
ax[0].annotate("", xy=(0.35, fracn[0.05][0]), xytext=(0.35, fracg[0.05][0]),
               arrowprops=dict(arrowstyle="<->", lw=0.9, color="#333333"))
ax[0].text(0.35, fracn[0.05][0] + 5.5,
           "+%.1f pp" % (fracn[0.05][0] - fracg[0.05][0]), fontsize=7.5,
           ha="center", va="bottom")
ax[0].set_title("(a) information collapse, 12,070 cadences")

ax[1].fill_between(x, lon, hin, color="#C44E52", alpha=0.18, lw=0, label="10--90 pct")
ax[1].plot(x, medn, "o-", ms=4.5, color="#C44E52", label=r"nightly $c_k$")
ax[1].fill_between(x, log, hig, color="#4C72B0", alpha=0.15, lw=0)
ax[1].plot(x, medg, "s--", ms=4, color="#4C72B0", mfc="white", label="one global offset")
ax[1].axhline(1.0, color="k", lw=0.7, ls=":")
ax[1].set_xticks(x); ax[1].set_xticklabels(LAB, rotation=45, ha="right")
ax[1].set_xlabel(r"$f_0L$")
ax[1].set_ylabel(r"$r_0^*$")
ax[1].set_ylim(0, 1.25)
ax[1].grid(alpha=0.25)
ax[1].legend(loc="lower right", framealpha=0.92)
ax[1].set_title(r"(b) flip threshold $r_0^*$: window vs zero points")

fig.tight_layout(rect=(0, 0.13, 1, 1))
fig.savefig(OUT, dpi=220)
print("wrote", OUT)
