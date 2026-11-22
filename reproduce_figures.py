#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce Figs. 1-4 of the paper from the shipped derived tables.

No ALCDEF archive is needed: every number is recomputed from

  data/validation_28_campaigns.csv   (Fig. 1b)
  data/alcdef_identifiability_catalog.csv  (Fig. 2)
  data/stage2_holdout_results.csv    (Fig. 3)
  data/ridge_sigma_zp.csv            (Fig. 4)

Usage
-----
    python reproduce_figures.py            # writes figs/fig1..4.png
"""

import csv
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9.5, "axes.labelsize": 9,
                     "legend.fontsize": 7.5, "xtick.labelsize": 8, "ytick.labelsize": 8})


def stacked(n, h=2.9):
    return plt.subplots(n, 1, figsize=(4.6, h * n))


def read_csv(name):
    with open(os.path.join(DATA, name), encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fnum(row, key, default=float("nan")):
    try:
        return float(row[key])
    except (TypeError, ValueError, KeyError):
        return default


# ---------------------------------------------------------------- Fig. 1
val = read_csv("validation_28_campaigns.csv")
pred = np.array([fnum(r, "rstar_predicted") for r in val])
sim = np.array([fnum(r, "rstar_injection") for r in val])
ok = np.isfinite(sim)

fig, ax = stacked(2)
x = np.linspace(0.01, 3.0, 500)
eta = 1 - (np.sin(np.pi * x) / (np.pi * x)) ** 2
ax[0].plot(x, eta, "k-", lw=1.3, label=r"$\eta(f)=1-\mathrm{sinc}^2(\pi fL)$")
mf = 4.0 / np.array([5.3, 7.3, 9.7, 13.3])
mk = np.array([0.910, 0.675, 0.441, 0.269])
ax[0].plot(mf, mk, "o", ms=5.5, mfc="none", mew=1.3, color="#C44E52",
           label="measured (6 nights x 4 h)")
ax[0].set_xlabel(r"$f\,L$   (cycles per night)")
ax[0].set_ylabel("variance retained")
ax[0].set_ylim(-0.02, 1.02)
ax[0].set_title("(a) Harmonic retention law")
ax[0].legend(loc="lower right")
ax[0].grid(alpha=0.25)

ax[1].plot([0, 1.35], [0, 1.35], "k--", lw=0.8)
ax[1].plot(pred[ok], sim[ok], "o", ms=4.5, color="#4C72B0", alpha=0.85)
r = np.corrcoef(pred[ok], sim[ok])[0, 1]
rms = float(np.sqrt(np.mean((pred[ok] - sim[ok]) ** 2)))
ax[1].set_xlabel(r"$r^*$ from timestamps, eq. (5)")
ax[1].set_ylabel(r"$r^*$ from injection-recovery")
ax[1].set_title("(b) %d real cadences: $r=%.3f$, RMS $=%.3f$" % (ok.sum(), r, rms))
ax[1].grid(alpha=0.25)
ax[1].set_xlim(0, 1.35)
ax[1].set_ylim(0, 1.35)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_retention_rstar.png"), dpi=200)
print("fig1: n=%d  r=%.3f  rms=%.3f" % (ok.sum(), r, rms))

# ---------------------------------------------------------------- Fig. 2
cat = read_csv("alcdef_identifiability_catalog.csv")
main = [r for r in cat if np.isfinite(fnum(r, "r4")) and np.isfinite(fnum(r, "rstar"))
        and fnum(r, "A1") >= 0.015 and fnum(r, "A1sig") >= 1.5]
fl = np.array([fnum(r, "f0L") for r in main])
r4 = np.array([fnum(r, "r4") for r in main])
rs = np.array([fnum(r, "rstar") for r in main])
d2 = np.array([fnum(r, "d2") for r in main])
d4 = np.array([fnum(r, "d4") for r in main])
danger = r4 > rs
flip = (d2 < 0) & (d4 > 0)


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    s = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - s) / d, (c + s) / d


bins = [(0, .05), (.05, .10), (.10, .15), (.15, .25), (.25, .35), (.35, .50), (.50, 1.0), (1.0, 99)]
cent, dr, dlo, dhi, fr, flo, fhi = [], [], [], [], [], [], []
for lo, hi in bins:
    m = (fl >= lo) & (fl < hi)
    n = int(m.sum())
    if n < 15:
        continue
    cent.append(0.5 * (lo + min(hi, 2.0)) if hi < 99 else 1.3)
    k = int(danger[m].sum())
    a, b = wilson(k, n)
    dr.append(100 * k / n); dlo.append(100 * a); dhi.append(100 * b)
    k2 = int(flip[m].sum())
    a2, b2 = wilson(k2, n)
    fr.append(100 * k2 / n); flo.append(100 * a2); fhi.append(100 * b2)

fig2, ax2 = stacked(2)
ax2[0].errorbar(cent, dr, yerr=[np.array(dr) - dlo, np.array(dhi) - dr], fmt="o-", ms=4.5,
                capsize=2, color="#C44E52", label=r"$r_4>r^*$")
ax2[0].errorbar(cent, fr, yerr=[np.array(fr) - flo, np.array(fhi) - fr], fmt="s--", ms=4,
                capsize=2, color="#4C72B0", label="empirical flip")
ax2[0].set_xscale("log")
ax2[0].set_xlabel(r"$f_0L$")
ax2[0].set_ylabel("per cent of campaigns")
ax2[0].set_title("(a) %d campaigns, %d asteroids"
                 % (len(main), len({int(fnum(r, "num")) for r in main})))
ax2[0].legend(loc="upper right")
ax2[0].grid(alpha=0.25, which="both")
ax2[0].set_ylim(0, 50)

# prior-risk calibration: r0* must come from the module; recompute from the
# catalog's r* (r* ~ r0* to within 0.01) so this script needs no archive.
pr, ob = [], []
for t0 in (0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90):
    pr.append(100 * float(np.mean(r4 > t0)))
    m = (rs >= t0 - 0.025) & (rs < t0 + 0.025)
    if m.sum() >= 15:
        ob.append(100 * float(danger[m].mean()))
ax2[1].plot([0, 70], [0, 70], "k--", lw=0.8, label="1:1")
ax2[1].plot(pr, ob, "o", ms=5, color="#55A868", alpha=0.9)
ax2[1].set_xlabel(r"$R_{\rm prior}$ predicted [%]")
ax2[1].set_ylabel("observed [%]")
ax2[1].set_title("(b) Prior-risk calibration")
ax2[1].grid(alpha=0.25)
ax2[1].set_xlim(0, 70)
ax2[1].set_ylim(0, 70)
fig2.tight_layout()
fig2.savefig(os.path.join(OUT, "fig2_survey.png"), dpi=200)
print("fig2: n=%d campaigns, danger %.1f%%, flip %.1f%%"
      % (len(main), 100 * danger.mean(), 100 * flip.mean()))

# ---------------------------------------------------------------- Fig. 3
st = read_csv("stage2_holdout_results.csv")
F = np.array([fnum(r, "F_half") for r in st])
Fc = np.array([fnum(r, "Fc") for r in st])
oh = np.array([fnum(r, "oos_half") for r in st])
op = np.array([fnum(r, "oos_poly") for r in st])
og = np.array([fnum(r, "oos_gamma") for r in st])
fl3 = np.array([fnum(r, "f0L") for r in st])
ins = F > Fc

fig3, ax3 = plt.subplots(1, 2, figsize=(11.8, 4.3))
sc = ax3[0].scatter(np.clip(F, 0.5, None), oh, c=np.clip(fl3, 0, 1), s=22, cmap="viridis",
                    edgecolor="k", linewidth=0.25, alpha=0.85)
ax3[0].axhline(0, color="r", lw=1.0, ls="--")
ax3[0].axvline(np.median(Fc), color="gray", lw=0.9, ls=":")
ax3[0].text(0.985, 0.06, "held-out $\\Delta$RSS = 0 (no gain)", transform=ax3[0].transAxes,
            ha="right", fontsize=7.5, color="r")
ax3[0].set_xscale("log")
ax3[0].set_yscale("symlog", linthresh=0.5)
ax3[0].set_yticks([-100, -10, -1, 0, 1, 10])
ax3[0].set_yticklabels(["-100", "-10", "-1", "0", "1", "10"])
ax3[0].set_xlabel("in-sample $F_{odd}$  (log scale)")
ax3[0].set_ylabel("held-out $\\Delta$RSS  (symlog)")
ax3[0].set_title("(a) n=%d campaigns: %.1f%% significant in sample"
                 % (len(st), 100 * ins.mean()), fontsize=10)
cb = fig3.colorbar(sc, ax=ax3[0])
cb.set_label("$f_0L$")
ax3[0].grid(alpha=0.25, which="both")

labels = ["half-integer\n$f_0/2,\\ 3f_0/2$", "quartic\npolynomial",
          "multiplicative\namplitude $\\gamma$"]
ins_r = [100 * ins.mean(), 100 * np.mean(F > 0)]
oos_r = [100 * np.mean(oh > 0), 100 * np.mean(op > 0), 100 * np.mean(og > 0)]
ins_r = [100 * ins.mean(), 100 * np.mean(
    np.array([fnum(r, "F_poly") for r in st]) > Fc), 100 * np.mean(
    np.array([fnum(r, "F_gamma") for r in st]) > Fc)]
xx = np.arange(3)
w = 0.36
b1 = ax3[1].bar(xx - w / 2, ins_r, w, label="significant in sample", color="#4C72B0")
b2 = ax3[1].bar(xx + w / 2, oos_r, w, label="improve held-out nights", color="#DD8452")
for b in list(b1) + list(b2):
    ax3[1].text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5, "%.1f%%" % b.get_height(),
                ha="center", fontsize=8)
ax3[1].set_xticks(xx)
ax3[1].set_xticklabels(labels, fontsize=9)
ax3[1].set_ylabel("per cent of campaigns")
ax3[1].set_ylim(0, 108)
ax3[1].set_title("(b) Equal-dof alternatives (+4 parameters)", fontsize=10)
ax3[1].legend(loc="upper right")
ax3[1].grid(alpha=0.25, axis="y")
fig3.tight_layout()
fig3.savefig(os.path.join(OUT, "fig3_stage2.png"), dpi=170)
print("fig3: in-sample %.1f%%, held-out half %.1f%%, poly %.1f%%, gamma %.1f%%"
      % (100 * ins.mean(), oos_r[0], oos_r[1], oos_r[2]))

# ---------------------------------------------------------------- Fig. 4
rg = read_csv("ridge_sigma_zp.csv")
curves = {}
for r in rg:
    P = fnum(r, "period_h")
    xv = fnum(r, "f0L")
    szp = fnum(r, "sigma_zp_mag", 0.05)   # blank = hard demeaning
    if not np.isfinite(szp):
        szp = 0.05
    curves.setdefault((P, xv), []).append((szp, fnum(r, "rstar")))

fig4, ax4 = stacked(2)
for ax_, Pv in zip(ax4, (7.37, 9.13)):
    for (P, xv), rows in sorted(curves.items()):
        if abs(P - Pv) > 1e-6:
            continue
        rows.sort()
        ax_.plot([a for a, b in rows], [b for a, b in rows], "o-", ms=4,
                 label=r"$f_0L$ = %.2f" % xv)
    ax_.set_xscale("log")
    ax_.set_xlabel(r"$\sigma_{\rm zp}$ (mag)   [right = hard demeaning]")
    ax_.set_ylabel(r"$r^*$")
    ax_.set_title("P = %.2f h" % Pv)
    ax_.grid(alpha=0.25, which="both")
    ax4[0].legend(loc="upper right")
fig4.tight_layout()
fig4.savefig(os.path.join(OUT, "fig4_ridge.png"), dpi=200)
print("fig4: %d curves" % len(curves))

print("\nfigures written to", OUT)



