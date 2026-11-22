# -*- coding: utf-8 -*-
"""v5 §4.2: four-arm calibration with the unified 4 -> 8 screen and the
absorption control 8 -> 12, on real cadences, with FIXED physical amplitudes.

Injection (identical physical ruler everywhere in v5):
    A1 = 0.05 or 0.10 mag,  A2 = 0.25,  A3 = 0.03,  A4 = 0.08
    true-P  arms: harmonics at  f, 2f, 3f, 4f
    true-2P arms: harmonics at  f/2, f, 3f/2, 2f
    noise = the campaign's own measured residual scatter
    arms B, C, D add apparition drift (amplitude modulation + phase drift)

Screens evaluated on the same simulated lightcurve (paired):
    S1 unabsorbed : H0={f,2f} (4)          vs H1=+{f/2,3f/2} (8)      dof=4
    S2 absorbed   : H0'={f,2f}+gamma12 (8) vs H1'=+{f/2,3f/2} (12)   dof=4
    controls for both: quartic {u..u^4} (4 extra), gamma12 (4 extra)
    accept = F_in > F_crit(4,dof)  and  held-out gain > 0  and  beats both

Sample: 200 cadences, one per asteroid, from the K>=5 physical selection.
Usage: python calib_v5_20261121.py [ncamp] [nrep] [A1] [suffix]
"""
import csv, os, sys, time
import numpy as np
from scipy.stats import f as fdist

H = r"C:\Users\jacky\.openclaw\workspace\ac_repo_plan"
ROOT = r"C:\Users\jacky\.openclaw\workspace\alcdef_data"
NCAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 200
NREP = int(sys.argv[2]) if len(sys.argv) > 2 else 8
A1MAG = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
SUF = sys.argv[4] if len(sys.argv) > 4 else "_a05"
NPZ = os.path.join(H, "calib_v8k%s.npz" % SUF)
SD_AMP, SD_PH = 0.15, 0.25


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
    if kind == "H0p":
        return np.column_stack([c1, s1, c2, s2, u * c1, u * s1, u * c2, u * s2])
    if kind == "H1p":
        return np.column_stack([c1, s1, c2, s2, u * c1, u * s1, u * c2, u * s2,
                                np.cos(w / 2), np.sin(w / 2),
                                np.cos(1.5 * w), np.sin(1.5 * w)])
    if kind == "C34":
        return np.column_stack([c1, s1, c2, s2,
                                np.cos(3 * w), np.sin(3 * w),
                                np.cos(4 * w), np.sin(4 * w)])
    if kind == "Cq":
        return np.column_stack([c1, s1, c2, s2, u, u ** 2, u ** 3, u ** 4])
    if kind == "Cg":
        return np.column_stack([c1, s1, c2, s2, u * c1, u * s1, u * c2, u * s2])
    raise ValueError(kind)


def fit_on(t, y, g, f, kind, u, idx):
    X = design(t, f, kind, u)[idx]
    yv = y[idx].copy()
    for gg in g:
        m = np.isin(idx, gg)
        if m.any():
            X[m] = X[m] - X[m].mean(axis=0)
            yv[m] = yv[m] - yv[m].mean()
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return None
    b, *_ = np.linalg.lstsq(X, yv, rcond=None)
    r = yv - X @ b
    return float(r @ r), b


def heldout(t, y, g, f, kind, u, b, k, idx):
    X = design(t, f, kind, u)
    m = np.isin(idx, g[k])
    Xh = X[m] - X[m].mean(axis=0)
    yh = y[m] - y[m].mean()
    r = yh - Xh @ b
    return float(r @ r)


def screen(t, y, g, f, kinds, u):
    """kinds = (base, alt) pairs; controls are Cq and Cg (4 extra parameters)."""
    allidx = np.arange(t.size)
    out = {}
    for kind in (kinds[0], kinds[1], "Cq", "Cg", "C34"):
        r = fit_on(t, y, g, f, kind, u, allidx)
        if r is None:
            return None
        oos = 0.0
        for k in range(len(g)):
            tr = allidx[~np.isin(allidx, g[k])]
            rr = fit_on(t, y, g, f, kind, u, tr)
            if rr is None:
                return None
            oos += heldout(t, y, g, f, kind, u, rr[1], k, allidx)
        out[kind] = (r[0], oos)
    return out


def verdict(out, dof_extra=4):
    dof = out.get("dof", None)
    pass


# ---- sample ----
feat = {}
for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"), encoding="utf-8")):
    if int(float(r["nights"])) < 5:
        continue
    feat.setdefault(int(r["num"]), []).append(r)
byast = list(feat.values())
rng_pick = np.random.default_rng(20261120)
order = np.arange(len(byast)); rng_pick.shuffle(order)

cad = []
for i in order:
    if len(cad) >= NCAMP:
        break
    r = sorted(byast[i], key=lambda z: (z["file"], int(z["app"])))[0]
    fp = os.path.join(ROOT, r["file"])
    if not os.path.exists(fp):
        continue
    jd, mag = parse(fp)
    if jd.size < 30:
        continue
    t_all = (jd - jd.min()) * 24.0
    o = np.argsort(t_all); s = t_all[o]
    cut = np.where(np.diff(s) / 24.0 > 90.0)[0] + 1
    parts = np.split(o, cut)
    k = int(r["app"])
    if k >= len(parts):
        continue
    ts = t_all[parts[k]]
    g = nights_of(ts)
    if len(g) < 5:
        continue
    cad.append(dict(t=ts, g=g, P=float(r["P"]), sig=max(float(r["sig"]), 0.02),
                    num=int(r["num"]), file=r["file"]))
print("cadences:", len(cad), " A1 = %.2f mag, reps %d" % (A1MAG, NREP), flush=True)

ARMS = ("A", "B", "C", "D", "E", "F", "G")
acc1 = {a: np.zeros(len(cad)) for a in ARMS}
rep1 = {a: np.zeros(len(cad)) for a in ARMS}
acc2 = {a: np.zeros(len(cad)) for a in ARMS}
ac4 = {a: np.zeros(len(cad)) for a in ARMS}
rep2 = {a: np.zeros(len(cad)) for a in ARMS}

SHAPE = [A1MAG, 0.25 * A1MAG / 0.05, 0.03 * A1MAG / 0.05, 0.08 * A1MAG / 0.05]
# arm F: a strict H0 waveform, A3 = A4 = 0 (only f and 2f present)
SHAPE_PURE = [A1MAG, 0.25 * A1MAG / 0.05, 0.0, 0.0]
rng = np.random.default_rng(7)
t0 = time.time()
for i, c in enumerate(cad):
    t = c["t"]; g = c["g"]; P = c["P"]; sig = c["sig"]
    T = t.max() - t.min(); u = (t - t.mean()) / max(T, 1e-9)
    fcat = 1.0 / P
    for _ in range(NREP):
        draw = [(rng.uniform(0, 2 * np.pi), rng.normal(0, SD_AMP), rng.normal(0, SD_PH))
                for _m in range(4)]
        noise = rng.normal(0, sig, t.size)
        for arm in ARMS:
            twoP = arm in ("A", "D")
            drift = arm in ("B", "C", "D", "G")
            shape = SHAPE_PURE if arm in ("F", "G") else SHAPE
            ftrue = fcat / 2.0 if twoP else fcat
            y = np.zeros_like(t)
            for m in range(4):
                ph, gm_, cm_ = draw[m]
                gm = gm_ if drift else 0.0
                cm = cm_ if drift else 0.0
                y += shape[m] * (1 + gm * u) * np.cos(2 * np.pi * (m + 1) * ftrue * t + ph + cm * u)
            y += noise
            y += rng.normal(0, sig, t.size)
            for tag, kinds, ac, rp in (("S1", ("H0", "H1"), acc1, rep1),
                                       ("S2", ("H0p", "H1p"), acc2, rep2)):
                out = screen(t, y, g, fcat, kinds, u)
                if out is None:
                    continue
                base, alt = kinds
                dof = max(t.size - 12 - len(g), 1)
                Fcrit = fdist.ppf(0.95, 4, dof)
                F = max((out[base][0] - out[alt][0]) / 4.0, 0.0) / max(out[alt][0] / dof, 1e-12)
                dh = out[base][1] - out[alt][1]
                dq = out[base][1] - out["Cq"][1]
                dg = out[base][1] - out["Cg"][1]
                d34 = out[base][1] - out["C34"][1]
                ok = (F > Fcrit) and (dh > 0) and (dh > dq) and (dh > dg)
                ok4 = ok and (dh > d34)
                rp[arm][i] += 1
                ac[arm][i] += int(ok)
                if tag == "S1":
                    ac4[arm][i] += int(ok4)
    if (i + 1) % 25 == 0:
        print("  %d/%d  %.1f min" % (i + 1, len(cad), (time.time() - t0) / 60), flush=True)

np.savez(NPZ,
         **{("acc1_" + a): acc1[a] for a in ARMS}, **{("rep1_" + a): rep1[a] for a in ARMS},
         **{("acc2_" + a): acc2[a] for a in ARMS}, **{("rep2_" + a): rep2[a] for a in ARMS},
         **{("acc4_" + a): ac4[a] for a in ARMS},
         nums=np.array([c["num"] for c in cad]), A1=A1MAG)

lab = dict(A="A true 2P, stationary      ", B="B true P, amplitude modul. ",
           C="C true P, four-harmonic drift", D="D true 2P + drift         ",
           E="E true P, stationary (null) ",
           F="F strict H0, A3=A4=0     ", G="G strict H0 + drift       ")
print("\nA1 = %.2f mag   (exact F_crit(4, dof) per campaign)" % A1MAG)
for tag, ac, rp in (("unabsorbed 4->8 ", acc1, rep1), ("absorbed 8->12  ", acc2, rep2),
                    ("unabs + beats 3f4f", ac4, rep1)):
    print(" screen %s" % tag)
    for a in ARMS:
        print("   %s  %5.1f%%  (%d/%d)" % (lab[a], 100 * ac[a].sum() / max(rp[a].sum(), 1),
                                           int(ac[a].sum()), int(rp[a].sum())))
print("->", NPZ, " %.1f min" % ((time.time() - t0) / 60))
