# -*- coding: utf-8 -*-
"""Audit: does the 7.2-hour nightly grouping merge independent photometric
systems inside one night?

For every campaign of the nested-test sample, each night (maximal run of
measurements separated by <7.2 h) is checked for the number of distinct

    lightcurve blocks (LCBLOCKID / block order),
    FILTER, MAGBAND (photometric passband),
    FACILITY, OBSERVERS,
    STANDARD (INTERNAL vs TRANSFORMED)

it contains. Campaigns are then split into "every night is single-system" and
"at least one night mixes systems", and the two groups are compared on the
residual scatter of Table 8 and on the four-condition pass rate.

Output: out_block_audit_20261121.csv, one row per campaign.
"""
import csv
import os
import sys
import time

import numpy as np
from scipy.stats import f as fdist

H = os.environ.get("AC_ROOT", os.path.dirname(os.path.abspath(__file__)))
ROOT = os.environ.get("ALCDEF_ROOT", os.path.join(H, "alcdef_data"))
OUT = os.path.join(H, "out_block_audit_20261121.csv")
KEYS = os.path.join(H, "out_lono_v5_c34.csv")


def parse_blocks(path):
    """Return jd, and per-point block id, filter, facility, observer, standard."""
    jd, blk, filt, fac, obs, std = [], [], [], [], [], []
    meta = {}
    bid = -1
    inmeta = False
    for line in open(path, "r", errors="ignore"):
        s = line.strip()
        if s == "STARTMETADATA":
            inmeta = True
            meta = {}
            bid += 1
            continue
        if s == "ENDMETADATA":
            inmeta = False
            continue
        if inmeta and "=" in s:
            k, v = s.split("=", 1)
            meta[k.strip()] = v.strip()
            continue
        if s.startswith("DATA="):
            p = s[5:].split("|")
            if len(p) < 2:
                continue
            try:
                j = float(p[0]); m = float(p[1])
            except ValueError:
                continue
            if not (j > 2.0e6 and -5 < m < 30):
                continue
            jd.append(j)
            blk.append(bid)
            filt.append((meta.get("FILTER", "?") + "/" + meta.get("MAGBAND", "?")))
            fac.append(meta.get("FACILITY", "?"))
            obs.append(meta.get("OBSERVERS", "?"))
            std.append(meta.get("STANDARD", "?"))
    return (np.array(jd), np.array(blk), np.array(filt), np.array(fac),
            np.array(obs), np.array(std))


def nights_of(t, gap=7.2):
    o = np.argsort(t)
    s = t[o]
    cut = np.where(np.diff(s) > gap)[0] + 1
    return [g for g in np.split(o, cut) if len(g) >= 3]


rows = [r for r in csv.DictReader(open(os.path.join(H, "out_full_lcdb_20261119.csv"),
                                      encoding="utf-8"))
        if int(float(r["nights"])) >= 5]
feat = {(r["file"], r["app"]): r for r in rows}
keys = [(r["file"], r["app"]) for r in csv.DictReader(open(KEYS, encoding="utf-8"))]
print("campaigns:", len(keys), flush=True)

t0 = time.time()
n = 0
with open(OUT, "w", newline="", encoding="utf-8") as fo:
    w = csv.writer(fo)
    w.writerow(["file", "app", "n_nights", "n_mixed_any", "n_mixed_block",
                "n_mixed_filter", "n_mixed_facility", "n_mixed_observer",
                "n_mixed_standard", "frac_mixed"])
    for base, appk in keys:
        p = os.path.join(ROOT, base)
        if not os.path.exists(p):
            continue
        jd, blk, filt, fac, obs, std = parse_blocks(p)
        if jd.size < 30:
            continue
        t = (jd - jd.min()) * 24.0
        o = np.argsort(t)
        s = t[o]
        cut = np.where(np.diff(s) / 24.0 > 90.0)[0] + 1
        parts = np.split(o, cut)
        k = int(appk)
        if k >= len(parts):
            continue
        idx = parts[k]
        ts = t[idx]
        g = nights_of(ts)
        if len(g) < 5:
            continue
        mn = mb = mf = mfa = mo = ms = 0
        for night in g:
            sel = idx[night]
            if len(set(blk[sel].tolist())) > 1:
                mb += 1; mn += 1
            if len(set(filt[sel].tolist())) > 1:
                mf += 1
            if len(set(fac[sel].tolist())) > 1:
                mfa += 1
            if len(set(obs[sel].tolist())) > 1:
                mo += 1
            if len(set(std[sel].tolist())) > 1:
                ms += 1
        any_mixed = sum(1 for night in g
                        if len(set(blk[idx[night]].tolist())) > 1
                        or len(set(filt[idx[night]].tolist())) > 1
                        or len(set(fac[idx[night]].tolist())) > 1
                        or len(set(std[idx[night]].tolist())) > 1)
        w.writerow([base, int(appk), len(g), any_mixed, mb, mf, mfa, mo, ms,
                    round(any_mixed / len(g), 4)])
        n += 1
        if n % 1000 == 0:
            print("  %d  %.1f min" % (n, (time.time() - t0) / 60), flush=True)

print("evaluated %d  (%.1f min) -> %s" % (n, (time.time() - t0) / 60, OUT), flush=True)
