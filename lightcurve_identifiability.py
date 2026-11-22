# -*- coding: utf-8 -*-
"""
lightcurve_identifiability
==========================

A priori identifiability diagnostics for ground-based asteroid
lightcurves whose photometric zero point is unknown on each night.

The module answers three questions before any period search is run:

1. r0_star  -- the geometric flip threshold on A4/A1 at a candidate period
               (computable from timestamps alone);
2. R_prior  -- the population probability that a asteroid with this cadence
               and this A4/A1 exceeds that threshold (P -> P/2 flip risk);
3. A1_min   -- the smallest fundamental amplitude for which an M = 4
               odd-harmonic test reaches 80% power at alpha = 0.05.

It also provides a Stage-2 disambiguation test that compares the
catalogued period against twice that period, together with the two
controls that the paper shows are necessary on real data:

* a flexibility-matched alternative (quartic polynomial), and
* a multiplicative amplitude slope gamma, which absorbs the seasonal
  amplitude evolution driven by the solar phase angle,
* plus a held-out-night cross-validation which is the only acceptance
  criterion that is not satisfied by generic flexibility.

References
----------
Cheung (2026), "Two identifiability regimes in ground-based asteroid
period determination", Astronomy and Computing (submitted).

Quickstart
----------
>>> from lightcurve_identifiability import blind_metrics, stage2_test
>>> m = blind_metrics(t, night_id, sigma=0.05, period_h=16.49)
>>> m['r0_star'], m['A1_min'], m['eta_f0']
>>> s = stage2_test(t, mag, night_id, period_h=16.49)
>>> s['F_odd'], s['heldout_drss_half'], s['verdict']
"""

from __future__ import annotations

import os
import math
from typing import Optional, Sequence

import numpy as np

__all__ = [
    "assign_nights",
    "gram_residuals",
    "blind_metrics",
    "population_flip_risk",
    "stage2_test",
    "minimum_identifiable_amplitude",
    "LAMBDA_CRIT",
]

#: Non-centrality required for 80% power at alpha = 0.05, dof = (4, 150).
LAMBDA_CRIT = 12.2


# ----------------------------------------------------------------------
# helper: design matrices
# ----------------------------------------------------------------------
def assign_nights(t: Sequence[float], gap_h: float = 7.2) -> np.ndarray:
    """Group observation times into nights.

    Parameters
    ----------
    t : array_like
        Times in hours, any zero point.
    gap_h : float
        A gap larger than this starts a new night.

    Returns
    -------
    night_id : ndarray of int
        Night index for each observation (0-based, in time order).
    """
    t = np.asarray(t, dtype=float)
    order = np.argsort(t)
    ident = np.empty(t.size, dtype=int)
    nid = 0
    prev = None
    for i in order:
        if prev is not None and (t[i] - prev) > gap_h:
            nid += 1
        ident[i] = nid
        prev = t[i]
    return ident


def _harmonic_design(t: np.ndarray, f0: float, max_order: int) -> np.ndarray:
    cols = []
    for m in range(1, max_order + 1):
        w = 2.0 * np.pi * m * f0 * t
        cols.append(np.cos(w))
        cols.append(np.sin(w))
    return np.column_stack(cols)


def _demean_within(X: np.ndarray, night_id: np.ndarray) -> np.ndarray:
    """Subtract the per-night mean of each column (nightly zero points)."""
    Xc = X.copy()
    for k in np.unique(night_id):
        m = night_id == k
        Xc[m] -= Xc[m].mean(axis=0)
    return Xc


# ----------------------------------------------------------------------
# 1. Gram-matrix diagnostics
# ----------------------------------------------------------------------
def gram_residuals(t: Sequence[float], night_id: Sequence[int],
                   period_h: float, orders: Sequence[int] = (1, 2, 3, 4)) -> dict:
    """Schur-complement residual fractions of the demeaned harmonic basis.

    Implements rho_m(S) = 1/2 tr[G_mm - G_mS G_SS^-1 G_Sm] (Eq. A.2 of the
    paper), where G is the Gram matrix of the *demeaned* harmonic basis.

    Returns
    -------
    dict with keys a1, a4, d3, r0_star, r_star, eta_f0, n_points, n_nights
    """
    t = np.asarray(t, dtype=float)
    night_id = np.asarray(night_id)
    f0 = 1.0 / float(period_h)
    n = t.size

    def block(m: int) -> np.ndarray:
        c = np.sqrt(2.0) * np.cos(2.0 * np.pi * m * f0 * t)
        s = np.sqrt(2.0) * np.sin(2.0 * np.pi * m * f0 * t)
        c = _demean_within(c[:, None], night_id)[:, 0]
        s = _demean_within(s[:, None], night_id)[:, 0]
        return np.array([[c @ c, c @ s], [s @ c, s @ s]]) / n

    G = {(a, b): block(a) if a == b else None for a in orders for b in orders}
    for a in orders:
        for b in orders:
            if G[(a, b)] is None:
                c1 = _demean_within((np.sqrt(2) * np.cos(2 * np.pi * a * f0 * t))[:, None], night_id)[:, 0]
                s1 = _demean_within((np.sqrt(2) * np.sin(2 * np.pi * a * f0 * t))[:, None], night_id)[:, 0]
                c2 = _demean_within((np.sqrt(2) * np.cos(2 * np.pi * b * f0 * t))[:, None], night_id)[:, 0]
                s2 = _demean_within((np.sqrt(2) * np.sin(2 * np.pi * b * f0 * t))[:, None], night_id)[:, 0]
                G[(a, b)] = np.array([[c1 @ c2, c1 @ s2], [s1 @ c2, s1 @ s2]]) / n

    def rho(m: int, model: Sequence[int]) -> float:
        p = 2 * len(model)
        GX = np.zeros((p, p))
        GXm = np.zeros((p, 2))
        for i, a in enumerate(model):
            for j, b in enumerate(model):
                GX[2 * i:2 * i + 2, 2 * j:2 * j + 2] = G[(a, b)]
            GXm[2 * i:2 * i + 2, :] = G[(a, m)]
        schur = G[(m, m)] - GXm.T @ np.linalg.pinv(GX) @ GXm
        return 0.5 * float(np.trace(schur))

    a1 = rho(1, (2, 4))
    a4 = rho(4, (1, 2))
    d3 = rho(3, (2, 4)) - rho(3, (1, 2))
    r0 = math.sqrt(a1 / a4) if a4 > 0 else float("nan")
    r3 = math.sqrt(max(a1 + d3 * 0.3 ** 2, 0.0) / a4) if a4 > 0 else float("nan")
    f0 = 1.0 / period_h
    x = math.pi * f0 * float(np.median([_night_length(t, night_id, k)
                                        for k in np.unique(night_id)]))
    eta = 1.0 - (math.sin(x) / x) ** 2 if x != 0 else 0.0
    return dict(a1=a1, a4=a4, d3=d3, r0_star=r0, r_star=r3, eta_f0=eta,
                n_points=int(n), n_nights=int(np.unique(night_id).size))


def _night_length(t: np.ndarray, night_id: np.ndarray, k: int) -> float:
    m = night_id == k
    return float(t[m].max() - t[m].min())


def minimum_identifiable_amplitude(r0_star: float, sigma: float, n_points: int,
                                   lam_crit: float = LAMBDA_CRIT) -> float:
    """A1_min = sigma * sqrt(2*lam_crit / (N * a1)), Eq. (12).

    ``r0_star`` enters through a1 = r0_star^2 * a4; since a4 ~ O(1) the
    approximation a1 ~ r0_star^2 / (1 + r0_star^2) is used only if the
    caller passes r0_star alone. Prefer :func:`blind_metrics`, which
    passes the exact a1.
    """
    a1 = r0_star ** 2 / (1.0 + r0_star ** 2)
    return sigma * math.sqrt(2.0 * lam_crit / max(n_points * a1, 1e-12))


def population_flip_risk(r0_star: float, population_r4: Sequence[float]) -> float:
    """R_prior = P_pop(A4/A1 > r0*), Eq. (7)."""
    if not np.isfinite(r0_star):
        return float("nan")
    arr = np.asarray(population_r4, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.mean(arr > r0_star))


def blind_metrics(t: Sequence[float], night_id: Sequence[int], sigma: float,
                  period_h: float, population_r4: Optional[Sequence[float]] = None,
                  amplitude_a1: Optional[float] = None) -> dict:
    """All pre-observing diagnostics for one target and one candidate period.

    Parameters
    ----------
    t, night_id : cadence (times in hours, night index per point)
    sigma : per-point photometric uncertainty (mag)
    period_h : candidate period (h)
    population_r4 : optional empirical A4/A1 sample for R_prior
    amplitude_a1 : optional measured fundamental amplitude (mag); if given,
                   the returned 'can_resolve' tells whether A1 > A1_min

    Returns
    -------
    dict with r0_star, r_star, eta_f0, a1, a4, d3, A1_min, R_prior,
    n_points, n_nights, f0L and (if amplitude given) can_resolve.
    """
    g = gram_residuals(t, night_id, period_h)
    n = g["n_points"]
    a1min = sigma * math.sqrt(2.0 * LAMBDA_CRIT / max(n * g["a1"], 1e-12))
    nights = [k for k in np.unique(np.asarray(night_id))]
    L = float(np.median([_night_length(np.asarray(t, float), np.asarray(night_id), k)
                         for k in nights]))
    out = dict(g)
    out["A1_min"] = a1min
    out["f0L"] = L / float(period_h)
    out["sigma"] = sigma
    if population_r4 is not None:
        out["R_prior"] = population_flip_risk(g["r0_star"], population_r4)
    if amplitude_a1 is not None:
        lam = (n / (2.0 * sigma ** 2)) * g["a1"] * amplitude_a1 ** 2
        out["lambda_odd"] = lam
        out["can_resolve"] = bool(amplitude_a1 >= a1min)
    return out


# ----------------------------------------------------------------------
# 2. Stage-2 disambiguation with the controls of the paper
# ----------------------------------------------------------------------
def _fit(X: np.ndarray, y: np.ndarray, night_id: np.ndarray):
    Xd = _demean_within(X, night_id)
    yd = y - y.mean()
    yd = _demean_within(yd[:, None], night_id)[:, 0]
    beta, *_ = np.linalg.lstsq(Xd, yd, rcond=1e-8)
    return float(np.sum((yd - Xd @ beta) ** 2)), beta


def _rss_holdout(X: np.ndarray, y: np.ndarray, night_id: np.ndarray,
                 beta: np.ndarray) -> float:
    r = y - X @ beta
    r = _demean_within(r[:, None], night_id)[:, 0]
    return float(r @ r)


def stage2_test(t: Sequence[float], mag: Sequence[float], night_id: Sequence[int],
                period_h: float, alpha: float = 0.05) -> dict:
    """Nested test of P_cat versus 2*P_cat with held-out validation.

    H0 : {f0, 2f0, 3f0, 4f0}                      (8 parameters)
    H1 : H0 + {f0/2, 3f0/2}                       (12 parameters)
    Alternatives with the same four extra parameters:
      * quartic polynomial in time,
      * multiplicative amplitude slope gamma,
        i.e. (t - tbar)/T times the even harmonics.

    The held-out statistic splits nights alternately into two halves,
    fits on one half and evaluates on the other, refitting only that
    half's nightly offsets. A genuine half-period signature must improve
    the held-out residual; generic flexibility does not.

    Returns
    -------
    dict with F_odd, F_crit, F_poly, F_gamma, heldout_drss_half,
    heldout_drss_poly, heldout_drss_gamma, verdict (one of
    'half-period candidate', 'unidentifiable / no gain',
    'not evaluable (rank deficient)').
    """
    t = np.asarray(t, dtype=float)
    mag = np.asarray(mag, dtype=float)
    night_id = np.asarray(night_id)
    y = mag - np.median(mag)
    f0 = 1.0 / float(period_h)
    span = max(float(np.ptp(t)), 1e-9)
    ts = 2.0 * (t - t.min()) / span - 1.0

    base = _harmonic_design(t, f0, 4)                       # f0..4f0
    design_half = _harmonic_design(t, f0 / 2.0, 3)           # f0/2, f0, 3f0/2
    half = np.column_stack([design_half[:, 0:2], design_half[:, 4:6]])
    ev = np.column_stack([base[:, 2:4], base[:, 6:8]])       # 2f0 and 4f0
    poly = np.column_stack([ts ** p for p in (1, 2, 3, 4)])
    gamma = ((t - t.mean()) / span)[:, None] * ev

    rss0, _ = _fit(base, y, night_id)
    dof = max(t.size - np.unique(night_id).size - 12, 1)
    from scipy.stats import f as _f  # local import: only needed here
    fcrit = float(_f.ppf(1.0 - alpha, 4, dof))

    def fstat(extra):
        r, _ = _fit(np.column_stack([base, extra]), y, night_id)
        return max(rss0 - r, 0.0) / 4.0 / (r / dof)

    out = dict(F_odd=fstat(half), F_poly=fstat(poly), F_gamma=fstat(gamma),
               F_crit=fcrit, dof=dof)

    # ---- held-out nights ----
    nights = np.unique(night_id)
    if nights.size < 6:
        out.update(heldout_drss_half=float("nan"), heldout_drss_poly=float("nan"),
                   heldout_drss_gamma=float("nan"),
                   verdict="not evaluable (fewer than 6 nights)")
        return out
    A = nights[::2]
    B = nights[1::2]
    ia = np.isin(night_id, A)
    ib = np.isin(night_id, B)
    if ia.sum() < 30 or ib.sum() < 30:
        out.update(heldout_drss_half=float("nan"), heldout_drss_poly=float("nan"),
                   heldout_drss_gamma=float("nan"),
                   verdict="not evaluable (halves too small)")
        return out

    _, b0 = _fit(base[ia], y[ia], night_id[ia])
    rss_b0 = _rss_holdout(base[ib], y[ib], night_id[ib], b0)
    if not np.isfinite(rss_b0) or math.sqrt(max(rss_b0, 0) / ib.sum()) > 0.6:
        out.update(heldout_drss_half=float("nan"), heldout_drss_poly=float("nan"),
                   heldout_drss_gamma=float("nan"),
                   verdict="not evaluable (rank deficient)")
        return out

    for name, extra in (("half", half), ("poly", poly), ("gamma", gamma)):
        _, bx = _fit(np.column_stack([base, extra])[ia], y[ia], night_id[ia])
        rss_b = _rss_holdout(np.column_stack([base, extra])[ib], y[ib], night_id[ib], bx)
        out["heldout_drss_" + name] = rss_b0 - rss_b

    if out["heldout_drss_half"] > 0 and out["heldout_drss_half"] > max(
            out["heldout_drss_poly"], out["heldout_drss_gamma"]):
        out["verdict"] = "half-period candidate"
    elif out["F_odd"] > fcrit:
        out["verdict"] = "unidentifiable / no gain"
    else:
        out["verdict"] = "no evidence for doubling"
    return out


# ----------------------------------------------------------------------
# 3. command-line demo
# ----------------------------------------------------------------------
def _demo() -> None:
    rng = np.random.default_rng(7)

    def make(P_true, night_h, n_night, n_pt, spacing, sigma=0.05):
        t, nights = [], []
        for k in range(n_night):
            start = k * spacing + rng.uniform(0, 1.5)
            ts = start + np.sort(rng.uniform(0, night_h, n_pt))
            t.extend(ts)
            nights.extend([k] * n_pt)
        t = np.array(t)
        nights = np.array(nights)
        f = 1.0 / P_true
        model = (0.12 * np.cos(2 * np.pi * f * t)
                 + 0.25 * np.cos(4 * np.pi * f * t + 0.7)
                 + 0.03 * np.cos(6 * np.pi * f * t)
                 + 0.08 * np.cos(8 * np.pi * f * t))
        return t, nights, model + rng.normal(0, sigma, t.size)

    for label, Ptrue, L, nnight, spacing in (("well covered", 15.37, 6.0, 12, 24.0),
                                            ("short nights", 31.40, 4.0, 8, 24.0)):
        Pcat = Ptrue / 2.0                     # catalogue holds half the truth
        t, nights, y = make(Ptrue, L, nnight, 30, spacing)
        m = blind_metrics(t, nights, sigma=0.05, period_h=Pcat, amplitude_a1=0.12)
        s = stage2_test(t, y, nights, period_h=Pcat)
        print("\n--- %s: catalogue P = %.2f h, true %.2f h, f0L = %.2f"
              % (label, Pcat, Ptrue, m["f0L"]))
        print("    r0* = %.3f   A1_min = %.3f mag   lambda_odd = %.1f"
              % (m["r0_star"], m["A1_min"], m["lambda_odd"]))
        print("    F_odd = %.1f (crit %.2f)   held-out: half %+.3f, poly %+.3f, gamma %+.3f"
              % (s["F_odd"], s["F_crit"], s["heldout_drss_half"],
                 s["heldout_drss_poly"], s["heldout_drss_gamma"]))
        print("    verdict: %s" % s["verdict"])


if __name__ == "__main__":
    _demo()
