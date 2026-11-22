#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Smoke tests for lightcurve_identifiability.

Run with:  python test_smoke.py

These are deliberately small and self-contained (no archive access). They
exist because two column-indexing bugs in the Stage-2 design matrix were
caught only by exactly this kind of synthetic check:

  * the "half-integer" block was silently built from {f0, 2f0} instead of
    {f0/2, 3f0/2};
  * the "even harmonics" block used for the amplitude-slope control was
    taken from the wrong columns.

Either mistake leaves F_odd looking reasonable while destroying the
interpretation of the test, so the layout of the design blocks is asserted
explicitly below rather than only through the final verdict.
"""

import os
import math
import unittest

import numpy as np

from lightcurve_identifiability import (
    assign_nights,
    blind_metrics,
    gram_residuals,
    population_flip_risk,
    stage2_test,
    LAMBDA_CRIT,
)


def synth(period_h=15.37, night_h=6.0, n_night=12, n_pt=30, spacing=24.0,
          sigma=0.05, seed=7, half_period_catalogue=True):
    """A lightcurve with harmonics 1..4 at the true period."""
    rng = np.random.default_rng(seed)
    t, nights = [], []
    for k in range(n_night):
        start = k * spacing + rng.uniform(0, 1.5)
        ts = start + np.sort(rng.uniform(0, night_h, n_pt))
        t.extend(ts)
        nights.extend([k] * n_pt)
    t = np.array(t)
    nights = np.array(nights)
    f = 1.0 / period_h
    model = (0.12 * np.cos(2 * np.pi * f * t)
             + 0.25 * np.cos(4 * np.pi * f * t + 0.7)
             + 0.03 * np.cos(6 * np.pi * f * t)
             + 0.08 * np.cos(8 * np.pi * f * t))
    return t, nights, model + rng.normal(0, sigma, t.size)


class TestNights(unittest.TestCase):
    def test_gap_splits_nights(self):
        t = np.array([0.0, 1.0, 22.0, 23.0, 50.0])
        nid = assign_nights(t, gap_h=7.2)
        self.assertEqual(nid.tolist(), [0, 0, 1, 1, 2])

    def test_dense_sampling_is_one_night(self):
        t = np.arange(0.0, 48.0, 0.08)
        self.assertEqual(len(set(assign_nights(t).tolist())), 1)


class TestGram(unittest.TestCase):
    def setUp(self):
        self.t, self.nights, _ = synth()

    def test_r0_matches_definition(self):
        g = gram_residuals(self.t, self.nights, 15.37)
        self.assertAlmostEqual(g["r0_star"], math.sqrt(g["a1"] / g["a4"]), places=10)

    def test_components_are_positive_and_bounded(self):
        g = gram_residuals(self.t, self.nights, 15.37)
        for key in ("a1", "a4"):
            self.assertGreater(g[key], 0.0)
            self.assertLess(g[key], 2.0)

    def test_retention_falls_with_coverage(self):
        """A short period means more cycles per night, hence less loss."""
        short = gram_residuals(self.t, self.nights, 5.0)["eta_f0"]
        long_ = gram_residuals(self.t, self.nights, 40.0)["eta_f0"]
        self.assertGreater(short, long_)


class TestDesignBlocks(unittest.TestCase):
    """Regression guards for the two indexing bugs described in the docstring."""

    def setUp(self):
        self.t, self.nights, _ = synth()
        self.f0 = 1.0 / 15.37

    def test_half_block_frequencies(self):
        from lightcurve_identifiability import _harmonic_design
        base = _harmonic_design(self.t, self.f0, 4)
        dh = _harmonic_design(self.t, self.f0 / 2.0, 3)
        half = np.column_stack([dh[:, 0:2], dh[:, 4:6]])
        ang = 2 * np.pi * self.t
        self.assertTrue(np.allclose(half[:, 0], np.cos(ang * self.f0 / 2)))
        self.assertTrue(np.allclose(half[:, 2], np.cos(ang * 1.5 * self.f0)))
        # and it must NOT coincide with the fundamental of the base block
        self.assertFalse(np.allclose(half[:, 0], base[:, 0]))

    def test_even_block_frequencies(self):
        from lightcurve_identifiability import _harmonic_design
        base = _harmonic_design(self.t, self.f0, 4)
        ev = np.column_stack([base[:, 2:4], base[:, 6:8]])
        ang = 2 * np.pi * self.t
        self.assertTrue(np.allclose(ev[:, 0], np.cos(ang * 2 * self.f0)))
        self.assertTrue(np.allclose(ev[:, 2], np.cos(ang * 4 * self.f0)))


class TestStage2(unittest.TestCase):
    def test_detects_a_true_doubled_period(self):
        t, nights, y = synth()
        s = stage2_test(t, y, nights, period_h=15.37 / 2)
        self.assertEqual(s["verdict"], "half-period candidate")
        self.assertGreater(s["F_odd"], s["F_crit"])
        self.assertGreater(s["heldout_drss_half"], 0)
        self.assertGreater(s["heldout_drss_half"],
                           max(s["heldout_drss_poly"], s["heldout_drss_gamma"]))

    def test_no_false_positive_when_period_is_correct(self):
        """Data generated at P must not be promoted to 2P."""
        t, nights, y = synth(half_period_catalogue=False)
        s = stage2_test(t, y, nights, period_h=15.37)
        self.assertNotEqual(s["verdict"], "half-period candidate")

    def test_short_nights_are_reported_as_not_evaluable(self):
        t, nights, y = synth(period_h=31.40, night_h=2.0, n_night=6)
        s = stage2_test(t, y, nights, period_h=31.40 / 2)
        self.assertTrue(s["verdict"].startswith("not evaluable"))

    def test_commensurate_period_is_flagged_not_trusted(self):
        """P = 8 h with 24 h sampling is degenerate; the module must not
        claim a detection there (the paper's rule: use incommensurate
        periods for synthetic work)."""
        t, nights, y = synth(period_h=16.0, night_h=4.0, n_night=8, spacing=24.0)
        s = stage2_test(t, y, nights, period_h=8.0)
        self.assertNotEqual(s["verdict"], "half-period candidate")


class TestBlindMetrics(unittest.TestCase):
    def test_keys_and_amplitude_flag(self):
        t, nights, _ = synth()
        m = blind_metrics(t, nights, sigma=0.05, period_h=15.37, amplitude_a1=0.12)
        for key in ("r0_star", "r_star", "eta_f0", "f0L", "A1_min", "lambda_odd",
                    "can_resolve"):
            self.assertIn(key, m)
        self.assertTrue(m["can_resolve"])
        small = blind_metrics(t, nights, sigma=0.05, period_h=15.37,
                              amplitude_a1=1e-4)
        self.assertFalse(small["can_resolve"])

    def test_priors(self):
        self.assertTrue(np.isclose(population_flip_risk(0.5, [0.2, 0.6, 0.9]), 2 / 3))
        self.assertTrue(np.isclose(population_flip_risk(0.0, [0.2, 0.6]), 1.0))
        self.assertEqual(LAMBDA_CRIT, 12.2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
