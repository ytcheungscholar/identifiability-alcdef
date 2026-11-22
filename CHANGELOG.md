# Changelog

## v1.2.0 — 22 November 2026

**Bug fixed (calibration only; no analytic result or population audit changes).**
In `reproduce_table7.py` (v1.1.1) and its sibling
`reproduce_table7_seventh_arm.py`, the measurement noise was drawn once per
campaign before the arm loop and then **again** inside the loop:

```python
noise = rng.normal(0, sig, t.size)      # once per campaign
for arm in ARMS:
    ...
    y += noise                           # intended
    y += rng.normal(0, sig, t.size)      # second draw — defect
```

Two consequences:

1. the effective noise variance was 2σ² instead of σ², so the injected signal had
   a signal-to-noise ratio a factor √2 lower than the paper described, and the
   power of the screen was systematically underestimated;
2. the second draw inside the arm loop meant the arms no longer shared one noise
   realisation, and adding an arm shifted the random stream of the others — the
   reason the v1.1.1 seventh arm differed from the six-arm run by up to 1.2
   percentage points.

Both scripts are replaced by a single `reproduce_table7.py` that generates all
seven arms in one pass with a single, shared noise draw.

**Table 7 (percentage of campaigns that pass), v1.2.0**

| arm | three 0.05 | three 0.10 | four 0.05 | four 0.10 |
|---|---|---|---|---|
| true 2P | 29.1 | 45.1 | 28.7 | 44.3 |
| true 2P + drift | 26.8 | 39.7 | 26.1 | 38.9 |
| true P | 5.0 | 5.5 | 1.2 | 1.0 |
| true P, amp-mod | 5.6 | 5.9 | 1.4 | 1.2 |
| true P, drift | 5.6 | 5.9 | 1.4 | 1.2 |
| strict H₀ | 1.7 | 1.7 | 1.6 | 1.6 |
| strict H₀ + drift | 2.1 | 2.2 | 2.0 | 2.2 |

Comparing with v1.1.1: power rises by 5.5–9.3 percentage points, the null arms
move by less than one point, and the cost of the fourth condition falls to
0.4/0.8 points (from 0.9/1.2).

**Knock-on changes**

* The four-condition nulls now bracket 1.0–2.2 % instead of 0.4–2.0 %; the
  widest interval is 1.2–3.0 % (strict H₀ + drift).
* The archive rate 3.35 % is still above every null *point estimate*, but the
  widest null interval now reaches it, so the manuscript no longer calls it an
  excess "clear of every null interval".
* `reproduce_table7_seventh_arm.py` retired; kept as
  `legacy_v111_reproduce_table7_seventh_arm.py`.
* Legacy artefacts renamed `legacy_v111_calib_v5b_a05/a10.npz` and
  `legacy_v111_calib_v8k_a05/a10.npz`. They reproduce the v1.1.1 table and must
  not be mixed with the v1.2.0 numbers.

**Packaging**

* `reproduce_block_split.py` added: it reproduces the two block-split rates of
  §4.1 (6.04 % for campaigns with a merged-block night, 2.76 % without) and the
  asteroid-cluster bootstrap quoted in §4.2.4 (5,744 clusters), which the v1.1.1
  release did not ship a script for. Percentiles of the bootstrap move by about
  0.05 between seeds; the rates do not.
* `out_zp_baseline_20261121.csv` added: it is read by `fig2.py` and
  `reproduce_table8.py` and was missing from the v1.1.1 archive.
* `requirements.txt` now lists `pandas` (used by `fig2.py`, `fig4.py`,
  `reproduce_table8.py`, `reproduce_closure.py`).
* Scripts read `AC_ROOT` and `ALCDEF_ROOT` from the environment instead of
  hard-coded Windows paths.
* Single flat layout, matching the GitHub repository root.

## v1.1.1 — 21 November 2026

Packaging only: tag created so that the Zenodo record matches the submitted
manuscript; no numerical change with respect to v1.1.0.

## v1.1.0 — 21 November 2026

Every population number recomputed on 12,070 cadences. Supersedes v1.0.0.

## v1.0.0 — 14 November 2026

First release.
