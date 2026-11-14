# Identifiability diagnostics for ground-based asteroid lightcurves

Companion code and derived data for

> **Two identifiability regimes in ground-based asteroid period determination:
> nightly zero points, the *P* → *P*/2 flip, and an information-collapse limit**
> Y. T. Cheung, submitted to *Astronomy and Computing*.

The package answers, *before any period search is run*, how much of a
rotation period a ground-based lightcurve can support. It is based on one
observation: with an unknown photometric zero point on every night, the
statistically correct estimator is equivalent to subtracting each night's
mean (Frisch–Waugh–Lovell), and that operation is a high-pass filter on the
lightcurve harmonics. The consequences are analytic and depend on the
timestamps, not on the fitting algorithm.

## Quickstart

```python
import numpy as np
from lightcurve_identifiability import blind_metrics, stage2_test, assign_nights

# t in hours, mag in magnitudes, sigma the per-point uncertainty
night_id = assign_nights(t)                    # groups into nights (7.2 h gap)
m = blind_metrics(t, night_id, sigma=0.05, period_h=16.49, amplitude_a1=0.06)
# -> r0_star, r_star, eta_f0, f0L, A1_min, lambda_odd, can_resolve

s = stage2_test(t, mag, night_id, period_h=16.49)
# -> F_odd, F_crit, heldout_drss_half / _poly / _gamma, verdict
```

`blind_metrics` returns

| key | meaning |
|---|---|
| `r0_star` | geometric flip threshold on $A_4/A_1$; if the object's $A_4/A_1$ exceeds it, an $M=2$ fit locks onto $P/2$ |
| `r_star` | same, including the third-harmonic correction |
| `eta_f0` | fraction of fundamental variance that survives the nightly projection |
| `f0L` | nightly phase coverage in cycles per night |
| `A1_min` | smallest fundamental amplitude for which the odd-harmonic test reaches 80 % power |
| `lambda_odd` | non-centrality of that test for the amplitude you pass in |
| `can_resolve` | `A1 >= A1_min` |

`stage2_test` returns a nested $F$-test of $P_\mathrm{cat}$ versus
$2P_\mathrm{cat}$ **plus the two controls the paper shows are necessary**:

* `heldout_drss_poly` — a quartic polynomial with the same four extra
  parameters (generic flexibility),
* `heldout_drss_gamma` — a multiplicative amplitude slope, the physical form
  of seasonal phase-angle evolution,
* `verdict` — `half-period candidate` only when the half-integer terms
  improve the held-out nights *and* beat both alternatives.

Run `python lightcurve_identifiability.py` for a two-case demonstration
(one well-covered campaign that is flagged, one short-night campaign that is
correctly reported as not evaluable).

## What is in this repository

```
lightcurve_identifiability.py      core module (~300 lines, numpy only)
test_smoke.py                      13 synthetic checks, run with python test_smoke.py
reproduce_figures.py               rebuilds Figs. 1-4 from the tables below
data/alcdef_identifiability_catalog.csv   6,164 campaigns, 4,172 asteroids
data/stage2_holdout_results.csv           562 evaluable campaigns, in-sample and held-out statistics
data/validation_28_campaigns.csv          28 real cadences, predicted vs injection-recovery r*
data/ridge_sigma_zp.csv                   softened threshold vs zero-point prior width
figs/                              output of reproduce_figures.py
```

`reproduce_figures.py` reproduces every figure with no archive access:
```
fig1: n=28  r=0.993  rms=0.026
fig2: n=1813 campaigns, danger 19.0%, flip 19.1%
fig3: in-sample 87.5%, held-out half 16.2%, poly 12.5%, gamma 29.7%
fig4: 6 curves
```

## Tests

```
python test_smoke.py        # 13 tests, < 1 s
```

The suite includes explicit regression guards on the layout of the Stage-2
design blocks (the half-integer columns must sit at $f_0/2$ and $3f_0/2$;
the amplitude-slope control must use the even harmonics). Two column-index
mistakes of exactly that kind were caught this way during development, and
either would have left `F_odd` looking plausible while destroying the
meaning of the test.

## The derived catalogue

`data/alcdef_identifiability_catalog.csv` is built in one pass over the
24,643 files of the ALCDEF archive matched to the LCDB period catalogue:
apparitions split at 90-day gaps, then physical cuts (≥30 measurements,
≥3 nights with ≥3 points, median night length 1–10 h, span ≤400 d).
**No amplitude cut is applied in the table** — the amplitude gate is applied
at analysis time, because the population fraction depends on it
(19.0 % for the full sample rising to 31.1 % when a strong second harmonic
is required). Columns:

`num, name, file, app, jd0, span_d, nights, npt, L, P, f0L, A1..A4, sig,
A1sig, A2sig, r4, rstar, d2, d4, U, amp_lcdb`

with `r4 = A4/A1`, `rstar` the flip threshold, `d2`/`d4` the M=2/M=4
residual differences at $f_0$ and $2f_0$, and `U` the LCDB quality code.

## Caveats worth knowing before you use it

1. **Never test at a period commensurate with the 24 h sampling cycle.**
   With $P = 8$ h exactly, the night-to-night phase pattern is periodic and
   the held-out halves become degenerate. The paper's synthetic work uses
   incommensurate periods for this reason, and the demonstration does too.
2. **The held-out test needs six nights.** With fewer, or with alternating
   nights that do not constrain the eight-parameter base, the verdict is
   `not evaluable`; in the archive that happens for 45 % of campaigns, which
   is itself a measurement of the identifiability limit.
3. **`R_prior` needs a population.** Pass `population_r4` (the `r4` column of
   the catalogue) to get a flip probability; without it the value is omitted.
4. **The threshold assumes one zero point per night.** Per-observer offsets,
   or nights split across two telescopes, must be passed as separate nights.

## Requirements

Python ≥ 3.9, `numpy`; `scipy` only for `stage2_test` (F-distribution
quantiles); `matplotlib` only for `reproduce_figures.py`.

## Data sources

* ALCDEF — Warner (2021), NASA PDS, doi:10.26033/b8cw-s522
* LCDB — Warner, Harris & Pravec (2009), Icarus 202, 134

## License

MIT — see `LICENSE`.

## Citation

See `CITATION.cff`. Please cite both the paper and this repository (Zenodo
DOI to be inserted on publication).

* Repository: <https://github.com/ytcheungscholar/identifiability-alcdef>
* Author ORCID: [0009-0006-7227-4993](https://orcid.org/0009-0006-7227-4993)
