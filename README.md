# Identifiability diagnostics for ground-based asteroid lightcurves

Companion code and derived data for

> **Two identifiability regimes in ground-based asteroid period determination:
> nightly zero points, the *P* → *P*/2 flip, and an information-collapse limit**
> Y. T. Cheung, *Astronomy and Computing* (submitted 21 November 2026).
> Preprint DOI: `10.5281/zenodo.XXXXXXX`

The package answers, *before any period search is run*, how much of a rotation
period a ground-based lightcurve can support. It rests on one observation: with
an unknown photometric zero point on every night, the statistically correct
estimator is equivalent to subtracting each night's mean
(Frisch–Waugh–Lovell), and that operation is a high-pass filter on the harmonic
content of the lightcurve. The consequences are analytic and depend on the
timestamps, not on the fitting algorithm.

**Version 1.1.0** (21 November 2026). Supersedes v1.0.0 (14 November 2026);
every population number changed in v1.1, so v1.0 tables must not be mixed with
these.

## Headline numbers

| quantity | value | source |
|---|---|---|
| cadences audited a priori | 12,070 (9,958 asteroids) | `data/v1.1/out_a1exact_20261119.csv` |
| physical sample before guards | 8,906 | `data/v1.1/out_full_lcdb_20261119.csv` |
| campaigns in the nested test | 6,657 | `data/v1.1/out_lono_v5_c34.csv` |
| in-sample rejection rate | 58.1 % | `reproduce_3f4f_control.py` |
| held-out improvement, subharmonic | 9.3 % | same |
| held-out improvement, quartic | 9.2 % | same |
| held-out improvement, slow modulation | 22.5 % | same |
| held-out improvement, higher harmonics {3f₀,4f₀} | **46.7 %** | same |
| three-condition screen, archive | 4.99 % (332/6,657) | same |
| four-condition screen, archive | **3.35 %** (223/6,657) | same |
| four-condition, archive, frequency re-fitted per fold | **3.00 %** (200/6,654) | `reproduce_closure.py` |
| four-condition nulls (strict / drift / stationary / amp-mod) | 1.9 / 0.8 / 1.1 / 1.4 % | `analyse_closure.py` |
| positive control $U=3$, three conditions | 55.5 % (1,264/2,277) | `reproduce_u3_control.py` |
| positive control $U=3$, four conditions | **53.4 %** (1,216/2,277) | `reproduce_closure.py` |
| fourth condition, in-sample rate | 59.2 % (3,942/6,656) | `reproduce_table6_harmonic.py` |
| seventh null arm: strict H₀ + drift, four conditions | 1.9 / 2.0 % (A₁ = 0.05 / 0.10) | `reproduce_table7_seventh_arm.py` |

The four-condition screen keeps 53.4 % of a known doubling while the archive
falls to 3.35 %: a 15.9× separation under one screen, one fold structure and one
exact per-campaign critical value. 3.35 % is above every cluster-bootstrap
interval of the four-condition nulls (highest upper bound 2.7 %); whether that
excess is astrophysical or a correlated-noise effect that the generator does not
reproduce is left open in the paper.

## Repository layout

```
code/     reproducing scripts (each one writes the CSV its table/figure needs)
data/v1.1 derived tables needed by code/ (no archive download required)
figs/     figures as shipped
data/     v1.0 tables, kept only for provenance -- superseded, do not use
```

## What regenerates what

| artefact | script | input | output |
|---|---|---|---|
| Table 2 (M=2 vs M=4, interaction grid) | `code/reproduce_table2.py` | synthetic cadences | `data/v1.1/out_table2_n200.csv` |
| Table 6 / Fig. 3a (held-out decomposition, four controls) | `code/reproduce_3f4f_control.py` | ALCDEF campaigns | `out_lono_v5_c34.csv` |
| Table 7 (six-arm calibration) | `code/reproduce_table7.py` | 200 real cadences | `calib_v5b_a05.npz`, `calib_v5b_a10.npz` |
| Table 7 bootstrap intervals | `code/analyse_closure.py` | the two npz files | printed intervals |
| Table 8 / Fig. 4 (calibration-prior kernel) | `code/reproduce_table8.py`, `code/fig4.py` | all 12,070 cadences | `out_table8_full_20261121.csv` |
| §4.2.2 frequency-refit pass | `code/reproduce_freqopt.py` | campaigns | `out_lono_v5_freqopt.csv` |
| four conditions at the re-fitted frequency | `code/reproduce_closure.py freqopt` | the above | `out_freqopt_c34.csv` |
| §4.2.4 positive control | `code/reproduce_u3_control.py` | $U=3$ campaigns | `out_u3_control_20261121.csv` |
| positive control, four conditions | `code/reproduce_closure.py u3` | the above | `out_u3_c34.csv` |
| Figs. 2–4 | `code/fig2.py`, `code/fig3.py`, `code/fig4.py` | shipped tables | `figs/*.png` |

Fig. 1 (a, b) is assembled from `data/v1.1/out_a1exact_20261119.csv` and
`data/v1.1/out_insample_v5.csv`; it is shipped as a PNG and is not regenerated
by a single script.

## Running

Every script sets two absolute paths at the top:

```python
H    = r"...\ac_repo_plan"     # where the derived CSVs live
ROOT = r"...\alcdef_data"      # the unpacked ALCDEF files
```

Point both at your checkout (`H` → `data/v1.1`, `ROOT` → the ALCDEF dump) and
run, e.g. `python code/reproduce_table8.py`. Requirements: Python ≥ 3.10,
`numpy`, `scipy`, `pandas`, `matplotlib`. The archive itself is public: ALCDEF
(Warner 2011) and the LCDB (Warner et al. 2009).

Runtime: Table 8 over 12,070 cadences takes minutes; the calibration arms
(six arms of the table plus the seventh drift null, each 200 cadences ×
8 realisations) take about half an hour per reference amplitude; the per-fold
frequency-refit pass takes about a quarter of an hour on a desktop CPU.

## Two conventions that matter

* **Exact critical values.** Every screen uses $F_{\rm crit}(4,\nu)$ with the
  campaign's own $\nu = N-8-K$, never a fixed 2.43. Using the fixed value
  shifts the archive rate from 4.99 % to 5.03 % and the four-condition count
  from 223 to 225 — small, but it is the difference between a reproducible and
  an approximate number.
* **Paired noise across arms.** In the calibration, the noise realisation and
  the harmonic phases are drawn once per campaign and shared by all arms, so
  adding an arm cannot move the numbers of the others. The six arms must be run
  together; running them separately breaks the pairing.

## Caveats the data will not let you forget

* The lowest coverage bin ($f_0L<0.05$) is a $0/0$ degeneracy: $r_0^*$ is large
  there for the wrong reason, and the recovery fraction is not a recovery
  measurement. Only the $U\ge3$, $\le30$ d subset gives a clean floor.
* The synthetic nulls carry white noise and harmonics up to the fourth only.
  They contain no within-night correlated noise, no non-Gaussian outliers and
  no $m\ge5$ content; the paper states the resulting excess as a signal whose
  origin is not yet established, not as a prevalence measurement.

## Citation

See `CITATION.cff`. If you use the derived catalogue or the diagnostics, please
cite both the paper and this repository.

## Licence

MIT (see `LICENSE`).
