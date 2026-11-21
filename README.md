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
| cadences audited a priori | 12,070 (9,958 asteroids) | `out_a1exact_20261119.csv` |
| physical sample before guards | 8,906 | `out_full_lcdb_20261119.csv` |
| campaigns in the nested test | 6,657 | `out_lono_v5_c34.csv` |
| in-sample rejection rate | 58.1 % | `reproduce_3f4f_control.py` |
| held-out improvement, subharmonic | 9.3 % | same |
| held-out improvement, quartic | 9.2 % | same |
| held-out improvement, slow modulation | 22.5 % | same |
| held-out improvement, higher harmonics {3f₀,4f₀} | **46.7 %** | same |
| three-condition screen, archive | 4.99 % (332/6,657) | same |
| four-condition screen, archive | **3.35 %** (223/6,657) | same |
| four-condition, archive, frequency re-fitted per fold | **3.00 %** (200/6,654) | `reproduce_closure.py` |
| four-condition, campaigns with no merged-block night | 2.76 % (151/5,464) | `reproduce_block_audit.py` |
| four-condition nulls (strict / drift / stationary / amp-mod) | 1.9 / 0.8 / 1.1 / 1.4 % | `analyse_closure.py` |
| positive control $U=3$, three conditions | 55.5 % (1,264/2,277) | `reproduce_u3_control.py` |
| positive control $U=3$, four conditions | **53.4 %** (1,216/2,277) | `reproduce_closure.py` |
| fourth condition, in-sample rate | 59.2 % (3,942/6,656) | `reproduce_table6_harmonic.py` |
| seventh null arm: strict H₀ + drift, four conditions | 1.9 / 2.0 % (A₁ = 0.05 / 0.10) | `reproduce_table7_seventh_arm.py` |

The four-condition screen keeps 53.4 % of a known doubling while the archive
falls to 3.35 %: a 15.9× separation under one screen, one fold structure and one
exact per-campaign critical value. 3.35 % is above every cluster-bootstrap
interval of the four-condition nulls (highest upper bound 2.7 %) but drops to
2.76 % once the nights that merge two photometric systems are removed, which is
the upper edge of those nulls. Whether the excess is astrophysical, a
correlated-noise effect, or photometric mixing that the generator does not
reproduce is left open in the paper (§4.1, §6h).

## Repository layout

Everything is at the top level — no subdirectories.

```
reproduce_*.py, analyse_closure.py, fig2.py, fig3.py, fig4.py
                    the paper's scripts; each writes the table it reproduces
out_*.csv, *.npz    the derived tables behind every number in the paper
fig*.png            Figures 1-4 as shipped
README.md           this file
CITATION.cff        citation metadata (v1.1.0)
LICENSE             MIT
```

Files that belong to the **v1.0 layout** and are kept only for provenance:
`lightcurve_identifiability.py`, `reproduce_figures.py`, `fill_placeholders.py`,
`_prior_bins*.py`, `_table7_recheck.py`, `alcdef_identifiability_catalog.csv`,
`ridge_sigma_zp.csv`, `stage2_holdout_results.csv`,
`validation_28_campaigns.csv`, `test_smoke.py`, `requirements.txt`. None of the
v1.1 scripts or tables depends on them.

## What regenerates what

| artefact | script | input | output |
|---|---|---|---|
| Table 2 (M=2 vs M=4, interaction grid) | `reproduce_table2.py` | synthetic cadences | `out_table2_n200.csv` |
| Table 6 / Fig. 3a (held-out decomposition, four controls) | `reproduce_3f4f_control.py` | ALCDEF campaigns | `out_lono_v5_c34.csv` |
| Table 6, fourth row (in-sample rate of {3f₀,4f₀}) | `reproduce_table6_harmonic.py` | ALCDEF campaigns | `out_insample_c34_20261121.csv` |
| Table 7 (six-arm calibration) | `reproduce_table7.py` | 200 real cadences | `calib_v5b_a05.npz`, `calib_v5b_a10.npz` |
| Table 7, seventh arm (strict H₀ + drift) | `reproduce_table7_seventh_arm.py` | the same 200 cadences | `calib_v8k_a05.npz`, `calib_v8k_a10.npz` |
| bootstrap intervals on those arms | `analyse_closure.py` | the npz files | printed intervals |
| Table 8 / Fig. 4 (calibration-prior kernel) | `reproduce_table8.py`, `fig4.py` | all 12,070 cadences | `out_table8_full_20261121.csv` |
| §4.1 block audit (nights that merge systems) | `reproduce_block_audit.py` | ALCDEF block metadata | `out_block_audit_20261121.csv` |
| §4.2.2 frequency-refit pass | `reproduce_freqopt.py` | campaigns | `out_lono_v5_freqopt.csv` |
| four conditions at the re-fitted frequency | `reproduce_closure.py freqopt` | the above | `out_freqopt_c34.csv` |
| §4.2.4 positive control | `reproduce_u3_control.py` | $U=3$ campaigns | `out_u3_control_20261121.csv` |
| positive control, four conditions | `reproduce_closure.py u3` | the above | `out_u3_c34.csv` |
| Figs. 2–4 | `fig2.py`, `fig3.py`, `fig4.py` | shipped tables | `fig*.png` |

Fig. 1 (a, b) is assembled from `out_a1exact_20261119.csv` and
`out_insample_v5.csv`; it is shipped as a PNG and is not regenerated by a single
script.

## Running

Each script sets two absolute paths at the top:

```python
H    = r"...\ac_repo_plan"     # where the derived CSVs live
ROOT = r"...\alcdef_data"      # the unpacked ALCDEF files
```

Point `H` at this checkout (the CSVs sit at the top level) and `ROOT` at an
ALCDEF dump, then run, e.g. `python reproduce_table8.py`. Requirements:
Python ≥ 3.10, `numpy`, `scipy`, `pandas`, `matplotlib`. The archive itself is
public: ALCDEF (Warner 2011) and the LCDB (Warner et al. 2009).

Runtime: Table 8 over 12,070 cadences takes minutes; the calibration arms (six
arms of the table plus the seventh drift null, each 200 cadences × 8
realisations) take about half an hour per reference amplitude; the per-fold
frequency-refit pass takes about a quarter of an hour on a desktop CPU.

## Two conventions that matter

* **Exact critical values.** Every screen uses $F_{\rm crit}(4,\nu)$ with the
  campaign's own $\nu = N-8-K$, never a fixed 2.43. Using the fixed value shifts
  the archive rate from 4.99 % to 5.03 % and the four-condition count from 223
  to 225 — small, but it is the difference between a reproducible and an
  approximate number.
* **Paired noise across arms.** In the calibration the noise realisation and the
  harmonic phases are drawn once per campaign and shared by all arms, so adding
  an arm cannot move the numbers of the others. The arms must be run together;
  running them separately breaks the pairing.

## Caveats the data will not let you forget

* The lowest coverage bin ($f_0L<0.05$) is a $0/0$ degeneracy: $r_0^*$ is large
  there for the wrong reason, and the recovery fraction is not a recovery
  measurement. Only the $U\ge3$, $\le30$ d subset gives a clean floor.
* The synthetic nulls carry white noise and harmonics up to the fourth only.
  They contain no within-night correlated noise, no non-Gaussian outliers and no
  $m\ge5$ content.
* 6.7 % of nights contain more than one ALCDEF lightcurve block, and campaigns
  with such a night pass the screen at twice the rate of the rest; the paper
  reports this as an open confound, not a measured artefact.

## Citation

See `CITATION.cff`. If you use the derived catalogue or the diagnostics, please
cite both the paper and this repository.

## Licence

MIT (see `LICENSE`).
