# Publishing v1.2.0

This folder is the release-ready **v1.2.0** tree (22 November 2026) for
`github.com/ytcheungscholar/identifiability-alcdef`. It is **flat**: every file
sits at the top level, exactly like the current GitHub repository.

## What changed in v1.2.0

- `reproduce_table7.py` now generates all **seven** calibration arms in one pass,
  with **one shared noise draw** per campaign/realisation. The v1.1.1 defect
  (noise added twice) is fixed.
- Table 7 / Fig. 3b were rerun: power 29.1/45.1 % (three conditions),
  28.7/44.3 % (four); nulls 1.0–2.2 %.
- New `calib_v8l_a05.npz`, `calib_v8l_a10.npz`.
- `reproduce_table7_seventh_arm.py` retired; kept as
  `legacy_v111_reproduce_table7_seventh_arm.py`.
- Old calibration snapshots renamed `legacy_v111_calib_v5b_*.npz` and
  `legacy_v111_calib_v8k_*.npz`. **Do not mix them with v1.2.0 numbers.**
- `out_zp_baseline_20261121.csv` added (read by `fig2.py` and
  `reproduce_table8.py`; it was missing from v1.1.1).
- `requirements.txt` now lists `pandas` (plus numpy, scipy, matplotlib).
- `CITATION.cff` is 1.2.0 / 2026-11-22.
- `README.md` and `CHANGELOG.md` updated.

## 1. Upload (no git needed)

Open <https://github.com/ytcheungscholar/identifiability-alcdef>, choose
**Add file → Upload files**, and drag the **contents** of this folder (not the
folder itself) into the page. Commit to `main`.

Then delete the five files that v1.2.0 supersedes, which the web uploader does
not remove automatically:

- `calib_v5b_a05.npz`
- `calib_v5b_a10.npz`
- `calib_v8k_a05.npz`
- `calib_v8k_a10.npz`
- `reproduce_table7_seventh_arm.py`

Their replacements are the `legacy_v111_*` files and `calib_v8l_*`.

## 2. Or with the portable git (OpenClaw)

This machine already has a portable git at:

```
C:\Users\jacky\AppData\Local\OpenClaw\deps\portable-git\mingw64\bin\git.exe
```

From a normal terminal (so that the credential prompt can appear):

```
cd <this folder>
git init -b main
git remote add origin https://github.com/ytcheungscholar/identifiability-alcdef.git
git add -A
git commit -m "v1.2.0: fix calibration noise defect; seven-arm Table 7 in one pass"
git push -u origin main
```

If the push asks for credentials, use a **fine-grained personal access token**
with `Contents: Read and write` on this repository as the password (GitHub no
longer accepts the account password). Do not paste the token into a chat.

## 3. Cut the release

GitHub → **Releases → Draft a new release**:

- Tag: **v1.2.0** (create the tag on publish)
- Title: `v1.2.0 — calibration noise fix and seven-arm Table 7`
- Notes: paste the top of `README.md` and the `CHANGELOG.md` v1.2.0 section.
- Publish.

The concept DOI is `10.5281/zenodo.24422013` and always resolves to the latest
archived version; v1.1.1 is `10.5281/zenodo.24422014`. If Zenodo is connected to
the repository, publishing the release mints a new version DOI for v1.2.0; the
paper cites only the concept DOI, so no manuscript change is needed.

## File map

```
README.md, CHANGELOG.md, CITATION.cff, LICENSE, requirements.txt
reproduce_*.py, analyse_closure.py, fig2.py, fig3.py, fig4.py
out_*.csv, calib_v8l_*.npz, legacy_v111_*.npz
fig*.png
```

`fill_placeholders.py`, `lightcurve_identifiability.py`, `reproduce_figures.py`,
`_prior_bins*.py`, `_table7_recheck.py`, `test_smoke.py` are v1.0/v1.1 provenance
files; they are kept for the record and are not needed by the v1.2.0 pipeline.
