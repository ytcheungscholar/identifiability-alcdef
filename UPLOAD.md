# How to get this to GitHub (3 steps, ~3 minutes)

This folder is the release-ready v1.1.0 tree for
`github.com/ytcheungscholar/identifiability-alcdef`.

## 1. Upload

**Easiest (no git needed):** open the repo on github.com → *Add file → Upload files*
→ drag the **contents** of this folder in (not the folder itself) → commit straight
to the default branch.
GitHub's web uploader takes up to 100 files / 25 MB per file, which this tree fits
(≈7.5 MB total, largest file 1.9 MB).

**Or with GitHub Desktop:** *File → Add local repository* → point at this folder →
*Commit to main* → *Publish/Push*.

**With git CLI (if you ever install it):**
```
cd <this folder>
git init -b main
git remote add origin https://github.com/ytcheungscholar/identifiability-alcdef.git
git add -A
git commit -m "v1.1.0: four-condition screen, higher-harmonic control, block audit"
git push -u origin main --force
```

## 2. Cut the release that Zenodo watches

GitHub → *Releases → Draft a new release*:

- Tag: **v1.1.0** (create the tag on publish)
- Title: `v1.1.0 — four-condition screen and ALCDEF block audit`
- Notes: paste the top section of `README.md` (headline numbers table)
- Publish.

If the Zenodo ↔ GitHub integration is enabled for this repo, publishing the
release mints the DOI immediately and shows it on the release page.

## 3. Send me the DOI

Give me the number (e.g. `10.5281/zenodo.12345678`) and I will:

- replace `10.5281/zenodo.XXXXXXX` in `main_v7.tex` (Data availability) **and** in
  `cover_letter_v8.tex`,
- recompile and verify `0 error / 0 undefined / no ??`,
- hand back the final PDF for submission.

## What is in here

```
README.md            v1.1.0, headline numbers, artefact -> script -> data map
CITATION.cff         1.1.0, 2026-11-21
code/                13 scripts (paper's Data availability names match exactly)
data/v1.1/           13 derived tables behind every number in the paper
figs/                Figures 1-4 as shipped
data/*.csv (root)    v1.0 tables, kept for provenance only -- superseded
```

Legacy files at the root (`reproduce_figures.py`, `fill_placeholders.py`,
`_prior_bins*.py`, `_table7_recheck.py`) belong to the v1.0 layout. Delete them
after the upload if you want a clean tree; nothing in `code/` or `data/v1.1/`
depends on them.
