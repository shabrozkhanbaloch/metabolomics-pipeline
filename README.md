# Metabolomics statistics pipeline (practice run)

This folder holds a small, working metabolomics data-analysis pipeline in Python.
It follows the guide at <https://codanics.com/metabolomics-complete-guide/>.

The pipeline was run on **simulated practice data**, not on real samples.
This file explains what is in the folder and how to run it. The results are
explained in [RESULTS.md](RESULTS.md).

## What is in this folder

| File | What it does |
|---|---|
| `mz_calculator.py` | Works out which m/z values a metabolite gives in the mass spectrometer, and matches a measured m/z back to a metabolite. |
| `make_example_data.py` | Creates the practice feature table (`data/feature_table.csv`). |
| `metabolomics_stats.py` | The main analysis: cleaning, normalization, PCA, statistical tests, a volcano plot and a peaks plot. |
| `data/feature_table.csv` | The practice data: 22 samples and 300 features. |
| `workflow_infographic.py` | Draws the one-page workflow picture (`results/workflow.png`). |
| `results/workflow.png` | The whole workflow in one picture, with the numbers of the last run. |
| `results/pca.png` | PCA plot. |
| `results/volcano.png` | Volcano plot. |
| `results/peaks.png` | Peaks plot: chromatogram view and mass spectrum view. |
| `results/results.csv` | Fold change, p-value and FDR for every tested feature. |
| `RESULTS.md` | Step-by-step explanation of every result and how to read the plots. |
| `requirements.txt` | The five Python packages needed, with the tested versions. |
| `LICENSE` | MIT license: anyone may use and change the code. |

Every script has comments that explain each step and why it is done.

## How to run it

Set up a Python environment once, from this folder:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Then run the four scripts in this order:

```bash
.venv/bin/python mz_calculator.py
.venv/bin/python make_example_data.py
.venv/bin/python metabolomics_stats.py
.venv/bin/python workflow_infographic.py
```

Tested with Python 3.14. `requirements.txt` installs the exact package
versions that were tested.

## Where this fits in a full metabolomics study

A full study has 12 steps. These scripts cover only the statistics part.

| Step | Covered here? |
|---|---|
| 1. Study design | No (planning) |
| 2–4. Sample collection, extraction, LC-MS or NMR measurement | No (lab work) |
| 5. Convert raw files to mzML | No (use ProteoWizard msConvert) |
| 6. Peak picking and alignment | No (use MZmine, MS-DIAL or XCMS) |
| 7. Filter features | **Yes** |
| 8. Drift and batch correction | No |
| 9. Normalize, log-transform, scale | **Yes** |
| 10. Statistics | **Yes** (PCA, Welch t-test, FDR) |
| 11. Identify metabolites | Only the simple m/z match in `mz_calculator.py` |
| 12. Pathway analysis | No (use MetaboAnalyst) |

## What the analysis does, and what it found

The analysis has nine steps: clean the feature table, normalize it, run a PCA
and a statistical test for every feature, then draw the plots. The practice
data hold 30 features that truly change, so the result can be checked.

- 300 features at the start, 254 left after cleaning.
- **18 features found as significant**, and all 18 are real (no false alarms).
- 12 of the 30 real changes were missed.

[RESULTS.md](RESULTS.md) has the full account: the workflow in one picture,
every step explained, all plots and tables, and why the 12 were missed.

## Limits of the pipeline

- It compares exactly two groups, named `Control` and `Treated`, and needs
  samples labelled `QC`.
- There is no correction for signal drift or batch effects.

## Using your own data

Save your feature table as `data/feature_table.csv` in this layout, then run
`metabolomics_stats.py`:

| sample | group | 181.0707_2.31 | 166.0863_5.12 | ... |
|---|---|---|---|---|
| C1 | Control | 152340 | 8812 | ... |
| T1 | Treated | 310221 | 9050 | ... |
| QC1 | QC | 228900 | 8930 | ... |

Do not run `make_example_data.py` after that, because it overwrites the file
with the practice data.
