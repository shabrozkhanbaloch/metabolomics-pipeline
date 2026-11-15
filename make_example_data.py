"""
make_example_data.py
====================
Creates a realistic but SIMULATED metabolomics feature table to practise on.

Why simulate?
    With made-up data we know the right answer in advance. We plant
    30 features that truly change between the groups, and later we check
    whether the statistics script finds them.

What the output looks like (data/feature_table.csv):
    - one row per sample
    - column "sample"  : sample name (C1..C8, T1..T8, QC1..QC6)
    - column "group"   : Control, Treated or QC
    - 300 feature columns, each named "<m/z>_<retention time>"
      and holding the peak intensity of that feature in that sample

This is the same shape of table that tools such as MZmine, MS-DIAL or XCMS
export after peak picking, so the statistics script will also work on real
data that is arranged this way.

Needs: numpy, pandas

Run it with:
    python make_example_data.py

Source of the method: https://codanics.com/metabolomics-complete-guide/
"""
from pathlib import Path

import numpy as np
import pandas as pd

# Where to save the table: a "data" folder next to this script.
OUT_FILE = Path(__file__).resolve().parent / "data" / "feature_table.csv"

# A fixed seed makes the "random" numbers the same on every run,
# so everybody who runs this script gets the identical table.
rng = np.random.default_rng(42)

# Size of the experiment:
#   300 features (metabolite signals)
#   8 samples in each biological group (Control and Treated)
#   6 pooled quality-control (QC) samples
n_features, n_per_group, n_qc = 300, 8, 6

# ---------------------------------------------------------------------------
# Step 1. Give every feature a typical intensity.
# Real LC-MS intensities are log-normal: a few very large peaks and many
# small ones. So we draw the baseline from a log-normal distribution.
# ---------------------------------------------------------------------------
base = rng.lognormal(mean=10, sigma=1.2, size=n_features)

# ---------------------------------------------------------------------------
# Step 2. Plant the "true" biology.
# `effect` is the fold change of each feature in Treated compared to Control.
# 1.0 means no change. We change only the first 30 features:
# ---------------------------------------------------------------------------
effect = np.ones(n_features)
effect[:15] = rng.uniform(2.0, 4.0, 15)     # features 1-15 go UP 2x to 4x
effect[15:30] = rng.uniform(0.25, 0.5, 15)  # features 16-30 go DOWN 2x to 4x
# The other 270 features stay at 1.0 (no real difference).


def samples(n, fold, noise=0.25):
    """
    Simulate `n` samples.

    Each intensity = baseline x fold change x random noise.
    `noise` is the spread of the random variation (biological + technical).
    The noise is multiplicative (log-normal), as in real MS data.
    """
    return base * fold * rng.lognormal(0, noise, size=(n, n_features))


# ---------------------------------------------------------------------------
# Step 3. Build the three groups of samples.
# ---------------------------------------------------------------------------
control = samples(n_per_group, 1.0)      # no effect applied
treated = samples(n_per_group, effect)   # planted effect applied

# A pooled QC is made in the lab by mixing a little of every sample.
# So its level is the average of Control (1) and Treated (effect).
# All QCs come from the same tube, so they differ only by small
# technical noise (0.08 instead of 0.25).
qc = samples(n_qc, (1 + effect) / 2, noise=0.08)

# ---------------------------------------------------------------------------
# Step 4. Add two problems that real data always has.
# ---------------------------------------------------------------------------
# Problem A: unequal sample amounts. Some samples are more concentrated
# than others (for example dilute urine, or a slightly different injection).
# Every feature in one sample is multiplied by the same factor (0.7 to 1.3).
# The statistics script must remove this with normalization.
dilution = rng.uniform(0.7, 1.3, size=(2 * n_per_group + n_qc, 1))
X = np.vstack([control, treated, qc]) * dilution

# Problem B: missing values. About 3% of the values are blanked out,
# as happens when a peak falls below the detection limit.
X[rng.random(X.shape) < 0.03] = np.nan

# ---------------------------------------------------------------------------
# Step 5. Add sample names, group labels and feature names.
# ---------------------------------------------------------------------------
names = (
    [f"C{i+1}" for i in range(n_per_group)]
    + [f"T{i+1}" for i in range(n_per_group)]
    + [f"QC{i+1}" for i in range(n_qc)]
)
groups = ["Control"] * n_per_group + ["Treated"] * n_per_group + ["QC"] * n_qc

# Each feature is identified by its m/z and its retention time (minutes).
# These are random here; they only serve as realistic-looking labels.
mz = np.round(rng.uniform(70, 900, n_features), 4)
rt = np.round(rng.uniform(0.5, 15, n_features), 2)
cols = [f"{m}_{r}" for m, r in zip(mz, rt)]   # feature name = m/z_retention time

# ---------------------------------------------------------------------------
# Step 6. Save as a CSV file.
# ---------------------------------------------------------------------------
df = pd.DataFrame(X, columns=cols)
df.insert(0, "group", groups)
df.insert(0, "sample", names)

OUT_FILE.parent.mkdir(exist_ok=True)   # create the data folder if needed
df.to_csv(OUT_FILE, index=False)
print(f"Saved {OUT_FILE}: {df.shape[0]} samples x {n_features} features")
