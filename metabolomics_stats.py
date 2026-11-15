"""
metabolomics_stats.py
=====================
Takes a feature table and produces: a PCA plot, a statistics table,
a volcano plot and a peaks plot.

Input  (data/feature_table.csv):
    columns "sample" and "group" (Control / Treated / QC),
    then one column per feature holding peak intensities.

Outputs (results/ folder):
    pca.png      - do the groups separate? are the QC samples tight?
    results.csv  - fold change, p-value and FDR for every feature
    volcano.png  - picture of which features changed
    peaks.png    - the peaks themselves: chromatogram and mass spectrum view

The nine steps:
    1. remove features with too many missing values
    2. fill the remaining gaps
    3. remove features that are unstable in the QC samples
    4. PQN normalization
    5. log transform and Pareto scaling
    6. PCA
    7. Welch t-test and FDR correction
    8. volcano plot
    9. peaks plot

Needs: numpy, pandas, scipy, scikit-learn, matplotlib

Run it with:
    python metabolomics_stats.py
(run make_example_data.py first if data/feature_table.csv does not exist)

Source of the method: https://codanics.com/metabolomics-complete-guide/
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")   # draw plots straight to files; no screen needed
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

# Folders, relative to where this script lives.
HERE = Path(__file__).resolve().parent
IN_FILE = HERE / "data" / "feature_table.csv"
OUT_DIR = HERE / "results"
OUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Load the table and split it into two parts:
#   meta = who each sample is (name and group)
#   X    = the numbers (rows = samples, columns = features)
# ---------------------------------------------------------------------------
df = pd.read_csv(IN_FILE)
meta, X = df[["sample", "group"]], df.drop(columns=["sample", "group"])
print("Start:", X.shape[1], "features")

# ---------------------------------------------------------------------------
# 1. Remove features missing in more than half of the samples.
# A feature that is absent most of the time cannot be tested reliably.
# X.isna().mean() gives, per feature, the fraction of samples with no value.
# ---------------------------------------------------------------------------
X = X.loc[:, X.isna().mean() <= 0.5]

# ---------------------------------------------------------------------------
# 2. Fill the remaining gaps (imputation).
# A missing value in MS data usually means "too low to detect", not "zero".
# So we fill each gap with a small number: 1/5 of the smallest value seen
# for that feature. This is a simple, widely used choice.
# ---------------------------------------------------------------------------
X = X.fillna(X.min() / 5)

# ---------------------------------------------------------------------------
# 3. Keep only features that are stable in the pooled QC samples.
# Every QC is the same mixture, so a good feature should give nearly the
# same intensity in all of them. We measure the spread with RSD
# (relative standard deviation = standard deviation / mean x 100).
# Features with RSD of 30% or more are too noisy and are dropped.
# ---------------------------------------------------------------------------
qc = X[meta.group == "QC"]
rsd = qc.std() / qc.mean() * 100
X = X.loc[:, rsd < 30]
print("After missing-value and QC-RSD filters:", X.shape[1], "features")

# ---------------------------------------------------------------------------
# 4. PQN normalization (Probabilistic Quotient Normalization).
# Problem: some samples are simply more concentrated than others, which
# raises ALL their features together. That is not biology.
# Solution:
#   a) build a reference profile = median QC intensity of each feature
#   b) divide each sample by the reference, feature by feature (quotients)
#   c) the median quotient of a sample is its dilution factor
#   d) divide the whole sample by that factor
# Using the median means a few truly changed features do not disturb it.
# ---------------------------------------------------------------------------
reference = X[meta.group == "QC"].median()
quotients = X.div(reference, axis=1)
X = X.div(quotients.median(axis=1), axis=0)

# ---------------------------------------------------------------------------
# 5. Log transform, then Pareto scaling.
# Log2: raw intensities are strongly skewed (a few huge peaks). Taking the
#   log makes them roughly bell-shaped, which the t-test expects, and turns
#   fold changes into simple differences (2x up = +1, 2x down = -1).
# Pareto scaling: centre each feature on zero and divide by the SQUARE ROOT
#   of its standard deviation. This stops large peaks from dominating the
#   PCA, without blowing up the noise of small peaks.
# L is used for the statistics, P only for the PCA.
# ---------------------------------------------------------------------------
L = np.log2(X)
P = (L - L.mean()) / np.sqrt(L.std())

# ---------------------------------------------------------------------------
# 6. PCA (Principal Component Analysis).
# PCA squeezes hundreds of features into two axes (PC1, PC2) that capture
# the largest differences between samples. It uses no group labels, so it
# is an honest first look at the data. We want to see:
#   - Control and Treated in separate clouds  -> there is a real effect
#   - QC samples packed tightly together      -> the measurement was stable
# ---------------------------------------------------------------------------
pca = PCA(n_components=2)
scores = pca.fit_transform(P)                 # position of each sample
var = pca.explained_variance_ratio_ * 100     # % of variation per axis

# Colour-blind-friendly colours for the three groups.
colours = {"Control": "#0072B2", "Treated": "#D55E00", "QC": "#999999"}
fig, ax = plt.subplots(figsize=(6, 5))
for g, c in colours.items():
    m = (meta.group == g).values              # True for samples of group g
    ax.scatter(scores[m, 0], scores[m, 1], c=c, label=g, s=60, edgecolor="white")
ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
ax.legend(frameon=False)
ax.set_title("PCA score plot")
fig.tight_layout()
fig.savefig(OUT_DIR / "pca.png", dpi=200)
print(f"PCA: PC1 explains {var[0]:.1f}% and PC2 explains {var[1]:.1f}% of the variation")

# ---------------------------------------------------------------------------
# 7. One statistical test per feature, then correct for multiple testing.
# Welch t-test: compares the mean of Treated against Control for one
#   feature. The Welch version does not assume equal spread in both groups.
#   QC samples are not part of the test.
# p-value: chance of seeing a difference this large if nothing really changed.
# The problem: with hundreds of tests, about 5% will have p < 0.05 by pure
#   luck. Benjamini-Hochberg FDR correction fixes this. The corrected value
#   (q) controls the share of false alarms among the features we call
#   significant.
# ---------------------------------------------------------------------------
c, t = L[meta.group == "Control"], L[meta.group == "Treated"]
p = stats.ttest_ind(t, c, equal_var=False).pvalue
q = stats.false_discovery_control(p, method="bh")   # Benjamini-Hochberg

# The data are already in log2 units, so the difference of the means IS
# the log2 fold change. +1 means doubled in Treated, -1 means halved.
log2fc = t.mean() - c.mean()

res = pd.DataFrame({"feature": L.columns, "log2FC": log2fc.values, "p": p, "q_FDR": q})

# A feature counts as significant only if it passes BOTH rules:
#   statistically reliable : q below 0.05
#   large enough to matter : at least 2-fold up or down (|log2FC| >= 1)
res["significant"] = (res.q_FDR < 0.05) & (res.log2FC.abs() >= 1)
res = res.sort_values("q_FDR")                       # strongest results first
res.to_csv(OUT_DIR / "results.csv", index=False)

n_up = int((res.significant & (res.log2FC > 0)).sum())
n_down = int((res.significant & (res.log2FC < 0)).sum())
print("Significant features (FDR < 0.05 and at least 2-fold):", int(res.significant.sum()))
print(f"  higher in Treated: {n_up}, lower in Treated: {n_down}")
print(res.head(5).to_string(index=False))

# ---------------------------------------------------------------------------
# 8. Volcano plot.
# One dot per feature.
#   x-axis: log2 fold change  (left = lower in Treated, right = higher)
#   y-axis: -log10 p-value    (higher = stronger evidence)
# Interesting features sit in the top-left and top-right corners.
#   orange = significantly higher in Treated
#   blue   = significantly lower in Treated
#   grey   = not significant
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 5))
col = np.where(~res.significant, "#BBBBBB", np.where(res.log2FC > 0, "#D55E00", "#0072B2"))
ax.scatter(res.log2FC, -np.log10(res.p), c=col, s=18)
ax.axvline(-1, ls="--", c="grey", lw=0.8)   # 2-fold down line
ax.axvline(1, ls="--", c="grey", lw=0.8)    # 2-fold up line
ax.set_xlabel("log2 fold change (Treated / Control)")
ax.set_ylabel("-log10 p-value")
ax.set_title("Volcano plot")
fig.tight_layout()
fig.savefig(OUT_DIR / "volcano.png", dpi=200)

# ---------------------------------------------------------------------------
# 9. Peaks plot.
# The plots above show statistics. This one shows the peaks themselves, in
# the two views an LC-MS instrument gives:
#   top    = chromatogram : signal against retention time
#   bottom = mass spectrum: signal against m/z
# Each feature name is "<m/z>_<retention time>", so we can split the name
# to know where every peak belongs on both axes.
# The heights are the mean PQN-normalized intensities of each group (X after
# step 4, before the log), so Control and Treated can be compared fairly.
# IMPORTANT: a feature table stores only ONE number per peak (its size), not
# its shape. So the chromatogram is REBUILT from the table: every peak is
# drawn as a bell curve of the same width. It is not the raw detector trace.
# ---------------------------------------------------------------------------
try:
    mz = np.array([float(name.split("_")[0]) for name in X.columns])
    rt = np.array([float(name.split("_")[1]) for name in X.columns])
except (ValueError, IndexError):
    mz = rt = None
    print("Peaks plot skipped: feature names are not in the form <m/z>_<retention time>")

if mz is not None:
    mean_c = X[meta.group == "Control"].mean().values   # peak heights, Control
    mean_t = X[meta.group == "Treated"].mean().values   # peak heights, Treated
    # True for the features that step 7 called significant (same order as X).
    sig = res.set_index("feature").significant.reindex(X.columns).values

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

    # --- Top: rebuilt chromatogram ---
    # A bell curve (Gaussian) at each retention time, as tall as the peak
    # intensity. Adding all the curves gives one trace per group.
    width = 0.04                                         # peak width (sigma), minutes
    grid = np.linspace(rt.min() - 0.5, rt.max() + 0.5, 4000)
    bells = np.exp(-0.5 * ((grid[:, None] - rt[None, :]) / width) ** 2)
    ax1.plot(grid, bells @ mean_c, c=colours["Control"], lw=1, label="Control (mean)")
    ax1.plot(grid, bells @ mean_t, c=colours["Treated"], lw=1, label="Treated (mean)")
    ax1.set_xlabel("Retention time (min)")
    ax1.set_ylabel("Intensity")
    ax1.set_title("Chromatogram rebuilt from the feature table")
    ax1.legend(frameon=False)

    # --- Bottom: mirror mass spectrum ---
    # One vertical stick per feature at its m/z. Control sticks point up and
    # Treated sticks point down, so the same peak can be compared directly:
    # a stick that is much longer on one side changed between the groups.
    ax2.vlines(mz, 0, mean_c, color=colours["Control"], lw=1, label="Control (mean)")
    ax2.vlines(mz, 0, -mean_t, color=colours["Treated"], lw=1, label="Treated (mean)")
    # Black dots mark the significant features, on the side where they are higher.
    tip = np.where(mean_t > mean_c, -mean_t, mean_c)
    ax2.scatter(mz[sig], tip[sig], c="black", s=14, zorder=3, label="Significant feature")
    ax2.axhline(0, c="black", lw=0.8)
    # Show the downward axis with positive numbers (it is still an intensity).
    ax2.yaxis.set_major_formatter(lambda y, _: f"{abs(y):g}")
    ax2.set_xlabel("m/z")
    ax2.set_ylabel("Intensity")
    ax2.set_title("Mass spectrum view: Control up, Treated down")
    # Leave empty space above the tallest stick for the legend.
    top = 1.3 * max(mean_c.max(), mean_t.max())
    ax2.set_ylim(-top, top)
    ax2.legend(frameon=False, ncol=3, loc="upper center")

    fig.tight_layout()
    fig.savefig(OUT_DIR / "peaks.png", dpi=200)

print(f"Saved the plots and results.csv in {OUT_DIR}")
