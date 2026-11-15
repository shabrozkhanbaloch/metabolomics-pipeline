"""
workflow_infographic.py
=======================
Draws a one-page picture of the whole analysis workflow, from the feature
table to the final plots. The numbers on the picture (samples, features,
significant features) are read from the data and the results, so the
picture always matches the last run.

Input:
    data/feature_table.csv   - the feature table
    results/results.csv      - made by metabolomics_stats.py
    results/pca.png, volcano.png, peaks.png

Output:
    results/workflow.png

Needs: numpy, pandas, matplotlib

Run it with (after metabolomics_stats.py):
    python workflow_infographic.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")   # draw straight to a file; no screen needed
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
RES = HERE / "results"

# ---------------------------------------------------------------------------
# Collect the numbers shown on the picture.
# The three filter lines below repeat steps 1-3 of metabolomics_stats.py,
# only to count how many features each step removes.
# ---------------------------------------------------------------------------
df = pd.read_csv(HERE / "data" / "feature_table.csv")
groups = df.group.value_counts()
X = df.drop(columns=["sample", "group"])
n_start = X.shape[1]
n_gaps = int(X.isna().sum().sum())
pct_gaps = 100 * n_gaps / X.size

X = X.loc[:, X.isna().mean() <= 0.5]            # step 1
n_after_missing = X.shape[1]
X = X.fillna(X.min() / 5)                       # step 2
qc = X[df.group == "QC"]
n_after_qc = int((qc.std() / qc.mean() * 100 < 30).sum())   # step 3

res = pd.read_csv(RES / "results.csv")
n_sig = int(res.significant.sum())
n_up = int((res.significant & (res.log2FC > 0)).sum())
n_down = n_sig - n_up

# ---------------------------------------------------------------------------
# Drawing helpers. The page is a 100 x 150 grid; y runs from top (150) down.
# ---------------------------------------------------------------------------
INK, MUTED, PAPER = "#1F2933", "#52606D", "#FFFFFF"
fig = plt.figure(figsize=(10, 15), facecolor="#F4F6F8")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 100)
ax.set_ylim(0, 150)
ax.axis("off")


def card(y_top, height, colour, stage, title):
    """A white rounded box with a coloured strip and a stage title."""
    ax.add_patch(FancyBboxPatch((6, y_top - height), 88, height,
                                boxstyle="round,pad=0,rounding_size=1.5",
                                fc=PAPER, ec="#D9DEE3", lw=1))
    ax.add_patch(FancyBboxPatch((6, y_top - height), 2.2, height,
                                boxstyle="round,pad=0,rounding_size=1.1",
                                fc=colour, ec="none"))
    ax.text(11, y_top - 3.2, stage, color=colour, fontsize=10, fontweight="bold", va="center")
    ax.text(11, y_top - 6.6, title, color=INK, fontsize=15, fontweight="bold", va="center")


def step(y, number, colour, text, note):
    """One numbered step: a coloured circle, what is done, and the result."""
    ax.scatter(13, y, s=420, color=colour, zorder=3)
    ax.text(13, y, str(number), color="white", fontsize=11, fontweight="bold",
            ha="center", va="center", zorder=4)
    ax.text(17, y + 1.0, text, color=INK, fontsize=11.5, va="center")
    ax.text(17, y - 1.6, note, color=MUTED, fontsize=9.5, va="center")


def badge(y, colour, big, small):
    """A large number on the right side of a card."""
    ax.text(84, y + 1.2, big, color=colour, fontsize=30, fontweight="bold", ha="center", va="center")
    ax.text(84, y - 2.2, small, color=MUTED, fontsize=9.5, ha="center", va="top", linespacing=1.3)


def arrow(y_from, y_to):
    """A down arrow between two cards."""
    ax.annotate("", xy=(50, y_to), xytext=(50, y_from),
                arrowprops=dict(arrowstyle="-|>", color="#9AA5B1", lw=2, mutation_scale=18))


# ---------------------------------------------------------------------------
# Title
# ---------------------------------------------------------------------------
ax.text(6, 145, "Metabolomics statistics workflow", color=INK, fontsize=24, fontweight="bold", va="center")
ax.text(6, 141, "From a feature table to the features that changed  |  simulated practice data",
        color=MUTED, fontsize=11.5, va="center")

# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------
c = "#52606D"
card(137, 15, c, "INPUT", "Feature table")
ax.text(11, 127.3, f"{len(df)} samples: {groups.get('Control', 0)} Control, "
        f"{groups.get('Treated', 0)} Treated, {groups.get('QC', 0)} QC",
        color=INK, fontsize=11.5, va="center")
ax.text(11, 124.4, f"One column per feature (m/z_retention time).  "
        f"{n_gaps} values missing ({pct_gaps:.1f}%).", color=MUTED, fontsize=9.5, va="center")
badge(130, c, str(n_start), "features")
arrow(121.5, 118.5)

# ---------------------------------------------------------------------------
# Clean
# ---------------------------------------------------------------------------
c = "#0072B2"
card(118, 28, c, "CLEAN", "Remove what cannot be trusted")
step(106.5, 1, c, "Remove features missing in more than half of the samples",
     f"{n_start - n_after_missing} removed")
step(100.0, 2, c, "Fill the remaining gaps with 1/5 of the feature's smallest value",
     f"{n_gaps} gaps filled")
step(93.5, 3, c, "Keep only features that are stable in the QC samples (RSD below 30%)",
     f"{n_after_missing - n_after_qc} removed")
badge(110, c, str(n_after_qc), "features kept")
arrow(89.5, 86.5)

# ---------------------------------------------------------------------------
# Normalize
# ---------------------------------------------------------------------------
c = "#009E73"
card(86, 21.5, c, "NORMALIZE", "Make samples comparable")
step(74.5, 4, c, "PQN normalization", "removes differences in sample concentration")
step(68.0, 5, c, "Log2 transform, then Pareto scaling",
     "log2 for the tests; Pareto scaling for the PCA only")
arrow(64.0, 61.0)

# ---------------------------------------------------------------------------
# Analyse
# ---------------------------------------------------------------------------
c = "#D55E00"
card(60.5, 24.5, c, "ANALYSE", "Find the features that changed")
step(49.0, 6, c, "PCA", "do the groups separate? are the QC samples tight?")
step(42.5, 7, c, "Welch t-test per feature, then Benjamini-Hochberg FDR",
     "significant = FDR below 0.05 and at least 2-fold change")
badge(53, c, str(n_sig), f"significant\n{n_up} higher, {n_down} lower\nin Treated")
arrow(35.5, 32.5)

# ---------------------------------------------------------------------------
# Output: the three plots as small pictures, plus the table
# ---------------------------------------------------------------------------
c = "#7B5EA7"
card(32, 29.5, c, "OUTPUT", "Plots and table (results folder)")
ax.scatter([78.5, 84.5], [26.4, 26.4], s=420, color=c, zorder=3)
for x, n in [(78.5, 8), (84.5, 9)]:
    ax.text(x, 26.4, str(n), color="white", fontsize=11, fontweight="bold", ha="center", va="center", zorder=4)
ax.text(75.5, 26.4, "steps", color=MUTED, fontsize=9.5, ha="right", va="center")

thumbs = [("pca.png", "pca.png\nPCA score plot"),
          ("volcano.png", "volcano.png\nvolcano plot (step 8)"),
          ("peaks.png", "peaks.png\npeaks plot (step 9)")]
for i, (name, label) in enumerate(thumbs):
    left = 0.11 + i * 0.21                       # figure fractions
    tax = fig.add_axes([left, 0.052, 0.19, 0.095])
    tax.imshow(mpimg.imread(RES / name))
    tax.set_xticks([]), tax.set_yticks([])
    for side in tax.spines.values():
        side.set_edgecolor("#D9DEE3")
    ax.text((left + 0.095) * 100, 5.9, label, color=MUTED, fontsize=8.5, ha="center", va="center",
            linespacing=1.3)

ax.add_patch(FancyBboxPatch((75.5, 8.3), 15.5, 13.2, boxstyle="round,pad=0,rounding_size=1",
                            fc="#F4F6F8", ec="#D9DEE3", lw=1))
ax.text(83.25, 17.2, "results.csv", color=INK, fontsize=11, fontweight="bold", ha="center", va="center")
ax.text(83.25, 12.6, "fold change,\np-value and FDR\nfor every feature", color=MUTED, fontsize=8.5,
        ha="center", va="center", linespacing=1.3)

fig.savefig(RES / "workflow.png", dpi=150, facecolor=fig.get_facecolor())
print(f"Saved {RES / 'workflow.png'}")
