"""Publication figure for CLEAN stage-1 (linear probe) training.

Parses clean/savedir_stage1/train.log and writes, next to this script:
  stage1_clean_training.{pdf,svg,png}  - 2x2 figure
  stage1_clean_metrics.csv             - per-epoch source data
"""
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
LOG = HERE.parent / "savedir_stage1" / "train.log"

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def parse(log):
    valid, train = {}, {}
    for line in log.read_text().splitlines():
        m = re.search(r"\| valid \| epoch (\d+) .*?\| loss ([\d.]+) \| valid_acc ([\d.]+)"
                      r" .*?valid_auc ([\d.]+) \| valid_bedroc ([\d.]+) \| valid_top3_acc ([\d.]+)", line)
        if m:
            e = int(m[1])
            valid[e] = dict(valid_loss=float(m[2]), top1=float(m[3]), auc=float(m[4]),
                            bedroc=float(m[5]), top3=float(m[6]))
            continue
        m = re.search(r"\| train \| epoch (\d+) \| scale [\d.]+ \| loss ([\d.]+)", line)
        if m:
            train[int(m[1])] = float(m[2])
    return [dict(epoch=e, train_loss=train[e], **valid[e]) for e in sorted(valid)]


rows = parse(LOG)
with open(HERE / "stage1_clean_metrics.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

ep = [r["epoch"] for r in rows]
col = lambda k: [r[k] for r in rows]
best = max(rows, key=lambda r: (r["bedroc"], -r["epoch"]))["epoch"]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})

fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.6), sharex=True)
(a, b), (c, d) = axes


def style(ax, ylabel, letter):
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.axvline(best, color=INK2, linestyle=(0, (3, 2)), linewidth=0.7, zorder=1)
    ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            color=INK, va="bottom")
    ax.set_xlim(0.5, ep[-1] + 0.5)
    ax.set_xticks([1, 5, 10, 15, 20, 24])


def line(ax, y, color=BLUE, label=None):
    ax.plot(ep, y, color=color, linewidth=1.4, marker="o", markersize=2.8,
            markeredgecolor="white", markeredgewidth=0.5, label=label, zorder=3)


line(a, col("bedroc"))
style(a, "Validation BEDROC (α = 80.5)", "a")
bv = rows[best - 1]["bedroc"]
a.set_ylim(top=bv + 0.007)
a.text(best, bv + 0.0035, f"best: epoch {best}, {bv:.3f}", ha="center", va="bottom",
       fontsize=6.5, color=INK, bbox=dict(facecolor="white", edgecolor="none", pad=0.6))

line(b, col("auc"))
style(b, "Validation ROC-AUC", "b")

line(c, col("top1"), BLUE, "Top-1")
line(c, col("top3"), ORANGE, "Top-3")
style(c, "Validation retrieval accuracy", "c")
c.legend(frameon=False, loc="lower right", fontsize=6.5, handlelength=1.6)
for k, lab, color in (("top1", "Top-1", BLUE), ("top3", "Top-3", ORANGE)):
    c.text(ep[-1] + 0.6, rows[-1][k], lab, color=INK, fontsize=6.5, va="center", clip_on=False)

line(d, col("train_loss"))
style(d, "Training loss (in-batch softmax)", "d")

for ax in (c, d):
    ax.set_xlabel("Epoch")

fig.tight_layout(w_pad=2.5, h_pad=1.2)
for ext in ("pdf", "svg"):
    fig.savefig(HERE / f"stage1_clean_training.{ext}")
fig.savefig(HERE / "stage1_clean_training.png", dpi=300)
print(f"{len(rows)} epochs, best epoch {best} (BEDROC {bv})")
