
import re
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
UPDATES_PER_EPOCH = 4161
STAGE1_BEST_EPOCH = 7
STAGE2_OFFSET = STAGE1_BEST_EPOCH * UPDATES_PER_EPOCH

BLUE = "#2878d8"
ORANGE = "#ee7040"


def read_log(path):
    """Extract training and validation metrics from a DrugCLIP log."""
    train = {}
    valid = {}

    for line in path.read_text(errors="replace").splitlines():
        if "| INFO | train_inner |" in line:
            update = re.search(r"num_updates=(\d+)", line)
            loss = re.search(r"\bloss=([0-9.]+)", line)

            if update and loss:
                # Keep the latest occurrence if a resumed run repeated an update.
                train[int(update.group(1))] = float(loss.group(1))

        elif "| INFO | valid |" in line:
            update = re.search(r"num_updates (\d+)", line)
            loss = re.search(r"\bloss ([0-9.]+)", line)
            auc = re.search(r"\bvalid_auc ([0-9.]+)", line)
            bedroc = re.search(r"\bvalid_bedroc ([0-9.]+)", line)

            if all((update, loss, auc, bedroc)):
                valid[int(update.group(1))] = {
                    "loss": float(loss.group(1)),
                    "auc": float(auc.group(1)),
                    "bedroc": float(bedroc.group(1)),
                }

    return train, valid


stage1_train, stage1_valid = read_log(HERE / "savedir_stage1/train.log")
stage2_train, stage2_valid = read_log(HERE / "savedir_stage2/train.log")

if not all((stage1_train, stage1_valid, stage2_train, stage2_valid)):
    raise RuntimeError("Missing training or validation metrics. Check the log paths.")

fig, axes = plt.subplots(2, 2, figsize=(13, 8))

fig.suptitle(
    "DrugCLIP fine-tuning on RBP pseudo-pockets "
    "(validation: 11,095 held-out pockets, UniProt-grouped split)",
    fontsize=13,
    fontweight="bold",
)

def plot_validation(ax, key, title):
    for data, offset, color, label in [
        (stage1_valid, 0, BLUE, "Stage 1: linear probe"),
        (stage2_valid, STAGE2_OFFSET, ORANGE, "Stage 2: last 4 layers"),
    ]:
        updates = sorted(data)
        ax.plot(


            [x + offset for x in updates],
            [data[x][key] for x in updates],
            "o-",
            color=color,
            linewidth=2.5,
            markersize=5,
            label=label,
        )

    ax.axvline(STAGE2_OFFSET, color="0.85", linewidth=1)
    ax.set_title(title, fontweight="bold")
    ax.grid(alpha=0.25)


plot_validation(axes[0, 0], "loss", "Validation loss (batch 128)")
plot_validation(axes[0, 1], "bedroc", "Validation BEDROC")
plot_validation(axes[1, 0], "auc", "Validation AUC")

ax = axes[1, 1]

for data, offset, color, label in [
    (stage1_train, 0, BLUE, "Stage 1: linear probe"),
    (stage2_train, STAGE2_OFFSET, ORANGE, "Stage 2: last 4 layers"),
]:
    updates = sorted(data)
    ax.plot(
        [x + offset for x in updates],
        [data[x] for x in updates],
        color=color,
        linewidth=1.5,
        label=label,
    )

ax.axvline(STAGE2_OFFSET, color="0.85", linewidth=1)
ax.set_title("Train loss (per 100 updates)", fontweight="bold")
ax.grid(alpha=0.25)

for ax in axes.flat:
    ax.set_xlabel("Cumulative training updates")
    ax.tick_params(labelsize=9)

axes[0, 0].set_ylabel("Validation loss")
axes[0, 1].set_ylabel("BEDROC")
axes[1, 0].set_ylabel("AUC")
axes[1, 1].set_ylabel("Training loss")

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center",
           bbox_to_anchor=(0.5, 0.94), ncol=2)

fig.tight_layout(rect=[0, 0, 1, 0.88])

output = HERE / "stage1_stage2_comparison"
fig.savefig(output.with_suffix(".png"), dpi=200)
fig.savefig(output.with_suffix(".pdf"))
plt.close(fig)

print(f"Stage 1: {len(stage1_valid)} validation points")
print(f"Stage 2: {len(stage2_valid)} validation points")
print(f"Saved: {output}.png")
print(f"Saved: {output}.pdf")
