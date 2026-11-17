"""Recompute a compact benchmark figure from the released aggregate counts."""
import argparse
from importlib.resources import files
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from armillary.statistics import failure_statistics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, default=Path("results/benchmark.pdf"))
args = parser.parse_args()
data = json.loads(files("armillary").joinpath("data/benchmark.json").read_text(encoding="utf-8"))
fig, ax = plt.subplots(figsize=(3.5, 2.8), layout="constrained")
for row, color in zip(data["pooled"], ("#4B6B9C", "#B16F49", "#537D68")):
    stats = failure_statistics(row["W1"] + row["W2"], row["N"], row["k"], data["rounds"])
    value = stats["equivalent_rate"]
    lo, hi = stats["equivalent_nominal_95"]
    ax.errorbar(row["d"], value, yerr=np.array([[value - lo], [hi - value]]),
                fmt="o", color=color, capsize=3, markersize=6)
ax.set(yscale="log", xticks=[12, 14, 16], xlabel="Code distance",
       ylabel=r"Logical error rate $P_{\mathrm{L}}$")
ax.set_title(r"$p_{\mathrm{phy}}=0.1\%$", fontsize=10)
ax.grid(axis="y", alpha=.2)
args.out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.out, dpi=240)
print(args.out)
