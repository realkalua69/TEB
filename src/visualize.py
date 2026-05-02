"""
Generate KDE density plot of inter-residue vector angles,
colored by preceding-residue size category.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

INPUT = "output/angles.tsv"
OUTPUT = "output/distribution.png"

df = pd.read_csv(INPUT, sep="\t")
df["angle_deg"] = ((df["angle_deg"] + 180) % 360) - 180

categories = ["Tiny", "Small", "Mid", "Large", "Bulky"]
palette = {
    "Tiny":  "#4e79a7",
    "Small": "#59a14f",
    "Mid":   "#f28e2b",
    "Large": "#e15759",
    "Bulky": "#76b7b2",
}

fig, ax = plt.subplots(figsize=(11, 6))
fig.patch.set_facecolor("#f0f0f0")
ax.set_facecolor("#f0f0f0")

x_range = np.linspace(-180, 180, 512)

for cat in categories:
    subset = df[df["size_cat"] == cat]["angle_deg"].values
    if len(subset) < 3:
        continue
    density = gaussian_kde(subset, bw_method=0.12)
    ax.plot(x_range, density(x_range), color=palette[cat], lw=2.0, label=cat)

for tick in range(-150, 180, 50):
    ax.axvline(tick, color="#999999", ls="--", lw=0.4, alpha=0.7)

ax.set_xlim(-180, 180)
ax.set_xticks(range(-150, 200, 50))
ax.set_xlabel("Inter-residue vector angle (\u00b0)", fontsize=12)
ax.set_ylabel("Density (a.u.)", fontsize=12)
ax.set_title(f"Helix Tripeptide Angle Distribution (N={len(df):,})", fontsize=13)
ax.legend(title="Size class", loc="upper right", framealpha=0.9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(OUTPUT, dpi=250, facecolor="#f0f0f0")
print(f"plot saved: {OUTPUT} ({len(df)} points)")
