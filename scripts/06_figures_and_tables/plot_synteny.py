import matplotlib.pyplot as plt
import numpy as np
import os

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Figures"
coords_file = "/mnt/d/Mitochondria/analysis/SRR9057080/synteny/srr9057080_vs_ref.coords"

ref_starts, ref_ends = [], []
qry_starts, qry_ends = [], []
identities = []

with open(coords_file) as f:
    for line in f:
        parts = line.strip().split('|')
        if len(parts) >= 6:
            try:
                s1, e1 = map(int, parts[0].split())
                s2, e2 = map(int, parts[1].split())
                idy = float(parts[3].strip())
                ref_starts.append(s1)
                ref_ends.append(e1)
                qry_starts.append(s2)
                qry_ends.append(e2)
                identities.append(idy)
            except Exception:
                continue

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Subplot A: Dot Plot of Syntenic Alignment
for s1, e1, s2, e2 in zip(ref_starts, ref_ends, qry_starts, qry_ends):
    color = '#1f77b4' if (e1 - s1) * (e2 - s2) > 0 else '#d62728' # Blue = forward collinear, Red = inverted
    ax1.plot([s1 / 1000.0, e1 / 1000.0], [s2 / 1000.0, e2 / 1000.0], color=color, lw=1.2, alpha=0.7)

ax1.set_title("A. High-Resolution Syntenic Dot Plot\n(NC_043914.1 Reference vs De Novo Assembly)", fontsize=11, fontweight='bold')
ax1.set_xlabel("Reference NC_043914.1 Coordinate (kb)", fontsize=10)
ax1.set_ylabel("Assembled Scaffolds Coordinate (kb)", fontsize=10)
ax1.grid(True, linestyle='--', alpha=0.5)

# Subplot B: Syntenic Identity Distribution
ax2.hist(identities, bins=25, color='#2ca02c', edgecolor='black', alpha=0.85)
ax2.axvline(x=np.mean(identities), color='red', linestyle='--', lw=1.5, label=f"Mean Identity: {np.mean(identities):.2f}%")
ax2.set_title("B. Sequence Identity Distribution Across 878 Syntenic Blocks", fontsize=11, fontweight='bold')
ax2.set_xlabel("Sequence Identity (%)", fontsize=10)
ax2.set_ylabel("Number of Syntenic Blocks", fontsize=10)
ax2.legend(loc='upper left')
ax2.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
fig4_out = os.path.join(out_dir, "Figure4_Synteny_and_Rearrangements.png")
plt.savefig(fig4_out, dpi=300)
plt.close()
print(f"Figure 4 saved to: {fig4_out}")
