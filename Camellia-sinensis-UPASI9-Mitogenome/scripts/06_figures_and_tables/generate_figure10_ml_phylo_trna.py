import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

exp_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis/Expanded_Phylogeny"
fig_dir = "/mnt/d/Mitochondria/Q1_Investigation/Figures"
out_png = os.path.join(fig_dir, "Figure10_Expanded_ML_Phylogeny_and_tRNA_Conservation.png")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": "#333333",
    "axes.linewidth": 0.9,
})

tree_masked = Phylo.read(os.path.join(exp_dir, "ML_12taxa_Masked_NT.treefile"), "newick")
tree_masked.root_with_outgroup("Arabidopsis_thaliana_NC037304")
tree_masked.ladderize(reverse=True)

DISPLAY_NAMES = {
    "C_sinensis_UPASI9_SRR9057080": "C. sinensis cv. UPASI-9 (SRR9057080, This Study)",
    "C_sinensis_RefSeq_NC043914": "C. sinensis RefSeq (NC_043914.1)",
    "C_sinensis_Wuyi_PZ267020": "C. sinensis cv. Wuyi Narcissus (PZ267020)",
    "C_sinensis_Xinyang10_PX398764": "C. sinensis cv. Xinyang10 (PX398764)",
    "C_sinensis_Yidan_PX759560": "C. sinensis cv. Yidan (PX759560)",
    "C_pitardii_PZ501402": "C. pitardii (PZ501402)",
    "C_drupifera_PZ378386": "C. drupifera (PZ378386)",
    "Actinidia_rufa_NC086781": "Actinidia rufa (NC_086781.1)",
    "Vaccinium_macrocarpon_NC023338": "Vaccinium macrocarpon (NC_023338.1)",
    "Carya_illinoinensis_PX519271": "Carya illinoinensis (PX519271.1; Chen et al., 2024)",
    "Fagus_sylvatica_MT446430": "Fagus sylvatica (MT446430.1)",
    "Arabidopsis_thaliana_NC037304": "Arabidopsis thaliana (NC_037304.1; Outgroup)",
}

ORDER_COLOR = {
    "C_sinensis_UPASI9_SRR9057080": "#C62828",
    "C_sinensis_RefSeq_NC043914": "#102C57",
    "C_sinensis_Wuyi_PZ267020": "#102C57",
    "C_sinensis_Xinyang10_PX398764": "#102C57",
    "C_sinensis_Yidan_PX759560": "#102C57",
    "C_pitardii_PZ501402": "#102C57",
    "C_drupifera_PZ378386": "#102C57",
    "Actinidia_rufa_NC086781": "#2E7D32",
    "Vaccinium_macrocarpon_NC023338": "#2E7D32",
    "Carya_illinoinensis_PX519271": "#B45309",
    "Fagus_sylvatica_MT446430": "#B45309",
    "Arabidopsis_thaliana_NC037304": "#455A64",
}

def get_cladogram_blended_coords(tree):
    terminals = tree.get_terminals()
    y_coords = {clade: i for i, clade in enumerate(terminals)}
    depths = tree.depths()
    unit_depths = tree.depths(unit_branch_lengths=True)
    max_d = max(depths.values()) or 1.0
    max_u = max(unit_depths.values()) or 1.0

    def set_y(clade):
        if clade.is_terminal():
            return y_coords[clade]
        ys = [set_y(c) for c in clade.clades]
        y_coords[clade] = sum(ys) / len(ys)
        return y_coords[clade]
    set_y(tree.root)

    # Blend 55% topological step + 45% sqrt(branch length) for crystal-clear branch separation
    x_coords = {}
    for clade in tree.find_clades():
        xd = (depths[clade] / max_d) ** 0.55
        xu = unit_depths[clade] / max_u
        x_coords[clade] = 0.50 * xu + 0.50 * xd
    return x_coords, y_coords

fig = plt.figure(figsize=(16.5, 10.2), dpi=300)
gs = fig.add_gridspec(2, 2, width_ratios=[1.42, 1.0], height_ratios=[1.05, 1.0], hspace=0.35, wspace=0.25)

# =====================================================================
# PANEL A: 12-Taxon IQ-TREE2 ML Phylogeny
# =====================================================================
axA = fig.add_subplot(gs[:, 0])
xc, yc = get_cladogram_blended_coords(tree_masked)

# Shaded order bands
axA.axhspan(-0.45, 6.45, color="#E8F0FE", alpha=0.52, zorder=0)
axA.axhspan(6.55, 8.45, color="#E8F5E9", alpha=0.52, zorder=0)
axA.axhspan(8.55, 10.45, color="#FEF3C7", alpha=0.58, zorder=0)
axA.axhspan(10.55, 11.45, color="#ECEFF1", alpha=0.65, zorder=0)

# Draw branches
for clade in tree_masked.find_clades(order="preorder"):
    x0, y0 = xc[clade], yc[clade]
    for child in clade.clades:
        x1, y1 = xc[child], yc[child]
        col = "#C62828" if child.name == "C_sinensis_UPASI9_SRR9057080" else "#263238"
        lw = 2.5 if child.name == "C_sinensis_UPASI9_SRR9057080" else 1.5
        axA.plot([x0, x0, x1], [y0, y1, y1], color=col, lw=lw, zorder=2)

# Label key internal node supports cleanly without collisions
for clade in tree_masked.get_nonterminals():
    if clade == tree_masked.root:
        continue
    sup = str(clade.name) if clade.name else ""
    if "/" in sup:
        sh, uf = sup.split("/")
        try:
            sh_v, uf_v = float(sh), float(uf)
            if uf_v >= 50 or sh_v >= 75:
                axA.text(xc[clade] - 0.015, yc[clade] - 0.14, f"{sh}/{uf}",
                         fontsize=7.3, fontweight="bold", color="#102C57",
                         ha="right", va="bottom",
                         bbox=dict(boxstyle="round,pad=0.14", fc="white", ec="#90A4AE", lw=0.55, alpha=0.94))
        except ValueError:
            pass

# Label tips
for tip in tree_masked.get_terminals():
    nm = tip.name
    lbl = DISPLAY_NAMES.get(nm, nm)
    col = ORDER_COLOR.get(nm, "#212529")
    fw = "bold" if nm in ("C_sinensis_UPASI9_SRR9057080", "Carya_illinoinensis_PX519271") else "normal"
    fs = 8.6 if nm in ("C_sinensis_UPASI9_SRR9057080", "Carya_illinoinensis_PX519271") else 8.2
    axA.scatter([xc[tip]], [yc[tip]], color=col, s=40 if fw == "bold" else 24, zorder=4)
    axA.text(xc[tip] + 0.022, yc[tip], lbl, fontsize=fs, fontweight=fw, color=col, va="center", ha="left")

# Order annotations cleanly positioned in dedicated right-hand strip
axA.text(2.14, 3.0, "ERICALES\nTheaceae\n(98/100)", fontsize=8.0, fontweight="bold", color="#102C57", ha="center", va="center", rotation=90)
axA.text(2.14, 7.5, "ERICALES\nOutgroups\n(99.8/100)", fontsize=7.6, fontweight="bold", color="#2E7D32", ha="center", va="center", rotation=90)
axA.text(2.14, 9.5, "FAGALES\n(97.5/99)", fontsize=7.8, fontweight="bold", color="#B45309", ha="center", va="center", rotation=90)
axA.text(2.14, 11.0, "OUTGROUP", fontsize=7.0, fontweight="bold", color="#455A64", ha="center", va="center", rotation=90)

# Place IQ-TREE2 ModelFinder Summary cleanly below y=11.7 so zero overlap with Arabidopsis!
info_txt = (
    "IQ-TREE2 ModelFinder & Multi-Matrix Concordance (5,832 bp / 1,944 aa; 8 Shared PCGs; 1,000 SH-aLRT / 1,000 UFBoot):\n"
    "• C-to-U Masked NT: TIM3+F+G4 (LogL = −11,164.69)  |  Unmasked NT: TIM3+F+I (LogL = −11,603.60)  →  +438.91 ΔLogL rescue\n"
    "• Translated AA (Li et al., 2024 protocol): HIVw+F+I (LogL = −7,325.49)  →  Ericales 98.6/100  |  Fagales (Carya+Fagus) 93.1/93"
)
axA.text(0.50, 0.015, info_txt, transform=axA.transAxes, fontsize=7.4, va="bottom", ha="center",
         bbox=dict(boxstyle="round,pad=0.45", fc="#F8FAFC", ec="#102C57", lw=0.95, alpha=0.98))

axA.set_xlim(-0.08, 2.26)
axA.set_ylim(-0.65, 13.65)
axA.invert_yaxis()
axA.set_yticks([])
axA.set_xticks([])
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.spines["left"].set_visible(False)
axA.set_xlabel("IQ-TREE2 Maximum Likelihood Topology (SH-aLRT / UFBoot Support Shown at Nodes)", fontsize=8.6, fontweight="bold", color="#102C57")
axA.set_title("A. 12-Taxon Multi-Order Angiosperm ML Phylogeny (C-to-U-Masked PCGs)", fontsize=10.6, fontweight="bold", color="#102C57", loc="left", pad=8)

# =====================================================================
# PANEL B: Pairwise Genetic Distance Collapse (Unmasked vs Masked)
# =====================================================================
axB = fig.add_subplot(gs[0, 1])
taxa_b = [
    ("C. sinensis Wuyi", 2.023, 0.017),
    ("C. sinensis RefSeq", 2.040, 0.034),
    ("C. pitardii", 2.040, 0.034),
    ("C. sinensis Xinyang10", 2.058, 0.051),
    ("C. drupifera", 2.058, 0.051),
    ("Actinidia rufa", 2.966, 1.063),
    ("Carya illinoinensis", 2.829, 1.646),
    ("Fagus sylvatica", 2.932, 1.698),
    ("Vaccinium macrocarpon", 3.035, 3.052),
    ("Arabidopsis thaliana", 4.853, 4.064),
]
labels_b = [t[0] for t in taxa_b]
unm_b = [t[1] for t in taxa_b]
msk_b = [t[2] for t in taxa_b]
y_pos = np.arange(len(labels_b))
h = 0.37

axB.barh(y_pos - h/2, unm_b, height=h, color="#D32F2F", alpha=0.86, label="Unmasked cDNA (LBA → Vaccinium 100/100)")
axB.barh(y_pos + h/2, msk_b, height=h, color="#1565C0", alpha=0.92, label="C-to-U Masked (Rescued → Camellia 98/100)")

for i, (u, m) in enumerate(zip(unm_b, msk_b)):
    axB.text(u + 0.07, i - h/2, f"{u:.2f}%", va="center", fontsize=6.9, color="#B71C1C", fontweight="bold")
    axB.text(m + 0.07, i + h/2, f"{m:.2f}%", va="center", fontsize=6.9, color="#0D47A1", fontweight="bold")

axB.set_yticks(y_pos)
axB.set_yticklabels(labels_b, fontsize=7.9)
axB.invert_yaxis()
axB.set_xlim(0, 5.85)
axB.set_xlabel("Pairwise Nucleotide p-Distance vs UPASI-9 (%)", fontsize=8.6, fontweight="bold", color="#102C57")
axB.set_title("B. RNA-Editing Long-Branch Artifact & Masking Rescue", fontsize=10.3, fontweight="bold", color="#102C57", loc="left", pad=6)
axB.legend(loc="upper right", fontsize=7.2, frameon=True, facecolor="white", edgecolor="#CFD8DC", framealpha=0.96)
axB.grid(axis="x", linestyle="--", alpha=0.35)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)

# =====================================================================
# PANEL C: Conserved tRNA T-Psi-C Loop (U-U-C-x-A-x2) & AT/GC Skew
# =====================================================================
axC = fig.add_subplot(gs[1, 1])
trna_data = [
    ("Carya illinoinensis\n(PX519271; 24 tRNAs)", 100.0, -0.0069, +0.0012, "#B45309"),
    ("Arabidopsis thaliana\n(NC_037304; 22 tRNAs)", 100.0, +0.0059, -0.0019, "#546E7A"),
    ("Actinidia rufa\n(NC_086781; 37 tRNAs)", 97.3, +0.0042, +0.0044, "#2E7D32"),
    ("Fagus sylvatica\n(MT446430; 19 tRNAs)", 94.7, -0.0044, -0.0091, "#D97706"),
    ("C. sinensis Yidan\n(PX759560; 33 tRNAs)", 93.9, +0.0015, +0.0061, "#1E40AF"),
    ("C. sinensis Xinyang10\n(PX398764; 29 tRNAs)", 93.1, +0.0031, -0.0004, "#1565C0"),
    ("C. sinensis RefSeq\n(NC_043914; 23 tRNAs)", 91.3, -0.0021, -0.0017, "#102C57"),
]
x_c = np.arange(len(trna_data))
pcts_c = [d[1] for d in trna_data]
cols_c = [d[4] for d in trna_data]
bars = axC.bar(x_c, pcts_c, width=0.68, color=cols_c, edgecolor="#222222", linewidth=0.6, alpha=0.92)

for bar, d in zip(bars, trna_data):
    h_val = bar.get_height()
    axC.text(bar.get_x() + bar.get_width()/2, h_val + 0.55, f"{h_val:.1f}%", ha="center", va="bottom", fontsize=7.6, fontweight="bold", color="#102C57")
    axC.text(bar.get_x() + bar.get_width()/2, 81.2, f"AT {d[2]:+.3f}\nGC {d[3]:+.3f}", ha="center", va="bottom", fontsize=5.8, color="white", fontweight="bold")

axC.set_ylim(80, 109.5)
axC.set_xticks(x_c)
axC.set_xticklabels([d[0] for d in trna_data], fontsize=6.9, rotation=20, ha="right")
axC.set_ylabel("tRNAs Conserving U-U-C-x-A-x2 Motif (%)", fontsize=8.4, fontweight="bold", color="#102C57")
axC.set_title("C. Conserved tRNA Ψ-Loop Motif (Chen et al., 2024) & Skew", fontsize=10.3, fontweight="bold", color="#102C57", loc="left", pad=6)
axC.grid(axis="y", linestyle="--", alpha=0.35)
axC.spines["top"].set_visible(False)
axC.spines["right"].set_visible(False)

up_trna_note = (
    "UPASI-9 Recovered MT tRNAs (Absent in NC_043914.1):\n"
    "• trnH-GUG (74 nt): UUCAAUC at pos 55 (100% conserved)\n"
    "• trnE-UUC (72 nt): UUCGAUU at pos 53 (100% conserved)"
)
axC.text(0.02, 0.97, up_trna_note, transform=axC.transAxes, fontsize=6.9, va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.35", fc="#FFF8E1", ec="#C62828", lw=0.85, alpha=0.96))

fig.suptitle("Expanded Multi-Order Angiosperm Maximum Likelihood Phylogeny (IQ-TREE2) and Conserved tRNA Ψ-Loop Architecture",
             fontsize=12.6, fontweight="bold", color="#102C57", y=0.985)

plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved refined Figure 10 to {out_png} ({os.path.getsize(out_png):,} bytes)")
