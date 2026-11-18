import os
import csv
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from collections import Counter
from Bio import SeqIO, PDB, Phylo

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Figures"
os.makedirs(out_dir, exist_ok=True)
sns.set_theme(style="ticks", font_scale=0.95)
plt.rcParams["font.family"] = "DejaVu Sans"

mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
clean_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"

# =========================================================================
# FIGURE 1: BUG-FIXED HIGH-RESOLUTION CIRCULAR MITOCHONDRIAL GENOME ATLAS
# =========================================================================
print("Generating Figure 1 (Fixed Exon-Level Circular Atlas)...")
rec_gb = next(SeqIO.parse(mt_gb, "genbank"))
ref_len = len(rec_gb.seq)
ref_seq_str = str(rec_gb.seq).upper()

# Extract individual EXON parts (NEVER span across trans-spliced gaps!)
exon_parts = []
gene_label_coords = {}
for feat in rec_gb.features:
    if feat.type in ("CDS", "rRNA"):
        name = feat.qualifiers.get("gene", feat.qualifiers.get("product", [feat.type]))[0]
        for p in feat.location.parts:
            s, e, st = int(p.start), int(p.end), (p.strand if p.strand is not None else 1)
            exon_parts.append((name, feat.type, s, e, st))
            if name not in gene_label_coords or (e - s) > (gene_label_coords[name][1] - gene_label_coords[name][0]):
                gene_label_coords[name] = (s, e, st)

# Load MTPT tracts (>=50bp)
mtpt_tracts = []
with open(os.path.join(clean_dir, "nc043914_vs_cp_mtpt.tsv")) as f:
    for row in csv.reader(f, delimiter="\t"):
        if int(row[3]) >= 50:
            mtpt_tracts.append((min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))))

# Load 2kb purified depth bins
with open(os.path.join(clean_dir, "nc043914_2kb_purified_depth_bins.json")) as f:
    depth_info = json.load(f)
bin_win = depth_info["window"]
bin_means = depth_info["means"]
bin_pos = [min(ref_len, i * bin_win + bin_win // 2) for i in range(len(bin_means))]

# Load canonical C-to-U CDS editing positions
edit_coords = []
with open("/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_variants_feature_annotated.csv") as f:
    for row in csv.DictReader(f):
        if row["region"] == "CDS" and ((row["strand"] == "+" and row["ref"] == "C" and row["alt"] == "T") or (row["strand"] == "-" and row["ref"] == "G" and row["alt"] == "A")):
            edit_coords.append(int(row["pos"]))

# GC 2kb sliding window
gc_pos, gc_vals = [], []
for i in range(0, ref_len - 2000, 2000):
    sub = ref_seq_str[i:i+2000]
    gc_pos.append(i + 1000)
    gc_vals.append((sub.count("G") + sub.count("C")) / 2000.0 * 100.0)

def gene_color(g, ftype):
    if ftype == "rRNA": return "#17becf"
    gl = g.lower()
    if gl.startswith("atp"): return "#1f77b4"
    if gl.startswith("nad"): return "#2ca02c"
    if gl.startswith("cox"): return "#d62728"
    if gl.startswith("cob"): return "#9467bd"
    if gl.startswith("ccm"): return "#ff7f0e"
    if gl.startswith(("rps", "rpl")): return "#e377c2"
    return "#8c564b"

fig = plt.figure(figsize=(11, 11))
ax = fig.add_subplot(111, polar=True)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)

# Backbone ring
th_all = np.linspace(0, 2 * np.pi, 1000)
ax.plot(th_all, [1.0] * len(th_all), color="black", lw=2.0)
ax.plot(th_all, [0.79] * len(th_all), color="#cccccc", lw=0.7)
ax.plot(th_all, [0.55] * len(th_all), color="#cccccc", lw=0.7)

# Draw Exon Parts ONLY (Forward r=0.90-0.97, Reverse r=0.81-0.88)
for name, ftype, s, e, st in exon_parts:
    t1 = (s / ref_len) * 2 * np.pi
    t2 = (e / ref_len) * 2 * np.pi
    # Make minimum visible angular width for tiny exons
    if abs(t2 - t1) < 0.008:
        mid = 0.5 * (t1 + t2)
        t1, t2 = mid - 0.004, mid + 0.004
    span = np.linspace(t1, t2, 12)
    c = gene_color(name, ftype)
    if st >= 0:
        ax.fill_between(span, 0.90, 0.97, color=c, alpha=0.95, linewidth=0)
    else:
        ax.fill_between(span, 0.81, 0.88, color=c, alpha=0.95, linewidth=0)

# Track: MTPT Tracts (gold) & Two Massive Repeats (purple) at r=0.73-0.78
for s, e in mtpt_tracts:
    t1, t2 = (s / ref_len) * 2 * np.pi, (e / ref_len) * 2 * np.pi
    if t2 - t1 < 0.008: t2 = t1 + 0.008
    ax.fill_between(np.linspace(t1, t2, 10), 0.73, 0.78, color="#e6ab02", alpha=0.95)
for s, e in [(1, 29795), (87015, 116814), (152144, 179092), (680512, 707441)]:
    t1, t2 = (s / ref_len) * 2 * np.pi, (e / ref_len) * 2 * np.pi
    ax.fill_between(np.linspace(t1, t2, 20), 0.70, 0.73, color="#6a3d9a", alpha=0.65)

# Track: Purified Competitive RNA-Seq Log10 Depth (r=0.56 to 0.69)
th_dep = (np.array(bin_pos) / ref_len) * 2 * np.pi
log_d = np.log10(np.array(bin_means) + 1.0)
norm_d = np.clip(log_d / 4.2, 0.0, 1.0)
r_dep = 0.56 + norm_d * 0.13
ax.fill_between(th_dep, 0.56, r_dep, color="#d62728", alpha=0.70)

# Track: GC Content (r=0.41 to 0.53)
th_gc = (np.array(gc_pos) / ref_len) * 2 * np.pi
norm_gc = np.clip((np.array(gc_vals) - 35.0) / 25.0, 0.0, 1.0)
ax.plot(th_gc, 0.41 + norm_gc * 0.12, color="#17becf", lw=1.0, alpha=0.9)

# Track: Canonical C-to-U Editing Ticks (r=0.35 to 0.39)
th_ed = (np.array(edit_coords) / ref_len) * 2 * np.pi
for te in th_ed:
    ax.plot([te, te], [0.35, 0.39], color="#2ca02c", lw=0.6, alpha=0.6)

# Outer gene labels cleanly de-overlapped around circle
placed_angles = []
for name, (s, e, st) in sorted(gene_label_coords.items(), key=lambda x: x[1][0]):
    if name.startswith("orf"): continue
    ang = ((0.5 * (s + e)) / ref_len) * 2 * np.pi
    if any(abs(ang - pa) < 0.095 for pa in placed_angles):
        continue
    placed_angles.append(ang)
    deg = np.degrees(ang)
    rot = 90 - deg if deg <= 180 else 270 - deg
    ha = "left" if deg <= 180 else "right"
    ax.text(ang, 1.04, name, fontsize=7.5, fontweight="bold", rotation=rot, rotation_mode="anchor", ha=ha, va="center")

# Clean 100 kb ticks (0 to 600 kb — avoids 700/0 overlap!)
tick_locs = np.arange(0, 700000, 100000)
ax.set_xticks((tick_locs / ref_len) * 2 * np.pi)
ax.set_xticklabels([f"{int(x/1000)} kb" for x in tick_locs], fontsize=9.5, fontweight="bold")
ax.tick_params(axis="x", pad=22)
ax.set_yticklabels([])
ax.set_ylim(0.15, 1.18)
ax.spines["polar"].set_visible(False)

handles = [
    mpatches.Patch(color="#1f77b4", label="ATP synthase (atp)"),
    mpatches.Patch(color="#2ca02c", label="NADH dehydr. (nad)"),
    mpatches.Patch(color="#d62728", label="Cytochrome oxidase (cox) / RNA-Seq Depth"),
    mpatches.Patch(color="#9467bd", label="Cytochrome b (cob)"),
    mpatches.Patch(color="#ff7f0e", label="Cyt c biogenesis (ccm) / MTPT Tracts"),
    mpatches.Patch(color="#e377c2", label="Ribosomal (rps/rpl)"),
    mpatches.Patch(color="#17becf", label="rRNAs & GC Content Track"),
    mpatches.Patch(color="#6a3d9a", label="Large Dispersed Repeats (29.8 / 27.0 kb)"),
]
ax.legend(handles=handles, loc="center", fontsize=7.5, frameon=True, framealpha=0.92)
plt.title("Camellia sinensis Reference Mitogenome (NC_043914.1, 707,441 bp)\nExon Architecture, MTPTs, Repeats & SRR9057080 Transcriptomic Tracks", fontsize=12, fontweight="bold", pad=24)
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure1_Complete_Circular_Mitochondrial_Genome.png"), dpi=300)
plt.close()
print("Saved Figure 1.")

# =========================================================================
# FIGURE 2: COMPETITIVE READ PARTITIONING & COVERAGE BREADTH DECONVOLUTION
# =========================================================================
print("Generating Figure 2 (Read Partition & Depth Deconvolution)...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))

# Panel A: Competitive Alignment Partition of 53.93M Clean Reads
parts = ["Nuclear / Other\nTranscripts\n(50.19M | 93.05%)", "Chloroplast DNA\n(NC_020019.1)\n(3.10M | 5.74%)", "Mitochondrial DNA\n(NC_043914.1)\n(0.65M | 1.21%)"]
vals_a = [50185808, 3097653, 652104]
cols_a = ["#8da0cb", "#66c2a5", "#fc8d62"]
bars_a = ax1.bar(parts, np.array(vals_a) / 1e6, color=cols_a, edgecolor="black", width=0.55)
ax1.set_yscale("log")
ax1.set_ylabel("Clean Paired-End Reads (Millions, Log10 Scale)", fontweight="bold")
ax1.set_title("A. Competitive Organellar Read Partitioning\n(53.93M Clean Reads vs Combined mtDNA + cpDNA Index)", fontweight="bold")
for b, v in zip(bars_a, vals_a):
    ax1.text(b.get_x() + b.get_width()/2, v/1e6 * 1.25, f"{v/1e6:.2f}M\n({v/53932552*100:.2f}%)", ha="center", va="bottom", fontsize=9, fontweight="bold")
ax1.set_ylim(0.2, 200)
ax1.grid(axis="y", linestyle="--", alpha=0.4)

# Panel B: Empirical Depth Deconvolution + Non-MTPT Breadth Thresholds
labels_b = [">=1x", ">=10x", ">=30x", ">=100x", "Mean Depth\n(Log10 x)"]
# Non-MTPT competitive breadths: 64.53%, 16.67%, 8.11%, 3.10%; Total Competitive: 65.01%, 17.65%, 9.18%, 4.14%
x = np.arange(4)
w = 0.35
naive_b = [65.09, 17.98, 9.42, 4.42]
pure_b = [64.53, 16.67, 8.11, 3.10]
b1 = ax2.bar(x - w/2, naive_b, w, label="Uncompetitive Single-Ref BAM (380.10x Mean)", color="#e78ac3", edgecolor="black")
b2 = ax2.bar(x + w/2, pure_b, w, label="Competitive Non-MTPT Core (28.52x Mean)", color="#1f78b4", edgecolor="black")
ax2.set_xticks(x)
ax2.set_xticklabels([">= 1x", ">= 10x", ">= 30x", ">= 100x"], fontweight="bold")
ax2.set_ylabel("Reference Mitogenome Breadth Covered (%)", fontweight="bold")
ax2.set_title("B. Mitogenome Coverage Breadth & Depth Recalculation\n(Eliminating 24,259.5x Artifact & 68.15% Plastid Cross-Mapping)", fontweight="bold")
ax2.set_ylim(0, 82)
for b in b1:
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 1.2, f"{b.get_height():.1f}%", ha="center", va="bottom", fontsize=8)
for b in b2:
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 1.2, f"{b.get_height():.1f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")
ax2.legend(loc="upper right", fontsize=8.5)
ax2.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure2_Read_Partition_and_Coverage_Breadth.png"), dpi=300)
plt.close()
print("Saved Figure 2.")

# =========================================================================
# FIGURE 4: VERIFIED SYNTENY (873 BLOCKS, 98.84% MEAN IDENTITY)
# =========================================================================
print("Generating Figure 4 (Verified NUCmer 873 Blocks)...")
coords_file = "/mnt/d/Mitochondria/analysis/SRR9057080/synteny/srr9057080_vs_ref.coords"
ref_s, ref_e, qry_s, qry_e, idys = [], [], [], [], []
with open(coords_file) as f:
    for line in f:
        parts = line.strip().split("|")
        if len(parts) >= 6:
            try:
                s1, e1 = map(int, parts[0].split())
                s2, e2 = map(int, parts[1].split())
                idy = float(parts[3].strip())
                ref_s.append(s1); ref_e.append(e1); qry_s.append(s2); qry_e.append(e2); idys.append(idy)
            except Exception:
                continue

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
for s1, e1, s2, e2 in zip(ref_s, ref_e, qry_s, qry_e):
    c = "#1f77b4" if (e1 - s1) * (e2 - s2) > 0 else "#d62728"
    ax1.plot([s1/1000.0, e1/1000.0], [s2/1000.0, e2/1000.0], color=c, lw=1.4, alpha=0.75)
ax1.set_title("A. Verified Syntenic Alignment Dot Plot (N = 873 Blocks)\n47.2% Forward Collinear (Blue) vs 52.8% Inverted (Red)", fontweight="bold", fontsize=10.5)
ax1.set_xlabel("Reference NC_043914.1 Coordinate (kb)", fontweight="bold")
ax1.set_ylabel("De Novo Scaffold Coordinate (kb)", fontweight="bold")
ax1.grid(True, linestyle="--", alpha=0.4)

ax2.hist(idys, bins=25, color="#2ca02c", edgecolor="black", alpha=0.85)
ax2.axvline(np.mean(idys), color="red", linestyle="--", lw=1.5, label=f"Mean Identity: {np.mean(idys):.2f}%")
ax2.axvline(np.median(idys), color="navy", linestyle=":", lw=1.5, label=f"Median Identity: {np.median(idys):.2f}%")
ax2.set_title(f"B. Sequence Identity Across All {len(idys)} Syntenic Blocks\n(60,735 bp Non-Redundant Reference Covered)", fontweight="bold", fontsize=10.5)
ax2.set_xlabel("Sequence Identity (%)", fontweight="bold")
ax2.set_ylabel("Number of Aligned Blocks", fontweight="bold")
ax2.legend(loc="upper left")
ax2.grid(True, linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure4_Synteny_and_Rearrangements.png"), dpi=300)
plt.close()
print("Saved Figure 4.")

# =========================================================================
# FIGURE 5: STRAND-SPECIFIC C-TO-U RNA EDITING & HETEROPLASMY (4 PANELS)
# =========================================================================
print("Generating Figure 5 (4-Panel Strand-Specific RNA Editing Proof)...")
ann_var = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_variants_feature_annotated.csv"
df_v = pd.read_csv(ann_var)
df_snp = df_v[df_v["is_snp"] == True].copy()
df_snp["pair"] = df_snp["ref"] + ">" + df_snp["alt"]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(13.5, 9.5))

# Panel A: VAF Spectrum with exact Ts/Tv annotations
sns.histplot(df_v["vaf"] * 100, bins=28, kde=True, color="#2ca02c", ax=ax1, edgecolor="black")
ax1.axvline(85, color="red", ls="--", lw=1.6, label="Fixed/Constitutive >=85% (N=496, Ts/Tv=6.64)")
ax1.axvline(10, color="orange", ls="--", lw=1.6, label="Intermediate VAF 10-85% (N=490, Ts/Tv=3.76)")
ax1.set_title("A. Curated Variant Allele Frequency Spectrum (N = 993)", fontweight="bold")
ax1.set_xlabel("Variant Allele Frequency (%)", fontweight="bold")
ax1.set_ylabel("Variant Count", fontweight="bold")
ax1.legend(loc="upper center", fontsize=8.5)
ax1.grid(True, ls="--", alpha=0.4)

# Panel B: Complete 12-Substitution Spectrum across 939 Curated SNPs
all_pairs = ["C>T", "G>A", "A>G", "T>C", "A>C", "C>A", "G>T", "T>G", "A>T", "T>A", "C>G", "G>C"]
pair_cnts = [int((df_snp["pair"] == p).sum()) for p in all_pairs]
bar_cols = ["#d62728" if p in ("C>T", "G>A") else ("#1f77b4" if p in ("A>G", "T>C") else "#999999") for p in all_pairs]
b_sub = ax2.bar(all_pairs, pair_cnts, color=bar_cols, edgecolor="black")
ax2.set_title("B. Single-Nucleotide Substitution Spectrum (939 SNPs, Ts/Tv = 4.91)\nC>T (340) + G>A (385) = 77.2% of All High-Confidence SNPs", fontweight="bold")
ax2.set_ylabel("SNP Count", fontweight="bold")
ax2.set_xlabel("Reference -> Variant Substitution Type", fontweight="bold")
for b in b_sub:
    ax2.text(b.get_x() + b.get_width()/2, b.get_height() + 5, f"{int(b.get_height())}", ha="center", va="bottom", fontsize=8, fontweight="bold")
ax2.set_ylim(0, 440)
ax2.grid(axis="y", ls="--", alpha=0.4)

# Panel C: Strand-Specific CDS Concordance!
df_cds = df_snp[df_snp["region"] == "CDS"]
p_ct = int(((df_cds["strand"] == "+") & (df_cds["pair"] == "C>T")).sum())
p_ga = int(((df_cds["strand"] == "+") & (df_cds["pair"] == "G>A")).sum())
m_ct = int(((df_cds["strand"] == "-") & (df_cds["pair"] == "C>T")).sum())
m_ga = int(((df_cds["strand"] == "-") & (df_cds["pair"] == "G>A")).sum())
x_c = np.arange(2)
ax3.bar(x_c - 0.18, [p_ct, m_ct], 0.35, label="Reference C>T Substitution", color="#e41a1c", edgecolor="black")
ax3.bar(x_c + 0.18, [p_ga, m_ga], 0.35, label="Reference G>A Substitution", color="#377eb8", edgecolor="black")
ax3.set_xticks(x_c)
ax3.set_xticklabels(["Forward (+)-Strand CDS\n(Sense C->U = C>T)", "Reverse (-)-Strand CDS\n(Sense C->U = G>A)"], fontweight="bold")
ax3.set_ylabel("CDS Variant Count", fontweight="bold")
ax3.set_title("C. Strand-Specific C-to-U RNA Editing Concordance\n(500 Concordant Edits vs Only 13 Discordant = 38.5x Enrichment)", fontweight="bold")
ax3.text(-0.18, p_ct + 6, str(p_ct), ha="center", fontweight="bold")
ax3.text(0.18, p_ga + 6, str(p_ga), ha="center", fontweight="bold")
ax3.text(0.82, m_ct + 6, str(m_ct), ha="center", fontweight="bold")
ax3.text(1.18, m_ga + 6, str(m_ga), ha="center", fontweight="bold")
ax3.set_ylim(0, 295)
ax3.legend(loc="upper center", fontsize=9)
ax3.grid(axis="y", ls="--", alpha=0.4)

# Panel D: Codon Position & Top Amino Acid Conversions (500 Canonical C-to-U CDS Edits)
labels_d = ["Pos 1\n(161; 32.2%)", "Pos 2\n(268; 53.6%)", "Pos 3\n(71; 14.2%)", "P->L\n(94)", "S->L\n(93)", "S->F\n(61)", "P->S\n(49)"]
vals_d = [161, 268, 71, 94, 93, 61, 49]
cols_d = ["#4daf4a", "#4daf4a", "#a6d854", "#984ea3", "#984ea3", "#984ea3", "#984ea3"]
b_d = ax4.bar(labels_d, vals_d, color=cols_d, edgecolor="black")
ax4.axvline(2.5, color="black", ls=":", lw=1.2)
ax4.set_ylabel("Canonical C-to-U Edit Count", fontweight="bold")
ax4.set_title("D. Codon Position Bias (85.8% 1st/2nd Pos) &\nTop Non-Synonymous Hydrophobic Rescues (83.6% NonSyn)", fontweight="bold")
for b in b_d:
    ax4.text(b.get_x() + b.get_width()/2, b.get_height() + 5, f"{int(b.get_height())}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
ax4.set_ylim(0, 310)
ax4.tick_params(axis="x", labelsize=8)
ax4.grid(axis="y", ls="--", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure5_Strand_Specific_RNA_Editing_and_Heteroplasmy.png"), dpi=300)
plt.close()
print("Saved Figure 5.")

# =========================================================================
# FIGURE 6: CODON USAGE, Ka/Ks EDITING INFLATION & PHYLOGENETIC TREES
# =========================================================================
print("Generating Figure 6 (Ka/Ks, RSCU & Phylogenetic Reconciliation)...")
df_kaks = pd.read_csv(os.path.join(clean_dir, "srr9057080_kaks_and_editing.csv")).sort_values(by="C_to_U_Sites", ascending=False)
df_rscu = pd.read_csv(os.path.join(clean_dir, "srr9057080_mito_rscu.csv"))

fig = plt.figure(figsize=(14, 9.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
ax_k = fig.add_subplot(gs[0, 0])
ax_r = fig.add_subplot(gs[0, 1])
ax_t1 = fig.add_subplot(gs[1, 0])
ax_t2 = fig.add_subplot(gs[1, 1])

# Panel A: Nd vs Sd & C->U Sites per Gene
xk = np.arange(len(df_kaks))
ax_k.bar(xk - 0.18, df_kaks["Nd"], 0.35, label="Non-Synonymous Diffs (Nd)", color="#d62728", edgecolor="black")
ax_k.bar(xk + 0.18, df_kaks["Sd"], 0.35, label="Synonymous Diffs (Sd)", color="#1f77b4", edgecolor="black")
ax_k.plot(xk, df_kaks["C_to_U_Sites"], color="black", marker="o", ms=4, lw=1.2, label="Canonical C->U Sites")
ax_k.set_xticks(xk)
ax_k.set_xticklabels(df_kaks["Gene"], rotation=45, ha="right", fontsize=8, fontweight="bold")
ax_k.set_ylabel("Substitutions vs DNA Reference", fontweight="bold")
ax_k.set_title("A. RNA-Editing-Driven Apparent Ka/Ks Inflation Across 19 MT Genes\n(Ks ~ 0 Proving DNA Conservation While C->U Editing Drives Nd)", fontweight="bold")
ax_k.legend(loc="upper right", fontsize=8)
ax_k.grid(axis="y", ls="--", alpha=0.4)

# Panel B: Preferred Codons RSCU (>1.0) vs Ending Nucleotide (A/U vs G/C)
df_pref = df_rscu[df_rscu["RSCU"] > 1.0].sort_values(by="RSCU", ascending=False)
cols_r = ["#2ca02c" if end else "#ff7f0e" for end in df_pref["Ends_with_AU"]]
ax_r.bar(range(len(df_pref)), df_pref["RSCU"], color=cols_r, edgecolor="black")
ax_r.set_xticks(range(len(df_pref)))
ax_r.set_xticklabels(df_pref["Codon"], rotation=60, fontsize=7.5, fontweight="bold")
ax_r.axhline(1.0, color="gray", ls="--")
ax_r.set_ylabel("RSCU Value", fontweight="bold")
ax_r.set_title("B. Codon Usage Bias Across 21 MT Genes (6,042 Codons, GC3s = 36.96%)\n26 of 28 Preferred Codons (92.9%, Green) Terminate in A/U", fontweight="bold")
ax_r.grid(axis="y", ls="--", alpha=0.4)

# Panel C & D: Clean custom rectangular phylogram drawer (no InnerX overlap, clean bootstrap badges!)
LABEL_MAP = {
    "C_drupifera_PZ378386": "C. drupifera (PZ378386)",
    "C_sinensis_Yidan_PX759560": "C. sinensis cv. Yidan (PX759560)",
    "C_sinensis_Wuyi_PZ267020": "C. sinensis cv. Wuyi (PZ267020)",
    "C_sinensis_RefSeq_NC043914": "C. sinensis RefSeq (NC_043914.1)",
    "C_sinensis_Xinyang10_PX398764": "C. sinensis cv. Xinyang10 (PX398764)",
    "C_pitardii_PZ501402": "C. pitardii (PZ501402)",
    "C_sinensis_UPASI9_cDNA_Raw": "C. sinensis UPASI-9 (Unmasked cDNA)",
    "C_sinensis_UPASI9_EditMasked": "C. sinensis UPASI-9 (Edit-Masked)",
    "Actinidia_rufa_NC086781": "Actinidia rufa (NC_086781.1)",
    "Vaccinium_macrocarpon_NC023338": "Vaccinium macrocarpon (NC_023338.1)",
}

def draw_clean_tree(tree, ax):
    # Extract confidence from Newick internal node names (e.g. 'Inner455.00' -> 55) if needed
    import re
    terms = tree.get_terminals()
    y_pos = {t: len(terms) - 1 - i for i, t in enumerate(terms)}
    x_pos = {}
    def get_x(clade, cur_x=0.0):
        bl = max(0.0, float(clade.branch_length or 0.0))
        nx = cur_x + bl
        x_pos[clade] = nx
        for ch in clade.clades:
            get_x(ch, nx)
    def get_y(clade):
        if clade.is_terminal():
            return y_pos[clade]
        ys = [get_y(ch) for ch in clade.clades]
        y_pos[clade] = 0.5 * (min(ys) + max(ys))
        return y_pos[clade]
    get_x(tree.root, 0.0)
    get_y(tree.root)
    max_x = max(x_pos.values())

    def plot_node(clade):
        cx, cy = x_pos[clade], y_pos[clade]
        if not clade.is_terminal():
            c_ys = [y_pos[ch] for ch in clade.clades]
            ax.plot([cx, cx], [min(c_ys), max(c_ys)], color="#222222", lw=1.5)
            # Parse bootstrap support
            bs = None
            if clade.confidence is not None:
                bs = int(round(float(clade.confidence)))
            elif clade.name and "Inner" in clade.name:
                m = re.search(r"Inner\d(\d{2,3})\.00", clade.name)
                if m: bs = int(m.group(1))
            if bs is not None and bs >= 40 and cx > 0.001:
                if abs(cx - 0.0294) < 0.0005:
                    ax.text(cx - max_x * 0.012, cy - 0.14, f"{bs}%", fontsize=7.0, fontweight="bold",
                            color="#1f78b4", ha="right", va="top",
                            bbox=dict(boxstyle="square,pad=0.06", fc="white", ec="none", alpha=0.88))
                else:
                    ax.text(cx - max_x * 0.015, cy + 0.14, f"{bs}%", fontsize=7.2, fontweight="bold",
                            color="#1f78b4", ha="right", va="bottom",
                            bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.88))
            for ch in clade.clades:
                hx = x_pos[ch]; hy = y_pos[ch]
                is_hi = "UPASI9" in (ch.name or "")
                ax.plot([cx, hx], [hy, hy], color="#d62728" if is_hi else "#222222", lw=2.2 if is_hi else 1.5)
                plot_node(ch)
        else:
            lbl = LABEL_MAP.get(clade.name, clade.name)
            is_hi = "UPASI-9" in lbl
            ax.text(cx + max_x * 0.02, cy, lbl, va="center", ha="left",
                    fontsize=8.2, fontweight="bold" if is_hi else "normal",
                    color="#d62728" if is_hi else "#111111")
    plot_node(tree.root)
    ax.set_xlim(-max_x * 0.03, max_x * 1.72)
    ax.set_ylim(-0.6, len(terms) - 0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.grid(axis="x", ls="--", alpha=0.35)

t_raw = Phylo.read(os.path.join(clean_dir, "phylogeny_unmasked_cdna.nwk"), "newick")
t_msk = Phylo.read(os.path.join(clean_dir, "phylogeny_masked_genomic.nwk"), "newick")
draw_clean_tree(t_raw, ax_t1)
ax_t1.set_title("C. Unmasked cDNA Phylogeny (8,089 bp, 10 Genes)\nRNA Editing Falsely Pushes UPASI-9 Outside Camellia (1.92% Dist)", fontweight="bold", fontsize=10)
ax_t1.set_xlabel("Substitutions / Site", fontweight="bold", fontsize=9)

draw_clean_tree(t_msk, ax_t2)
ax_t2.set_title("D. RNA-Editing-Masked Genomic Phylogeny (500 Bootstraps)\nUPASI-9 Correctly Nests Inside Camellia sinensis (0.025% Dist)", fontweight="bold", fontsize=10)
ax_t2.set_xlabel("Substitutions / Site", fontweight="bold", fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure6_Codon_Usage_KaKs_and_Phylogeny.png"), dpi=300)
plt.close()
print("Saved Figure 6.")

# =========================================================================
# FIGURE 7: 28-ACCESSION PAN-COHORT & DROUGHT STRESS DYNAMICS
# =========================================================================
print("Generating Figure 7 (28-Accession Cohort & 2x2 Drought Stress)...")
df_coh = pd.read_csv("/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables/Table_S8_Cohort_28_PanMitochondrial_Summary.csv").sort_values(by="Mapped_to_NC043914", ascending=False)
df_str = pd.read_csv(os.path.join(clean_dir, "drought_stress_mito_expression_prjna539837.csv")).sort_values(by="UPASI9_Log2FC", ascending=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
colors_c = ["#d62728" if r == "SRR9057080" else ("#ff7f0e" if bp == "PRJNA539837" else "#1f77b4") for r, bp in zip(df_coh["Run"], df_coh["BioProject"])]
ax1.bar(range(len(df_coh)), df_coh["Mapped_to_NC043914"] / 1e6, color=colors_c, edgecolor="black")
ax1.set_xticks(range(len(df_coh)))
ax1.set_xticklabels(df_coh["Run"], rotation=90, fontsize=7.5, fontweight="bold")
ax1.set_ylabel("Mapped Organellar Reads (Millions)", fontweight="bold")
ax1.set_title("A. Mitochondrial Yield Across 28 Camellia RNA-Seq Accessions\n(PRJNA539837 Deep Quartet in Red/Orange; 4,266 Non-MTPT SNPs)", fontweight="bold", fontsize=10.5)
ax1.grid(axis="y", ls="--", alpha=0.4)

# Panel B: Top & Bottom Stress-Responsive Mitochondrial Genes (UPASI-9 vs TRF-1)
sub_s = pd.concat([df_str.head(7), df_str.tail(5)])
xs = np.arange(len(sub_s))
ax2.barh(xs - 0.18, sub_s["UPASI9_Log2FC"], 0.35, label="UPASI-9 (Stress / Ctrl)", color="#d62728", edgecolor="black")
ax2.barh(xs + 0.18, sub_s["TRF1_Log2FC"], 0.35, label="TRF-1 (Stress / Ctrl)", color="#2ca02c", edgecolor="black")
ax2.set_yticks(xs)
ax2.set_yticklabels(sub_s["Gene"], fontweight="bold", fontsize=9)
ax2.axvline(0, color="black", lw=1.0)
ax2.set_xlabel("Normalized Expression Log2 Fold Change (Water Stress / Control)", fontweight="bold")
ax2.set_title("B. Water-Stress Mitochondrial Transcriptional Response\n(Coordinated Induction of atp9, rps7, nad4l, cox2, nad6)", fontweight="bold", fontsize=10.5)
ax2.legend(loc="lower left", fontsize=8.2, framealpha=0.95)
ax2.set_xlim(-1.05, 1.75)
ax2.grid(axis="x", ls="--", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure7_Pan_Cohort_and_Drought_Stress_Dynamics.png"), dpi=300)
plt.close()
print("Saved Figure 7.")

# =========================================================================
# FIGURE 8: COMPARATIVE CAMELLIA SIZES + MISSING REFSEQ GENE RECOVERY
# =========================================================================
print("Generating Figure 8 (Comparative Sizes & Recovered Missing Genes Heatmap)...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw={"width_ratios": [1.05, 1.0]})

comp_data = [
    ("C. drupifera (2 Chr)", 978.6),
    ("C. sinensis cv. Yidan", 931.6),
    ("C. oleifera (7 Chr)", 916.3),
    ("C. sinensis cv. Wuyi (2 Chr)", 910.2),
    ("C. pitardii", 895.3),
    ("C. sinensis cv. Xinyang10 (2 Chr)", 845.1),
    ("C. sinensis RefSeq (NC_043914)", 707.4),
    ("UPASI-9 Purified MT+MTPT RNA", 526.7),
    ("UPASI-9 Legacy GetOrganelle", 258.7),
]
df_cp = pd.DataFrame(comp_data, columns=["Name", "kb"])
cols_cp = ["#d62728" if "Purified" in n else ("#ff9896" if "Legacy" in n else ("#1f77b4" if "NC_043914" in n else "#440154")) for n in df_cp["Name"]]
bars_cp = ax1.barh(df_cp["Name"][::-1], df_cp["kb"][::-1], color=cols_cp[::-1], edgecolor="black", alpha=0.90)
ax1.set_xlabel("Genome / Transcribed Assembly Length (kb)", fontweight="bold")
ax1.set_title("A. Comparative Mitochondrial Sizes in Genus Camellia", fontweight="bold")
ax1.set_xlim(0, 1180)
for b in bars_cp:
    ax1.text(b.get_width() + 12, b.get_y() + b.get_height()/2, f"{b.get_width():.1f} kb", va="center", fontsize=8.5, fontweight="bold")
ax1.grid(axis="x", ls="--", alpha=0.4)

# Panel B: Core Genes Missing from NC_043914.1 RefSeq but Intact in Xinyang10 & Recovered in UPASI-9
genes_rec = ["cox1 (527 aa)", "cox3 (265 aa)", "atp8 (159 aa)", "nad7 (323 aa*)", "nad9 (190 aa)", "rps3 (563 aa)", "rps4 (276 aa)", "rpl16 (144 aa*)", "trnH-GUG", "trnE-UUC", "atp1 / cob / nad6 (21 RefSeq)"]
# Matrix: [NC_043914_RefSeq, Xinyang10_PX398764, UPASI9_SRR9057080] (1=Present/Recovered, 0=Missing/Unannotated)
mat = np.array([
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [0, 1, 1],
    [1, 1, 1],
])
sns.heatmap(mat, annot=[["Absent" if v==0 else "Intact" for v in r] for r in mat], fmt="", cmap="YlGn", cbar=False,
            xticklabels=["NC_043914.1\n(707 kb)", "Xinyang10\n(845 kb)", "UPASI-9\n(SRR9057080)"],
            yticklabels=genes_rec, linewidths=0.8, linecolor="black", ax=ax2)
ax2.tick_params(axis="x", labelsize=9)
ax2.set_title("B. Conserved Genes Missing in NC_043914.1\n(Previously Misclassified as 'Novel Class D ORFs')", fontweight="bold", fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure8_Comparative_Camellia_Mitogenomes_and_Gene_Recovery.png"), dpi=300)
plt.close()
print("Saved Figure 8.")

# =========================================================================
# FIGURE 9: GPU ESMFOLD WT vs EDITED ATP1 + RECOVERED COX1 3D PROTEOMICS
# =========================================================================
print("Generating Figure 9 (WT vs Edited ATP1 Structural Rescue + COX1)...")
wt_pdb = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NC_043914_WT_atp1_predicted.pdb"
ed_pdb = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NODE_17_length_1781_cov_62.065632_1_predicted.pdb"
cox1_pdb = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NODE_15_COX1_predicted.pdb"

def load_ca(pdb_p):
    parser = PDB.PDBParser(QUIET=True)
    st = parser.get_structure("S", pdb_p)
    cas = [a for a in st.get_atoms() if a.get_name() == "CA"]
    pl = [a.get_bfactor()*100.0 if a.get_bfactor() <= 1.0 else a.get_bfactor() for a in cas]
    return cas, np.array(pl)

wt_ca, wt_pl = load_ca(wt_pdb)
ed_ca, ed_pl = load_ca(ed_pdb)
_, cx_pl = load_ca(cox1_pdb)

sup = PDB.Superimposer()
sup.set_atoms(wt_ca, ed_ca)
sup.apply(ed_ca)
wt_xyz = np.array([a.get_coord() for a in wt_ca])
ed_xyz = np.array([a.get_coord() for a in ed_ca])

df_mut = pd.read_csv("/mnt/d/Mitochondria/Q1_Investigation/Protein_Structure/atp1_wt_vs_edited_structural_comparison.csv")

fig = plt.figure(figsize=(14, 9.2))
gs = fig.add_gridspec(2, 2)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0], projection="3d")
ax4 = fig.add_subplot(gs[1, 1])

# Panel A: Per-Residue pLDDT Rescue Trace (WT vs Edited)
res_idx = np.arange(1, 510)
ax1.plot(res_idx, wt_pl, color="#ff7f0e", lw=1.2, alpha=0.85, label=f"Unedited DNA WT ATP1 (Mean pLDDT={wt_pl.mean():.2f})")
ax1.plot(res_idx, ed_pl, color="#1f77b4", lw=1.5, label=f"RNA-Edited UPASI-9 ATP1 (Mean pLDDT={ed_pl.mean():.2f})")
for m in df_mut["Residue"]:
    ax1.axvline(m, color="red", ls=":", lw=0.9, alpha=0.6)
ax1.annotate("P472L (+22.0)\nP497L (+24.0)", xy=(472, 76), xytext=(360, 62),
             arrowprops=dict(arrowstyle="->", color="red"), fontsize=8.5, fontweight="bold", color="#b2182b")
ax1.axhline(90, color="green", ls="--", lw=1.0)
ax1.set_title("A. ESMFold C-Terminal Alpha-Helical Stability Rescue via RNA Editing", fontweight="bold")
ax1.set_xlabel("ATP1 Amino Acid Residue", fontweight="bold")
ax1.set_ylabel("pLDDT Confidence Score", fontweight="bold")
ax1.set_ylim(55, 102)
ax1.legend(loc="lower left", fontsize=8.5)
ax1.grid(True, ls="--", alpha=0.4)

# Panel B: Hydropathy Gain & Local pLDDT Improvement at 8 Sites
AA3_TO_1 = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}
xm = np.arange(len(df_mut))
lbls_m = [f"{AA3_TO_1.get(r['WT_AA'], r['WT_AA'][0])}{r['Residue']}{AA3_TO_1.get(r['Edited_AA'], r['Edited_AA'][0])}" for _, r in df_mut.iterrows()]
ax2.bar(xm - 0.18, df_mut["Delta_pLDDT"], 0.35, color="#2ca02c", edgecolor="black", label="Local pLDDT Gain (Points)")
ax2.bar(xm + 0.18, df_mut["Delta_KD"], 0.35, color="#9467bd", edgecolor="black", label="Hydropathy Shift (Delta KD)")
ax2.set_xticks(xm)
ax2.set_xticklabels(lbls_m, rotation=30, fontweight="bold", fontsize=8.5)
ax2.axhline(0, color="black", lw=0.8)
ax2.set_title("B. Hydrophobic Restoration (+25.9 KD Total) & Local pLDDT Gain", fontweight="bold")
ax2.set_ylabel("Magnitude of Change (Edited - WT)", fontweight="bold")
ax2.legend(loc="upper left", fontsize=8.5)
ax2.grid(axis="y", ls="--", alpha=0.4)

# Panel C: 3D Superposed C-alpha Backbone of WT (Orange) vs Edited (Blue) + Mutated Sites (Red Stars)
ax3.plot(wt_xyz[:, 0], wt_xyz[:, 1], wt_xyz[:, 2], color="#ff7f0e", lw=1.4, alpha=0.75, label="WT Reference Backbone")
ax3.plot(ed_xyz[:, 0], ed_xyz[:, 1], ed_xyz[:, 2], color="#1f77b4", lw=1.6, alpha=0.90, label="RNA-Edited UPASI-9")
m_idxs = df_mut["Residue"].values - 1
ax3.scatter(ed_xyz[m_idxs, 0], ed_xyz[m_idxs, 1], ed_xyz[m_idxs, 2], color="red", s=55, edgecolor="black", zorder=10, label="8 Edited Loci (347-497)")
ax3.set_title(f"C. 3D Backbone Superposition (Global RMSD = {sup.rms:.2f} \u00c5)", fontweight="bold")
ax3.set_xlabel("X (\u00c5)"); ax3.set_ylabel("Y (\u00c5)"); ax3.set_zlabel("Z (\u00c5)")
ax3.legend(loc="upper left", fontsize=7.5)

# Panel D: Newly Recovered COX1 (527 aa) Structural Confidence Profile
ax4.plot(np.arange(1, len(cx_pl) + 1), cx_pl, color="#8c564b", lw=1.4)
ax4.fill_between(np.arange(1, len(cx_pl) + 1), 0, cx_pl, color="#8c564b", alpha=0.18)
ax4.axhline(90, color="green", ls="--", lw=1.0, label="Crystal-Grade >= 90 (78.7% of residues)")
ax4.set_title(f"D. ESMFold 3D Confidence of Newly Recovered COX1 (527 aa)\n(Absent in NC_043914.1; Mean pLDDT = {cx_pl.mean():.2f}, Median = {np.median(cx_pl):.1f})", fontweight="bold")
ax4.set_xlabel("COX1 Amino Acid Residue", fontweight="bold")
ax4.set_ylabel("pLDDT Score", fontweight="bold")
ax4.set_ylim(30, 103)
ax4.legend(loc="lower right", fontsize=8.5)
ax4.grid(True, ls="--", alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, "Figure9_ATP1_WT_vs_Edited_3D_Structural_Remodeling.png"), dpi=300)
plt.close()
print("Saved Figure 9.")
print("All primary Q1 figures v2 generated cleanly!")
