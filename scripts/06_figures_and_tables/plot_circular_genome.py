import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
import os

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Figures"
os.makedirs(out_dir, exist_ok=True)

# 1. Load reference genome annotations & lengths
ref_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
ref_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"

for r in SeqIO.parse(ref_fasta, "fasta"):
    ref_len = len(r.seq)
    ref_seq = str(r.seq)
    break

# Extract CDS positions
cds_features = []
for record in SeqIO.parse(ref_gb, "genbank"):
    for feat in record.features:
        if feat.type == "CDS":
            gene = feat.qualifiers.get("gene", ["CDS"])[0]
            start = int(feat.location.start)
            end = int(feat.location.end)
            strand = feat.location.strand
            cds_features.append((gene, start, end, strand))

print(f"Total Reference CDS Extracted: {len(cds_features)}")

# Calculate GC content in 2 kb windows
window_size = 2000
step = 1000
positions = []
gc_content = []
for i in range(0, ref_len - window_size, step):
    sub = ref_seq[i:i+window_size]
    gc = (sub.count('G') + sub.count('C') + sub.count('g') + sub.count('c')) / window_size * 100.0
    positions.append(i + window_size // 2)
    gc_content.append(gc)

# Set up Circular Plot (Polar Projection)
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, polar=True)
ax.set_theta_direction(-1) # Clockwise
ax.set_theta_offset(np.pi / 2.0) # Start from 12 o'clock

# Radius layers:
# Outer track (r = 1.00 - 1.05): Scale & Backbone
# Forward CDS track (r = 0.90 - 0.96)
# Reverse CDS track (r = 0.82 - 0.88)
# GC Content track (r = 0.65 - 0.78)
# Assembled Core Coverage (r = 0.50 - 0.60)

# Backbone
theta_all = np.linspace(0, 2*np.pi, 1000)
ax.plot(theta_all, [1.0]*len(theta_all), color='black', lw=2)

# Gene functional color mapping
def get_gene_color(name):
    if name.startswith('atp'): return '#1f77b4' # Blue
    elif name.startswith('nad'): return '#2ca02c' # Green
    elif name.startswith('cox'): return '#d62728' # Red
    elif name.startswith('cob'): return '#9467bd' # Purple
    elif name.startswith('ccm'): return '#ff7f0e' # Orange
    elif name.startswith('rps') or name.startswith('rpl'): return '#e377c2' # Pink
    elif name.startswith('mat'): return '#8c564b' # Brown
    else: return '#7f7f7f' # Grey

# Draw Genes
for gene, start, end, strand in cds_features:
    t_start = (start / ref_len) * 2 * np.pi
    t_end = (end / ref_len) * 2 * np.pi
    t_span = np.linspace(t_start, t_end, 20)
    color = get_gene_color(gene)
    
    strand_val = strand if strand is not None else 1
    if strand_val >= 0:
        ax.fill_between(t_span, 0.90, 0.96, color=color, alpha=0.9)
    else:
        ax.fill_between(t_span, 0.82, 0.88, color=color, alpha=0.9)

# Draw GC Content
theta_gc = (np.array(positions) / ref_len) * 2 * np.pi
norm_gc = (np.array(gc_content) - 40.0) / 20.0 # Map 40-60% to 0-1
r_gc = 0.65 + norm_gc * 0.12
ax.plot(theta_gc, r_gc, color='#17becf', lw=1.2, alpha=0.85)

# Ticks and labels around the circle (every 100 kb)
tick_locs = np.arange(0, ref_len, 100000)
tick_thetas = (tick_locs / ref_len) * 2 * np.pi
tick_labels = [f"{int(x/1000)} kb" for x in tick_locs]
ax.set_xticks(tick_thetas)
ax.set_xticklabels(tick_labels, fontsize=9, fontweight='bold')
ax.set_yticklabels([])
ax.set_ylim(0.40, 1.15)
ax.spines['polar'].set_visible(False)

# Add Title & Legend
plt.title("Camellia sinensis Complete Mitochondrial Genome Architecture (707,441 bp)\nReference NC_043914.1 & Expressed Core Gene Organization", fontsize=12, fontweight='bold', pad=25)

# Custom legend
custom_lines = [
    plt.Rectangle((0,0), 1, 1, color='#1f77b4'),
    plt.Rectangle((0,0), 1, 1, color='#2ca02c'),
    plt.Rectangle((0,0), 1, 1, color='#d62728'),
    plt.Rectangle((0,0), 1, 1, color='#9467bd'),
    plt.Rectangle((0,0), 1, 1, color='#ff7f0e'),
    plt.Rectangle((0,0), 1, 1, color='#e377c2'),
    plt.Rectangle((0,0), 1, 1, color='#17becf')
]
ax.legend(custom_lines, ['ATP synthase (atp)', 'NADH dehydrogenase (nad)', 'Cytochrome c oxidase (cox)', 'Cytochrome b (cob)', 'Cytochrome c biogenesis (ccm)', 'Ribosomal proteins (rps/rpl)', 'GC Content Track'],
          loc='center', bbox_to_anchor=(0.5, 0.5), fontsize=8.5, frameon=True)

plt.tight_layout()
fig_out = os.path.join(out_dir, "Figure1_Complete_Circular_Mitochondrial_Genome.png")
plt.savefig(fig_out, dpi=300)
plt.close()
print(f"Figure 1 saved to: {fig_out}")
