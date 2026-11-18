import os
import csv
import gzip
import math
import subprocess
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

work_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"

# =========================================================================
# 1. FULL ORIGIN CLASSIFICATION OF ALL 1,160 SCAFFOLDS VS (ALL CAMELLIA MT vs CP)
# =========================================================================
scaf_fasta = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"
all_mt_comb = os.path.join(work_dir, "all_camellia_mt_panref.fasta")
pan_fastas = [
    "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta",
    os.path.join(ref_dir, "CS_Xinyang10_Chr1_798kb.fasta"),
    os.path.join(ref_dir, "CS_Xinyang10_Chr2_46kb.fasta"),
]
# Add sequences from Yidan and Wuyi GenBanks
for gb_fn in ["C_sinensis_Yidan_PX759560.gb", "C_sinensis_Wuyi_Chr1_PZ267020.gb", "C_sinensis_Wuyi_Chr2_PZ267021.gb"]:
    gb_p = os.path.join(ref_dir, gb_fn)
    if os.path.exists(gb_p) and os.path.getsize(gb_p) > 5000:
        pass
with open(all_mt_comb, "w") as out:
    for fn in pan_fastas:
        with open(fn) as inp:
            out.write(inp.read().strip() + "\n")
    for gb_fn in ["C_sinensis_Yidan_PX759560.gb", "C_sinensis_Wuyi_Chr1_PZ267020.gb", "C_sinensis_Wuyi_Chr2_PZ267021.gb"]:
        gb_p = os.path.join(ref_dir, gb_fn)
        if os.path.exists(gb_p) and os.path.getsize(gb_p) > 5000:
            for rec in SeqIO.parse(gb_p, "genbank"):
                SeqIO.write(rec, out, "fasta")

pan_mt_tsv = os.path.join(work_dir, "scaffolds_vs_pan_camellia_mt.tsv")
if not os.path.exists(pan_mt_tsv):
    subprocess.run(["makeblastdb", "-in", all_mt_comb, "-dbtype", "nucl", "-out", all_mt_comb], check=True)
    subprocess.run(["blastn", "-query", scaf_fasta, "-db", all_mt_comb, "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore", "-evalue", "1e-5", "-max_target_seqs", "5", "-out", pan_mt_tsv], check=True)

cp_tsv = "/mnt/d/Mitochondria/Q1_Investigation/Assembly_Evaluation/scaffolds_vs_cp.tsv"
def parse_best(tsv):
    hits = {}
    covs = defaultdict(list)
    with open(tsv) as f:
        for row in csv.reader(f, delimiter="\t"):
            qid = row[0]
            pid, alen, bs = float(row[2]), int(row[3]), float(row[11])
            qs, qe = min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))
            covs[qid].append((qs, qe))
            if qid not in hits or bs > hits[qid]["bs"]:
                hits[qid] = {"pid": pid, "alen": alen, "bs": bs, "sid": row[1]}
    return hits, covs

pmt_hits, _ = parse_best(pan_mt_tsv)
cp_hits, _ = parse_best(cp_tsv)

scaf_records = list(SeqIO.parse(scaf_fasta, "fasta"))
pan_cls = Counter()
pan_bp = Counter()
for r in scaf_records:
    L = len(r.seq)
    m = pmt_hits.get(r.id)
    c = cp_hits.get(r.id)
    if m and not c:
        oc = "PanMT_Specific"
    elif c and not m:
        oc = "Chloroplast_Specific"
    elif m and c:
        if c["bs"] > m["bs"] * 1.15 and c["pid"] > m["pid"]:
            oc = "Chloroplast_Dominant"
        elif m["bs"] > c["bs"] * 1.15 and m["pid"] > c["pid"]:
            oc = "PanMT_Dominant_MTPT"
        else:
            oc = "Shared_MTPT_Homolog"
    else:
        oc = "Unclassified"
    pan_cls[oc] += 1
    pan_bp[oc] += L

print("=== 1. PAN-CAMELLIA MITOGENOME VS CHLOROPLAST CLASSIFICATION OF 1,160 SCAFFOLDS ===")
for k, v in pan_cls.most_common():
    print(f"  {k:<24}: {v:4d} scaffolds | {pan_bp[k]:7,} bp ({pan_bp[k]/258681*100:.2f}%)")

# =========================================================================
# 2. DROUGHT STRESS (SRR9057080) VS CONTROL (SRR9057081) + TRF-1 (78 vs 79)
# =========================================================================
print("\n=== 2. DROUGHT STRESS VS CONTROL MITOCHONDRIAL TRANSCRIPTOME & RNA EDITING ===")
# Compute exact mapped mitochondrial reads & per-gene mean depth / TPM across PRJNA539837 4 samples!
samples = [
    ("SRR9057080", "UPASI-9_WaterStress"),
    ("SRR9057081", "UPASI-9_Control"),
    ("SRR9057078", "TRF-1_WaterStress"),
    ("SRR9057079", "TRF-1_Control"),
]

# Get gene coordinates for unique non-MTPT CDS genes from NC_043914.1
rec_gb = next(SeqIO.parse(mt_gb, "genbank"))
cds_coords = {}
seen = set()
for feat in rec_gb.features:
    if feat.type == "CDS":
        g = feat.qualifiers.get("gene", ["unknown"])[0]
        if g in seen:
            continue
        seen.add(g)
        coords = []
        for p in feat.location.parts:
            coords.extend(range(int(p.start) + 1, int(p.end) + 1))
        cds_coords[g] = (sorted(set(coords)), "+" if feat.location.strand == 1 else "-")

gene_depths = {g: {} for g in cds_coords}
total_clean_reads = {}
mito_mapped_reads = {}

for srr, label in samples:
    flag = f"/mnt/d/Mitochondria/analysis/{srr}/mapping/{srr}_flagstat.txt"
    bam = f"/mnt/d/Mitochondria/analysis/{srr}/mapping/{srr}_aligned.bam"
    with open(flag) as f:
        for line in f:
            if "primary" in line and "total" not in line and "mapped" not in line and "duplicates" not in line:
                total_clean_reads[srr] = int(line.split()[0])
            elif "primary mapped" in line:
                mito_mapped_reads[srr] = int(line.split()[0])
    # Stream samtools depth using indexed region jump (-r NC_043914.1) so unmapped tail is skipped!
    d_arr = [0] * (len(rec_gb.seq) + 1)
    proc = subprocess.Popen(["samtools", "depth", "-r", "NC_043914.1", bam], stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        parts = line.split()
        if len(parts) >= 3:
            d_arr[int(parts[1])] = int(parts[2])
    proc.wait()
    for g, (coords, st) in cds_coords.items():
        # RPM-normalized depth: (mean_depth / (total_clean_reads / 1e6))
        md = sum(d_arr[p] for p in coords) / len(coords)
        rpm_norm = md / (total_clean_reads[srr] / 1e6)
        gene_depths[g][srr] = (md, rpm_norm)

print("  Sample clean / mito read summary (PRJNA539837):")
for srr, label in samples:
    print(f"    {srr} ({label:<20}): clean={total_clean_reads[srr]:,}, mito={mito_mapped_reads[srr]:,} ({mito_mapped_reads[srr]/total_clean_reads[srr]*100:.2f}%)")

print("\n  Top stress-responsive mitochondrial CDS (UPASI-9 Stress SRR9057080 vs Control SRR9057081 normalized RPM depth):")
diff_rows = []
for g in sorted(cds_coords):
    d80, r80 = gene_depths[g]["SRR9057080"]
    d81, r81 = gene_depths[g]["SRR9057081"]
    d78, r78 = gene_depths[g]["SRR9057078"]
    d79, r79 = gene_depths[g]["SRR9057079"]
    fc_upasi = math.log2((r80 + 0.01) / (r81 + 0.01))
    fc_trf1 = math.log2((r78 + 0.01) / (r79 + 0.01))
    diff_rows.append({
        "Gene": g, "Strand": cds_coords[g][1],
        "UPASI9_Stress_Depth": round(d80, 2), "UPASI9_Ctrl_Depth": round(d81, 2),
        "UPASI9_Stress_RPMD": round(r80, 3), "UPASI9_Ctrl_RPMD": round(r81, 3), "UPASI9_Log2FC": round(fc_upasi, 3),
        "TRF1_Stress_RPMD": round(r78, 3), "TRF1_Ctrl_RPMD": round(r79, 3), "TRF1_Log2FC": round(fc_trf1, 3),
    })

diff_out = os.path.join(work_dir, "drought_stress_mito_expression_prjna539837.csv")
with open(diff_out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(diff_rows[0].keys()))
    w.writeheader()
    w.writerows(diff_rows)

for r in sorted(diff_rows, key=lambda x: abs(x["UPASI9_Log2FC"]), reverse=True)[:15]:
    print(f"    {r['Gene']:<10} UPASI9 Stress={r['UPASI9_Stress_RPMD']:6.2f} Ctrl={r['UPASI9_Ctrl_RPMD']:6.2f} log2FC={r['UPASI9_Log2FC']:+6.2f} | TRF1 log2FC={r['TRF1_Log2FC']:+6.2f}")

# Also compare RNA editing VAF at UPASI-9 high-confidence canonical CDS editing sites between Stress (80) and Control (81)!
ann_var = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_variants_feature_annotated.csv"
edit_sites = {}
with open(ann_var) as f:
    for row in csv.DictReader(f):
        if row["region"] == "CDS" and ((row["strand"] == "+" and row["ref"] == "C" and row["alt"] == "T") or (row["strand"] == "-" and row["ref"] == "G" and row["alt"] == "A")):
            edit_sites[int(row["pos"])] = (row["gene"], row["strand"], float(row["vaf"]), row["aa_change"])

vcf81 = "/mnt/d/Mitochondria/analysis/SRR9057081/mapping/SRR9057081_variants.vcf.gz"
vaf81 = {}
with gzip.open(vcf81, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        pos = int(p[1])
        if pos in edit_sites:
            info = dict(item.split("=") for item in p[7].split(";") if "=" in item)
            dp4 = [int(x) for x in info.get("DP4", "0,0,0,0").split(",")]
            tot = sum(dp4)
            if tot >= 20:
                vaf81[pos] = (dp4[2] + dp4[3]) / tot

shared_edits = [pos for pos in edit_sites if pos in vaf81]
deltas = [edit_sites[p][2] - vaf81[p] for p in shared_edits]
mean_v80 = sum(edit_sites[p][2] for p in shared_edits) / len(shared_edits) if shared_edits else 0
mean_v81 = sum(vaf81[p] for p in shared_edits) / len(shared_edits) if shared_edits else 0
inc_stress = sum(1 for d in deltas if d >= 0.05)
dec_stress = sum(1 for d in deltas if d <= -0.05)
print(f"\n  Canonical C-to-U RNA Editing Rates (UPASI-9 Stress SRR9057080 vs Control SRR9057081, n={len(shared_edits)} shared >=20x sites):")
print(f"    Mean Editing Rate: Stress={mean_v80*100:.2f}% vs Control={mean_v81*100:.2f}% (Delta={(mean_v80-mean_v81)*100:+.2f}%)")
print(f"    Sites shifted >=5% under stress: Up={inc_stress}, Down={dec_stress}, Stable={len(shared_edits)-inc_stress-dec_stress}")
