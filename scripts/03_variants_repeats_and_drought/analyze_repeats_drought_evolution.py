import os
import csv
import gzip
import math
from collections import defaultdict, Counter
from Bio import SeqIO

gb_file = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
rec = next(SeqIO.parse(gb_file, "genbank"))
seq_len = len(rec.seq)

repeats_tsv = "/mnt/d/Mitochondria/Q1_Investigation/Repeats_SV/nc043914_self_repeats.tsv"
repeat_bins = {"20-50bp": [], "51-100bp": [], "101-500bp": [], "501-1000bp": [], ">1000bp": []}
all_repeats = []

if os.path.exists(repeats_tsv):
    with open(repeats_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                pident = float(parts[0])
                length = int(parts[1])
                s1, e1 = int(parts[2]), int(parts[3])
                s2, e2 = int(parts[4]), int(parts[5])
                if length == seq_len:  # self full match
                    continue
                all_repeats.append((length, pident, s1, e1, s2, e2))
                if 20 <= length <= 50:
                    repeat_bins["20-50bp"].append(length)
                elif 51 <= length <= 100:
                    repeat_bins["51-100bp"].append(length)
                elif 101 <= length <= 500:
                    repeat_bins["101-500bp"].append(length)
                elif 501 <= length <= 1000:
                    repeat_bins["501-1000bp"].append(length)
                elif length > 1000:
                    repeat_bins[">1000bp"].append(length)

print("="*70)
print(f"1. REPEAT DISTRIBUTION IN CAMELLIA SINENSIS MITOGENOME (NC_043914.1, {seq_len:,} bp)")
print("="*70)
print(f"Total self-repeat pairs detected: {len(all_repeats)}")
total_rep_bp = sum(sum(lens) for lens in repeat_bins.values())
for b, lens in repeat_bins.items():
    cnt = len(lens)
    sbp = sum(lens)
    pct_cnt = (cnt / len(all_repeats) * 100) if all_repeats else 0
    pct_bp = (sbp / total_rep_bp * 100) if total_rep_bp else 0
    avg = (sbp / cnt) if cnt else 0
    print(f"  Bin {b:<12}: {cnt:4d} repeats ({pct_cnt:5.1f}%) | Total bp = {sbp:7,d} ({pct_bp:5.1f}%) | Avg = {avg:6.1f} bp")

# Now check large repeats (>1kb)
large_reps = [r for r in all_repeats if r[0] > 1000]
print("\nLarge repeats (>1,000 bp):")
for r in large_reps:
    print(f"  Length: {r[0]:,} bp | Identity: {r[1]:.2f}% | Region 1: {r[2]:,}-{r[3]:,} | Region 2: {r[4]:,}-{r[5]:,}")

# 2. Check RNA editing shifts between Drought Stress (SRR9057080) and Control (SRR9057081)
ann_var = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_variants_feature_annotated.csv"
edit_sites = {}
with open(ann_var) as f:
    for row in csv.DictReader(f):
        if row["region"] == "CDS" and ((row["strand"] == "+" and row["ref"] == "C" and row["alt"] == "T") or (row["strand"] == "-" and row["ref"] == "G" and row["alt"] == "A")):
            edit_sites[int(row["pos"])] = (row["gene"], row["strand"], float(row["vaf"]), row["aa_change"], row["ref"], row["alt"])

vcf81 = "/mnt/d/Mitochondria/analysis/SRR9057081/mapping/SRR9057081_variants.vcf.gz"
vaf81 = {}
if os.path.exists(vcf81):
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

shared = [pos for pos in edit_sites if pos in vaf81]
print("\n" + "="*70)
print(f"2. RNA EDITING DYNAMICS UNDER DROUGHT STRESS (SRR9057080 Stress vs SRR9057081 Control)")
print("="*70)
print(f"Shared high-confidence editing sites (>=20x in both): n = {len(shared)}")
if shared:
    v80_list = [edit_sites[p][2] for p in shared]
    v81_list = [vaf81[p] for p in shared]
    mean_80 = sum(v80_list) / len(v80_list)
    mean_81 = sum(v81_list) / len(v81_list)
    deltas = [edit_sites[p][2] - vaf81[p] for p in shared]
    up = [p for p in shared if (edit_sites[p][2] - vaf81[p]) >= 0.05]
    down = [p for p in shared if (edit_sites[p][2] - vaf81[p]) <= -0.05]
    stable = [p for p in shared if abs(edit_sites[p][2] - vaf81[p]) < 0.05]
    print(f"  Mean Editing Rate: Drought Stress = {mean_80*100:.2f}% | Control = {mean_81*100:.2f}% | Delta = {(mean_80-mean_81)*100:+.2f}%")
    print(f"  Sites significantly Up-regulated (>= +5% VAF): {len(up)} sites ({len(up)/len(shared)*100:.1f}%)")
    print(f"  Sites significantly Down-regulated (<= -5% VAF): {len(down)} sites ({len(down)/len(shared)*100:.1f}%)")
    print(f"  Sites stable (< 5% shift): {len(stable)} sites ({len(stable)/len(shared)*100:.1f}%)")

    print("\n  Top Up-regulated Editing Sites Under Drought Stress in Core OXPHOS Subunits:")
    sorted_up = sorted(up, key=lambda p: edit_sites[p][2] - vaf81[p], reverse=True)
    for p in sorted_up[:10]:
        g, st, v80, aach, ref, alt = edit_sites[p]
        v81 = vaf81[p]
        print(f"    Pos {p:6d} | Gene: {g:<8} ({st}) | Stress VAF: {v80*100:5.1f}% | Ctrl VAF: {v81*100:5.1f}% | Shift: {(v80-v81)*100:+5.1f}% | Effect: {aach}")

    print("\n  Top Down-regulated Editing Sites Under Drought Stress:")
    sorted_down = sorted(down, key=lambda p: edit_sites[p][2] - vaf81[p])
    for p in sorted_down[:10]:
        g, st, v80, aach, ref, alt = edit_sites[p]
        v81 = vaf81[p]
        print(f"    Pos {p:6d} | Gene: {g:<8} ({st}) | Stress VAF: {v80*100:5.1f}% | Ctrl VAF: {v81*100:5.1f}% | Shift: {(v80-v81)*100:+5.1f}% | Effect: {aach}")

# 3. Check NODE_2_ORF2 chimeric candidate
print("\n" + "="*70)
print("3. CHIMERIC ORF CANDIDATE (NODE_2_ORF2) CONTEXT & REPEAT JUNCTIONS")
print("="*70)
orf_cur_file = "/mnt/d/Mitochondria/Q1_Investigation/ORF_Validation/orf_curation_summary.csv"
if os.path.exists(orf_cur_file):
    with open(orf_cur_file) as f:
        for r in csv.DictReader(f):
            if r.get("orf_id") == "NODE_2_ORF2" or "NODE_2" in r.get("orf_id", ""):
                print(r)
