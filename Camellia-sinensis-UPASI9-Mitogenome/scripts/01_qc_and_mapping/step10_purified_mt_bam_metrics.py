import os
import csv
import subprocess
from collections import Counter
from Bio import SeqIO

work_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
comb_bam = os.path.join(work_dir, "SRR9057080_vs_mt_cp_sorted.bam")
mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
mtpt_tsv = os.path.join(work_dir, "nc043914_vs_cp_mtpt.tsv")

# 1. Exact depth on NC_043914.1 AFTER competitive cpDNA removal!
mtpt_pos = set()
with open(mtpt_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        if int(row[3]) >= 50:
            qs, qe = min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))
            for p in range(qs, qe + 1):
                mtpt_pos.add(p)

ref_len = 707441
depths = []
non_mtpt_depths = []
c1, c10, c30, c100 = 0, 0, 0, 0
nc1, nc10, nc30, nc100 = 0, 0, 0, 0

proc = subprocess.Popen(["samtools", "depth", "-a", "-r", "NC_043914.1", comb_bam], stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    parts = line.split()
    if len(parts) >= 3:
        pos = int(parts[1])
        d = int(parts[2])
        depths.append(d)
        if d >= 1: c1 += 1
        if d >= 10: c10 += 1
        if d >= 30: c30 += 1
        if d >= 100: c100 += 1
        if pos not in mtpt_pos:
            non_mtpt_depths.append(d)
            if d >= 1: nc1 += 1
            if d >= 10: nc10 += 1
            if d >= 30: nc30 += 1
            if d >= 100: nc100 += 1
proc.wait()

depths.sort()
non_mtpt_depths.sort()
n_non = len(non_mtpt_depths)
print("=== PURIFIED COMPETITIVE MT COVERAGE ON NC_043914.1 (AFTER cpDNA REMOVAL) ===")
print(f"  All 707,441 bp: Mapped bases={sum(depths):,} bp | Mean={sum(depths)/len(depths):.2f}x | Median={depths[len(depths)//2]}x | Max={depths[-1]:,}x")
print(f"  Breadth (>=1x, >=10x, >=30x, >=100x): {c1/ref_len*100:.2f}%, {c10/ref_len*100:.2f}%, {c30/ref_len*100:.2f}%, {c100/ref_len*100:.2f}%")
print(f"  Non-MTPT Mitochondrial Core ({n_non:,} bp): Mapped bases={sum(non_mtpt_depths):,} bp | Mean={sum(non_mtpt_depths)/n_non:.2f}x | Median={non_mtpt_depths[n_non//2]}x | Max={non_mtpt_depths[-1]:,}x")
print(f"  Non-MTPT Breadth (>=1x, >=10x, >=30x, >=100x): {nc1/n_non*100:.2f}%, {nc10/n_non*100:.2f}%, {nc30/n_non*100:.2f}%, {nc100/n_non*100:.2f}%")
