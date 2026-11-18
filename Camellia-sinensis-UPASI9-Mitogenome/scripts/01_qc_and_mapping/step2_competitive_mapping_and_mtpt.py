import os
import subprocess
import csv
from collections import Counter, defaultdict
from Bio import SeqIO

work_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
os.makedirs(work_dir, exist_ok=True)

mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
cp_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_chloroplast_NC_020019.1.fasta"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
fq1 = "/mnt/d/Mitochondria/analysis/SRR9057080/qc/SRR9057080_clean_1.fq.gz"
fq2 = "/mnt/d/Mitochondria/analysis/SRR9057080/qc/SRR9057080_clean_2.fq.gz"

# 1. Identify exact MTPT (plastid-derived) tracts inside NC_043914.1
mtpt_tsv = os.path.join(work_dir, "nc043914_vs_cp_mtpt.tsv")
cp_db = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/ref_db/cs_cp_ref"
subprocess.run([
    "blastn", "-query", mt_fasta, "-db", cp_db,
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
    "-evalue", "1e-5", "-out", mtpt_tsv
], check=True)

mtpt_ivs = []
mtpt_hits = []
with open(mtpt_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        qs, qe = min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))
        alen = int(row[3])
        pident = float(row[2])
        if alen >= 50:
            mtpt_ivs.append((qs, qe))
            mtpt_hits.append((qs, qe, alen, pident, int(row[8]), int(row[9])))

def merge_iv(ivs):
    if not ivs: return []
    ivs.sort()
    m = [list(ivs[0])]
    for s, e in ivs[1:]:
        if s <= m[-1][1] + 1:
            m[-1][1] = max(m[-1][1], e)
        else:
            m.append([s, e])
    return m

merged_mtpt = merge_iv(mtpt_ivs)
mtpt_bp = sum(e - s + 1 for s, e in merged_mtpt)
print(f"=== MTPT REGIONS IN NC_043914.1 (BLASTn vs NC_020019.1, len >= 50 bp) ===")
print(f"  Merged MTPT tracts: {len(merged_mtpt)}, Total MTPT bases: {mtpt_bp:,} bp ({mtpt_bp/707441*100:.2f}% of mitogenome)")
for qs, qe, alen, pid, ss, se in sorted(mtpt_hits, key=lambda x: x[2], reverse=True)[:15]:
    print(f"    MT:{qs}-{qe} ({alen} bp, {pid:.1f}%) <-> CP:{ss}-{se}")

# 2. Build combined MT + CP reference and run competitive BWA-MEM
comb_ref = os.path.join(work_dir, "combined_mt_cp.fasta")
if not os.path.exists(comb_ref):
    with open(comb_ref, "w") as out:
        for fn in [mt_fasta, cp_fasta]:
            with open(fn) as inp:
                out.write(inp.read().strip() + "\n")
    subprocess.run(["bwa", "index", comb_ref], check=True)

comb_bam = os.path.join(work_dir, "SRR9057080_vs_mt_cp_sorted.bam")
if not os.path.exists(comb_bam):
    print("Running competitive BWA-MEM against combined mtDNA + cpDNA reference...")
    cmd = f"bwa mem -t 24 {comb_ref} {fq1} {fq2} | samtools view -@ 8 -bS - | samtools sort -@ 8 -o {comb_bam}"
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", comb_bam], check=True)

# 3. Parse competitive BAM stats
print("=== COMPETITIVE ALIGNMENT FLAGSTAT & IDXSTATS ===")
idx = subprocess.check_output(["samtools", "idxstats", comb_bam], text=True)
print(idx)
flag = subprocess.check_output(["samtools", "flagstat", comb_bam], text=True)
print(flag)
