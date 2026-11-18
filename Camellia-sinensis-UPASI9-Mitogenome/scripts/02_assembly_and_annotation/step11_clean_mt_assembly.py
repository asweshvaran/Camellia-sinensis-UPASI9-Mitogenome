import os
import csv
import subprocess
from collections import Counter
from Bio import SeqIO

work_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
os.makedirs(work_dir, exist_ok=True)
mt_ref = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
xy1 = "/mnt/d/Mitochondria/reference/complete_mitogenomes/CS_Xinyang10_Chr1_798kb.fasta"
xy2 = "/mnt/d/Mitochondria/reference/complete_mitogenomes/CS_Xinyang10_Chr2_46kb.fasta"
cp_ref = "/mnt/d/Mitochondria/reference/Camellia_sinensis_chloroplast_NC_020019.1.fasta"
fq1 = "/mnt/d/Mitochondria/analysis/SRR9057080/qc/SRR9057080_clean_1.fq.gz"
fq2 = "/mnt/d/Mitochondria/analysis/SRR9057080/qc/SRR9057080_clean_2.fq.gz"

# Build unified Xinyang10 + cpDNA competitive index (so cox1, cox3, atp8, nad7, nad9, rps3, rps4, rpl16, trnH, trnE are included!)
xy_cp_ref = os.path.join(work_dir, "xinyang10_plus_cp.fasta")
if not os.path.exists(xy_cp_ref + ".bwt"):
    with open(xy_cp_ref, "w") as out:
        for fn in [xy1, xy2, cp_ref]:
            with open(fn) as inp:
                out.write(inp.read().strip() + "\n")
    subprocess.run(["bwa", "index", xy_cp_ref], check=True)

# Fast competitive mapping keeping ONLY mapped reads (-F 4) so BAM is tiny (~150 MB) and sorts in ~15 seconds!
xy_cp_bam = os.path.join(work_dir, "SRR9057080_vs_xy_cp_mapped_only.bam")
if not os.path.exists(xy_cp_bam):
    print("Running fast competitive BWA-MEM against Xinyang10 (845 kb) + cpDNA (157 kb), keeping mapped reads (-F 4)...")
    cmd = f"bwa mem -t 24 {xy_cp_ref} {fq1} {fq2} | samtools view -@ 8 -F 4 -bS - | samtools sort -@ 8 -o {xy_cp_bam}"
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(["samtools", "index", xy_cp_bam], check=True)

idx = subprocess.check_output(["samtools", "idxstats", xy_cp_bam], text=True)
print("=== COMPETITIVE MAPPING AGAINST XINYANG10 + cpDNA ===")
print(idx)

# Extract read pairs hitting Xinyang10 Chr1/Chr2 (PX398764.1 / PX398765.1) competitively (NOT cpDNA!)
r1_pur = os.path.join(work_dir, "purified_xy_mt_R1.fq.gz")
r2_pur = os.path.join(work_dir, "purified_xy_mt_R2.fq.gz")
if not os.path.exists(r1_pur):
    print("Extracting purified mitochondrial read pairs (PX398764.1 & PX398765.1)...")
    # Filter primary alignments (-F 2308) on PX398764.1 & PX398765.1 -> sort by read name -> paired fastq
    cmd_ext = f"samtools view -@ 8 -u -F 2308 {xy_cp_bam} PX398764.1 PX398765.1 | samtools sort -n -@ 8 -m 512M - | samtools fastq -@ 4 -1 {r1_pur} -2 {r2_pur} -0 /dev/null -s /dev/null -n"
    subprocess.run(cmd_ext, shell=True, check=True)

# Run SPAdes --rna on purified mitochondrial reads!
spades_out = os.path.join(work_dir, "spades_rna_purified_mt")
tx_fasta = os.path.join(spades_out, "transcripts.fasta")
if not os.path.exists(tx_fasta):
    print("Running SPAdes --rna on purified mitochondrial reads...")
    subprocess.run(["spades.py", "--rna", "-t", "16", "-m", "10", "-1", r1_pur, "-2", r2_pur, "-o", spades_out], check=True)

# Evaluate purified assembly & BLAST against unified MT vs CP!
recs = list(SeqIO.parse(tx_fasta, "fasta"))
lens = sorted([len(r.seq) for r in recs], reverse=True)
tot_bp = sum(lens)
cs = 0; n50, l50 = 0, 0
for i, l in enumerate(lens):
    cs += l
    if cs >= tot_bp / 2.0:
        n50, l50 = l, i + 1; break
gc = sum(r.seq.upper().count("G") + r.seq.upper().count("C") for r in recs) / tot_bp * 100 if tot_bp else 0
print("\n=== PURIFIED MITOCHONDRIAL rnaSPAdes ASSEMBLY METRICS (SRR9057080) ===")
print(f"  Transcripts / Contigs : {len(recs):,}")
print(f"  Total Length          : {tot_bp:,} bp")
print(f"  Longest Transcript    : {lens[0]:,} bp")
print(f"  N50 / L50             : {n50:,} bp / L50={l50}")
print(f"  GC Content            : {gc:.2f}%")

# Origin check against Xinyang10 vs cpDNA
xy_db = os.path.join(work_dir, "Xinyang10_comb.fasta")
cp_db = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/ref_db/cs_cp_ref"
p_mt_tsv = os.path.join(work_dir, "purified_xy_vs_mt.tsv")
p_cp_tsv = os.path.join(work_dir, "purified_xy_vs_cp.tsv")
subprocess.run(["blastn", "-query", tx_fasta, "-db", xy_db, "-outfmt", "6 qseqid pident length qlen bitscore", "-evalue", "1e-5", "-max_target_seqs", "2", "-out", p_mt_tsv], check=True)
subprocess.run(["blastn", "-query", tx_fasta, "-db", cp_db, "-outfmt", "6 qseqid pident length qlen bitscore", "-evalue", "1e-5", "-max_target_seqs", "2", "-out", p_cp_tsv], check=True)

mt_h = {}
with open(p_mt_tsv) as f:
    for r in csv.reader(f, delimiter="\t"):
        if r[0] not in mt_h or float(r[4]) > float(mt_h[r[0]][4]): mt_h[r[0]] = r
cp_h = {}
with open(p_cp_tsv) as f:
    for r in csv.reader(f, delimiter="\t"):
        if r[0] not in cp_h or float(r[4]) > float(cp_h[r[0]][4]): cp_h[r[0]] = r

pure_mt, shared, cp_only, unclass = 0, 0, 0, 0
pure_bp, shared_bp, cp_bp, un_bp = 0, 0, 0, 0
for r in recs:
    m, c = mt_h.get(r.id), cp_h.get(r.id)
    L = len(r.seq)
    if m and not c: pure_mt += 1; pure_bp += L
    elif m and c:
        if float(c[4]) > float(m[4]) * 1.15: cp_only += 1; cp_bp += L
        else: shared += 1; shared_bp += L
    elif c and not m: cp_only += 1; cp_bp += L
    else: unclass += 1; un_bp += L

print(f"  Mitochondrial-Specific: {pure_mt:4d} contigs | {pure_bp:7,} bp ({pure_bp/tot_bp*100:.1f}%)")
print(f"  MTPT / Shared Homolog : {shared:4d} contigs | {shared_bp:7,} bp ({shared_bp/tot_bp*100:.1f}%)")
print(f"  Chloroplast-Dominant  : {cp_only:4d} contigs | {cp_bp:7,} bp ({cp_bp/tot_bp*100:.1f}%)")
print(f"  Unclassified          : {unclass:4d} contigs | {un_bp:7,} bp ({un_bp/tot_bp*100:.1f}%)")

# Run Prodigal + BLASTp vs Unified DB on purified transcripts to verify full-length genes recovered cleanly!
pur_faa = os.path.join(work_dir, "purified_mt_proteins.faa")
subprocess.run(["prodigal", "-i", tx_fasta, "-a", pur_faa, "-p", "meta", "-q"], check=True)
uni_db = os.path.join(work_dir, "unified_camellia_mt_cp_proteins.faa")
pur_blastp = os.path.join(work_dir, "purified_proteins_vs_unified.tsv")
subprocess.run(["blastp", "-query", pur_faa, "-db", uni_db, "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "1e-5", "-max_target_seqs", "1", "-out", pur_blastp], check=True)

mt_found = set()
cp_found = set()
with open(pur_blastp) as f:
    for r in csv.reader(f, delimiter="\t"):
        if float(r[2]) >= 80.0:
            parts = r[1].split("|")
            if parts[0] == "MT": mt_found.add(parts[2].lower())
            else: cp_found.add(parts[2].lower())

print(f"  Unique Mitochondrial CDS recovered ({len(mt_found)}): {sorted(mt_found)}")
print(f"  Co-assembled Chloroplast CDS ({len(cp_found)}): {sorted(cp_found)}")
