import os
import csv
import subprocess
from Bio import SeqIO

work_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
os.makedirs(work_dir, exist_ok=True)
scaf_fasta = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
xy_chr1 = "/mnt/d/Mitochondria/reference/complete_mitogenomes/CS_Xinyang10_Chr1_798kb.fasta"
xy_chr2 = "/mnt/d/Mitochondria/reference/complete_mitogenomes/CS_Xinyang10_Chr2_46kb.fasta"
xy_comb = os.path.join(work_dir, "Xinyang10_comb.fasta")

with open(xy_comb, "w") as out:
    for fn in [xy_chr1, xy_chr2]:
        with open(fn) as inp:
            out.write(inp.read().strip() + "\n")

subprocess.run(["makeblastdb", "-in", xy_comb, "-dbtype", "nucl", "-out", xy_comb], check=True)
out_tsv = os.path.join(work_dir, "scaffolds_vs_xinyang10.tsv")
subprocess.run(["blastn", "-query", scaf_fasta, "-db", xy_comb, "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore", "-evalue", "1e-5", "-max_target_seqs", "3", "-out", out_tsv], check=True)

hits = {}
with open(out_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        qid = row[0]
        if qid not in hits or float(row[11]) > float(hits[qid][11]):
            hits[qid] = row

print("=== CHECK KEY SCAFFOLDS AGAINST XINYANG10 (845 kb MITOGENOME) ===")
for prefix in ["NODE_1_", "NODE_2_", "NODE_3_", "NODE_4_", "NODE_5_", "NODE_6_", "NODE_7_", "NODE_8_", "NODE_9_", "NODE_10_", "NODE_14_", "NODE_15_", "NODE_17_", "NODE_61_", "NODE_89_"]:
    matched = [(k, v) for k, v in hits.items() if k.startswith(prefix)]
    if matched:
        k, v = matched[0]
        print(f"  {k[:32]:<32} -> {v[1]} | pid={v[2]}% len={v[3]} (qlen={v[4]})")
    else:
        print(f"  {prefix:<32} -> NO HIT IN XINYANG10")

# Also extract nuc sequence of NODE_10, NODE_15, NODE_8, NODE_61, NODE_89 length & GC
print("\n=== ALSO CHECK NC_043914.1 UNANNOTATED ORFs / tblastn of NODE_15 (527aa) & NODE_10 (563aa) ===")
faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
test_faa = os.path.join(work_dir, "class_cd_test.faa")
test_recs = [r for r in SeqIO.parse(faa_path, "fasta") if any(r.id.startswith(p) for p in ["NODE_2_length_3671_cov_47.563937_2", "NODE_8_", "NODE_10_", "NODE_15_", "NODE_61_", "NODE_89_"])]
SeqIO.write(test_recs, test_faa, "fasta")
mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
print("--- tblastn vs NC_043914.1 ---")
res1 = subprocess.check_output(["tblastn", "-query", test_faa, "-subject", mt_fasta, "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "1e-3"], text=True)
print(res1 or "  (none)")
print("--- tblastn vs Xinyang10 ---")
res2 = subprocess.check_output(["tblastn", "-query", test_faa, "-subject", xy_comb, "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "1e-3"], text=True)
print(res2 or "  (none)")
