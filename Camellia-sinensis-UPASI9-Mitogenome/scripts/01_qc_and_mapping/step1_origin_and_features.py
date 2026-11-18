import os
import csv
import subprocess
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np

out_dir = "/mnt/d/Mitochondria/Q1_Investigation"
mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
cp_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_chloroplast_NC_020019.1.fasta"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
scaf_fasta = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
path_fasta = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/embplant_mt.K105.scaffolds.graph1.1.path_sequence.fasta"
bam_path = "/mnt/d/Mitochondria/analysis/SRR9057080/mapping/SRR9057080_aligned.bam"
var_csv = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_high_confidence_variants.csv"

os.makedirs("/mnt/d/Mitochondria/Q1_Investigation/Annotation", exist_ok=True)

# 1. Check/make BLAST DBs for mtDNA and cpDNA
cp_db = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/ref_db/cs_cp_ref"
mt_db = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/ref_db/cs_mt_ref"
if not os.path.exists(cp_db + ".nhr"):
    subprocess.run(["makeblastdb", "-in", cp_fasta, "-dbtype", "nucl", "-out", cp_db], check=True)

# BLAST 1,160 scaffolds against mtDNA and cpDNA
blast_mt = "/mnt/d/Mitochondria/Q1_Investigation/Assembly_Evaluation/scaffolds_vs_mt.tsv"
blast_cp = "/mnt/d/Mitochondria/Q1_Investigation/Assembly_Evaluation/scaffolds_vs_cp.tsv"
outfmt = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
subprocess.run(["blastn", "-query", scaf_fasta, "-db", mt_db, "-outfmt", outfmt, "-evalue", "1e-5", "-max_target_seqs", "5", "-out", blast_mt], check=True)
subprocess.run(["blastn", "-query", scaf_fasta, "-db", cp_db, "-outfmt", outfmt, "-evalue", "1e-5", "-max_target_seqs", "5", "-out", blast_cp], check=True)

def best_blast(tsv_path):
    hits = {}
    cov_intervals = defaultdict(list)
    with open(tsv_path) as f:
        for row in csv.reader(f, delimiter="\t"):
            qid, sid, pident, alen, qlen, slen, qs, qe, ss, se, ev, bs = row
            pident, bs = float(pident), float(bs)
            qs, qe = min(int(qs), int(qe)), max(int(qs), int(qe))
            cov_intervals[qid].append((qs, qe))
            if qid not in hits or bs > hits[qid]["bs"]:
                hits[qid] = {"pident": pident, "alen": int(alen), "qlen": int(qlen), "bs": bs, "ss": int(ss), "se": int(se)}
    # merge intervals for query coverage
    qcov = {}
    for qid, ivs in cov_intervals.items():
        ivs.sort()
        merged = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        qcov[qid] = sum(e - s + 1 for s, e in merged)
    return hits, qcov

mt_hits, mt_qcov = best_blast(blast_mt)
cp_hits, cp_qcov = best_blast(blast_cp)

records = list(SeqIO.parse(scaf_fasta, "fasta"))
class_counts = Counter()
class_bp = Counter()
scaf_origins = {}

origin_tsv = "/mnt/d/Mitochondria/Q1_Investigation/Assembly_Evaluation/scaffold_origin_classification.tsv"
with open(origin_tsv, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Scaffold_ID", "Length_bp", "Origin_Class", "MT_Bitscore", "MT_Pident", "MT_Cov_pct", "CP_Bitscore", "CP_Pident", "CP_Cov_pct"])
    for r in records:
        qid = r.id
        L = len(r.seq)
        m = mt_hits.get(qid)
        c = cp_hits.get(qid)
        m_bs = m["bs"] if m else 0.0
        c_bs = c["bs"] if c else 0.0
        m_pid = m["pident"] if m else 0.0
        c_pid = c["pident"] if c else 0.0
        m_cov = (mt_qcov.get(qid, 0) / L) * 100.0
        c_cov = (cp_qcov.get(qid, 0) / L) * 100.0

        if m and not c:
            oc = "Mitochondrial_Specific"
        elif c and not m:
            oc = "Chloroplast_Contamination"
        elif m and c:
            # Both hit: distinguish MTPT / shared vs dominant organelle
            if c_bs > m_bs * 1.15 and c_pid > m_pid:
                oc = "Chloroplast_Dominant"
            elif m_bs > c_bs * 1.15 and m_pid > c_pid:
                oc = "Mitochondrial_Dominant_MTPT"
            else:
                oc = "Shared_Homolog_MTPT"
        else:
            oc = "Unclassified_Nuclear_or_Novel"
        class_counts[oc] += 1
        class_bp[oc] += L
        scaf_origins[qid] = oc
        w.writerow([qid, L, oc, m_bs, m_pid, f"{m_cov:.1f}", c_bs, c_pid, f"{c_cov:.1f}"])

print("=== SCAFFOLD ORIGIN CLASSIFICATION (1,160 scaffolds, 258,681 bp) ===")
for k, v in class_counts.most_common():
    print(f"  {k:<32}: {v:4d} scaffolds | {class_bp[k]:7,} bp ({class_bp[k]/258681*100:.2f}%)")

print("\n=== CHECK KEY TOP SCAFFOLDS & tRNA/ORF SCAFFOLDS ===")
for key in ["NODE_1_", "NODE_2_", "NODE_3_", "NODE_4_", "NODE_5_", "NODE_6_", "NODE_7_", "NODE_8_", "NODE_9_", "NODE_10_", "NODE_14_", "NODE_15_", "NODE_17_", "NODE_23_", "NODE_61_", "NODE_89_"]:
    for r in records:
        if r.id.startswith(key):
            m = mt_hits.get(r.id)
            c = cp_hits.get(r.id)
            print(f"  {r.id[:32]:<32} len={len(r.seq):4d} | Class={scaf_origins[r.id]:<28} | MT=({m['pident'] if m else 0:.1f}%, {mt_qcov.get(r.id,0)}bp) | CP=({c['pident'] if c else 0:.1f}%, {cp_qcov.get(r.id,0)}bp)")

# Check path_sequence.fasta BLAST
if os.path.exists(path_fasta):
    subprocess.run(["blastn", "-query", path_fasta, "-db", mt_db, "-outfmt", outfmt, "-evalue", "1e-5", "-max_target_seqs", "2"], check=False)
    print("\n--- path_sequence.fasta vs cpDNA ---")
    subprocess.run(["blastn", "-query", path_fasta, "-db", cp_db, "-outfmt", outfmt, "-evalue", "1e-5", "-max_target_seqs", "2"], check=False)
