import os
import csv
import subprocess
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"
out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
os.makedirs(out_dir, exist_ok=True)

gb_files = [
    ("NC_043914_RefSeq", "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"),
    ("Xinyang10_Chr1", os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb")),
    ("Xinyang10_Chr2", os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb")),
    ("Yidan", os.path.join(ref_dir, "C_sinensis_Yidan_PX759560.gb")),
    ("Wuyi_Chr1", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr1_PZ267020.gb")),
    ("Wuyi_Chr2", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb")),
    ("C_pitardii", os.path.join(ref_dir, "C_pitardii_PZ501402.gb")),
    ("C_drupifera_Chr1", os.path.join(ref_dir, "C_drupifera_Chr1_PZ378386.gb")),
    ("C_drupifera_Chr2", os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb")),
    ("cpDNA_NC_020019", "/mnt/d/Mitochondria/reference/Camellia_sinensis_chloroplast_NC_020019.1.gb"),
]

unified_faa = os.path.join(out_dir, "unified_camellia_mt_cp_proteins.faa")
recs = []
gene_catalog = defaultdict(set)

for label, gbp in gb_files:
    if not os.path.exists(gbp) or os.path.getsize(gbp) < 1000:
        print(f"Skipping {label}: missing/small")
        continue
    cnt_cds = 0
    for rec in SeqIO.parse(gbp, "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":
                g = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["unknown"]))[0]
                p = feat.qualifiers.get("product", ["unknown"])[0]
                t = feat.qualifiers.get("translation", [None])[0]
                if t:
                    cnt_cds += 1
                    gene_catalog[label].add(g)
                    org_tag = "CP" if label.startswith("cpDNA") else "MT"
                    recs.append(SeqRecord(Seq(t), id=f"{org_tag}|{label}|{g}", description=p))
    print(f"  Parsed {label:<20}: {cnt_cds:3d} CDS ({len(gene_catalog[label])} unique gene names)")

SeqIO.write(recs, unified_faa, "fasta")
subprocess.run(["makeblastdb", "-in", unified_faa, "-dbtype", "prot", "-out", unified_faa], check=True)

# Compare genes across NC_043914 vs Xinyang10 vs Yidan
mt_xy = gene_catalog["Xinyang10_Chr1"] | gene_catalog["Xinyang10_Chr2"]
mt_ref = gene_catalog["NC_043914_RefSeq"]
print("\n=== GENES PRESENT IN XINYANG10 MITOGENOME BUT UNANNOTATED/MISSING IN NC_043914.1 ===")
print("  Missing from NC_043914:", sorted(mt_xy - mt_ref))
print("  Present only in NC_043914:", sorted(mt_ref - mt_xy))

# BLASTp all 849 SRR9057080 Prodigal ORFs against unified MT+CP protein DB!
faa_849 = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
out_blast = os.path.join(out_dir, "orfs849_vs_unified_mt_cp.tsv")
subprocess.run([
    "blastp", "-query", faa_849, "-db", unified_faa,
    "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle",
    "-evalue", "1e-5", "-max_target_seqs", "3", "-out", out_blast
], check=True)

best_hits = {}
with open(out_blast) as f:
    for row in csv.reader(f, delimiter="\t"):
        qid = row[0]
        bs = float(row[11])
        if qid not in best_hits or bs > float(best_hits[qid][11]):
            best_hits[qid] = row

print("\n=== RE-ANNOTATION OF PREVIOUS 'CLASS C/D' & NOTABLE ORFs USING UNIFIED CAM-MT + CP DB ===")
targets = [
    "NODE_1_length_7316_cov_129.055748_3", "NODE_1_length_7316_cov_129.055748_4",
    "NODE_2_length_3671_cov_47.563937_2",
    "NODE_7_length_3049_cov_278.377038_4", "NODE_7_length_3049_cov_278.377038_5",
    "NODE_8_length_2913_cov_26.061610_2", "NODE_8_length_2913_cov_26.061610_3", "NODE_8_length_2913_cov_26.061610_4",
    "NODE_9_length_2846_cov_208.218533_2", "NODE_9_length_2846_cov_208.218533_3",
    "NODE_10_length_2350_cov_28.737194_1", "NODE_10_length_2350_cov_28.737194_2",
    "NODE_13_length_2010_cov_3.946982_2", "NODE_15_length_1916_cov_19.127554_1",
    "NODE_21_length_1417_cov_2.653963_2", "NODE_24_length_1279_cov_5.706985_1",
    "NODE_28_length_1243_cov_3.372583_1", "NODE_51_length_753_cov_6.608025_1",
]
for tid in targets:
    h = best_hits.get(tid)
    if h:
        print(f"  {tid:<38} -> {h[1]:<35} pid={h[2]}% len={h[3]}/{h[12]}aa (subj {h[13]}aa) E={h[10]}")
    else:
        print(f"  {tid:<38} -> NO HIT IN MT OR CP DB")

# Count unique MT vs CP proteins recovered in SRR9057080 assembly!
mt_genes_recovered = set()
cp_genes_recovered = set()
for qid, h in best_hits.items():
    pid = float(h[2])
    if pid >= 80.0:
        parts = h[1].split("|")
        if parts[0] == "MT":
            mt_genes_recovered.add(parts[2])
        elif parts[0] == "CP":
            cp_genes_recovered.add(parts[2])

print(f"\n=== TOTAL UNIQUE GENES RECOVERED IN SRR9057080 DE NOVO ASSEMBLY (>=80% ID) ===")
print(f"  True Mitochondrial CDS genes recovered ({len(mt_genes_recovered)}): {sorted(mt_genes_recovered)}")
print(f"  Co-assembled Chloroplast CDS genes ({len(cp_genes_recovered)}): {sorted(cp_genes_recovered)}")
