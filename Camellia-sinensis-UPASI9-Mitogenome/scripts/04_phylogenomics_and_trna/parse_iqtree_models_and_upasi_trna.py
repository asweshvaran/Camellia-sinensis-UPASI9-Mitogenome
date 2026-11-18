import os
import re
import csv
from Bio import SeqIO

exp_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis/Expanded_Phylogeny"
s13_path = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables/Table_S13_Phylogenetic_Models_and_Support.csv"

models = {}
for label in ["ML_12taxa_Unmasked_NT", "ML_12taxa_Masked_NT", "ML_12taxa_Masked_AA"]:
    iqf = os.path.join(exp_dir, label + ".iqtree")
    best_m = "NA"
    with open(iqf) as f:
        for line in f:
            if "Best-fit model" in line and ":" in line:
                best_m = line.split(":", 1)[1].strip()
                break
    models[label] = best_m
    print(f"  {label:<24} -> Best-fit model (BIC): {best_m}")

# Update Table S13 with exact best-fit models
rows = list(csv.DictReader(open(s13_path)))
for r in rows:
    r["Best_Fit_Model_BIC"] = models.get(r["Analysis_Matrix"], r["Best_Fit_Model_BIC"])
with open(s13_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print("Updated Table_S13 with exact ModelFinder BIC models.")

# Check UPASI-9 de novo recovered trnH-GUG (NODE_61) and trnE-UUC (NODE_89) for U-U-C-x-A-x2 motif
scaf_fa = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
seqs = {r.id: r.seq for r in SeqIO.parse(scaf_fa, "fasta")}
print("\nChecking UPASI-9 de novo recovered mitochondrial tRNAs:")
for name, scid, beg, end in [
    ("UPASI-9_trnH-GUG", "NODE_61_length_551_cov_1.560538", 6, 79),
    ("UPASI-9_trnE-UUC", "NODE_89_length_336_cov_0.584416", 297, 226),
]:
    s = seqs[scid]
    if beg < end:
        tseq = str(s[beg-1:end]).upper().replace("T", "U")
    else:
        tseq = str(s[end-1:beg].reverse_complement()).upper().replace("T", "U")
    matches = list(re.finditer(r"UUC[ACGU]A[ACGU]{2}", tseq))
    print(f"  {name:<18} | Len={len(tseq)} nt | Seq={tseq}")
    print(f"     UUCxAx2 motif: {[ (m.group(0), m.start()+1) for m in matches ]}")
