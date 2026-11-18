import numpy as np
from Bio import SeqIO

pdb_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NODE_17_length_1781_cov_62.065632_1_predicted.pdb"
ref_faa = "/mnt/d/Mitochondria/reference/ref_mitochondrial_proteins.faa"

# 1. Parse pLDDT from PDB B-factors
ca_plddt = []
res_names = []
with open(pdb_path) as f:
    for line in f:
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            b_factor = float(line[60:66].strip())
            # Convert 0-1 scale to 0-100 if needed
            score = b_factor * 100.0 if b_factor <= 1.0 else b_factor
            ca_plddt.append(score)
            res_names.append(line[17:20].strip())

ca_plddt = np.array(ca_plddt)
mean_plddt = np.mean(ca_plddt)
median_plddt = np.median(ca_plddt)
pct_vhigh = np.sum(ca_plddt >= 90) / len(ca_plddt) * 100 # Very high confidence
pct_high = np.sum((ca_plddt >= 70) & (ca_plddt < 90)) / len(ca_plddt) * 100 # Confident
pct_low = np.sum((ca_plddt >= 50) & (ca_plddt < 70)) / len(ca_plddt) * 100 # Low
pct_vlow = np.sum(ca_plddt < 50) / len(ca_plddt) * 100 # Very low

ref_atp1 = None
for r in SeqIO.parse(ref_faa, "fasta"):
    if r.id.startswith("atp1"):
        ref_atp1 = str(r.seq)
        break

query_atp1 = None
faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
for r in SeqIO.parse(faa_path, "fasta"):
    if r.id.startswith("NODE_17_length_1781_cov_62.065632_1"):
        query_atp1 = str(r.seq).rstrip("*")
        break

diffs = []
min_len = min(len(ref_atp1), len(query_atp1))
for i in range(min_len):
    if ref_atp1[i] != query_atp1[i]:
        diffs.append((i + 1, ref_atp1[i], query_atp1[i], ca_plddt[i] if i < len(ca_plddt) else 0.0))

report = f"""==========================================================
ATP1 (ATP SYNTHASE ALPHA SUBUNIT) STRUCTURAL BIOINFORMATICS
==========================================================
Gene Name                     : atp1 (F1-ATP synthase regulatory alpha subunit)
Reference Protein Length      : {len(ref_atp1)} amino acids (NC_043914.1)
Assembled Protein Length      : {len(query_atp1)} amino acids (SRR9057080, UPASI-9)
Sequence Identity to Ref      : {(min_len - len(diffs)) / min_len * 100:.2f}%
Total Amino Acid Substitutions: {len(diffs)}
----------------------------------------------------------
ESMFOLD 3D PREDICTION CONFIDENCE (Standard 0-100 pLDDT Scale):
  - Mean Global pLDDT         : {mean_plddt:.2f} / 100
  - Median pLDDT              : {median_plddt:.2f} / 100
  - Very High Conf (>= 90)    : {pct_vhigh:.1f}% of residues (Crystal-grade)
  - Confident (70 - 90)       : {pct_high:.1f}% of residues
  - Low Confidence (50 - 70)  : {pct_low:.1f}% of residues
  - Disordered / Very Low (<50): {pct_vlow:.1f}% of residues
----------------------------------------------------------
Identified Amino Acid Substitutions (Putative RNA Editing Signatures):
"""
for pos, r_aa, q_aa, pld in diffs:
    report += f"  Residue {pos:<4}: Ref {r_aa} -> Var {q_aa} | Local pLDDT: {pld:.1f}\n"

report += "==========================================================\n"
print(report)

out_report = "/mnt/d/Mitochondria/Q1_Investigation/Protein_Structure/atp1_structural_metrics.txt"
with open(out_report, "w") as f:
    f.write(report)
