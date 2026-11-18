import os
import csv
import numpy as np
from Bio import PDB
from Bio.Seq import Seq

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Protein_Structure"
os.makedirs(out_dir, exist_ok=True)

wt_pdb = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NC_043914_WT_atp1_predicted.pdb"
ed_pdb = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures/NODE_17_length_1781_cov_62.065632_1_predicted.pdb"

parser = PDB.PDBParser(QUIET=True)
s_wt = parser.get_structure("WT", wt_pdb)
s_ed = parser.get_structure("ED", ed_pdb)

wt_ca = [atom for atom in s_wt.get_atoms() if atom.get_name() == "CA"]
ed_ca = [atom for atom in s_ed.get_atoms() if atom.get_name() == "CA"]

assert len(wt_ca) == len(ed_ca) == 509, f"Length mismatch: {len(wt_ca)} vs {len(ed_ca)}"

# Global C-alpha RMSD
sup = PDB.Superimposer()
sup.set_atoms(wt_ca, ed_ca)
global_rmsd = sup.rms
sup.apply(s_ed.get_atoms())

# C-terminal domain (residues 340..509 where all 8 substitutions reside) RMSD
sup_ct = PDB.Superimposer()
sup_ct.set_atoms(wt_ca[339:], ed_ca[339:])
ct_rmsd = sup_ct.rms

# Per-residue C-alpha distance & B-factor (pLDDT) & Kyte-Doolittle hydrophobicity
kd = {'ALA':1.8,'ARG':-4.5,'ASN':-3.5,'ASP':-3.5,'CYS':2.5,'GLN':-3.5,'GLU':-3.5,'GLY':-0.4,'HIS':-3.2,'ILE':4.5,'LEU':3.8,'LYS':-3.9,'MET':1.9,'PHE':2.8,'PRO':-1.6,'SER':-0.8,'THR':-0.7,'TRP':-0.9,'TYR':-1.3,'VAL':4.2}

wt_plddt = [a.get_bfactor() * 100.0 if a.get_bfactor() <= 1.0 else a.get_bfactor() for a in wt_ca]
ed_plddt = [a.get_bfactor() * 100.0 if a.get_bfactor() <= 1.0 else a.get_bfactor() for a in ed_ca]

print("=== 3D STRUCTURAL SUPERPOSITION: WILD-TYPE vs RNA-EDITED ATP1 (509 aa) ===")
print(f"  WT Mean Global pLDDT     : {np.mean(wt_plddt):.2f} / 100 (Median: {np.median(wt_plddt):.2f})")
print(f"  Edited Mean Global pLDDT : {np.mean(ed_plddt):.2f} / 100 (Median: {np.median(ed_plddt):.2f})")
print(f"  Global C-alpha RMSD      : {global_rmsd:.3f} A (over all 509 residues)")
print(f"  C-term Domain RMSD (340-509): {ct_rmsd:.3f} A")

mut_positions = [347, 355, 393, 406, 431, 446, 472, 497]
wt_coords = np.array([a.get_coord() for a in wt_ca])
ed_coords = np.array([a.get_coord() for a in ed_ca])

rows = []
print("\n  Substituted / RNA-Edited Residue Comparison:")
tot_kd_wt, tot_kd_ed = 0.0, 0.0
for pos in mut_positions:
    idx = pos - 1
    w_res = wt_ca[idx].get_parent().get_resname()
    e_res = ed_ca[idx].get_parent().get_resname()
    w_p = wt_plddt[idx]
    e_p = ed_plddt[idx]
    ca_shift = float(np.linalg.norm(wt_coords[idx] - ed_coords[idx]))
    # Packing neighbors within 8.0 A of CA
    w_nbrs = int(np.sum(np.linalg.norm(wt_coords - wt_coords[idx], axis=1) <= 8.0)) - 1
    e_nbrs = int(np.sum(np.linalg.norm(ed_coords - ed_coords[idx], axis=1) <= 8.0)) - 1
    h_wt, h_ed = kd.get(w_res, 0), kd.get(e_res, 0)
    tot_kd_wt += h_wt; tot_kd_ed += h_ed
    print(f"    Res {pos}: {w_res}->{e_res} | pLDDT {w_p:.1f}->{e_p:.1f} ({e_p-w_p:+.1f}) | KD {h_wt:+.1f}->{h_ed:+.1f} ({h_ed-h_wt:+.1f}) | CA_shift={ca_shift:.2f}A | Nbrs8A {w_nbrs}->{e_nbrs}")
    rows.append({
        "Residue": pos, "WT_AA": w_res, "Edited_AA": e_res,
        "WT_pLDDT": round(w_p, 2), "Edited_pLDDT": round(e_p, 2), "Delta_pLDDT": round(e_p - w_p, 2),
        "WT_Hydropathy_KD": h_wt, "Edited_Hydropathy_KD": h_ed, "Delta_KD": round(h_ed - h_wt, 2),
        "CA_Shift_A": round(ca_shift, 3), "Contacts_8A_WT": w_nbrs, "Contacts_8A_Edited": e_nbrs
    })

print(f"\n  Cumulative Hydropathy across 8 edited positions: WT={tot_kd_wt:+.1f} -> Edited={tot_kd_ed:+.1f} (Delta={tot_kd_ed-tot_kd_wt:+.1f} kcal/mol scale)")

out_csv = os.path.join(out_dir, "atp1_wt_vs_edited_structural_comparison.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"Saved detailed WT vs Edited structural comparison to {out_csv}")
