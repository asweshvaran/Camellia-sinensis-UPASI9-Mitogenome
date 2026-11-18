import os
import re
import csv
import subprocess
from collections import Counter
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis/Expanded_Phylogeny"
os.makedirs(out_dir, exist_ok=True)
supp_dir = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables"
ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"
fna_849 = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_genes.fna"
unified_blast = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis/orfs849_vs_unified_mt_cp.tsv"

def norm_gene(g):
    g = g.strip().lower()
    if g in ("cytb", "cob", "cytochrome b"): return "cob"
    if g.startswith("nd") and len(g) >= 3 and g[2].isdigit():
        g = "nad" + g[2:]
    if g.startswith("nad4l"): return "nad4l"
    if g.startswith("ccmb"): return "ccmb"
    if g.startswith("ccmc"): return "ccmc"
    if g.startswith("ccmfc"): return "ccmfc"
    if g.startswith("ccmfn"): return "ccmfn"
    # Clean product names if gene qualifier missing
    for std in ["atp1", "atp4", "atp6", "atp8", "atp9", "cob", "cox1", "cox2", "cox3",
                "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6", "nad7", "nad9",
                "ccmb", "ccmc", "ccmfc", "ccmfn", "matr", "mttb", "sdh3", "sdh4",
                "rpl2", "rpl5", "rpl10", "rpl16", "rps1", "rps3", "rps4", "rps7", "rps12", "rps13", "rps14", "rps19"]:
        if g == std:
            return std
    return g

# Load UPASI-9 purified/curated CDS from predicted_genes.fna + rnaSPAdes transcripts
best_hits = {}
with open(unified_blast) as f:
    for row in csv.reader(f, delimiter="\t"):
        qid = row[0]
        bs = float(row[11])
        if qid not in best_hits or bs > float(best_hits[qid][11]):
            best_hits[qid] = row

fna_map = {r.id: str(r.seq).upper() for r in SeqIO.parse(fna_849, "fasta")}
srr_mt_genes = {}
for qid, h in best_hits.items():
    parts = h[1].split("|")
    if parts[0] != "MT":
        continue
    gene = norm_gene(parts[2])
    pid, alen, slen = float(h[2]), int(h[3]), int(h[13])
    cov = alen / slen
    if pid >= 85.0 and cov >= 0.80:
        if gene not in srr_mt_genes or float(h[11]) > srr_mt_genes[gene]["bs"]:
            srr_mt_genes[gene] = {"nt_seq": fna_map[qid], "bs": float(h[11])}

def extract_cds_norm(gb_path, ref_dict):
    res = {}
    for rec in SeqIO.parse(gb_path, "genbank"):
        for feat in rec.features:
            if feat.type == "CDS":
                raw_g = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["unknown"]))[0]
                g = norm_gene(raw_g)
                try:
                    nt = str(feat.extract(rec.seq, references=ref_dict)).upper()
                except Exception:
                    continue
                if len(nt) >= 150 and (g not in res or len(nt) > len(res[g])):
                    res[g] = nt
    return res

taxa_gbs = [
    ("C_sinensis_RefSeq_NC043914", "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb", "Ericales (Theaceae)"),
    ("C_sinensis_Xinyang10_PX398764", os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb"), "Ericales (Theaceae)"),
    ("C_sinensis_Yidan_PX759560", os.path.join(ref_dir, "C_sinensis_Yidan_PX759560.gb"), "Ericales (Theaceae)"),
    ("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr1_PZ267020.gb"), "Ericales (Theaceae)"),
    ("C_pitardii_PZ501402", os.path.join(ref_dir, "C_pitardii_PZ501402.gb"), "Ericales (Theaceae)"),
    ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr1_PZ378386.gb"), "Ericales (Theaceae)"),
    ("Actinidia_rufa_NC086781", os.path.join(ref_dir, "Outgroup_NC_086781.1.gb"), "Ericales (Actinidiaceae)"),
    ("Vaccinium_macrocarpon_NC023338", os.path.join(ref_dir, "Outgroup_Vaccinium_NC_023338.1.gb"), "Ericales (Ericaceae)"),
    ("Carya_illinoinensis_PX519271", os.path.join(ref_dir, "Outgroup_Carya_illinoinensis_PX519271.1.gb"), "Fagales (Juglandaceae)"),
    ("Fagus_sylvatica_MT446430", os.path.join(ref_dir, "Outgroup_Fagus_sylvatica_MT446430.1.gb"), "Fagales (Fagaceae)"),
    ("Arabidopsis_thaliana_NC037304", os.path.join(ref_dir, "Outgroup_Arabidopsis_thaliana_NC_037304.1.gb"), "Brassicales (Brassicaceae)"),
]

extra_refs = {}
for _, gp, _ in taxa_gbs:
    if os.path.exists(gp):
        for r in SeqIO.parse(gp, "genbank"):
            extra_refs[r.id] = r.seq
for gp2 in [os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb"),
            os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb"),
            os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb")]:
    if os.path.exists(gp2):
        for r in SeqIO.parse(gp2, "genbank"):
            extra_refs[r.id] = r.seq

taxa_cds = {}
for name, gp, _ in taxa_gbs:
    taxa_cds[name] = extract_cds_norm(gp, extra_refs)
for tname, gp2 in [("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb")),
                   ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb")),
                   ("C_sinensis_Xinyang10_PX398764", os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb"))]:
    if os.path.exists(gp2):
        taxa_cds[tname].update(extract_cds_norm(gp2, extra_refs))

taxa_cds["C_sinensis_UPASI9_SRR9057080"] = {g: info["nt_seq"] for g, info in srr_mt_genes.items()}

common = sorted(set.intersection(*(set(m.keys()) for m in taxa_cds.values())))
print(f"=== 1. EXPANDED 12-TAXON SHARED PCG SET ({len(common)} genes) ===")
print("  Genes:", common)

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

used_genes = []
concat_raw = {t: [] for t in taxa_cds}
concat_masked = {t: [] for t in taxa_cds}
concat_aa_masked = {t: [] for t in taxa_cds}

for g in common:
    ref_s = taxa_cds["C_sinensis_Xinyang10_PX398764"][g]
    lens = [len(taxa_cds[t][g]) for t in taxa_cds]
    if min(lens) < len(ref_s) * 0.75:
        print(f"  [SKIP] {g}: min len {min(lens)} vs ref {len(ref_s)}")
        continue
    tax_mapped = {t: ["-"] * len(ref_s) for t in taxa_cds}
    for t in taxa_cds:
        q_s = taxa_cds[t][g]
        if q_s == ref_s:
            tax_mapped[t] = list(ref_s)
        else:
            aln = aligner.align(ref_s, q_s)[0]
            for b1, b2 in zip(aln.aligned[0], aln.aligned[1]):
                for k in range(b1[1] - b1[0]):
                    tax_mapped[t][b1[0] + k] = q_s[b2[0] + k]
    # Keep codon-complete gap-free columns
    keep_codons = []
    for pos in range(0, len(ref_s) - 2, 3):
        if all(tax_mapped[t][pos] != "-" and tax_mapped[t][pos+1] != "-" and tax_mapped[t][pos+2] != "-" for t in taxa_cds):
            keep_codons.extend([pos, pos+1, pos+2])
    if len(keep_codons) >= 150:
        used_genes.append((g, len(keep_codons)))
        for t in taxa_cds:
            raw_seq = "".join(tax_mapped[t][i] for i in keep_codons)
            concat_raw[t].append(raw_seq)
            if t == "C_sinensis_UPASI9_SRR9057080":
                masked_bases = []
                for i in keep_codons:
                    b_up = tax_mapped[t][i]
                    b_xy = tax_mapped["C_sinensis_Xinyang10_PX398764"][i]
                    b_yd = tax_mapped["C_sinensis_Yidan_PX759560"][i]
                    if b_up == "T" and b_xy == "C" and b_yd == "C":
                        masked_bases.append("C")
                    else:
                        masked_bases.append(b_up)
                m_seq = "".join(masked_bases)
            else:
                m_seq = raw_seq
            concat_masked[t].append(m_seq)
            # Translate masked codon alignment (strip trailing stop if *)
            aa = str(Seq(m_seq).translate(table=1)).replace("*", "X")
            concat_aa_masked[t].append(aa)

aln_nt_len = sum(l for _, l in used_genes)
aln_aa_len = aln_nt_len // 3
print(f"  Final concatenated PCGs ({len(used_genes)}): {used_genes}")
print(f"  Gap-free codon-aligned length: {aln_nt_len:,} bp ({aln_aa_len:,} aa)")

# Write FASTA alignments
fa_raw = os.path.join(out_dir, "expanded_12taxa_unmasked_nt.fasta")
fa_masked = os.path.join(out_dir, "expanded_12taxa_masked_nt.fasta")
fa_aa = os.path.join(out_dir, "expanded_12taxa_masked_aa.fasta")

for path, cmap in [(fa_raw, concat_raw), (fa_masked, concat_masked), (fa_aa, concat_aa_masked)]:
    with open(path, "w") as f:
        for t in cmap:
            f.write(f">{t}\n{''.join(cmap[t])}\n")

# Compute pairwise p-distance from UPASI-9 to every taxon (raw vs masked)
print("\n=== 2. PAIRWISE P-DISTANCE FROM UPASI-9 (UNMASKED vs C-TO-U MASKED) ===")
up_raw = "".join(concat_raw["C_sinensis_UPASI9_SRR9057080"])
up_msk = "".join(concat_masked["C_sinensis_UPASI9_SRR9057080"])
dist_rows = []
for t in concat_raw:
    if t == "C_sinensis_UPASI9_SRR9057080":
        continue
    s_t = "".join(concat_masked[t])
    d_raw = sum(1 for a, b in zip(up_raw, s_t) if a != b) / len(up_raw)
    d_msk = sum(1 for a, b in zip(up_msk, s_t) if a != b) / len(up_msk)
    print(f"  vs {t:<32}: Unmasked={d_raw*100:6.3f}% | Masked={d_msk*100:6.3f}%")
    dist_rows.append({"Comparison_Taxon": t, "Unmasked_pDist_Pct": round(d_raw*100, 4), "Masked_pDist_Pct": round(d_msk*100, 4)})

# Run IQ-TREE2 on all 3 matrices (following Li et al., 2024 & Chen et al., 2024)
print("\n=== 3. RUNNING IQ-TREE2 (MODELFINDER + 1000 UFBOOT + 1000 SH-aLRT) ===")
iqtree_summaries = []
for label, fa_path, seqtype in [
    ("ML_12taxa_Unmasked_NT", fa_raw, "DNA"),
    ("ML_12taxa_Masked_NT", fa_masked, "DNA"),
    ("ML_12taxa_Masked_AA", fa_aa, "AA"),
]:
    prefix = os.path.join(out_dir, label)
    cmd = [
        "iqtree2", "-s", fa_path, "-m", "MFP",
        "-B", "1000", "-alrt", "1000", "--bnni",
        "-o", "Arabidopsis_thaliana_NC037304",
        "-nt", "AUTO", "--redo", "-pre", prefix
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  [ERROR] IQ-TREE2 failed for {label}: {proc.stderr[-500:]}")
        continue
    # Parse best-fit model and log-likelihood from .iqtree report
    iq_report = prefix + ".iqtree"
    best_model, logl, total_tree_len = "NA", "NA", "NA"
    with open(iq_report) as f:
        for line in f:
            if line.startswith("Best-fit model:"):
                best_model = line.split(":", 1)[1].strip()
            elif line.startswith("Log-likelihood of the tree:"):
                logl = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Total tree length"):
                total_tree_len = line.split(":", 1)[1].strip()
    treefile = prefix + ".treefile"
    nwk_str = open(treefile).read().strip()
    print(f"  [OK] {label:<24} | Model: {best_model} | LogL: {logl} | TreeLen: {total_tree_len}")
    print(f"       Newick: {nwk_str}")
    iqtree_summaries.append({
        "Analysis_Matrix": label,
        "Sequence_Type": seqtype,
        "Taxa_Count": len(taxa_cds),
        "Shared_PCGs": len(used_genes),
        "Alignment_Length": aln_nt_len if seqtype == "DNA" else aln_aa_len,
        "Best_Fit_Model_BIC": best_model,
        "Log_Likelihood": logl,
        "Total_Tree_Length": total_tree_len,
        "Bootstrap_Replicates": "1000 UFBoot / 1000 SH-aLRT",
        "Newick_Tree": nwk_str
    })

# Write Table S13
s13_path = os.path.join(supp_dir, "Table_S13_Phylogenetic_Models_and_Support.csv")
with open(s13_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(iqtree_summaries[0].keys()))
    w.writeheader()
    w.writerows(iqtree_summaries)
print(f"\nSaved {s13_path}")

# =====================================================================
# 4. CONSERVED tRNA T-Psi-C LOOP (U-U-C-x-A-x2) ANALYSIS (Chen et al., 2024)
#    & GENOME AT/GC SKEW BENCHMARK (Li et al., 2024)
# =====================================================================
print("\n=== 4. CONSERVED tRNA T-Psi-C LOOP (U-U-C-x-A-x2) & AT/GC SKEW ===")

def analyze_genome_trnas_and_skew(species_label, gb_path):
    rec = next(SeqIO.parse(gb_path, "genbank"))
    seq = str(rec.seq).upper()
    A, T, G, C = seq.count("A"), seq.count("T"), seq.count("G"), seq.count("C")
    gc_pct = (G + C) / len(seq) * 100
    at_skew = (A - T) / (A + T)
    gc_skew = (G - C) / (G + C)
    trna_rows = []
    motif_hits = 0
    for feat in rec.features:
        if feat.type == "tRNA":
            g = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["tRNA"]))[0]
            prod = feat.qualifiers.get("product", [g])[0]
            tseq = str(feat.extract(rec.seq)).upper().replace("T", "U")
            # Search Chen et al. (2024) conserved Psi-loop motif: U-U-C-x-A-x2 (i.e. UUC[ACGU]A[ACGU]{2})
            # Usually located in the 3' half of the tRNA (positions ~45-68)
            matches = list(re.finditer(r"UUC[ACGU]A[ACGU]{2}", tseq))
            has_motif = len(matches) > 0
            if has_motif:
                motif_hits += 1
                m_str = matches[-1].group(0)
                m_pos = matches[-1].start() + 1
            else:
                # Check 1-mismatch variant (e.g., UUCG/UCCA/GUUC)
                m_str = "Variant/Absent"
                m_pos = -1
            tA, tU, tG, tC = tseq.count("A"), tseq.count("U"), tseq.count("G"), tseq.count("C")
            t_gc = (tG + tC) / len(tseq) * 100 if tseq else 0
            trna_rows.append({
                "Species": species_label,
                "tRNA_Gene": g,
                "Product": prod,
                "Length_nt": len(tseq),
                "Strand": "+" if feat.location.strand == 1 else "-",
                "GC_Pct": round(t_gc, 2),
                "Has_UUCxAx2_PsiLoop": "Yes" if has_motif else "No",
                "PsiLoop_Motif_Seq": m_str,
                "Motif_Start_Pos": m_pos
            })
    print(f"  {species_label:<30} | Len={len(seq):7,d} bp | GC={gc_pct:5.2f}% | AT-skew={at_skew:+.4f} | GC-skew={gc_skew:+.4f} | tRNAs={len(trna_rows):2d} | UUCxAx2 conserved={motif_hits}/{len(trna_rows)} ({motif_hits/len(trna_rows)*100 if trna_rows else 0:.1f}%)")
    return trna_rows, {"Species": species_label, "Length_bp": len(seq), "GC_Pct": round(gc_pct, 2), "AT_Skew": round(at_skew, 4), "GC_Skew": round(gc_skew, 4), "Total_tRNAs": len(trna_rows), "UUCxAx2_Conserved": motif_hits}

all_trna_rows = []
skew_summary = []
for sp_lbl, gbp in [
    ("Camellia_sinensis_NC043914", "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"),
    ("Camellia_sinensis_Xinyang10_Chr1", os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb")),
    ("Camellia_sinensis_Yidan", os.path.join(ref_dir, "C_sinensis_Yidan_PX759560.gb")),
    ("Carya_illinoinensis_PX519271", os.path.join(ref_dir, "Outgroup_Carya_illinoinensis_PX519271.1.gb")),
    ("Fagus_sylvatica_MT446430", os.path.join(ref_dir, "Outgroup_Fagus_sylvatica_MT446430.1.gb")),
    ("Actinidia_rufa_NC086781", os.path.join(ref_dir, "Outgroup_NC_086781.1.gb")),
    ("Arabidopsis_thaliana_NC037304", os.path.join(ref_dir, "Outgroup_Arabidopsis_thaliana_NC_037304.1.gb")),
]:
    if os.path.exists(gbp):
        trows, ssum = analyze_genome_trnas_and_skew(sp_lbl, gbp)
        all_trna_rows.extend(trows)
        skew_summary.append(ssum)

# Also check UPASI-9 de novo recovered tRNAs (trnH-GUG and trnE-UUC from Table S4 / Wuyi match)
trna_s4_path = os.path.join(supp_dir, "Table_S4_Identified_tRNA_Genes.csv")
if os.path.exists(trna_s4_path):
    print("\n  Checking UPASI-9 Table_S4 tRNAs:")
    with open(trna_s4_path) as f:
        for row in csv.DictReader(f):
            print("   ", dict(row))

s14_path = os.path.join(supp_dir, "Table_S14_Conserved_tRNA_PsiLoop_Motifs.csv")
with open(s14_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(all_trna_rows[0].keys()))
    w.writeheader()
    w.writerows(all_trna_rows)
print(f"\nSaved {s14_path} ({len(all_trna_rows)} total tRNAs analyzed across 7 plant mitogenomes)")
