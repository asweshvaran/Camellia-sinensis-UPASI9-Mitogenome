import os
import csv
import shutil
import subprocess
from collections import defaultdict
from Bio import SeqIO, PDB

supp_dir = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables"
clean_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
ann_dir = "/mnt/d/Mitochondria/Q1_Investigation/Annotation"
prot_dir = "/mnt/d/Mitochondria/Q1_Investigation/Protein_Structure"
var_dir = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants"
eval_dir = "/mnt/d/Mitochondria/Q1_Investigation/Assembly_Evaluation"
script_arch = "/mnt/d/Mitochondria/Q1_Investigation/Analysis_Scripts"
leg_arch = os.path.join(script_arch, "Legacy_v1")
os.makedirs(script_arch, exist_ok=True)
os.makedirs(leg_arch, exist_ok=True)

# 0. Clean stray .py/.sh scripts out of Supplementary_Tables/ into Analysis_Scripts/Legacy_v1/
moved = 0
for fn in os.listdir(supp_dir):
    if fn.endswith((".py", ".sh")):
        shutil.move(os.path.join(supp_dir, fn), os.path.join(leg_arch, fn))
        moved += 1
print(f"Archived {moved} legacy scripts from Supplementary_Tables/ into Analysis_Scripts/Legacy_v1/")

# 1. Reconcile Table_S2: Side-by-Side Legacy GetOrganelle vs Purified rnaSPAdes Assembly
tx_fasta = os.path.join(clean_dir, "spades_rna_purified_mt", "transcripts.fasta")
pur_recs = list(SeqIO.parse(tx_fasta, "fasta"))
pur_lens = sorted([len(r.seq) for r in pur_recs], reverse=True)
tot_bp = sum(pur_lens)
cs = 0; p_n50, p_l50 = 0, 0
for i, l in enumerate(pur_lens):
    cs += l
    if cs >= tot_bp / 2.0: p_n50, p_l50 = l, i + 1; break
p_gc = sum(r.seq.upper().count("G") + r.seq.upper().count("C") for r in pur_recs) / tot_bp * 100.0
recs500 = [r for r in pur_recs if len(r.seq) >= 500]
l500 = sorted([len(r.seq) for r in recs500], reverse=True)
bp500 = sum(l500)
cs = 0; n50_500, l50_500 = 0, 0
for i, l in enumerate(l500):
    cs += l
    if cs >= bp500 / 2.0: n50_500, l50_500 = l, i + 1; break

rows_s2 = [
    ["Metric", "Legacy GetOrganelle (extended_spades)", "GetOrganelle K105 Selected Graph", "Purified Xinyang10-Guided rnaSPAdes (All)", "Purified Xinyang10-Guided rnaSPAdes (>=500 bp)"],
    ["Input / Filtering Strategy", "Unfiltered GetOrganelle embplant_mt recruitment (53.93M reads)", "K105 embplant_mt path selection", "Competitive BWA-MEM vs Xinyang10 (PX398764-65) + cpDNA (NC_020019.1)", "Post-assembly length filter (>=500 bp) on purified transcripts"],
    ["Total Contigs / Scaffolds", "1,160", "2 isolated nodes", "2,295", "217"],
    ["Total Assembled Length", "258,681 bp (36.57% of NC_043914.1)", "3,333 bp (1,574 + 1,759 bp)", "854,700 bp", "286,843 bp"],
    ["Assembly N50 / L50", "193 bp / L50=202", "1,759 bp / L50=1", f"{p_n50:,} bp / L50={p_l50}", f"{n50_500:,} bp / L50={l50_500}"],
    ["Longest Sequence", "7,316 bp (NODE_1 = Chloroplast DNA!)", "1,759 bp (Plastid rRNA)", f"{pur_lens[0]:,} bp (NODE_1 = MT rpl2-rps19-rpl16)", f"{l500[0]:,} bp (NODE_1 = MT rpl2-rps19-rpl16)"],
    ["GC Content (%)", "49.82%", "54.25% (rRNA)", f"{p_gc:.2f}%", f"{sum(r.seq.upper().count('G')+r.seq.upper().count('C') for r in recs500)/bp500*100:.2f}%"],
    ["Mitochondrial-Specific Sequence", "73,801 bp (106 scaffolds; 28.53%)", "0 bp (0.0%)", "466,552 bp (941 contigs; 54.6%)", "237,395 bp (186 contigs; 82.8%)"],
    ["MTPT / Shared Organellar Homologs", "116,989 bp (729 scaffolds; 45.23%)", "3,333 bp (plastid rRNA)", "60,192 bp (65 contigs; 7.0%)", "48,907 bp (30 contigs; 17.0%)"],
    ["Pure Chloroplast Contamination", "60,955 bp (313 scaffolds; 23.56%)", "100% plastid-biased rRNA (24,259.5x artifact)", "0 bp (0 contigs; 0.0%)", "0 bp (0 contigs; 0.0%)"],
    ["Unclassified Sequence", "6,936 bp (12 scaffolds; 2.68%)", "0 bp", "327,956 bp (1,289 short contigs <500bp)", "541 bp (1 low-cov poly-AT repeat; 0.2%)"],
    ["GFA Graph Nodes / Links", "1,535 nodes / 1,963 links (815 hubs)", "2 nodes / 0 links", "4,222 segments / 114 links (2,295 paths)", "Subgraph >=500 bp"],
    ["Normalized Canonical MT CDS Recovered", "30 MT (+ 16 CP CDS)", "0 MT CDS", "37 MT CDS (+ 3 MTPT co-assembled)", "37 MT CDS (full-length polycistronic)"],
]
with open(os.path.join(supp_dir, "Table_S2_De_Novo_Assembly_Metrics.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows_s2)
print("Updated Table_S2_De_Novo_Assembly_Metrics.csv")

# 2. Reconcile Table_S3: Full Reconciled Inventory of All Recovered Mitochondrial Genes + Reference Coverage
# Merge NC_043914_feature_depth_profile.csv (primary depth/breadth) + purified_proteins_vs_unified.tsv + RNA editing counts
feat_cov = {}
with open(os.path.join(ann_dir, "nc043914_feature_depth_profile.csv")) as f:
    for r in csv.DictReader(f):
        g = r["Gene"]
        # Keep highest-coverage primary locus if duplicated (e.g., rps13)
        if g not in feat_cov or float(r["Mean_Depth"]) > float(feat_cov[g]["Mean_Depth"]):
            feat_cov[g] = r

alias = {"nd1":"nad1", "nd2":"nad2", "nd3":"nad3", "nd4":"nad4", "nd4l":"nad4l", "nd5":"nad5", "nd6":"nad6", "nd7":"nad7", "nd9":"nad9", "mttb":"mttB", "matr":"matR", "ccmb":"ccmB", "ccmc":"ccmC", "ccmfc":"ccmFc", "ccmfn":"ccmFn"}
pur_mt_hits = {}
with open(os.path.join(clean_dir, "purified_proteins_vs_unified.tsv")) as f:
    for r in csv.reader(f, delimiter="\t"):
        pid, aln, ql, sl = float(r[2]), int(r[3]), int(r[4]), int(r[5])
        if pid >= 80.0 and r[1].startswith("MT|"):
            src, g = r[1].split("|")[1], alias.get(r[1].split("|")[2].lower(), r[1].split("|")[2].lower())
            cov = min(100.0, aln / sl * 100.0)
            if g not in pur_mt_hits or (cov, pid) > (pur_mt_hits[g][0], pur_mt_hits[g][1]):
                pur_mt_hits[g] = (cov, pid, ql, sl, r[0], src)

# Count canonical CDS edits per gene from curated_variants_feature_annotated.csv
cds_edits = defaultdict(int)
with open(os.path.join(var_dir, "curated_variants_feature_annotated.csv")) as f:
    for r in csv.DictReader(f):
        if r["region"] == "CDS" and ((r["strand"] == "+" and r["ref"] == "C" and r["alt"] == "T") or (r["strand"] == "-" and r["ref"] == "G" and r["alt"] == "A")):
            cds_edits[r["gene"]] += 1

all_genes_s3 = sorted(set(pur_mt_hits.keys()) | {g for g, r in feat_cov.items() if r["Type"] in ("CDS", "rRNA")})
rows_s3 = [[
    "Gene", "Feature_Type", "Annotated_in_NC043914.1", "Pan_Camellia_Status",
    "Purified_rnaSPAdes_Contig", "Query_AA_Len", "Subject_AA_Len", "Protein_Coverage_pct", "Protein_Identity_pct",
    "NC043914_Mean_Depth", "NC043914_Breadth_10x_pct", "Canonical_CtoU_Edits_vs_NC043914", "Reconciliation_Notes"
]]
for g in sorted(all_genes_s3, key=lambda x: (0 if x in pur_mt_hits else 1, x)):
    fc = feat_cov.get(g)
    ph = pur_mt_hits.get(g)
    ftype = fc["Type"] if fc else "CDS"
    in_ref = "Yes" if fc else "No (Missing/Unannotated)"
    pan_st = "Conserved Core Mitogenome" if (ph or g in ("rps16", "rrnL", "rrnS")) else "RefSeq-Specific Hypothetical ORF"
    if ph:
        cov, pid, ql, sl, qid, src = ph
        contig_short = qid.rsplit("_", 1)[0]
        cov_s, pid_s, ql_s, sl_s = f"{cov:.1f}%", f"{pid:.2f}%", str(ql), str(sl)
    else:
        contig_short = "Reference-Mapped Only" if g in ("rps16", "orf100", "orf115b") else "rRNA BLASTn (NODE_3 / NODE_11)"
        cov_s, pid_s, ql_s, sl_s = ("100.0%", "99.47%", "NA", "3581 nt") if g == "rrnL" else (("100.0%", "100.0%", "NA", "1117 nt") if g == "rrnS" else ("NA", "NA", "NA", "NA"))
    mdepth = fc["Mean_Depth"] if fc else "N/A (Unannotated in RefSeq)"
    br10 = fc["Breadth_10x_pct"] if fc else "N/A"
    ed = str(cds_edits.get(g, 0)) if fc and ftype == "CDS" else ("N/A (rRNA)" if ftype == "rRNA" else "See Ka/Ks vs Xinyang10")
    if g in ("atp8", "cox1", "cox3", "nad7", "nad9", "rpl2", "rpl16", "rps3", "rps4", "rps19"):
        note = "Missing in NC_043914.1 RefSeq; resolved via Xinyang10/Wuyi pan-MT DB"
    elif g == "rrnS":
        note = "100% MTPT Plastid-Co-Mapped (24,259.5x artifact source)"
    elif g == "rrnL":
        note = "Core Mitochondrial 26S rRNA (0% MTPT; verified NODE_3 99.47%)"
    elif ph:
        note = "Core Mitochondrial Locus"
    elif g == "rps16":
        note = "Reference-mapped core locus (divergent/pseudogenized in Xinyang10)"
    else:
        note = "RefSeq-annotated expressed hypothetical ORF"
    rows_s3.append([g, ftype, in_ref, pan_st, contig_short, ql_s, sl_s, cov_s, pid_s, mdepth, br10, ed, note])

with open(os.path.join(supp_dir, "Table_S3_Annotated_Core_Mitochondrial_Genes.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows_s3)
print(f"Updated Table_S3_Annotated_Core_Mitochondrial_Genes.csv ({len(rows_s3)-1} loci)")

# 3. Reconcile Table_S4: Identified tRNA Genes with Organellar Origin Verification (Legacy + Wuyi Pan-MT Check)
scaf_orig = {}
with open(os.path.join(eval_dir, "scaffold_origin_classification.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        scaf_orig[r["Scaffold_ID"].split("_length_")[0]] = (r["Scaffold_ID"], r["Origin_Class"], r["MT_Pident"], r["CP_Pident"])

trna_rows = []
with open(os.path.join(supp_dir, "Table_S4_Identified_tRNA_Genes.csv")) as f:
    trna_orig = list(csv.DictReader(f))
rows_s4 = [["tRNA_Isotype", "Anticodon", "Legacy_Scaffold", "Full_Scaffold_ID", "Begin_bp", "End_bp", "Length", "Intron", "Infernal_Score", "Organellar_Origin", "Pan_Camellia_MT_Confirmation"]]
for r in trna_orig:
    sc = r.get("Scaffold") or r.get("Legacy_Scaffold")
    info = scaf_orig.get(sc, (sc, "Unknown", "0", "0"))
    if sc in ("NODE_61", "NODE_89"):
        pan_note = "100% Exact Match to C. sinensis cv. Wuyi (PZ267020.1) & Xinyang10 Complete Mitogenomes (Absent in NC_043914.1)"
        org_lbl = "True Mitochondrial (Missing in RefSeq)"
    else:
        pan_note = f"Legacy scaffold origin: {info[1]} (cpDNA id={info[3]}%, mtDNA id={info[2]}%)"
        org_lbl = info[1]
    rows_s4.append([r["tRNA_Isotype"], r["Anticodon"], sc, info[0], r["Begin_bp"], r["End_bp"], r["Length"], r["Intron"], r["Infernal_Score"], org_lbl, pan_note])

with open(os.path.join(supp_dir, "Table_S4_Identified_tRNA_Genes.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows_s4)
print("Updated Table_S4_Identified_tRNA_Genes.csv")

# 4. Create Table_S9: Curated Variants with Feature, Strand Concordance, Codon Position & AA Impact
shutil.copy2(os.path.join(var_dir, "curated_variants_feature_annotated.csv"), os.path.join(supp_dir, "Table_S9_RNA_Editing_Sites_and_Strand_Bias.csv"))
print("Created Table_S9_RNA_Editing_Sites_and_Strand_Bias.csv")

# 5. Create Table_S10 (Merged Ka/Ks + Water Stress Response) and Table_S10b (RSCU)
str_map = {}
with open(os.path.join(clean_dir, "drought_stress_mito_expression_prjna539837.csv")) as f:
    for r in csv.DictReader(f):
        str_map[r["Gene"]] = r
kaks_map = {}
with open(os.path.join(clean_dir, "srr9057080_kaks_and_editing.csv")) as f:
    for r in csv.DictReader(f):
        kaks_map[r["Gene"]] = r

all_k = sorted(set(str_map.keys()) | set(kaks_map.keys()))
rows_s10 = [[
    "Gene", "Aligned_Codons", "Sd", "Nd", "C_to_U_Sites", "Ka", "Ks", "Apparent_Ka_Ks",
    "UPASI9_Stress_RPMD", "UPASI9_Ctrl_RPMD", "UPASI9_Log2FC", "TRF1_Stress_RPMD", "TRF1_Ctrl_RPMD", "TRF1_Log2FC"
]]
for g in all_k:
    km = kaks_map.get(g, {})
    sm = str_map.get(g, {})
    rows_s10.append([
        g, km.get("Aligned_Codons", "NA"), km.get("Sd", "NA"), km.get("Nd", "NA"), km.get("C_to_U_Sites", "NA"),
        km.get("Ka", "NA"), km.get("Ks", "NA"), km.get("Ka_Ks", "NA"),
        sm.get("UPASI9_Stress_RPMD", "NA"), sm.get("UPASI9_Ctrl_RPMD", "NA"), sm.get("UPASI9_Log2FC", "NA"),
        sm.get("TRF1_Stress_RPMD", "NA"), sm.get("TRF1_Ctrl_RPMD", "NA"), sm.get("TRF1_Log2FC", "NA")
    ])
with open(os.path.join(supp_dir, "Table_S10_KaKs_RSCU_and_Stress_Expression.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows_s10)
shutil.copy2(os.path.join(clean_dir, "srr9057080_mito_rscu.csv"), os.path.join(supp_dir, "Table_S10b_RSCU_Codon_Usage.csv"))
print("Created Table_S10 & Table_S10b")

# 6. Create Table_S11: ATP1 WT vs Edited Site-Specific & Global Structural Summary
AA3_TO_1 = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"}
rows_s11 = [["Model_or_Site", "Mutation_Code", "WT_pLDDT", "Edited_pLDDT", "Delta_pLDDT", "WT_KD_Hydropathy", "Edited_KD_Hydropathy", "Delta_KD", "Calpha_Displacement_A", "WT_Neighbor_Contacts_8A", "Edited_Neighbor_Contacts_8A"]]
with open(os.path.join(prot_dir, "atp1_wt_vs_edited_structural_comparison.csv")) as f:
    for r in csv.DictReader(f):
        code = f"{AA3_TO_1.get(r['WT_AA'], r['WT_AA'][0])}{r['Residue']}{AA3_TO_1.get(r['Edited_AA'], r['Edited_AA'][0])}"
        rows_s11.append([
            f"ATP1 Site {r['Residue']}", code, r["WT_pLDDT"], r["Edited_pLDDT"], r["Delta_pLDDT"],
            r["WT_Hydropathy_KD"], r["Edited_Hydropathy_KD"], r["Delta_KD"], r["CA_Shift_A"], r["Contacts_8A_WT"], r["Contacts_8A_Edited"]
        ])
rows_s11.append(["ATP1 Global (509 aa)", "8 C-terminal edits (347-497)", "91.18 (mean)", "95.02 (mean)", "+3.84", "Global WT", "+25.9 cumulative gain", "+25.9", "1.251 A (C-term 1.953 A)", "91.7% >=90 pLDDT (Edited)", "ESMFold v1"])
rows_s11.append(["COX1 De Novo (527 aa)", "Newly recovered (Absent in NC_043914.1)", "N/A", "90.75 (mean, median 99.0)", "N/A", "N/A", "N/A", "N/A", "N/A", "78.7% >=90 pLDDT", "ESMFold v1"])
with open(os.path.join(supp_dir, "Table_S11_ATP1_Structural_Comparison.csv"), "w", newline="") as f:
    csv.writer(f).writerows(rows_s11)
print("Created Table_S11_ATP1_Structural_Comparison.csv")

# 7. Archive canonical v2/Q1 scratch scripts cleanly into Analysis_Scripts/
scratch_dir = "/mnt/c/Users/user/.gemini/antigravity/scratch"
copied = 0
for fn in sorted(os.listdir(scratch_dir)):
    if fn.endswith((".py", ".sh")):
        shutil.copy2(os.path.join(scratch_dir, fn), os.path.join(script_arch, fn))
        copied += 1
with open(os.path.join(script_arch, "README.md"), "w") as f:
    f.write("# Canonical Reproducible Pipeline Scripts (SRR9057080 Camellia sinensis Mitochondrial Investigation)\n\n"
            "- `step1_origin_and_features.py`: Competitive BLASTn origin deconvolution of legacy 1,160 scaffolds.\n"
            "- `step2_competitive_mapping_and_mtpt.py`: MTPT detection & competitive dual-organelle BWA-MEM alignment (53.93M reads).\n"
            "- `step3_variant_feature_and_repeats.py`: Strand-specific C-to-U concordance, codon position bias, and dispersed repeat scan.\n"
            "- `step4_deep_annotation_reconciliation.py`: Exon-level feature depth profiling (`nc043914_feature_depth_profile.csv`), rRNA BLAST, and SSR identification.\n"
            "- `step5_gfa_and_cohort_meta.py`: SPAdes & GetOrganelle GFA topology parsing.\n"
            "- `fetch_ncbi_references_and_sra.py` & `fetch_outgroup_and_check_gbs.py`: Windows NCBI E-utilities / SRA RunInfo metadata acquisition.\n"
            "- `step7_unified_pan_camellia_annotation.py` & `check_rem_orfs.py`: Pan-Camellia MT + CP BLASTp resolution of all legacy Class D ORFs.\n"
            "- `step8_drought_stress_and_evolution.py`: 2x2 factorial drought-stress expression & editing buffering analysis (PRJNA539837).\n"
            "- `step9_rscu_kaks_phylogeny.py` & `step9c_complete_phylogeny.py`: 21-gene RSCU, 19-gene Ka/Ks inflation proof, and 9-taxon unmasked vs edit-masked trees.\n"
            "- `step10_purified_mt_bam_metrics.py`: Empirical coverage breadth & non-MTPT core depth (`28.52x` Non-MTPT / `118.88x` total).\n"
            "- `step11_clean_mt_assembly.py`: Xinyang10 + cpDNA competitive read purification & `rnaSPAdes` de novo transcript assembly (`37` MT CDS).\n"
            "- `fold_wt_atp1.py` & `step12_structural_wt_vs_edited.py`: GPU ESMFold WT vs Edited ATP1 & COX1 3D structural modeling.\n"
            "- `step13_reconcile_28cohort_metadata.py`: 28-accession RNA-Seq cohort metadata & Non-MTPT SNP reconciliation.\n"
            "- `build_verified_bibliography.py` & `fix_verified_bibliography.py`: Crossref/PubMed API DOI verification (`references.bib` / `references.ris`).\n"
            "- `generate_q1_figures_v2.py` & `generate_figure3_assembly.py`: Publication-grade Figures 1–9.\n"
            "- `reconcile_supp_tables.py`: Master supplementary tables S1–S11 reconciliation.\n"
            "- `compile_complete_manuscript_v3.py` & audit/reviewer report generators: Final DOCX deliverables.\n")
print(f"Archived {copied} canonical scripts + README.md to {script_arch}")
