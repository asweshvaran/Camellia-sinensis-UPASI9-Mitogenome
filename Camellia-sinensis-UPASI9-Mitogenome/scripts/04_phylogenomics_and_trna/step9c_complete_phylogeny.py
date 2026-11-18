import os
import csv
import random
from collections import Counter
from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"
fna_849 = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_genes.fna"
unified_blast = os.path.join(out_dir, "orfs849_vs_unified_mt_cp.tsv")

def norm_gene(g):
    g = g.strip().lower()
    if g in ("cytb", "cob"): return "cob"
    if g.startswith("nd") and len(g) >= 3 and g[2].isdigit():
        g = "nad" + g[2:]
    if g.startswith("nad4l"): return "nad4l"
    if g.startswith("ccmb"): return "ccmb"
    if g.startswith("ccmc"): return "ccmc"
    if g.startswith("ccmfc"): return "ccmfc"
    if g.startswith("ccmfn"): return "ccmfn"
    return g

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
                if g not in res or len(nt) > len(res[g]):
                    res[g] = nt
    return res

taxa_gbs = [
    ("C_sinensis_RefSeq_NC043914", "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"),
    ("C_sinensis_Xinyang10_PX398764", os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb")),
    ("C_sinensis_Yidan_PX759560", os.path.join(ref_dir, "C_sinensis_Yidan_PX759560.gb")),
    ("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr1_PZ267020.gb")),
    ("C_pitardii_PZ501402", os.path.join(ref_dir, "C_pitardii_PZ501402.gb")),
    ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr1_PZ378386.gb")),
    ("Actinidia_rufa_NC086781", os.path.join(ref_dir, "Outgroup_NC_086781.1.gb")),
    ("Vaccinium_macrocarpon_NC023338", os.path.join(ref_dir, "Outgroup_Vaccinium_NC_023338.1.gb")),
]

extra_refs = {}
for _, gp in taxa_gbs:
    if os.path.exists(gp):
        for r in SeqIO.parse(gp, "genbank"): extra_refs[r.id] = r.seq
for gp2 in [os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb"), os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb"), os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb")]:
    if os.path.exists(gp2):
        for r in SeqIO.parse(gp2, "genbank"): extra_refs[r.id] = r.seq

taxa_cds = {}
for name, gp in taxa_gbs:
    taxa_cds[name] = extract_cds_norm(gp, extra_refs)
for tname, gp2 in [("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb")),
                   ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb")),
                   ("C_sinensis_Xinyang10_PX398764", os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb"))]:
    if os.path.exists(gp2):
        taxa_cds[tname].update(extract_cds_norm(gp2, extra_refs))

taxa_cds["C_sinensis_UPASI9_cDNA_raw"] = {g: info["nt_seq"] for g, info in srr_mt_genes.items()}

# Find genes shared across all 9 taxa (or 8 excluding NC043914 gaps? Wait - see common across all 9!)
common = sorted(set.intersection(*(set(m.keys()) for m in taxa_cds.values())))
print(f"Core genes shared across ALL {len(taxa_cds)} taxa ({len(common)}): {common}")

# Build alignment relative to Xinyang10
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1

used_genes = []
concat_raw = {t: [] for t in taxa_cds}
# Also construct RNA-editing-masked UPASI-9 (revert C->T differences against Camellia DNA consensus back to C so only true non-C->T genomic differences remain!)
concat_masked = {t: [] for t in taxa_cds}

for g in common:
    ref_s = taxa_cds["C_sinensis_Xinyang10_PX398764"][g]
    lens = [len(taxa_cds[t][g]) for t in taxa_cds]
    if min(lens) < len(ref_s) * 0.80:
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
    keep_cols = [i for i in range(len(ref_s)) if all(tax_mapped[t][i] != "-" for t in taxa_cds)]
    if len(keep_cols) >= 150:
        used_genes.append((g, len(keep_cols)))
        for t in taxa_cds:
            concat_raw[t].append("".join(tax_mapped[t][i] for i in keep_cols))
            if t == "C_sinensis_UPASI9_cDNA_raw":
                # Mask C->T RNA editing: if all DNA C. sinensis refs have 'C' and UPASI9 has 'T', revert to 'C'
                masked_bases = []
                for i in keep_cols:
                    b_up = tax_mapped[t][i]
                    b_xy = tax_mapped["C_sinensis_Xinyang10_PX398764"][i]
                    b_yd = tax_mapped["C_sinensis_Yidan_PX759560"][i]
                    if b_up == "T" and b_xy == "C" and b_yd == "C":
                        masked_bases.append("C")
                    else:
                        masked_bases.append(b_up)
                concat_masked[t].append("".join(masked_bases))
            else:
                concat_masked[t].append("".join(tax_mapped[t][i] for i in keep_cols))

print(f"Concatenated genes ({len(used_genes)}): {used_genes}")
aln_len = sum(l for _, l in used_genes)
print(f"Total gap-free concatenated alignment length: {aln_len:,} bp")

def build_bootstrapped_tree(concat_map, upasi_label, out_nwk):
    recs = [SeqRecord(Seq("".join(concat_map[t])), id=(upasi_label if "UPASI9" in t else t)) for t in sorted(concat_map)]
    msa = MultipleSeqAlignment(recs)
    calc = DistanceCalculator("identity")
    dm = calc.get_distance(msa)
    constructor = DistanceTreeConstructor(calc, "nj")
    tree = constructor.build_tree(msa)

    random.seed(42)
    n_taxa = len(msa)
    seqs_str = [str(r.seq) for r in msa]
    taxon_ids = [r.id for r in msa]
    def get_splits(tr):
        all_l = set(c.name for c in tr.get_terminals())
        sp = set()
        for clade in tr.get_nonterminals():
            sub = frozenset(c.name for c in clade.get_terminals())
            if 1 < len(sub) < len(all_l):
                comp = frozenset(all_l - sub)
                sp.add(min(sub, comp, key=lambda s: tuple(sorted(s))))
        return sp
    base_splits = get_splits(tree)
    sp_cnt = Counter()
    for _ in range(500):
        cols = [random.randrange(aln_len) for _ in range(aln_len)]
        brecs = [SeqRecord(Seq("".join(seqs_str[ti][c] for c in cols)), id=taxon_ids[ti]) for ti in range(n_taxa)]
        bt = constructor.build_tree(MultipleSeqAlignment(brecs))
        for sp in get_splits(bt):
            if sp in base_splits:
                sp_cnt[sp] += 1
    all_l = set(taxon_ids)
    for clade in tree.get_nonterminals():
        sub = frozenset(c.name for c in clade.get_terminals())
        if 1 < len(sub) < len(all_l):
            comp = frozenset(all_l - sub)
            key = min(sub, comp, key=lambda s: tuple(sorted(s)))
            clade.confidence = int(round(sp_cnt[key] / 500.0 * 100))
    for t in tree.get_terminals():
        if "Vaccinium" in t.name:
            tree.root_with_outgroup(t)
            break
    Phylo.write(tree, out_nwk, "newick")
    return dm, tree

print("\n=== TREE A: UNMASKED cDNA (Demonstrating RNA-Editing Branch Artifact) ===")
dm_raw, t_raw = build_bootstrapped_tree(concat_raw, "C_sinensis_UPASI9_cDNA_Raw", os.path.join(out_dir, "phylogeny_unmasked_cdna.nwk"))
print(dm_raw)
Phylo.draw_ascii(t_raw)

print("\n=== TREE B: RNA-EDITING-MASKED GENOMIC PHYLOGENY ===")
dm_msk, t_msk = build_bootstrapped_tree(concat_masked, "C_sinensis_UPASI9_EditMasked", os.path.join(out_dir, "phylogeny_masked_genomic.nwk"))
print(dm_msk)
Phylo.draw_ascii(t_msk)
