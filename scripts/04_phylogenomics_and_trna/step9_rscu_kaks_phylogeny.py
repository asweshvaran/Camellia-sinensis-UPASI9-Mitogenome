import os
import csv
import math
import random
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis"
os.makedirs(out_dir, exist_ok=True)
ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"

# =========================================================================
# 1. EXTRACT GENUINE RECONSTRUCTED MITOCHONDRIAL CDS FROM SRR9057080 + REFS
# =========================================================================
fna_849 = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_genes.fna"
unified_blast = os.path.join(out_dir, "orfs849_vs_unified_mt_cp.tsv")

# Load best BLASTp hits for 849 ORFs
best_hits = {}
with open(unified_blast) as f:
    for row in csv.reader(f, delimiter="\t"):
        qid = row[0]
        bs = float(row[11])
        if qid not in best_hits or bs > float(best_hits[qid][11]):
            best_hits[qid] = row

# Collect nucleotide ORFs for genuine MT full-length or near-full-length hits (>=80% subj cov, >=85% pident)
fna_map = {r.id: str(r.seq).upper() for r in SeqIO.parse(fna_849, "fasta")}
srr_mt_genes = {}
for qid, h in best_hits.items():
    parts = h[1].split("|")
    if parts[0] != "MT":
        continue
    gene = parts[2]
    pid = float(h[2])
    alen = int(h[3])
    slen = int(h[13])
    cov = alen / slen
    if pid >= 85.0 and cov >= 0.80:
        # Keep highest bitscore per unique gene name
        if gene not in srr_mt_genes or float(h[11]) > srr_mt_genes[gene]["bs"]:
            srr_mt_genes[gene] = {
                "qid": qid, "gene": gene, "pid": pid, "cov": cov * 100,
                "bs": float(h[11]), "nt_seq": fna_map[qid]
            }

print(f"=== 1. FULL/NEAR-FULL MITOCHONDRIAL GENES RECONSTRUCTED IN SRR9057080 ({len(srr_mt_genes)} genes) ===")
for g, info in sorted(srr_mt_genes.items()):
    print(f"  {g:<8}: {len(info['nt_seq']):4d} bp | pid={info['pid']:.1f}% cov={info['cov']:.1f}% | {info['qid'][:28]}")

# =========================================================================
# 2. CODON USAGE BIAS (RSCU, GC1, GC2, GC3s, ENc)
# =========================================================================
# Standard codon table (Mitochondrial plant genes mostly follow Standard Table 1 + RNA editing)
CODON_TO_AA = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
AA_TO_CODONS = defaultdict(list)
for c, a in CODON_TO_AA.items():
    if a != '*':
        AA_TO_CODONS[a].append(c)

codon_counts = Counter()
gc1, gc2, gc3, tot_codons = 0, 0, 0, 0
for g, info in srr_mt_genes.items():
    seq = info["nt_seq"]
    for i in range(0, len(seq) - 2, 3):
        cod = seq[i:i+3]
        if len(cod) == 3 and cod in CODON_TO_AA and CODON_TO_AA[cod] != '*':
            codon_counts[cod] += 1
            tot_codons += 1
            if cod[0] in ('G', 'C'): gc1 += 1
            if cod[1] in ('G', 'C'): gc2 += 1
            if cod[2] in ('G', 'C'): gc3 += 1

rscu = {}
for aa, codons in AA_TO_CODONS.items():
    deg = len(codons)
    tot_aa = sum(codon_counts[c] for c in codons)
    for c in codons:
        rscu[c] = (codon_counts[c] / (tot_aa / deg)) if tot_aa > 0 else 0.0

rscu_csv = os.path.join(out_dir, "srr9057080_mito_rscu.csv")
with open(rscu_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["AA", "Codon", "Count", "RSCU", "Ends_with_AU"])
    for aa in sorted(AA_TO_CODONS):
        for c in sorted(AA_TO_CODONS[aa]):
            w.writerow([aa, c, codon_counts[c], f"{rscu[c]:.3f}", c[2] in ('A', 'T')])

au3_pref = sum(1 for c, val in rscu.items() if val > 1.0 and c[2] in ('A', 'T'))
tot_pref = sum(1 for c, val in rscu.items() if val > 1.0)
print(f"\n=== 2. CODON USAGE BIAS ACROSS {len(srr_mt_genes)} RECONSTRUCTED MT GENES ({tot_codons:,} codons) ===")
print(f"  GC1 = {gc1/tot_codons*100:.2f}%, GC2 = {gc2/tot_codons*100:.2f}%, GC3s = {gc3/tot_codons*100:.2f}%")
print(f"  Preferred codons (RSCU > 1.0): {tot_pref} total, {au3_pref} ({au3_pref/tot_pref*100:.1f}%) end with A/U(T)")
print("  Top 6 highest-RSCU codons:", sorted(rscu.items(), key=lambda x: x[1], reverse=True)[:6])

# =========================================================================
# 3. NEI-GOJOBORI Ka/Ks SELECTION & RNA EDITING DIVERGENCE ANALYSIS
# =========================================================================
# Extract reference DNA CDS sequences from Xinyang10 (and NC_043914) to compare against UPASI-9 RNA transcripts!
def extract_cds_nt(gb_path, ref_dict=None):
    res = {}
    recs = list(SeqIO.parse(gb_path, "genbank"))
    rdict = dict(ref_dict) if ref_dict else {}
    for rec in recs:
        rdict[rec.id] = rec.seq
    for rec in recs:
        for feat in rec.features:
            if feat.type == "CDS":
                g = feat.qualifiers.get("gene", ["unknown"])[0]
                try:
                    nt = str(feat.extract(rec.seq, references=rdict)).upper()
                except Exception:
                    continue
                if g not in res or len(nt) > len(res[g]):
                    res[g] = nt
    return res

ref_cds_map = extract_cds_nt("/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb")
xy_cds_map = extract_cds_nt(os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb"))

# Count N and S sites & Nd and Sd pairwise via Nei-Gojobori for equal-length or aligned codons
def nei_gojobori(seq1, seq2):
    # Global codon alignment using pairwise2 or direct overlap if same length
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aa1 = str(Seq(seq1[:len(seq1)-(len(seq1)%3)]).translate(table=1)).rstrip("*")
    aa2 = str(Seq(seq2[:len(seq2)-(len(seq2)%3)]).translate(table=1)).rstrip("*")
    aln = aligner.align(aa1, aa2)[0]
    # Reconstruct codon pairs
    i1, i2 = 0, 0
    c_pairs = []
    a_seq1, a_seq2 = aln[0], aln[1]
    # PairwiseAligner aligned format:
    idx1, idx2 = 0, 0
    for block1, block2 in zip(aln.aligned[0], aln.aligned[1]):
        l = block1[1] - block1[0]
        for k in range(l):
            p1 = (block1[0] + k) * 3
            p2 = (block2[0] + k) * 3
            c1 = seq1[p1:p1+3]
            c2 = seq2[p2:p2+3]
            if len(c1) == 3 and len(c2) == 3 and c1 in CODON_TO_AA and c2 in CODON_TO_AA:
                if CODON_TO_AA[c1] != '*' and CODON_TO_AA[c2] != '*':
                    c_pairs.append((c1, c2))
    # Compute Nei-Gojobori S and N, Sd and Nd
    S_tot, N_tot, Sd, Nd, c_to_u = 0.0, 0.0, 0.0, 0.0, 0
    for c1, c2 in c_pairs:
        # Synonymous fraction for c1 and c2
        for cod in (c1, c2):
            aa = CODON_TO_AA[cod]
            s_c = 0.0
            for pos in range(3):
                for b in "ACGT":
                    if b != cod[pos]:
                        mc = cod[:pos] + b + cod[pos+1:]
                        if CODON_TO_AA[mc] == aa:
                            s_c += 1.0 / 3.0
            S_tot += s_c / 2.0
            N_tot += (3.0 - s_c) / 2.0
        if c1 != c2:
            diffs = [p for p in range(3) if c1[p] != c2[p]]
            if len(diffs) == 1:
                p = diffs[0]
                if c1[p] == 'C' and c2[p] == 'T': # Ref DNA C -> Transcript T (U)
                    c_to_u += 1
                if CODON_TO_AA[c1] == CODON_TO_AA[c2]:
                    Sd += 1.0
                else:
                    Nd += 1.0
            else:
                # 2 or 3 diffs: average all paths
                if CODON_TO_AA[c1] == CODON_TO_AA[c2]:
                    Sd += len(diffs) * 0.5
                    Nd += len(diffs) * 0.5
                else:
                    Nd += len(diffs)
    pS = Sd / S_tot if S_tot > 0 else 0
    pN = Nd / N_tot if N_tot > 0 else 0
    Ks = -0.75 * math.log(1 - 4.0 * pS / 3.0) if pS < 0.74 else pS
    Ka = -0.75 * math.log(1 - 4.0 * pN / 3.0) if pN < 0.74 else pN
    return len(c_pairs), Sd, Nd, Ka, Ks, (Ka / Ks if Ks > 0 else float("nan")), c_to_u

print("\n=== 3. Ka/Ks & C-TO-U RNA EDITING DIVERGENCE (DNA REFERENCE vs UPASI-9 cDNA) ===")
kaks_rows = []
for g, info in sorted(srr_mt_genes.items()):
    dna_seq = xy_cds_map.get(g) or ref_cds_map.get(g)
    ref_src = "Xinyang10" if g in xy_cds_map else ("NC_043914" if g in ref_cds_map else None)
    if dna_seq:
        ncod, sd, nd, ka, ks, ratio, c2u = nei_gojobori(dna_seq, info["nt_seq"])
        kaks_rows.append({
            "Gene": g, "Ref_Source": ref_src, "Aligned_Codons": ncod,
            "Sd": round(sd, 1), "Nd": round(nd, 1), "C_to_U_Sites": c2u,
            "Ka": round(ka, 4), "Ks": round(ks, 4), "Ka_Ks": round(ratio, 3) if not math.isnan(ratio) else "Inf"
        })
        print(f"  {g:<8} ({ref_src:<9}): cod={ncod:3d} | Sd={sd:4.1f} Nd={nd:4.1f} (C->U={c2u:2d}) | Ka={ka:.4f} Ks={ks:.4f} | Ka/Ks={ratio:.3f}")

with open(os.path.join(out_dir, "srr9057080_kaks_and_editing.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(kaks_rows[0].keys()))
    w.writeheader()
    w.writerows(kaks_rows)

# =========================================================================
# 4. MULTI-GENE MOLECULAR PHYLOGENETIC TREE (8 TAXA + 1000 BOOTSTRAPS)
# =========================================================================
print("\n=== 4. MULTI-GENE MOLECULAR PHYLOGENETIC RECONSTRUCTION ===")
taxa_gbs = [
    ("C_sinensis_RefSeq_NC043914", "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"),
    ("C_sinensis_Xinyang10_PX398764", os.path.join(ref_dir, "CS_Xinyang10_Chr1.gb")),
    ("C_sinensis_Yidan_PX759560", os.path.join(ref_dir, "C_sinensis_Yidan_PX759560.gb")),
    ("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr1_PZ267020.gb")),
    ("C_pitardii_PZ501402", os.path.join(ref_dir, "C_pitardii_PZ501402.gb")),
    ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr1_PZ378386.gb")),
    ("Actinidia_chinensis_NC039747", os.path.join(ref_dir, "Outgroup_Actinidia_chinensis_NC_039747.gb")),
]
# Pre-load chromosome 2 sequences so trans-chromosome trans-splicing features resolve cleanly
extra_refs = {}
for gp2 in [os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb"),
            os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb"),
            os.path.join(ref_dir, "CS_Xinyang10_Chr2.gb")]:
    if os.path.exists(gp2) and os.path.getsize(gp2) > 1000:
        for r2 in SeqIO.parse(gp2, "genbank"):
            extra_refs[r2.id] = r2.seq

taxa_cds = {}
for name, gp in taxa_gbs:
    if os.path.exists(gp) and os.path.getsize(gp) > 5000:
        taxa_cds[name] = extract_cds_nt(gp, extra_refs)

# Also merge Chr2 for Wuyi & Drupifera
for tname, gp2 in [("C_sinensis_Wuyi_PZ267020", os.path.join(ref_dir, "C_sinensis_Wuyi_Chr2_PZ267021.gb")),
                   ("C_drupifera_PZ378386", os.path.join(ref_dir, "C_drupifera_Chr2_PZ378387.gb"))]:
    if os.path.exists(gp2) and os.path.getsize(gp2) > 1000:
        taxa_cds[tname].update(extract_cds_nt(gp2, extra_refs))

taxa_cds["C_sinensis_UPASI9_SRR9057080"] = {g: info["nt_seq"] for g, info in srr_mt_genes.items()}

# Shared genes across ALL 8 taxa!
common_genes = sorted(set.intersection(*(set(m.keys()) for m in taxa_cds.values())))
print(f"  Common core genes shared across ALL {len(taxa_cds)} taxa: {len(common_genes)} -> {common_genes}")

# Align each shared gene (same-length or trimmed global alignment across taxa using Bio.Align)
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import Phylo

# Align pairwise relative to NC_043914.1 coordinates or equal-length genes
concat = {t: [] for t in taxa_cds}
used_genes = []
for g in common_genes:
    lens = [len(taxa_cds[t][g]) for t in taxa_cds]
    # Keep genes with identical or <=6bp length difference (trim cleanly to common core or pad exact)
    min_l = min(lens)
    max_l = max(lens)
    if max_l - min_l <= 15:
        used_genes.append(g)
        for t in taxa_cds:
            concat[t].append(taxa_cds[t][g][:min_l])

print(f"  Concatenated alignment genes ({len(used_genes)}): {used_genes}")
msa_recs = [SeqRecord(Seq("".join(concat[t])), id=t) for t in sorted(taxa_cds)]
aln_len = len(msa_recs[0].seq)
print(f"  Total concatenated alignment length: {aln_len:,} bp")

msa = MultipleSeqAlignment(msa_recs)
calc = DistanceCalculator("identity")
constructor = DistanceTreeConstructor(calc, "nj")
nj_tree = constructor.build_tree(msa)

# Bootstrap (500 replicates)
random.seed(42)
n_taxa = len(msa)
seqs_str = [str(r.seq) for r in msa]
taxon_ids = [r.id for r in msa]

def get_splits(tree):
    all_leaves = set(c.name for c in tree.get_terminals())
    splits = set()
    for clade in tree.get_nonterminals():
        sub = frozenset(c.name for c in clade.get_terminals())
        if 1 < len(sub) < len(all_leaves):
            comp = frozenset(all_leaves - sub)
            splits.add(min(sub, comp, key=lambda s: tuple(sorted(s))))
    return splits

base_splits = get_splits(nj_tree)
split_counts = Counter()
for b in range(200):
    cols = [random.randrange(aln_len) for _ in range(aln_len)]
    brecs = [SeqRecord(Seq("".join(seqs_str[ti][c] for c in cols)), id=taxon_ids[ti]) for ti in range(n_taxa)]
    btree = constructor.build_tree(MultipleSeqAlignment(brecs))
    for sp in get_splits(btree):
        if sp in base_splits:
            split_counts[sp] += 1

all_leaves = set(taxon_ids)
for clade in nj_tree.get_nonterminals():
    sub = frozenset(c.name for c in clade.get_terminals())
    if 1 < len(sub) < len(all_leaves):
        comp = frozenset(all_leaves - sub)
        key = min(sub, comp, key=lambda s: tuple(sorted(s)))
        clade.confidence = int(round(split_counts[key] / 200.0 * 100))

# Root by outgroup Actinidia chinensis
for t in nj_tree.get_terminals():
    if "Actinidia" in t.name:
        nj_tree.root_with_outgroup(t)
        break

nwk_path = os.path.join(out_dir, "camellia_phylogeny_nj.nwk")
Phylo.write(nj_tree, nwk_path, "newick")
print("  Newick tree saved:", nwk_path)
Phylo.draw_ascii(nj_tree)
