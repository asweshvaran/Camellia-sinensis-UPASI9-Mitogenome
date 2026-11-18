import os
import csv
import subprocess
import re
from collections import Counter, defaultdict
from Bio import SeqIO

work_dir = "/mnt/d/Mitochondria/Q1_Investigation"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
bam_path = "/mnt/d/Mitochondria/analysis/SRR9057080/mapping/SRR9057080_aligned.bam"
origin_tsv = os.path.join(work_dir, "Assembly_Evaluation/scaffold_origin_classification.tsv")
orf_tsv = os.path.join(work_dir, "ORF_Validation/orf_curated_classification.tsv")
mtpt_tsv = os.path.join(work_dir, "Clean_Reanalysis/nc043914_vs_cp_mtpt.tsv")

# 1. Reconcile Scaffold Origins with Class A/B/C/D/E ORFs and tRNAs!
scaf_origin = {}
with open(origin_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        scaf_origin[row["Scaffold_ID"]] = row

print("=== 1. CROSS-CHECK CLASS A/B/C/D/E ORFs AGAINST SCAFFOLD ORIGIN ===")
by_cls_orig = Counter()
with open(orf_tsv) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        oid = row["ORF_ID"]
        # ORF_ID is NODE_X_length_Y_cov_Z_N
        scaf_id = "_".join(oid.split("_")[:-1])
        oc = scaf_origin.get(scaf_id, {}).get("Origin_Class", "Unknown")
        cls = row["Class"]
        by_cls_orig[(cls, oc)] += 1
        if cls in ("A", "C", "D"):
            print(f"  [{cls}] {oid:<42} len={row['Length_aa']:>3}aa target={row['Homology_Target']:<6} | Origin={oc}")

print("\n  Summary Class x Origin:")
for (c, o), cnt in sorted(by_cls_orig.items()):
    print(f"    Class {c} x {o:<30}: {cnt}")

# 2. Check barrnap --kingdom euk / mito and direct rRNA BLAST against NC_043914.1 rRNAs
print("\n=== 2. FIX rRNA ANNOTATION (barrnap + direct NC_043914.1 rRNA BLAST) ===")
record = next(SeqIO.parse(mt_gb, "genbank"))
mt_rrna_fasta = os.path.join(work_dir, "Annotation/nc043914_rrnas.fasta")
rrna_recs = []
for feat in record.features:
    if feat.type == "rRNA":
        gname = feat.qualifiers.get("gene", feat.qualifiers.get("product", ["rRNA"]))[0]
        seq = feat.extract(record.seq)
        from Bio.SeqRecord import SeqRecord
        rrna_recs.append(SeqRecord(seq, id=gname, description=""))
SeqIO.write(rrna_recs, mt_rrna_fasta, "fasta")
scaf_fasta = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
rrna_blast = os.path.join(work_dir, "Annotation/scaffolds_vs_mt_rrna.tsv")
subprocess.run(["makeblastdb", "-in", mt_rrna_fasta, "-dbtype", "nucl", "-out", mt_rrna_fasta], check=True)
subprocess.run(["blastn", "-query", scaf_fasta, "-db", mt_rrna_fasta, "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore", "-evalue", "1e-5", "-out", rrna_blast], check=True)
with open(rrna_blast) as f:
    rows = list(csv.reader(f, delimiter="\t"))
print(f"  Scaffolds hitting NC_043914.1 annotated rrnS/rrnL: {len(rows)} alignments")
for row in rows[:12]:
    print(f"    {row[0][:30]:<30} -> {row[1]} | pid={row[2]}% len={row[3]} (qlen={row[4]}, slen={row[5]}) | {scaf_origin.get(row[0],{}).get('Origin_Class','?')}")

# 3. Gene-by-gene depth profile across NC_043914.1 (all 32 CDS, 2 rRNA, 23 tRNA)
print("\n=== 3. GENE-BY-GENE DEPTH, BREADTH & TRANSCRIPT EXPRESSION ON NC_043914.1 ===")
# Load MTPT mask
mtpt_pos = set()
with open(mtpt_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        if int(row[3]) >= 50:
            qs, qe = min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))
            for p in range(qs, qe + 1):
                mtpt_pos.add(p)

# Load exact base depth array (1-based index -> depth_arr[pos])
ref_len = len(record.seq)
depth_arr = [0] * (ref_len + 1)
proc = subprocess.Popen(["samtools", "depth", "-a", bam_path], stdout=subprocess.PIPE, text=True)
for line in proc.stdout:
    parts = line.split()
    if len(parts) >= 3:
        depth_arr[int(parts[1])] = int(parts[2])
proc.wait()

# Feature depth summary (handle duplicate gene copies separately with index)
feat_rows = []
gene_seen = Counter()
for feat in record.features:
    if feat.type in ("CDS", "rRNA", "tRNA"):
        name = feat.qualifiers.get("gene", feat.qualifiers.get("product", [feat.type]))[0]
        gene_seen[(feat.type, name)] += 1
        label = name if gene_seen[(feat.type, name)] == 1 else f"{name}_copy{gene_seen[(feat.type, name)]}"
        coords = []
        for part in feat.location.parts:
            coords.extend(range(int(part.start) + 1, int(part.end) + 1))
        coords = sorted(set(coords))
        L = len(coords)
        vals = [depth_arr[p] for p in coords]
        mean_d = sum(vals) / L if L else 0
        svals = sorted(vals)
        med_d = svals[L // 2] if L else 0
        max_d = max(vals) if vals else 0
        b1 = sum(1 for v in vals if v >= 1) / L * 100
        b10 = sum(1 for v in vals if v >= 10) / L * 100
        b30 = sum(1 for v in vals if v >= 30) / L * 100
        mtpt_frac = sum(1 for p in coords if p in mtpt_pos) / L * 100
        strand = "+" if feat.location.strand == 1 else "-"
        feat_rows.append({
            "Feature_ID": label, "Gene": name, "Type": feat.type, "Strand": strand,
            "Start": int(feat.location.start) + 1, "End": int(feat.location.end),
            "Exonic_Len_bp": L, "Mean_Depth": round(mean_d, 2), "Median_Depth": med_d,
            "Max_Depth": max_d, "Breadth_1x_pct": round(b1, 2), "Breadth_10x_pct": round(b10, 2),
            "Breadth_30x_pct": round(b30, 2), "MTPT_Overlap_pct": round(mtpt_frac, 1)
        })

out_gene_cov = os.path.join(work_dir, "Annotation/nc043914_feature_depth_profile.csv")
with open(out_gene_cov, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(feat_rows[0].keys()))
    w.writeheader()
    w.writerows(feat_rows)

print(f"  Saved {len(feat_rows)} annotated features to {out_gene_cov}")
print("  Top expressed CDS genes (>= 10x mean depth):")
for r in sorted([x for x in feat_rows if x["Type"] == "CDS"], key=lambda x: x["Mean_Depth"], reverse=True):
    print(f"    {r['Feature_ID']:<12} ({r['Strand']}) len={r['Exonic_Len_bp']:4d}bp | mean={r['Mean_Depth']:7.1f}x med={r['Median_Depth']:5d}x | >=10x:{r['Breadth_10x_pct']:5.1f}% | MTPT:{r['MTPT_Overlap_pct']:4.1f}%")
print("  rRNA & tRNA features:")
for r in [x for x in feat_rows if x["Type"] in ("rRNA", "tRNA")]:
    if r["Mean_Depth"] >= 5 or r["Type"] == "rRNA":
        print(f"    [{r['Type']}] {r['Feature_ID']:<14} ({r['Strand']}) len={r['Exonic_Len_bp']:4d}bp | mean={r['Mean_Depth']:8.1f}x | >=10x:{r['Breadth_10x_pct']:5.1f}% | MTPT:{r['MTPT_Overlap_pct']:4.1f}%")

# 4. SSR (Microsatellite) scan on NC_043914.1
print("\n=== 4. SSR / MICROSATELLITE SCAN ON NC_043914.1 ===")
seq_str = str(record.seq).upper()
ssr_counts = Counter()
# Standard thresholds: mono>=10, di>=6, tri>=5, tetra>=4, penta>=3, hexa>=3
patterns = [
    ("Mononucleotide (>=10x)", r"([ACGT])\1{9,}"),
    ("Dinucleotide (>=6x)", r"((?!AA|CC|GG|TT)[ACGT]{2})\1{5,}"),
    ("Trinucleotide (>=5x)", r"((?!([ACGT])\2\2)[ACGT]{3})\1{4,}"),
    ("Tetranucleotide (>=4x)", r"([ACGT]{4})\1{3,}"),
    ("Pentanucleotide (>=3x)", r"([ACGT]{5})\1{2,}"),
    ("Hexanucleotide (>=3x)", r"([ACGT]{6})\1{2,}"),
]
occupied = [False] * len(seq_str)
motif_top = Counter()
for label, pat in patterns:
    cnt = 0
    for m in re.finditer(pat, seq_str):
        s, e = m.start(), m.end()
        if any(occupied[s:e]):
            continue
        for k in range(s, e):
            occupied[k] = True
        cnt += 1
        motif_top[m.group(1)] += 1
    ssr_counts[label] = cnt
print(f"  Total non-overlapping SSRs: {sum(ssr_counts.values())}")
for k, v in ssr_counts.items():
    print(f"    {k:<24}: {v}")
print("  Top 8 SSR motifs:", motif_top.most_common(8))
