import os
import csv
import subprocess
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

work_dir = "/mnt/d/Mitochondria/Q1_Investigation"
mt_gb = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.gb"
mt_fasta = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
var_csv = os.path.join(work_dir, "Heteroplasmy_Variants/curated_high_confidence_variants.csv")

record = next(SeqIO.parse(mt_gb, "genbank"))
ref_seq = record.seq
ref_len = len(ref_seq)
print(f"Loaded GenBank {record.id}: {ref_len:,} bp")

# 1. Extract all annotated features (CDS, rRNA, tRNA, exon, intron)
cds_list = []
rna_list = []
gene_features = []

for feat in record.features:
    if feat.type in ("CDS", "rRNA", "tRNA"):
        name = feat.qualifiers.get("gene", feat.qualifiers.get("product", [feat.type]))[0]
        strand = feat.location.strand # 1 or -1
        parts = [(int(p.start) + 1, int(p.end), p.strand) for p in feat.location.parts]
        item = {
            "type": feat.type,
            "name": name,
            "strand": strand,
            "start": int(feat.location.start) + 1,
            "end": int(feat.location.end),
            "parts": parts,
            "feat": feat
        }
        gene_features.append(item)
        if feat.type == "CDS":
            cds_list.append(item)
        else:
            rna_list.append(item)

print(f"Annotated features in NC_043914.1: CDS={len(cds_list)}, rRNA={sum(1 for x in rna_list if x['type']=='rRNA')}, tRNA={sum(1 for x in rna_list if x['type']=='tRNA')}")
for x in rna_list:
    if x["type"] == "rRNA":
        print(f"  rRNA {x['name']}: {x['start']}-{x['end']} (len={x['end']-x['start']+1}, strand={x['strand']})")

# Build lookup of 1-based ref pos -> coding info (CDS gene, strand, codon_pos, aa_ref, codon_ref)
pos_to_cds = {}
pos_to_rna = {}

for item in rna_list:
    for s, e, st in item["parts"]:
        for p in range(s, e + 1):
            pos_to_rna[p] = (item["type"], item["name"], item["strand"])

for item in cds_list:
    # Build ordered list of coding coordinates matching spliced transcript 5'->3'
    strand = item["strand"]
    coords = []
    # Note: feat.location.parts are already ordered 5'->3' in GenBank CompoundLocation? Verify!
    # In BioPython, CompoundLocation parts for -1 strand are stored 5'->3' along transcript!
    # Let's extract carefully using feat.extract and tracking coordinates explicitly:
    ordered_coords = []
    for part in item["feat"].location.parts:
        ps = int(part.start) + 1
        pe = int(part.end)
        pst = part.strand
        if pst == 1:
            ordered_coords.extend(range(ps, pe + 1))
        else:
            ordered_coords.extend(range(pe, ps - 1, -1))
    # Verify reconstructed CDS matches feat.extract(ref_seq)
    rebuilt = "".join(str(ref_seq[c - 1]) if strand == 1 else str(Seq(str(ref_seq[c - 1])).complement()) for c in ordered_coords)
    extracted = str(item["feat"].extract(ref_seq))
    # Wait: trans-spliced or mixed-strand parts might exist! Use per-part strand!
    ordered_bases = []
    ordered_info = []
    for part in item["feat"].location.parts:
        ps = int(part.start) + 1
        pe = int(part.end)
        pst = part.strand
        rng = range(ps, pe + 1) if pst == 1 else range(pe, ps - 1, -1)
        for c in rng:
            b = str(ref_seq[c - 1]).upper()
            if pst == -1:
                b = str(Seq(b).complement())
            ordered_bases.append(b)
            ordered_info.append((c, pst))
    rebuilt2 = "".join(ordered_bases)
    if rebuilt2 != extracted:
        print(f"  WARN mismatch on {item['name']}: len {len(rebuilt2)} vs {len(extracted)}")
    codon_start = int(item["feat"].qualifiers.get("codon_start", ["1"])[0]) - 1
    for idx in range(codon_start, len(ordered_info)):
        c_pos, pst = ordered_info[idx]
        aa_idx = (idx - codon_start) // 3
        c_frame = (idx - codon_start) % 3 # 0=1st, 1=2nd, 2=3rd
        triplet_start = codon_start + aa_idx * 3
        if triplet_start + 3 <= len(rebuilt2):
            codon = rebuilt2[triplet_start:triplet_start + 3]
        else:
            codon = "NNN"
        pos_to_cds[c_pos] = {
            "gene": item["name"],
            "strand": pst,
            "codon_idx": aa_idx + 1,
            "codon_pos": c_frame + 1,
            "codon": codon
        }

# 2. Annotate all 993 curated variants!
sub_by_region = Counter()
cds_strand_subs = Counter()
codon_pos_edits = Counter()
aa_changes = Counter()
annotated_vars = []

with open(var_csv) as f:
    for row in csv.DictReader(f):
        pos = int(row["pos"])
        ref, alt = row["ref"], row["alt"]
        vaf = float(row["vaf"])
        cat = row["category"]
        is_snp = row["is_snp"] == "True"
        pair = f"{ref}>{alt}" if is_snp else "INDEL"

        if pos in pos_to_cds:
            cinfo = pos_to_cds[pos]
            region = "CDS"
            gene = cinfo["gene"]
            strand = "+" if cinfo["strand"] == 1 else "-"
            cp = cinfo["codon_pos"]
            codon_ref = cinfo["codon"]
            aa_change_str = "N/A"
            if is_snp and codon_ref != "NNN":
                # transcript nucleotide at mutated base
                alt_tx = alt if cinfo["strand"] == 1 else str(Seq(alt).complement())
                codon_alt_list = list(codon_ref)
                codon_alt_list[cp - 1] = alt_tx
                codon_alt = "".join(codon_alt_list)
                aa_r = str(Seq(codon_ref).translate(table=1))
                aa_a = str(Seq(codon_alt).translate(table=1))
                aa_change_str = f"{aa_r}->{aa_a} ({codon_ref}->{codon_alt})"
                cds_strand_subs[(strand, pair)] += 1
                if (strand == "+" and pair == "C>T") or (strand == "-" and pair == "G>A"):
                    codon_pos_edits[cp] += 1
                    aa_changes[(aa_r, aa_a, "NonSyn" if aa_r != aa_a else "Syn")] += 1
        elif pos in pos_to_rna:
            rtype, rname, rst = pos_to_rna[pos]
            region = rtype
            gene = rname
            strand = "+" if rst == 1 else "-"
            cp = 0
            aa_change_str = "N/A"
        else:
            region = "Intergenic_or_Intron"
            gene = "Intergenic/Intron"
            strand = "."
            cp = 0
            aa_change_str = "N/A"

        sub_by_region[(region, pair)] += 1
        annotated_vars.append({
            **row, "region": region, "gene": gene, "strand": strand, "codon_pos": cp, "aa_change": aa_change_str
        })

out_ann_var = os.path.join(work_dir, "Heteroplasmy_Variants/curated_variants_feature_annotated.csv")
with open(out_ann_var, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(annotated_vars[0].keys()))
    w.writeheader()
    w.writerows(annotated_vars)

print("\n=== VARIANT REGION DISTRIBUTION (N=993) ===")
reg_tot = Counter(v["region"] for v in annotated_vars)
for k, v in reg_tot.most_common():
    print(f"  {k:<22}: {v:4d} ({v/len(annotated_vars)*100:.1f}%)")

print("\n=== STRAND-SPECIFIC CDS SUBSTITUTION CONCORDANCE (EXACT RNA-EDITING TEST) ===")
for (st, pair), cnt in cds_strand_subs.most_common(12):
    print(f"  CDS strand ({st}) {pair:<6}: {cnt:4d}")

print("\n=== CANONICAL STRAND-CONCORDANT C-to-U EDITS IN CDS BY CODON POSITION ===")
tot_ed = sum(codon_pos_edits.values())
for cp in [1, 2, 3]:
    cnt = codon_pos_edits[cp]
    print(f"  Codon pos {cp}: {cnt:4d} ({cnt/tot_ed*100 if tot_ed else 0:.1f}%)")

syn_cnt = sum(cnt for (r, a, st), cnt in aa_changes.items() if st == "Syn")
nsyn_cnt = sum(cnt for (r, a, st), cnt in aa_changes.items() if st == "NonSyn")
print(f"  Canonical C-to-U CDS changes: Non-synonymous={nsyn_cnt} ({nsyn_cnt/(syn_cnt+nsyn_cnt)*100 if (syn_cnt+nsyn_cnt) else 0:.1f}%), Synonymous={syn_cnt}")
print("  Top AA transitions from canonical C-to-U edits:")
for (r, a, st), cnt in aa_changes.most_common(10):
    print(f"    {r}->{a} ({st}): {cnt}")

# 3. Repeat & SSR Architecture of NC_043914.1
print("\n=== DISPERSED REPEATS & SSR LANDSCAPE IN NC_043914.1 ===")
rep_tsv = os.path.join(work_dir, "Repeats_SV/nc043914_self_repeats.tsv")
os.makedirs(os.path.dirname(rep_tsv), exist_ok=True)
mt_db = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/ref_db/cs_mt_ref"
subprocess.run([
    "blastn", "-query", mt_fasta, "-db", mt_db,
    "-outfmt", "6 pident length qstart qend sstart send evalue",
    "-evalue", "1e-5", "-word_size", "7", "-out", rep_tsv
], check=True)

fwd_reps, inv_reps = [], []
with open(rep_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        pid = float(row[0])
        alen = int(row[1])
        qs, qe, ss, se = int(row[2]), int(row[3]), int(row[4]), int(row[5])
        if qs == ss and qe == se:
            continue # self-hit
        if qs > ss:
            continue # deduplicate symmetric pair (A,B) vs (B,A)
        if alen >= 30 and pid >= 80.0:
            if ss < se:
                fwd_reps.append((alen, pid, qs, qe, ss, se))
            else:
                inv_reps.append((alen, pid, qs, qe, ss, se))

all_reps = fwd_reps + inv_reps
print(f"  Dispersed repeats (>=30 bp, >=80% id): Total={len(all_reps)} (Forward={len(fwd_reps)}, Inverted/Palindromic={len(inv_reps)})")
print(f"    30-99 bp  : {sum(1 for r in all_reps if 30 <= r[0] < 100)}")
print(f"    100-499 bp: {sum(1 for r in all_reps if 100 <= r[0] < 500)}")
print(f"    >=500 bp  : {sum(1 for r in all_reps if r[0] >= 500)} (max={max((r[0] for r in all_reps), default=0):,} bp)")
for r in sorted(all_reps, reverse=True)[:5]:
    print(f"    Large Repeat: len={r[0]:,} bp, pid={r[1]:.2f}%, {r[2]}-{r[3]} <-> {r[4]}-{r[5]}")
