import os
import glob
import csv
import json
import gzip
from collections import defaultdict, Counter

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables"
meta_csv = "/mnt/d/Mitochondria/Q1_Investigation/Comparative_Genomics/cohort_28_sra_metadata.csv"

meta = {}
with open(meta_csv, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        meta[row["Run"]] = row

mtpt_pos = set()
mtpt_tsv = "/mnt/d/Mitochondria/Q1_Investigation/Clean_Reanalysis/nc043914_vs_cp_mtpt.tsv"
with open(mtpt_tsv) as f:
    for row in csv.reader(f, delimiter="\t"):
        if int(row[3]) >= 50:
            qs, qe = min(int(row[6]), int(row[7])), max(int(row[6]), int(row[7]))
            for p in range(qs, qe + 1):
                mtpt_pos.add(p)

rows = []
non_mtpt_sites = defaultdict(list)
mtpt_sites = defaultdict(list)

for d in sorted(glob.glob("/mnt/d/Mitochondria/analysis/SRR*")):
    srr = os.path.basename(d)
    fjson = os.path.join(d, "qc", f"{srr}_fastp.json")
    flag = os.path.join(d, "mapping", f"{srr}_flagstat.txt")
    vcf = os.path.join(d, "mapping", f"{srr}_variants.vcf.gz")

    clean_reads, q30, dup, gc = 0, 0.0, 0.0, 0.0
    if os.path.exists(fjson):
        with open(fjson) as f:
            j = json.load(f)
            af = j["summary"]["after_filtering"]
            clean_reads = af["total_reads"]
            q30 = af["q30_rate"] * 100
            dup = j["duplication"]["rate"] * 100
            gc = af["gc_content"] * 100

    tot, mapped, pct = 0, 0, 0.0
    if os.path.exists(flag):
        with open(flag) as f:
            for line in f:
                if "primary" in line and "total" not in line and "mapped" not in line and "duplicates" not in line:
                    tot = int(line.split()[0])
                elif "primary mapped" in line:
                    mapped = int(line.split()[0])
                    pct = float(line.split("(")[1].split("%")[0])

    # Count high-quality variants (QUAL >= 20, DP >= 10) partitioned into Non-MTPT vs MTPT!
    hq_non_mtpt, hq_mtpt, ts, tv = 0, 0, 0, 0
    if os.path.exists(vcf):
        with gzip.open(vcf, "rt") as f:
            for line in f:
                if line.startswith("#"): continue
                p = line.strip().split("\t")
                pos = int(p[1])
                ref, alt = p[3], p[4].split(",")[0]
                qual = float(p[5]) if p[5] != "." else 0.0
                if qual < 20: continue
                info = dict(item.split("=") for item in p[7].split(";") if "=" in item)
                dp = int(info.get("DP", 0))
                mq = float(info.get("MQ", 0))
                if dp >= 10 and mq >= 30 and len(ref) == 1 and len(alt) == 1:
                    if pos in mtpt_pos:
                        hq_mtpt += 1
                        mtpt_sites[pos].append(srr)
                    else:
                        hq_non_mtpt += 1
                        non_mtpt_sites[pos].append(srr)
                        pair = f"{ref}>{alt}"
                        if pair in ("A>G", "G>A", "C>T", "T>C"): ts += 1
                        else: tv += 1

    m = meta.get(srr, {})
    rows.append({
        "Run": srr,
        "BioProject": m.get("BioProject", "?"),
        "Strategy": m.get("LibraryStrategy", "?"),
        "Cultivar_Sample": m.get("Title", "?")[:50],
        "Clean_Reads": clean_reads,
        "Q30_pct": round(q30, 2),
        "Dup_pct": round(dup, 2),
        "Mapped_to_NC043914": mapped,
        "Mapped_pct": round(pct, 2),
        "HQ_NonMTPT_SNPs": hq_non_mtpt,
        "HQ_MTPT_SNPs": hq_mtpt,
        "NonMTPT_TsTv": round(ts / tv, 2) if tv > 0 else "N/A"
    })

out_csv = os.path.join(out_dir, "Table_S8_Cohort_28_PanMitochondrial_Summary.csv")
with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print(f"Saved 28-accession summary to {out_csv}")
print(f"Non-MTPT HQ SNPs across 28 accessions: {len(non_mtpt_sites):,} unique positions")
print(f"  Shared >=50% (>=14 accessions): {sum(1 for v in non_mtpt_sites.values() if len(v) >= 14):,}")
print(f"  Cultivar-specific (<=2 accessions): {sum(1 for v in non_mtpt_sites.values() if len(v) <= 2):,}")
print(f"MTPT-zone HQ SNPs across 28 accessions: {len(mtpt_sites):,} unique positions")
for r in sorted(rows, key=lambda x: x["Mapped_to_NC043914"], reverse=True)[:8]:
    print(f"  {r['Run']} ({r['BioProject']}) | Mito={r['Mapped_to_NC043914']:,} ({r['Mapped_pct']}%) | NonMTPT_SNPs={r['HQ_NonMTPT_SNPs']} (Ts/Tv={r['NonMTPT_TsTv']}) | {r['Cultivar_Sample'][:35]}")
