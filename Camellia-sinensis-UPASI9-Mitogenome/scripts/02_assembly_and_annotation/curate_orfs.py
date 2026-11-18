import csv
from Bio import SeqIO
from collections import Counter

faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
blast_tsv = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/blastp_vs_ref.tsv"

# 1. Parse blast hits
blast_hits = {}
with open(blast_tsv) as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        qid, sid, pident, length, mismatch, gapopen, qs, qe, ss, se, evalue, bitscore, qlen, slen = row
        pident = float(pident)
        evalue = float(evalue)
        bitscore = float(bitscore)
        qlen = int(qlen)
        slen = int(slen)
        length = int(length)
        cov = (length / slen) * 100.0
        
        if qid not in blast_hits or bitscore > blast_hits[qid]['bitscore']:
            blast_hits[qid] = {
                'sid': sid,
                'gene': sid.split('|')[0],
                'pident': pident,
                'cov': cov,
                'evalue': evalue,
                'bitscore': bitscore,
                'qlen': qlen,
                'slen': slen,
                'length': length
            }

# 2. Iterate through all 849 predicted ORFs
classification = {
    'A': [], # Conserved mitochondrial genes (>=85% identity, >=70% cov)
    'B': [], # Divergent mitochondrial genes (homology, but lower identity or fragmented)
    'C': [], # Potential chimeric ORFs (chimeric termini or partial overlap with mitochondrial core)
    'D': [], # Hypothetical protein candidates (length >= 100 aa, full ORF with ATG, no blast hit)
    'E': []  # Unsupported or likely spurious ORFs (length < 100 aa, edge start, or random intergenic noise)
}

all_records = list(SeqIO.parse(faa_path, "fasta"))

for r in all_records:
    qid = r.id
    seq_len = len(r.seq)
    desc = r.description
    
    is_edge = 'start_type=Edge' in desc or 'partial=1' in desc
    has_atg = 'start_type=ATG' in desc
    
    if qid in blast_hits and blast_hits[qid]['evalue'] <= 1e-5:
        hit = blast_hits[qid]
        if hit['pident'] >= 90.0 and hit['cov'] >= 75.0:
            classification['A'].append((r, hit))
        elif hit['pident'] >= 75.0:
            # Check if chimeric candidate: covers part of reference gene with extra flanking length >= 50 aa
            if seq_len > hit['length'] + 50 and is_edge == False:
                classification['C'].append((r, hit))
            else:
                classification['B'].append((r, hit))
        else:
            classification['B'].append((r, hit))
    else:
        # Unmatched ORFs
        if seq_len >= 100 and not is_edge and has_atg:
            classification['D'].append((r, None))
        else:
            classification['E'].append((r, None))

total = len(all_records)
print("=================================================================")
print("EVIDENCE-BASED CLASSIFICATION OF 849 PREDICTED ORFS (SRR9057080)")
print("=================================================================")
print(f"Total Predicted ORFs in Raw Prodigal Output : {total}")
print(f"Class A: Conserved Core Mitochondrial Genes  : {len(classification['A'])} ({len(classification['A'])/total*100:.2f}%)")
print(f"Class B: Divergent / Fragmented Homologues   : {len(classification['B'])} ({len(classification['B'])/total*100:.2f}%)")
print(f"Class C: Candidate Chimeric ORFs (CMS-like)  : {len(classification['C'])} ({len(classification['C'])/total*100:.2f}%)")
print(f"Class D: Hypothetical Robust ORFs (>=100 aa) : {len(classification['D'])} ({len(classification['D'])/total*100:.2f}%)")
print(f"Class E: Likely Spurious ORFs / Noise (<100aa): {len(classification['E'])} ({len(classification['E'])/total*100:.2f}%)")
print("=================================================================\n")

# Save detailed table to file
out_tsv = "/mnt/d/Mitochondria/Q1_Investigation/ORF_Validation/orf_curated_classification.tsv"
with open(out_tsv, "w", newline="") as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(["ORF_ID", "Length_aa", "Class", "Description", "Homology_Target", "Pident", "Coverage_Subject", "Evalue"])
    for cls_name, items in classification.items():
        for r, hit in items:
            target = hit['gene'] if hit else "None"
            pident = f"{hit['pident']:.1f}" if hit else "N/A"
            cov = f"{hit['cov']:.1f}" if hit else "N/A"
            evalue = f"{hit['evalue']:.2e}" if hit else "N/A"
            cls_desc = {
                'A': "Conserved Core Mitochondrial Gene",
                'B': "Divergent/Fragmented Homologue",
                'C': "Candidate Chimeric ORF",
                'D': "Hypothetical Protein Candidate (>=100aa)",
                'E': "Spurious ORF / Non-coding Noise (<100aa)"
            }[cls_name]
            writer.writerow([r.id, len(r.seq), cls_name, cls_desc, target, pident, cov, evalue])

print(f"Curated ORF classification table written to {out_tsv}")
