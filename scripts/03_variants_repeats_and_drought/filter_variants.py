import gzip
import csv

vcf_path = "/mnt/d/Mitochondria/analysis/SRR9057080/mapping/SRR9057080_variants.vcf.gz"

raw_total = 0
filtered_pass = 0
failed_mq = 0
failed_dp = 0
failed_qual = 0
failed_strand = 0
failed_alt_depth = 0

pass_snps = 0
pass_indels = 0
pass_ts = 0
pass_tv = 0

pass_fixed = 0 # VAF >= 0.85
pass_heteroplasmic = 0 # 0.10 <= VAF < 0.85
pass_low_freq = 0 # VAF < 0.10

curated_variants = []

with gzip.open(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        raw_total += 1
        parts = line.strip().split("\t")
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4].split(",")[0]
        qual = float(parts[5]) if parts[5] != "." else 0.0
        info_str = parts[7]
        fmt = parts[8].split(":")
        sample = parts[9].split(":")
        fmt_dict = dict(zip(fmt, sample))
        
        info = dict(item.split("=") for item in info_str.split(";") if "=" in item)
        
        # 1. Quality filter
        if qual < 20.0:
            failed_qual += 1
            continue
            
        # 2. Depth filter
        dp = int(info.get("DP", 0))
        if dp < 20:
            failed_dp += 1
            continue
            
        # 3. Mapping quality filter
        mq = float(info.get("MQ", 0))
        if mq < 30.0:
            failed_mq += 1
            continue
            
        # 4. Strand bias filter using DP4 (ref_fwd, ref_rev, alt_fwd, alt_rev)
        dp4 = [int(x) for x in info.get("DP4", "0,0,0,0").split(",")]
        alt_fwd = dp4[2]
        alt_rev = dp4[3]
        alt_tot = alt_fwd + alt_rev
        ref_tot = dp4[0] + dp4[1]
        
        if alt_fwd < 2 or alt_rev < 2:
            failed_strand += 1
            continue
            
        if alt_tot < 4:
            failed_alt_depth += 1
            continue
            
        # Variant Passed Strict Curation!
        filtered_pass += 1
        vaf = alt_tot / (ref_tot + alt_tot) if (ref_tot + alt_tot) > 0 else 0.0
        
        is_snp = len(ref) == 1 and len(alt) == 1
        if is_snp:
            pass_snps += 1
            pair = f"{ref}>{alt}"
            if pair in ['A>G', 'G>A', 'C>T', 'T>C']:
                pass_ts += 1
            else:
                pass_tv += 1
        else:
            pass_indels += 1
            
        if vaf >= 0.85:
            pass_fixed += 1
            cat = "Fixed_Homoplasmic"
        elif vaf >= 0.10:
            pass_heteroplasmic += 1
            cat = "Candidate_Heteroplasmy"
        else:
            pass_low_freq += 1
            cat = "Low_Frequency_Subclonal"
            
        curated_variants.append({
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "qual": qual,
            "dp": dp,
            "mq": mq,
            "alt_depth": alt_tot,
            "vaf": vaf,
            "category": cat,
            "is_snp": is_snp
        })

ts_tv = pass_ts / pass_tv if pass_tv > 0 else 0

print("=================================================================")
print("RIGOROUS VARIANT & HETEROPLASMY AUDIT (SRR9057080)")
print("=================================================================")
print(f"Total Raw Unfiltered Variants Reported       : {raw_total:,}")
print(f"Filtered out due to Low Quality (QUAL < 20)  : {failed_qual:,}")
print(f"Filtered out due to Low Depth (DP < 20)      : {failed_dp:,}")
print(f"Filtered out due to Low MapQ (MQ < 30 / NUMT): {failed_mq:,}")
print(f"Filtered out due to Strand Bias (Alt F/R < 2): {failed_strand:,}")
print(f"Filtered out due to Insufficient Alt Depth   : {failed_alt_depth:,}")
print("-----------------------------------------------------------------")
print(f"High-Confidence Curated Variants (PASSED)    : {filtered_pass:,} ({filtered_pass/raw_total*100:.1f}%)")
print(f"  - Curated SNPs                             : {pass_snps:,} ({pass_snps/filtered_pass*100:.1f}%)")
print(f"  - Curated Indels                           : {pass_indels:,} ({pass_indels/filtered_pass*100:.1f}%)")
print(f"  - Transitions (Ts) : {pass_ts:,} | Transversions (Tv) : {pass_tv:,} | Ts/Tv : {ts_tv:.2f}")
print("-----------------------------------------------------------------")
print("CURATED ALLELIC FREQUENCY ARCHITECTURE:")
print(f"  - Homoplasmic / Fixed Variants (VAF >= 85%) : {pass_fixed:,} ({pass_fixed/filtered_pass*100:.1f}%)")
print(f"  - Candidate Heteroplasmic Sites (10% <= VAF < 85%): {pass_heteroplasmic:,} ({pass_heteroplasmic/filtered_pass*100:.1f}%)")
print(f"  - Low-Frequency Subclonal Variants (VAF < 10%): {pass_low_freq:,} ({pass_low_freq/filtered_pass*100:.1f}%)")
print("=================================================================")

# Save curated variants to file
out_csv = "/mnt/d/Mitochondria/Q1_Investigation/Heteroplasmy_Variants/curated_high_confidence_variants.csv"
with open(out_csv, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["pos", "ref", "alt", "qual", "dp", "mq", "alt_depth", "vaf", "category", "is_snp"])
    writer.writeheader()
    writer.writerows(curated_variants)
print(f"Saved {len(curated_variants)} curated variants to {out_csv}")
