import subprocess
import os

bam_path = "/mnt/d/Mitochondria/analysis/SRR9057080/mapping/SRR9057080_aligned.bam"
ref_path = "/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"

# 1. Total reference length
with open(ref_path) as f:
    ref_len = sum(len(line.strip()) for line in f if not line.startswith('>'))

# 2. samtools depth calculation
cmd = f"samtools depth -a {bam_path}"
proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)

depths = []
covered_1x = 0
covered_10x = 0
covered_30x = 0
covered_100x = 0

for line in proc.stdout:
    parts = line.strip().split()
    if len(parts) >= 3:
        d = int(parts[2])
        depths.append(d)
        if d >= 1: covered_1x += 1
        if d >= 10: covered_10x += 1
        if d >= 30: covered_30x += 1
        if d >= 100: covered_100x += 1

proc.wait()

total_pos = len(depths)
if total_pos > 0:
    mean_depth = sum(depths) / total_pos
    depths.sort()
    median_depth = depths[total_pos // 2]
    total_mapped_bases = sum(depths)
else:
    mean_depth = 0
    median_depth = 0
    total_mapped_bases = 0

print("=====================================================")
print("EXACT BAM COVERAGE METRICS (SRR9057080)")
print("=====================================================")
print(f"Reference Sequence Length      : {ref_len:,} bp")
print(f"Total Mapped Bases to Mt       : {total_mapped_bases:,} bp")
print(f"True Mean Sequencing Depth     : {mean_depth:.2f}x")
print(f"True Median Sequencing Depth   : {median_depth:.2f}x")
print(f"Breadth of Coverage >= 1x      : {covered_1x:,} bp ({covered_1x / ref_len * 100:.2f}%)")
print(f"Breadth of Coverage >= 10x     : {covered_10x:,} bp ({covered_10x / ref_len * 100:.2f}%)")
print(f"Breadth of Coverage >= 30x     : {covered_30x:,} bp ({covered_30x / ref_len * 100:.2f}%)")
print(f"Breadth of Coverage >= 100x    : {covered_100x:,} bp ({covered_100x / ref_len * 100:.2f}%)")
print("=====================================================")
