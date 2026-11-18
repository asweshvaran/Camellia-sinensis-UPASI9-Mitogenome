import json

fastp_json = "/mnt/d/Mitochondria/analysis/SRR9057080/qc/SRR9057080_fastp.json"
with open(fastp_json) as f:
    qc = json.load(f)

summary = qc['summary']
filtering = qc['filtering_result']

print("==================================================")
print("RAW SEQUENCING AND QC AUDIT: SRR9057080")
print("==================================================")
print(f"Total Raw Reads (Before QC)    : {summary['before_filtering']['total_reads']:,}")
print(f"Total Raw Bases (Before QC)    : {summary['before_filtering']['total_bases']:,} bp")
print(f"Raw Q20 Rate                   : {summary['before_filtering']['q20_rate']*100:.2f}%")
print(f"Raw Q30 Rate                   : {summary['before_filtering']['q30_rate']*100:.2f}%")
print(f"Raw GC Content                 : {summary['before_filtering']['gc_content']*100:.2f}%")
print("--------------------------------------------------")
print(f"Clean Reads (After QC)         : {summary['after_filtering']['total_reads']:,}")
print(f"Clean Bases (After QC)         : {summary['after_filtering']['total_bases']:,} bp")
print(f"Clean Q20 Rate                 : {summary['after_filtering']['q20_rate']*100:.2f}%")
print(f"Clean Q30 Rate                 : {summary['after_filtering']['q30_rate']*100:.2f}%")
print(f"Clean GC Content               : {summary['after_filtering']['gc_content']*100:.2f}%")
print("--------------------------------------------------")
print(f"Passed Filter Reads            : {filtering['passed_filter_reads']:,}")
print(f"Low Quality Reads Filtered     : {filtering['low_quality_reads']:,}")
print(f"Too Many N Reads Filtered      : {filtering['too_many_N_reads']:,}")
print(f"Adapter Trimmed Reads          : {filtering['adapter_trimmed_reads']:,}")
print(f"Duplication Rate               : {qc['duplication']['rate']*100:.2f}%")
print(f"Insert Size Peak               : {qc['insert_size']['peak']} bp")
print("==================================================")
