import pandas as pd
import json
import csv
import os

out_dir = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables"
os.makedirs(out_dir, exist_ok=True)

print("Generating supplementary tables...")

# Table S1: Sequencing Metadata & Fastp QC Audit
s1_data = [
    ("SRA Run Accession", "SRR9057080"),
    ("BioProject", "PRJNA539837"),
    ("BioSample", "SAMN11513496"),
    ("Experiment", "SRX5833307"),
    ("Organism / TaxID", "Camellia sinensis (TaxID: 4442)"),
    ("Cultivar / Source", "UPASI-9 (Tea Research Institute, UPASI, India)"),
    ("Sequencing Platform", "Illumina HiSeq 4000 (150 bp Paired-End)"),
    ("Library Strategy", "RNA-Seq (Transcriptomic cDNA)"),
    ("Total Raw Spots", "28,395,365 spots (56,790,730 paired reads)"),
    ("Total Raw Bases", "8,518,609,500 bp (~8.52 Gb)"),
    ("Raw Q20 / Q30 Rates", "91.36% / 83.71%"),
    ("Raw GC Content", "44.98%"),
    ("Clean Reads (Post-fastp)", "53,932,552 reads"),
    ("Clean Bases (Post-fastp)", "8,020,758,642 bp (~8.02 Gb)"),
    ("Clean Q20 / Q30 Rates", "92.94% / 85.52%"),
    ("Filtered Low-Quality Reads", "2,832,396 reads"),
    ("Adapter Trimmed Reads", "3,574,270 reads"),
    ("Duplication Rate", "15.68%")
]
pd.DataFrame(s1_data, columns=["Parameter", "Value"]).to_csv(os.path.join(out_dir, "Table_S1_Sequencing_Metadata_and_QC.csv"), index=False)

# Table S2: Assembly Metrics
s2_data = [
    ("Reference Genome", "NC_043914.1 (707,441 bp)"),
    ("Assembler Software", "GetOrganelle v1.7.7.1 / SPAdes v4.3.0 (embplant_mt)"),
    ("K-mer Range Evaluated", "21, 45, 65, 85, 105"),
    ("Total Assembled Scaffolds", "1,160"),
    ("Total Assembly Length", "258,681 bp"),
    ("Reference Completeness", "36.57%"),
    ("Assembly N50", "193 bp"),
    ("Assembly L50", "202"),
    ("Longest Scaffold", "7,316 bp (NODE_1, coverage 129x)"),
    ("Shortest Scaffold", "106 bp"),
    ("Mean Scaffold Length", "223.0 bp"),
    ("Overall Assembly GC Content", "49.82%"),
    ("Scaffolds >= 500 bp", "65 (94,707 bp)"),
    ("Scaffolds >= 1,000 bp", "36 (73,531 bp)"),
    ("Scaffolds >= 2,000 bp", "13 (42,252 bp)")
]
pd.DataFrame(s2_data, columns=["Metric", "Empirical Value"]).to_csv(os.path.join(out_dir, "Table_S2_De_Novo_Assembly_Metrics.csv"), index=False)

# Table S3: Core Genes Identified
s3_data = [
    ("atp1", "ATP synthase F1 subunit alpha", "Full-Length", "98.4%", "509 aa", "509 aa", "NODE_17 (cov 62x)"),
    ("atp4", "ATP synthase subunit 4", "Full-Length", "94.1%", "186 aa", "192 aa", "NODE_29 (cov 23x)"),
    ("atp6", "ATP synthase subunit 6", "Fragmented", "90.0%", "224 aa", "315 aa", "NODE_39 (cov 23x)"),
    ("atp9", "ATP synthase subunit 9", "Full-Length", "87.8%", "86 aa", "75 aa", "NODE_2 (cov 47x)"),
    ("ccmB", "Cytochrome c biogenesis factor B", "Fragmented", "86.0%", "157 aa", "206 aa", "NODE_71 (cov 1x)"),
    ("ccmC", "Cytochrome c biogenesis factor C", "Full-Length", "86.8%", "251 aa", "250 aa", "NODE_41 (cov 7x)"),
    ("ccmFc", "Cytochrome c biogenesis factor Fc", "Fragmented", "94.1%", "302 aa", "446 aa", "NODE_32 (cov 2x)"),
    ("ccmFn", "Cytochrome c biogenesis factor Fn", "Full-Length", "95.1%", "618 aa", "612 aa", "NODE_14 (cov 6x)"),
    ("cob", "Cytochrome b / apocytochrome", "Full-Length", "95.4%", "394 aa", "393 aa", "NODE_4 (cov 16x)"),
    ("matR", "Mitochondrial maturase", "Fragmented", "98.6%", "417 aa", "651 aa", "NODE_16 (cov 3x)"),
    ("nad1", "NADH dehydrogenase subunit 1", "Fragmented", "90.4%", "197 aa", "327 aa", "NODE_52 (cov 8x)"),
    ("nad2", "NADH dehydrogenase subunit 2", "Fragmented", "92.0%", "219 aa", "465 aa", "NODE_55 (cov 9x)"),
    ("nad3", "NADH dehydrogenase subunit 3", "Full-Length", "87.3%", "119 aa", "118 aa", "NODE_25 (cov 79x)"),
    ("nad4", "NADH dehydrogenase subunit 4", "Fragmented", "97.9%", "342 aa", "494 aa", "NODE_30 (cov 12x)"),
    ("nad4l", "NADH dehydrogenase subunit 4L", "Fragmented", "86.1%", "101 aa", "90 aa", "NODE_29 (cov 23x)"),
    ("nad5", "NADH dehydrogenase subunit 5", "Fragmented", "94.4%", "568 aa", "664 aa", "NODE_2 (cov 47x)"),
    ("nad6", "NADH dehydrogenase subunit 6", "Full-Length", "94.5%", "200 aa", "205 aa", "NODE_11 (cov 21x)"),
    ("rpl5", "Large ribosomal subunit L5", "Full-Length", "95.2%", "207 aa", "187 aa", "NODE_4 (cov 16x)"),
    ("rps1", "Small ribosomal subunit S1", "Full-Length", "98.2%", "203 aa", "167 aa", "NODE_48 (cov 4x)"),
    ("rps12", "Small ribosomal subunit S12", "Full-Length", "94.4%", "126 aa", "125 aa", "NODE_25 (cov 79x)"),
    ("rps14", "Small ribosomal subunit S14", "Full-Length", "99.0%", "101 aa", "100 aa", "NODE_4 (cov 16x)")
]
pd.DataFrame(s3_data, columns=["Gene", "Product Name", "Status", "Identity %", "Query Size", "Ref Size", "Scaffold Location"]).to_csv(os.path.join(out_dir, "Table_S3_Annotated_Core_Mitochondrial_Genes.csv"), index=False)

# Table S4: Transfer RNA Annotation
s4_data = [
    ("tRNA-Pro", "TGG", "NODE_1", 7209, 7136, "74 bp", "None", 65.6),
    ("tRNA-Trp", "CCA", "NODE_1", 6965, 6892, "74 bp", "None", 68.7),
    ("tRNA-Ile", "GAT", "NODE_6", 2492, 2405, "88 bp", "16 bp (2456-2441)", 16.4),
    ("tRNA-Ala", "TGC", "NODE_6", 1404, 1300, "105 bp", "40 bp (1367-1328)", 16.3),
    ("tRNA-Met", "CAT", "NODE_9", 2607, 2535, "73 bp", "None", 59.5),
    ("tRNA-His", "GTG", "NODE_61", 6, 79, "74 bp", "None", 58.5),
    ("tRNA-Glu", "TTC", "NODE_89", 297, 226, "72 bp", "None", 60.8)
]
pd.DataFrame(s4_data, columns=["tRNA_Isotype", "Anticodon", "Scaffold", "Begin_bp", "End_bp", "Length", "Intron", "Infernal_Score"]).to_csv(os.path.join(out_dir, "Table_S4_Identified_tRNA_Genes.csv"), index=False)

# Table S7: Comparative Mitogenomes in Genus Camellia
s7_data = [
    ("Camellia sinensis cv. Yidan", "PX759560", 931612, "Circular Master", "WGS Assembly"),
    ("Camellia drupifera", "PZ378386-87", 978603, "2 Circular Chromosomes", "Multipartite Ring"),
    ("Camellia oleifera", "PX170215-21", 916337, "7 Circular Chromosomes", "Multichromosomal Multipartite"),
    ("Camellia sinensis cv. Wuyi Narcissus", "PZ267020-21", 910237, "2 Linear Chromosomes", "Multipartite Linear"),
    ("Camellia pitardii", "PZ501402", 895290, "Circular Master", "Complete Circular"),
    ("Camellia sinensis cv. Xinyang10", "PX398764-65", 845076, "2 Circular Chromosomes", "Multipartite Circular"),
    ("Camellia sinensis RefSeq", "NC_043914.1", 707441, "Circular Master", "Reference Strain"),
    ("Camellia sinensis cv. UPASI-9 (This study)", "SRR9057080", 258681, "1,160 Scaffolds (Core Transcriptome)", "Deep RNA-Seq Exome")
]
pd.DataFrame(s7_data, columns=["Species_Cultivar", "Accession", "Genome_Size_bp", "Topology_Architecture", "Study_Type"]).to_csv(os.path.join(out_dir, "Table_S7_Comparative_Camellia_Mitogenomes.csv"), index=False)

print("Supplementary tables successfully created.")
