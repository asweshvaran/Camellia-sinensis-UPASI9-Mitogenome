import os
import csv

out_csv = "/mnt/d/Mitochondria/Q1_Investigation/Supplementary_Tables/Table_S12_Repeat_Distribution_and_Expansion.csv"

rows = [
    {"Repeat_Class": "Short Repeats (20-50 bp)", "Size_Range_bp": "20-50", "Repeat_Count": 333, "Percentage_of_Total_Repeats": 43.0, "Total_Sequence_bp": 13387, "Percentage_of_Repeat_Sequence": 8.0, "Average_Length_bp": 40.2, "Recombination_and_Evolutionary_Mechanism": "Microhomology-mediated end joining (MMEJ), replication slippage, and genome inflation (Wang et al., 2024; Christensen, 2013)"},
    {"Repeat_Class": "Short Repeats (51-100 bp)", "Size_Range_bp": "51-100", "Repeat_Count": 342, "Percentage_of_Total_Repeats": 44.1, "Total_Sequence_bp": 23815, "Percentage_of_Repeat_Sequence": 14.1, "Average_Length_bp": 69.6, "Recombination_and_Evolutionary_Mechanism": "Break-induced replication (BIR) hot spots, short-repeat expansion, syntenic breakpoint scrambling"},
    {"Repeat_Class": "Subtotal: All Short Repeats (<100 bp)", "Size_Range_bp": "20-100", "Repeat_Count": 675, "Percentage_of_Total_Repeats": 87.1, "Total_Sequence_bp": 37202, "Percentage_of_Repeat_Sequence": 22.1, "Average_Length_bp": 55.1, "Recombination_and_Evolutionary_Mechanism": "Primary driver of mitogenome size variation and interspecific hybridization fossil record in genus Camellia"},
    {"Repeat_Class": "Intermediate Repeats (101-500 bp)", "Size_Range_bp": "101-500", "Repeat_Count": 96, "Percentage_of_Total_Repeats": 12.4, "Total_Sequence_bp": 17618, "Percentage_of_Repeat_Sequence": 10.5, "Average_Length_bp": 183.5, "Recombination_and_Evolutionary_Mechanism": "Low-frequency non-homologous recombination, generation of candidate chimeric open reading frames (e.g. NODE_2_ORF2)"},
    {"Repeat_Class": "Intermediate Repeats (501-1000 bp)", "Size_Range_bp": "501-1000", "Repeat_Count": 0, "Percentage_of_Total_Repeats": 0.0, "Total_Sequence_bp": 0, "Percentage_of_Repeat_Sequence": 0.0, "Average_Length_bp": 0.0, "Recombination_and_Evolutionary_Mechanism": "Absent in NC_043914.1 reference architecture"},
    {"Repeat_Class": "Large Repeat Pair 1 (Direct)", "Size_Range_bp": "29807", "Repeat_Count": 2, "Percentage_of_Total_Repeats": 0.26, "Total_Sequence_bp": 59614, "Percentage_of_Repeat_Sequence": 35.4, "Average_Length_bp": 29807.0, "Recombination_and_Evolutionary_Mechanism": "High-frequency reversible homologous recombination generating two autonomous subgenomic circular isoforms"},
    {"Repeat_Class": "Large Repeat Pair 2 (Inverted)", "Size_Range_bp": "26965", "Repeat_Count": 2, "Percentage_of_Total_Repeats": 0.26, "Total_Sequence_bp": 53930, "Percentage_of_Repeat_Sequence": 32.0, "Average_Length_bp": 26965.0, "Recombination_and_Evolutionary_Mechanism": "High-frequency genomic inversion flipping orientation of intervening segment and maintaining multipartite dynamic network"},
    {"Repeat_Class": "Subtotal: Large Repeats (>1000 bp)", "Size_Range_bp": ">1000", "Repeat_Count": 4, "Percentage_of_Total_Repeats": 0.52, "Total_Sequence_bp": 113544, "Percentage_of_Repeat_Sequence": 67.4, "Average_Length_bp": 28386.0, "Recombination_and_Evolutionary_Mechanism": "Structural isomerization, multipartite chromosome equilibrium, and subgenomic branching in assembly graphs"},
    {"Repeat_Class": "Total Mitogenome Repeats", "Size_Range_bp": "20-29807", "Repeat_Count": 775, "Percentage_of_Total_Repeats": 100.0, "Total_Sequence_bp": 168364, "Percentage_of_Repeat_Sequence": 100.0, "Average_Length_bp": 217.2, "Recombination_and_Evolutionary_Mechanism": "Total repeat content comprises 23.80% of the 707,441 bp master circular mitogenome"}
]

with open(out_csv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

print("Generated Table_S12_Repeat_Distribution_and_Expansion.csv successfully!")
