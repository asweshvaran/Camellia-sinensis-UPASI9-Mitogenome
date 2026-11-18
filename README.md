# *Camellia sinensis* cv. UPASI-9 Mitochondrial Transcriptome, Phylogenomics & Structural Proteomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![BioProject](https://img.shields.io/badge/NCBI%20BioProject-PRJNA539837-green)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA539837)
[![Primary SRA](https://img.shields.io/badge/NCBI%20SRA-SRR9057080-orange)](https://www.ncbi.nlm.nih.gov/sra/SRR9057080)

Reproducible bioinformatics pipelines, phylogenetic alignments, ESMFold 3D protein models, 300-DPI publication figures, and supplementary datasets accompanying the manuscript:

> **Competitive Organellar Deconvolution, Short-Repeat Evolutionary Dynamics, Multi-Order Maximum Likelihood Phylogenomics, and Drought-Stress Bioenergetic Remodeling of the Tea Plant (*Camellia sinensis* cv. UPASI-9) Mitochondrial Transcriptome**
>
> **Asweshvaran R.**<sup>1</sup>, **Nalina M**<sup>2</sup>, and **Chandrashekara KN**<sup>3,*</sup>
>
> <sup>1</sup> Shriram Research Institute, Bengaluru 560048, Karnataka, India  
> <sup>2</sup> Research and Development Division, Capsber Agriscience, Doddaballapur, Bangalore Rural 561203, Karnataka, India  
> <sup>3</sup> Plant Physiology and Biotechnology Division, UPASI Tea Research Foundation, Tea Research Institute, Valparai 642127, Coimbatore District, Tamil Nadu, India

---

## Key Scientific Highlights

1. **Competitive Dual-Organelle Deconvolution (`SRR9057080`, 53.93M clean reads, 8.02 Gb):** Eliminated a 24,259.5× legacy coverage artifact and resolved 1,395,409 cross-mapping chloroplast reads (68.15%) across 10 MTPT tracts (10,641 bp), establishing an authentic non-MTPT core depth of **28.52×** (118.88× total competitive depth).
2. **Purified *De Novo* Transcriptome Assembly (`rnaSPAdes`):** Reduced pure chloroplast contamination from **23.56% to 0.0%**, expanded mitochondrial-specific sequence **6.3-fold to 466.6 kb** (N50 = 1,730 bp in ≥500 bp core), recovered **37 canonical mitochondrial PCGs**, and identified candidate chimeric CMS locus **`NODE_2_ORF2`** (128 aa, *nad5* transmembrane fusion).
3. **Short-Repeat Evolutionary Dynamics:** Stratification of 775 dispersed repeats in `NC_043914.1` confirmed **87.1% short-repeat dominance (<100 bp; 675 pairs)**, supporting microhomology-mediated double-strand break repair (MMEJ/BIR) and hybridization-driven mitogenome expansion across *Camellia*.
4. **Strand-Specific C-to-U RNA Editing & 12-Taxon IQ-TREE2 ML Phylogenomics:** Proved **38.5-fold strand-specific C-to-U editing polarity** across 993 curated variants. Masking C-to-U editing produced a **+438.91 log-likelihood rescue** (`TIM3+F+G4` NT; `HIVw+F+I` AA; 1,000 SH-aLRT / 1,000 UFBoot), cleanly resolving **Order Ericales (99.8/100)** from **Order Fagales (*Carya illinoinensis* + *Fagus sylvatica*; 97.5/99)** and confirming **91.3%–100.0% conservation of the tRNA `U-U-C-x-A-x2` TΨC-loop motif** across 187 plant tRNAs.
5. **Drought Bioenergetics & GPU ESMFold Structural Rescue:** Water deficit induced *atp9* (**+71.3%** to 3,825.5× in UPASI-9; **+162.5%** in TRF-1) while selectively repressing Complex II (*sdh3*, log₂FC = −0.681) to mitigate ROS. ESMFold v1 proved eight C-terminal C-to-U edits in ATP1 increase hydropathy by **+25.9 KD** and rescue folding stability (**pLDDT 91.18 → 95.02**; C-terminal RMSD 1.953 Å).

---

## Repository Structure

```text
Camellia-sinensis-UPASI9-Mitogenome/
├── README.md
├── LICENSE
├── environment.yml
├── scripts/
│   ├── 01_qc_and_mapping/               # SRA retrieval, fastp QC, competitive BWA-MEM & MTPT deconvolution
│   ├── 02_assembly_and_annotation/      # Purified rnaSPAdes, GFA audit, 5-tier ORF curation, MUMmer4 synteny
│   ├── 03_variants_repeats_and_drought/ # Strand-specific editing, short-repeat binning, PRJNA539837 drought cohort
│   ├── 04_phylogenomics_and_trna/       # 12-taxon IQ-TREE2 ML (NT + AA), Ka/Ks, RSCU, 187-tRNA Ψ-loop scanner
│   ├── 05_structural_modeling/          # GPU ESMFold v1 modeling & Kabsch Cα superposition (ATP1, COX1, NODE_2)
│   └── 06_figures_and_tables/           # Publication 300-DPI figure & Supplementary Table S1–S14 generators
├── supplementary_tables/                # Supplementary Tables S1–S14 (.csv / .tsv)
├── phylogeny/                           # 12-taxon codon alignments (.fasta), IQ-TREE2 logs, & Newick .treefile
├── structures/                          # ESMFold PDB models (ATP1 WT/Edited, COX1, NODE_2_ORF2) & metrics
├── figures/                             # High-resolution 300-DPI Figures 1–10 (.png)
└── references/                          # Crossref-verified BibTeX, Zotero RIS, & MDPI ACS JSON (N = 28)
```

---

## Data Sources & Reproducibility

Raw sequencing reads are publicly archived in NCBI SRA (not stored in GitHub due to file size limits):
- **Factorial Drought Cohort (`PRJNA539837`):** `SRR9057080` (UPASI-9 Stress), `SRR9057081` (UPASI-9 Control), `SRR9057078` (TRF-1 Stress), `SRR9057079` (TRF-1 Control).
- **Comparative RNA-Seq Panel:** 24 *Camellia* accessions from `PRJNA303133` and `PRJNA520786` (listed in `supplementary_tables/Table_S8_Cohort_28_PanMitochondrial_Summary.csv`).
- **Complete Organellar References (NCBI GenBank):** `NC_043914.1`, `NC_020019.1`, `PX398764.1`–`PX398765.1`, `PX759560.1`, `PZ267020.1`–`PZ267021.1`, `PZ501402.1`, `PZ378386.1`–`PZ378387.1`, `NC_086781.1`, `NC_023338.1`, `PX519271.1`, `MT446430.1`, `NC_037304.1`.

### Environment Setup

```bash
mamba env create -f environment.yml
mamba activate camellia-mito-q1
```

---

## Citation

If you use these pipelines, alignments, structural models, or curated tables, please cite:

> Asweshvaran, R.; Nalina, M.; Chandrashekara, K.N. Competitive Organellar Deconvolution, Short-Repeat Evolutionary Dynamics, Multi-Order Maximum Likelihood Phylogenomics, and Drought-Stress Bioenergetic Remodeling of the Tea Plant (*Camellia sinensis* cv. UPASI-9) Mitochondrial Transcriptome. *Diversity* **2026** (Submitted).
