#!/bin/bash
set -e
DIR="/mnt/d/Mitochondria/analysis/SRR9057080/ncrna"
ASSEMBLY="/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
mkdir -p "$DIR"
cd "$DIR"

echo "=== Running tRNAscan-SE ==="
tRNAscan-SE -O -o trna_organellar.txt -f trna_structures.txt -m trna_stats.txt "$ASSEMBLY"

echo "=== Running Barrnap for rRNA ==="
barrnap --kingdom bac "$ASSEMBLY" > rrna_bac.gff
barrnap --kingdom euk "$ASSEMBLY" > rrna_euk.gff

echo "=== Completed ncRNA ==="
