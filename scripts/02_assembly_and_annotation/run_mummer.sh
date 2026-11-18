#!/bin/bash
set -e
DIR="/mnt/d/Mitochondria/analysis/SRR9057080/synteny"
REF="/mnt/d/Mitochondria/reference/Camellia_sinensis_mitochondrion_NC_043914.1.fasta"
QRY="/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/scaffolds.fasta"
mkdir -p "$DIR"
cd "$DIR"

nucmer --prefix=srr9057080_vs_ref "$REF" "$QRY"
show-coords -r -c -l srr9057080_vs_ref.delta > srr9057080_vs_ref.coords
echo "Synteny alignment complete. Total aligned segments:"
wc -l srr9057080_vs_ref.coords
