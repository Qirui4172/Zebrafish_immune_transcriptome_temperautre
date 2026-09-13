#!/usr/bin/env bash
set -euo pipefail

GENOME_INDEX="/home/qirui/spleen_miRNA/genomeIndex/Danio_rerio.GRCz11.dna.primary_assembly_modified"
GENOME_FASTA="${GENOME_INDEX}.fa"
HAIRPIN="dre_hairpin.fa"
MATURE="dre_mature_modified.fa"
QUANTIFIER_HAIRPIN="dre_hairpin+novel.fa"
QUANTIFIER_MATURE="dre_mature+novel.fa"

Groups=("24ct" "24lps" "28ct" "28lps" "32ct" "32lps")

for group in "${Groups[@]}"; do
	# Mapping
	mapper.pl "mapper${group}.txt" -e -d -h -i -j -l 18 -m -n -o 16 \
		-p "$GENOME_INDEX" -s "${group}.pool.fa" -t "${group}.pool.arf"

	# miRDeep2
	miRDeep2.pl "${group}.pool.fa" "$GENOME_FASTA" "${group}.pool.arf" \
		"$MATURE" none "$HAIRPIN" -t Zebrafish -p -T 16 2 >> report.log

	# Quantify
	quantifier.pl -p "$QUANTIFIER_HAIRPIN" -m "$QUANTIFIER_MATURE" \
		-P -r "${group}.pool.fa" -t dre
done

