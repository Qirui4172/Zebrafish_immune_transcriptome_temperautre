#!/bin/sh


Groups=("24ct" "24lps" "28ct" "28lps" "32ct" "32lps")

for group in ${Groups[@]};do
	# Mapping
	mapper.pl mapper${group}.txt -e -d -h -i -j -l 18 -m -n -o 16 -p /home/qirui/spleen_miRNA/genomeIndex/Danio_rerio.GRCz11.dna.primary_assembly_modified -s ${group}.pool.fa -t ${group}.pool.arf

	# miRDeep2
	miRDeep2.pl ${group}.pool.fa Danio_rerio.GRCz11.dna.primary_assembly_modified.fa ${group}.pool.arf dre_mature_modified.fa none dre_hairpin.fa -t Zebrafish -p -T 16 2 >> report.log

	# Quantify
	quantifier.pl -p /genomeIndex/dre_hairpin+novel.fa -m dre_mature+novel.fa -P -r ${group}.pool.fa -t dre
done

