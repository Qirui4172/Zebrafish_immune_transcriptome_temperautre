#!/bin/sh

Samples=("t24ct2" "t24ct3" "t24ct4" "t24ct5" "t24ct8" "t24lps1" "t24lps2" "t24lps3" "t24lps4" "t24lps8" "t28ct1" "t28ct2" "t28ct3" "t28ct5" "t28ct7" "t28lps1" "t28lps2" "t28lps3" "t28lps4" "t28lps8" "t32ct2" "t32ct3" "t32ct6" "t32ct8" "t32ct10" "t32lps1" "t32lps2" "t32lps3" "t32lps4" "t32lps5")

for sample in ${Samples[@]};do
	cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAA -q 20 --quality-base=33 -m 20 --trim-n -o ${sample}.tmp1.fastq ${sample}_S.*_R1_001.fastq

	cutadapt -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -q 20 --quality-base=33 -m 20 --trim-n -o ${sample}.clean.fastq ${sample}.tmp.fastq
done

