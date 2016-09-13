#!/bin/bash

# GENOME ANALYZER 0.1 - Calls the VCFmerger.sh script for the entire genome

# Add DOWNLOAD option

WINDOW=10000

gpraw="/home/roger/Documents/2_GenomicsData/1000GP/Chromosomes"
gpdat="/home/roger/Documents/2_GenomicsData/1000GP/"
alraw="/home/roger/Documents/2_GenomicsData/Alns/Chromosomes"
aldat="/home/roger/Documents/2_GenomicsData/Alns/"

for i in `seq 1 22` X Y; do
	gpfile=ls -d *chr$i.*vcf.gz
	alnfile=chr$i_aln.mfa

	bcftools view -Oz -S ^InbredIndividuals.txt $gpfile > chr$i.vcf.gz # Filter out inbred individuals
	tabix -p chr$i.vcf.gz
	
	grep -v "score" chr$i.mfa > $alnfile # Remove 'score' lines from the MFA file
	MFAtoVCF.py -s chr$i.fa -q Chimp -c $i $alnfile # Convert MFA to VCF. VCF.GZ file is automatically tabixed.

	VCFmerger.sh $MASK chr$i.vcf.gz $chr$i_aln.vcf.gz -w $WINDOW -db Genomics$i # Merges and analyzes variation data
