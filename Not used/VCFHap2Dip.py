#python
# a script that makes a haploid vcf diploid. Expects a whole genome VCF

import gzip
import csv
import argparse
import sys

parser = argparse.ArgumentParser(description="A script that makes a haploid vcf diploid. Expects a whole genome VCF.")
parser.add_argument("-v", "--vcf", action="store", required=True, help="Input VCF file. Should be a multisample, whole genome VCF from Shore with haploid calls.")
parser.add_argument("-o", "--out", action="store", required=True, help="Output filename")
parser.add_argument("-g", "--gzip", action="store_true", required=False, help="Set if the VCF is gzipped.")

args = parser.parse_args()

vcf_in = args.vcf
out_name = args.out

if args.gzip:
    opener = gzip.open
else:
    opener = open

with opener(vcf_in, 'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    vcf_out = csv.writer(open(out_name, 'w'), delimiter='\t', lineterminator="\n")
    for row in tsvin:
        if any('##' in strings for strings in row):
            vcf_out.writerow(row)
            continue
        if any('#CHROM' in strings for strings in row):
            vcf_out.writerow(row)
            continue
        chrom,pos,id,ref,alt,qual,filter,info,format=row[0:9]
        out = [chrom,pos,id,ref,alt,qual,filter,info,format]
        haplotypes = row[9:]
        for hap in haplotypes:
            call = hap.split(":")[0] # get the 0
            fixed_hap = call+"|"+call+":"+hap.split(":")[1]+":"+hap.split(":")[2]
            out.append(fixed_hap)
        vcf_out.writerow(out)