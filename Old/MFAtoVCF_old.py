import argparse

#######################
## PARSING ARGUMENTS ##
#######################

parser = argparse.ArgumentParser(description='Obtain a VCF file from a MFA using SNP-sites')
parser.add_argument('-i','--input',type=str,required=True,help="Input file, in MFA format")
parser.add_argument('-o','--output',type=str,default="output.vcf",help="Output file, in VCF format")
parser.add_argument('-s','--start',type=int,default=10000,help="Window size for metrics calculation (default 10000)")
parser.add_argument('-c','--chromosome',default="window",help="Step size for metrics calculation (default 0)")

args = parser.parse_args()

#####################################
## CREATING VCF FILES BY SEQ PAIRS ##
#####################################

with open(args.input) as file:
	for 



