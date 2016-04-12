#!/usr/bin/perl

my $gpfile = $ARGV[0]; # 1000 GP file
my $alnfile = $ARGV[1] # Alignment file
my $merge = $ARGV[2] # Result of the merger

open (GPFILE, '<', $gpfile) or die "Cannot open input file";; # Opens the 1000GP file (stored in $gpfile)
open (ALNFILE, '<', $alnfile) or die "Cannot open input file";; # Opens the alignment file (stored in $alnfile)
open(STDOUT, '>', './output.vcf') or die "Cannot create output file"; # Redirects output to a FASTA file in the same folder

while (<GPFILE>) {
	my $line = $_;
	if ($line =~ /^#/){ # Looks for a match that starts with '>' followed by any character (identifer). That excludes white and sequence lines.
