#!/usr/bin/perl

my $gpfile = $ARGV[0]; # 1000 GP file
my $alnfile = $ARGV[1] # Alignment file
my $merge = $ARGV[2] # Result of the merger

open (GPFILE, '<', $gpfile) or die "Cannot open input file";; # Opens the 1000GP file (stored in $gpfile)
open (ALNFILE, '<', $alnfile) or die "Cannot open input file";; # Opens the alignment file (stored in $alnfile)

while <GPFILE> 