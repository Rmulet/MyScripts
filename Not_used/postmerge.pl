#!/usr/bin/perl

# Postmerge.pl, Version 0.8. This script is no longer in use because the same can be achieved with the experimental feature bcftools --missing-to-ref. It is faster and also overcomes one downside of this script: the variants not found in the mouse genotype are also replaced with 0/0 (that is, we assume they are monomorphic).

use strict;

my $merger = $ARGV[0]; # Merger human-mouse
my $output = $ARGV[1]; # Output in VCF format
my $individuals = -1;
my $counter = 0;

open (my $file, '<', $merger) or die "Cannot open input $merger: $!"; # Opens the 1000GP file (stored in $gpfile)
open(STDOUT, '>', $output) or die "Cannot create output file"; # Redirects output to a file in the same folder

while (<$file>) {
	my $line = $_;
	if ($line =~ m|^#|){ # Looks for a match that starts with '#' (headers)
		print $line;
	}
	elsif ($line =~ m|(.+?)AN=0\tGT\t((\./\.\t)+)?(\./\.$)|){ # Unknown genotype (./.) in mouse
		if ($individuals == -1) { # Count the number of individuals from the 1000 GP
			my @matches = ($2 =~ m|\./\.|g);
			$individuals = @matches;
		}
		print $1 . "AN=4522\t" . "GT\t" . "0/0\t" x $individuals . "./." . "\n";
		$counter++
	}
	elsif ($line =~ m|(.+?)AN=2(.+)(1/1$)|){ # Alternative allele (1/1) in mouse
		if ($individuals == -1) { # Count the number of individuals from the 1000 GP
			my @matches = ($2 =~ m|\./\.|g);
			$individuals = @matches;
		}	
		print $1 . "AN=4524;AC=2\t" . "GT\t" . "0/0\t" x $individuals . "1/1" . "\n";
		$counter++
	}
	else {
		print $line;
	}
}

close($file);
close(STDOUT);

print(sprintf ("Replaced %d alleles in a total of %d individuals",$counter,$individuals));
