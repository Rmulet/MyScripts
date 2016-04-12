#!/usr/bin/perl

use strict;

my $merger = $ARGV[0]; # Merger human-mouse
my $individuals = -1;
my $counter = 0;

open (my $file, '<', $merger) or die "Cannot open input $merger: $!"; # Opens the 1000GP file (stored in $gpfile)
#open(STDOUT, '>', './output.vcf') or die "Cannot create output file"; # Redirects output to a FASTA file in the same folder

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

sprintf ("Replaced %d alleles in a total of %d individuals",$counter,$individuals)