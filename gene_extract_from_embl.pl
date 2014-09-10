#!/usr/bin/perl 
# gene_extract_from_embl.pl
use strict; use warnings;

die "usage: gene_extract_from_embl.pl <embl file> <search_term> > <output.fa>\n" unless @ARGV == 2;

my $search = $ARGV[1];  

my $seq = $ARGV[0];

#extract genome from genbank file
my $genome = ();

#use IN for input to open file
open (my $in, "<$ARGV[0]") or die "error reading $ARGV[0]: $!";
#while loops through only input
while (my $line = <$in>) {
	if ($line =~ /^SQ\s+/) {
		$line = <$in>;
		do {
			$line =~ s/\d+//g;
			$line =~ s/\s//g;
			$line =~ s/\n//g;
			$genome .= substr($line, 0, length($line));
			$line = <$in>;
		} until $line =~ /\/\//;
	}
}
close $in;

open ($in, "<$ARGV[0]") or die "error reading $ARGV[0]: $!";
#search for gene
while (my $line = <$in>) {
	if ($line =~ /^FT\s{3}gene/) {
		while ($line = <$in>) {
			if ($line =~ /\/gene="\S*$search\S*"/) {
				my ($name) = $line =~ m/"(.+)"/;
			$line = <$in>;
			if ($line =~ /^FT\s{3}CDS\s+/) {
				my ($beg, $end) = $line =~ /(\d+)\.\.(\d+)/;
				my $tru_seq = ();			
				#for CDS on reverse frame
				if ($line =~ /complement/) {
					my $seq = substr($genome, $beg - 1, ($end - $beg) + 1);
					#reverse the DNA
					$tru_seq = reverse $seq;
					#complement the DNA
					$tru_seq =~ tr/ATCGatcg/TAGCtagc/;	
				} 
				#for CDS on forward frame 
				else { 
				$tru_seq = substr($genome, $beg -1, ($end - $beg) + 1);
				}
				print ">$name\t$seq\n$tru_seq\n";
			}
			}
		}	
	}
}
close $in; 
