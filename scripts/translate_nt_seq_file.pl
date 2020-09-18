#! /usr/bin/env perl

# PROGRAM: Translate a file of (un-aligned) nucleotide sequences, print to protein file

# Copyright (C) 2016 Chase W. Nelson

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# DATE CREATED: January 18, 2016
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: Special Volunteer, Division of Cancer Epidemiology & Genetics, National
#     Cancer Institute, National Institutes of Health, Rockville, MD 20850, USA
# AFFILIATION3: BigPlant Consortium, Center for Genomics and Systems Biology, New York 
#     University, New York, NY 10003, USA

# Takes in as an argument:
#	[0] one nucleotide FASTA file containing the original (un-aligned) nucleotide sequences.

use strict;
#use warnings;
use Data::Dumper;

my $fasta_file_name = $ARGV[0];
unless($fasta_file_name =~ /.fa/) { die "\n\n# Nucleotide FASTA file must contain .fa " . 
	"or .fasta extension. TERMINATED\n\n"; }

# Read in the group of sequences from the fasta file
my $seq = '';
my @seqs_arr;
my $header = '';
my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file_name") or die "Could not open file $fasta_file_name\n";

print "\nRecording coding sequence data for $fasta_file_name...\n";

while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$seq_num ++;
		} else {
			push(@seqs_arr,$seq);
			push(@headers_arr,$header);
			$header = $_;
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			$last_seq_length = $this_seq_length;
			$seq = '';
			
			
		}
	} else {
		$seq .= $_;
	}
}

close IN_FASTA;

push(@seqs_arr,$seq);
push(@headers_arr,$header);

print "\n";

# Translate sequences and print to file
my $translated_filename;
if($fasta_file_name =~ '.fasta') {
	$translated_filename = $` . "_aa.fasta"; # translated
} elsif($fasta_file_name =~ '.fa') {
	$translated_filename = $` . "_aa.fa"; # translated
} else {
	$translated_filename = "translated.fa";
}

open(OUT, ">>$translated_filename");

for(my $i=0; $i<@seqs_arr; $i++) {
	my $curr_seq = $seqs_arr[$i];
	my $curr_length = length($curr_seq);
	my $curr_header = $headers_arr[$i];
	
	if ($curr_length < 3) {
		die "\n\n## WARNING: sequence $curr_header is less than 3 nucleotides.\n\n";
	}
	
	if(($curr_length % 3) > 0) {
		print "\n\n## WARNING: $curr_header length is not a multiple of 3 (complete codons).\n\n";
	}
	
	my $amino_acid_seq;
	
	# Construct (translate) amino acid sequence
	for(my $index = 0; $index < $curr_length; $index+=3) {
		my $codon = substr($curr_seq, $index, 3);
		$amino_acid_seq .= &get_amino_acid($codon);
	}
	
	my $amino_acid_seq_length = length($amino_acid_seq);
	
	my $new_header = $curr_header;
	
	print OUT "$new_header\n";
	
	for(my $j=0; $j<($amino_acid_seq_length); $j+=60) {
		if($j>($amino_acid_seq_length - 60)) {
			my $line = substr($amino_acid_seq, $j);
			print OUT "$line\n";
			last;
		} else {
			my $line = substr($amino_acid_seq, $j, 60);
			print OUT "$line\n";
		}
	}
	
}

close OUT;


#########################################################################################
#########################################################################################
##################################  SUBROUTINES  ########################################
#########################################################################################
#########################################################################################

#########################################################################################
# Get the amino acid (single-letter code) encoded by a given DNA or RNA codon
sub get_amino_acid {
	my ($codon) = @_;
	$codon = uc($codon); # uc returns uppercase
	$codon =~ tr/U/T/;
	my $amino_acid;
	
	# Establish genetic code for use with synonymous sites; DNA or RNA
	my %code = (
		"AAA" => "K",
		"AAC" => "N",
		"AAG" => "K",
		"AAT" => "N",
		"AAU" => "N",
		"ACA" => "T",
		"ACC" => "T",
		"ACG" => "T",
		"ACT" => "T",
		"ACU" => "T",
		"AGA" => "R",
		"AGC" => "S",
		"AGG" => "R",
		"AGT" => "S",
		"AGU" => "S",
		"ATA" => "I",
		"ATC" => "I",
		"ATG" => "M",
		"ATT" => "I",
		"AUA" => "I",
		"AUC" => "I",
		"AUG" => "M",
		"AUU" => "I",
		"CAA" => "Q",
		"CAC" => "H",
		"CAG" => "Q",
		"CAT" => "H",
		"CAU" => "H",
		"CCA" => "P",
		"CCC" => "P",
		"CCG" => "P",
		"CCT" => "P",
		"CCU" => "P",
		"CGA" => "R",
		"CGC" => "R",
		"CGG" => "R",
		"CGT" => "R",
		"CGU" => "R",
		"CTA" => "L",
		"CTC" => "L",
		"CTG" => "L",
		"CTT" => "L",
		"CUA" => "L",
		"CUC" => "L",
		"CUG" => "L",
		"CUU" => "L",
		"GAA" => "E",
		"GAC" => "D",
		"GAG" => "E",
		"GAT" => "D",
		"GAU" => "D",
		"GCA" => "A",
		"GCC" => "A",
		"GCG" => "A",
		"GCT" => "A",
		"GCU" => "A",
		"GGA" => "G",
		"GGC" => "G",
		"GGG" => "G",
		"GGT" => "G",
		"GGU" => "G",
		"GTA" => "V",
		"GTC" => "V",
		"GTG" => "V",
		"GTT" => "V",
		"GUA" => "V",
		"GUC" => "V",
		"GUG" => "V",
		"GUU" => "V",
		"TAA" => "*",
		"TAC" => "Y",
		"TAG" => "*",
		"TAT" => "Y",
		"UAA" => "*",
		"UAC" => "Y",
		"UAG" => "*",
		"UAU" => "Y",
		"TCA" => "S",
		"TCC" => "S",
		"TCG" => "S",
		"TCT" => "S",
		"UCA" => "S",
		"UCC" => "S",
		"UCG" => "S",
		"UCU" => "S",
		"TGA" => "*",
		"TGC" => "C",
		"TGG" => "W",
		"TGT" => "C",
		"UGA" => "*",
		"UGC" => "C",
		"UGG" => "W",
		"UGU" => "C",
		"TTA" => "L",
		"TTC" => "F",
		"TTG" => "L",
		"TTT" => "F",
		"UUA" => "L",
		"UUC" => "F",
		"UUG" => "L",
		"UUU" => "F",
	);
	
	$amino_acid = $code{$codon};
	
	if($amino_acid eq '') {
		$amino_acid = 'X';
	}
	
	return $amino_acid;
}

