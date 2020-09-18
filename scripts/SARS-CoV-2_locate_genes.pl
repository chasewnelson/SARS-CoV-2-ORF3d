#! /usr/bin/perl

# PROGRAM: Locate gene start and stop sites by finding the first sequences beginning
# and ending with the pre-specified sequences, hardcoded below, taken from the reference
# sequence: https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3

#########################################################################################
## LICENSE
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# AUTHOR: Chase W. Nelson
# Copyright (C) 2019 Chase W. Nelson
# DATE CREATED: April 2019

# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com

# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA

# CITATION1: OLGenie, https://github.com/chasewnelson/OLGenie
# CITATION2: Nelson CW, Ardern Z, Wei X. OLGenie: detecting natural selection to identify functional overlapping genes. In preparation.
# CITATION2: Nelson CW, Moncla LH, Hughes AL (2015) SNPGenie: estimating evolutionary 
#	parameters to detect natural selection using pooled next-generation sequencing data. 
#	Bioinformatics 31(22):3709-11, doi: 10.1093/bioinformatics/btv449.

use strict;
use Data::Dumper;

# Get the time
my $time1 = time;
my $local_time1 = localtime;

STDOUT->autoflush(1);

my $fasta_file = $ARGV[0];

unless(-f "$fasta_file") {
	die "\n### DIE: no fasta file.\n\n";
}

# Read in the group of sequences from the fasta file
my %header2sequence;
my $seq = '';
my $header = '';
my @headers_arr;
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file") or die "Could not open file $fasta_file\n";

while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$header =~ s/^>//; # get rid of FASTA header indicator
			$header =~ tr/\|/_/; # convert pipes (|) to underscores (_), following IQTree convention
			$header =~ s/\s.*$//; # trim anything including and after the first whitespace
			
			$seq_num ++;
		} else {
			$seq = uc($seq);
			$seq =~ tr/U/T/;
			
			if($header =~ /[^\w^\-^\.]/) {
				die "\n### TAXA NAMES CONTAIN INAPPROPRIATE CHARACTERS IN FASTA FILE: $header.\n" .
					"### Only alphanumeric characters (a-z, A-Z, 0-9), underscores (_), dashes (-), and periods (.) may be used. SCRIPT TERMINATED.\n\n";
			}
			
			if(exists $header2sequence{$header}) {
				die "\n\n### DIE: Each sequence in the FASTA multiple sequence alignment must have a unique header (up to first whitespace) for identification. TERMINATED.\n\n";
			} else {
				$header2sequence{$header} = $seq;
				push(@headers_arr, $header);
			}
			
			$header = $_;
			$header =~ s/^>//; # get rid of FASTA header indicator
			$header =~ tr/\|/_/; # convert pipes (|) to underscores (_), following IQTree convention
			$header =~ s/\s.*$//; # trim anything including and after the first whitespace
			
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\n### DIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
				$seq = '';
			}
		}
	} else {
		$seq .= $_;
	}
}

close IN_FASTA;

$seq = uc($seq);
$seq =~ tr/U/T/;

if($header =~ /^[^\w^\-^\.]/) {
	die "\n### TAXA NAMES CONTAIN INAPPROPRIATE CHARACTERS IN FASTA FILE: $header.\n" .
		"### Only alphanumeric characters (a-z, A-Z, 0-9), underscores (_), dashes (-), and periods (.) may be used. SCRIPT TERMINATED.\n\n";
}
$header2sequence{$header} = $seq;
push(@headers_arr, $header);


##########################################################################################
# Hard-code gene beginnings and endings
my %gene_to_termini;

# Hard-code

# ORF1ab: full gene
$gene_to_termini{'ORF1ab_1'}->{'start'} = 'ATGGAGAGCCTTGTCCCTGGTTTCAA';
$gene_to_termini{'ORF1ab_1'}->{'end'} = 'TCAGCTGATGCACAATCGTTTTTAAAC';

$gene_to_termini{'ORF1ab_2'}->{'start'} = 'CGGGTTTGCGGTGTAAGTGCAGCCCG';
$gene_to_termini{'ORF1ab_2'}->{'end'} = 'TAGTGATGTTCTTGTTAACAACTAA';

# ORF1ab_1: NOL (Wuhan-Hu-1 266-13465)
$gene_to_termini{'ORF1ab_1_NOL'}->{'start'} = 'ATGGAGAGCCTTGTCCCTGGTTTCAAC';
$gene_to_termini{'ORF1ab_1_NOL'}->{'end'} = 'AGCTGATGCACAATCGTTTTTA';

# ORF1ab: OL (Wuhan-Hu-1 13469-13483)
$gene_to_termini{'ORF1ab_1_OL_ORF1ab_2'}->{'start'} = 'GGGTTTGCGGTGTAA'; # whole thing
$gene_to_termini{'ORF1ab_1_OL_ORF1ab_2'}->{'end'} = 'GGGTTTGCGGTGTAA';

# ORF1ab_2: NOL (Wuhan-Hu-1 266-13465)
$gene_to_termini{'ORF1ab_2_NOL'}->{'start'} = 'GCAGCCCGTCTTACACCGTGCGGCACAGGC';#
$gene_to_termini{'ORF1ab_2_NOL'}->{'end'} = 'ATGTTCTTGTTAACAACTAA';

# ORF1a: full gene
$gene_to_termini{'ORF1a'}->{'start'} = 'ATGGAGAGCCTTGTCCCTGGTTTCAACGA';
$gene_to_termini{'ORF1a'}->{'end'} = 'CGGGTTTGCGGTGTAA';

# Nonstructural proteins: mature peptides
$gene_to_termini{'nsp1'}->{'start'} = 'ATGGAGAGCCTTGTCCCTGGTTTCAACGAGA';
$gene_to_termini{'nsp1'}->{'end'} = 'CGTGAACTCATGCGTGAGCTTAACGGAGGG';

$gene_to_termini{'nsp2'}->{'start'} = 'GCATACACTCGCTATGTCGATAACAACTTCT';
$gene_to_termini{'nsp2'}->{'end'} = 'ACAATACCTTCACACTCAAAGGCGGT';

$gene_to_termini{'nsp3'}->{'start'} = 'GCACCAACAAAG';
$gene_to_termini{'nsp3'}->{'end'} = 'CAAAGATAGCACTTAAGGGTGGT';

$gene_to_termini{'nsp4'}->{'start'} = 'AAAATTGTTAAT';
$gene_to_termini{'nsp4'}->{'end'} = 'CAAACCTCTATCACCTCAGCTGTTTTGCAG';

$gene_to_termini{'nsp5'}->{'start'} = 'AGTGGTTTTAGAAAAATGGCATTCCCATC';
$gene_to_termini{'nsp5'}->{'end'} = 'CAATGCTCAGGTGTTACTTTCCAA';

$gene_to_termini{'nsp6'}->{'start'} = 'AGTGCAGTGAAAAGAACAATCAAGGGTAC';
$gene_to_termini{'nsp6'}->{'end'} = 'CTTGTATCAAAGTAGCCACTGTACAG';

$gene_to_termini{'nsp7'}->{'start'} = 'TCTAAAATGTCAGATGTAAAGTGCACATCA';
$gene_to_termini{'nsp7'}->{'end'} = 'TGCTGGACAACAGGGCAACCTTACAA';

$gene_to_termini{'nsp8'}->{'start'} = 'GCTATAGCCTCAGAGTTTAGTTCCCTTCCA';
$gene_to_termini{'nsp8'}->{'end'} = 'CAATTCTGCTGTCAAATTACAG';

$gene_to_termini{'nsp9'}->{'start'} = 'AATAATGAGCTTAGTCCTGTTGCACTACGA';
$gene_to_termini{'nsp9'}->{'end'} = 'AGTTTAGCTGCCACAGTACGTCTACAA';

$gene_to_termini{'nsp10'}->{'start'} = 'GCTGGTAATGCAACAGAAGTGCCTGCCAA';
$gene_to_termini{'nsp10'}->{'end'} = 'GTGATCAACTCCGCGAACCCATGCTTCAG';

$gene_to_termini{'nsp11'}->{'start'} = 'TCAGCTGATGCACAATCGTTTTTAAAC';
$gene_to_termini{'nsp11'}->{'end'} = 'CGGGTTTGCGGTG';

$gene_to_termini{'nsp12_1'}->{'start'} = 'TCAGCTGATGCACAATCGTTTTTA';
$gene_to_termini{'nsp12_1'}->{'end'} = 'TCAGCTGATGCACAATCGTTTTTAAAC';

$gene_to_termini{'nsp12_2'}->{'start'} = 'CGGGTTTGCGGTGTAAGTGCAGCCC';
$gene_to_termini{'nsp12_2'}->{'end'} = 'TACACACCGCATACAGTCTTACAG';

$gene_to_termini{'nsp13'}->{'start'} = 'GCTGTTGGGGCTTGTGTTCTTTGCAATT';
$gene_to_termini{'nsp13'}->{'end'} = 'AATTCCACGTAGGAATGTGGCAACTTTACAA';

$gene_to_termini{'nsp14'}->{'start'} = 'GCTGAAAATGTAACAGGACTCTTTAAA';
$gene_to_termini{'nsp14'}->{'end'} = 'CTGGAACACTTTTACAAGACTTCAG';

$gene_to_termini{'nsp15'}->{'start'} = 'AGTTTAGAAAATGTGGCTTTTAATGTT';
$gene_to_termini{'nsp15'}->{'end'} = 'TAGAAACATTTTACCCAAAATTACAA';

$gene_to_termini{'nsp16'}->{'start'} = 'TCTAGTCAAGCGTGGCAACCGGGTGT';
$gene_to_termini{'nsp16'}->{'end'} = 'CTAGTGATGTTCTTGTTAACAACTAA';

# Structural and accessory
$gene_to_termini{'S'}->{'start'} = 'ATGTTTGTTTTTCTTGTT';
$gene_to_termini{'S'}->{'end'} = 'GTCAAATTACATTACACATAA';

$gene_to_termini{'ORF3a'}->{'start'} = 'ATGGATTTGTTTATGAGAATCTTCAC';
$gene_to_termini{'ORF3a'}->{'end'} = 'GACTACTAGCGTGCCTTTGTAA';

# ORF3a: NOL1 (Wuhan-Hu-1 25393-25521)
$gene_to_termini{'ORF3a_NOL1'}->{'start'} = 'ATGGATTTGTTTATGAGAATCTTCAC';
$gene_to_termini{'ORF3a_NOL1'}->{'end'} = 'GATACAAGCCTCACTCCCTTTC';

# ORF3a: OL (Wuhan-Hu-1 25522-25698)
$gene_to_termini{'ORF3a_OL_ORF3d'}->{'start'} = 'GGATGGCTTATTGTTGGCGTTGCA';
$gene_to_termini{'ORF3a_OL_ORF3d'}->{'end'} = 'CTCGTTGCTGCTGGCCTTGAA';

# ORF3a: NOL2 (Wuhan-Hu-1 25699-26220)
$gene_to_termini{'ORF3a_NOL2'}->{'start'} = 'GCCCCTTTTCTCTATCTTTATGCTTT';
$gene_to_termini{'ORF3a_NOL2'}->{'end'} = 'GACTACTAGCGTGCCTTTGTAA';

$gene_to_termini{'ORF3d'}->{'start'} = 'ATGGCTTATTGTTGGCGTTGCACTTC';
$gene_to_termini{'ORF3d'}->{'end'} = 'TTTGCTCGTTGCTGCTGGCCTTGA';

$gene_to_termini{'E'}->{'start'} = 'ATGTACTCATTCGTTTCGGAAGAGACAGGT';
$gene_to_termini{'E'}->{'end'} = 'GTTCCTGATCTTCTGGTCTAA';

$gene_to_termini{'M'}->{'start'} = 'ATGGCAGATTCCAACGGTACTATTACCGT';
$gene_to_termini{'M'}->{'end'} = 'GCTTTGCTTGTACAGTAA';

# M_OL_MB: TTATGGCCAGTAACTTTAGCTTGTTTTGTGCTT .. TCCATGTGGTCATTCAATCCAGAAACTAAC

$gene_to_termini{'ORF6'}->{'start'} = 'ATGTTTCATCTCGTTGACTTTCAGG';
$gene_to_termini{'ORF6'}->{'end'} = 'AGCAACCAATGGAGATTGATTAA';

$gene_to_termini{'ORF7a'}->{'start'} = 'ATGAAAATTATTCTTTTCTTGGCA';
$gene_to_termini{'ORF7a'}->{'end'} = 'CAAAAGAAAGACAGAATGA';

# ORF7a: NOL (Wuhan-Hu-1 27394-27753)
$gene_to_termini{'ORF7a_NOL'}->{'start'} = 'ATGAAAATTATTCTTTTCTTGGCA';
$gene_to_termini{'ORF7a_NOL'}->{'end'} = 'CACTCAAAAGAAAGACA';

$gene_to_termini{'ORF7b'}->{'start'} = 'ATGATTGAACTTTCATTAATTGAC';
$gene_to_termini{'ORF7b'}->{'end'} = 'GAAACTTGTCACGCCTAA';

# ORF7b: NOL (Wuhan-Hu-1: 27762-27887)
$gene_to_termini{'ORF7b_NOL'}->{'start'} = 'GAACTTTCATTAATTGACTTCT';
$gene_to_termini{'ORF7b_NOL'}->{'end'} = 'GAAACTTGTCACGCCTAA';

# ORF8 starts identically to ORF8'; probably the same
$gene_to_termini{'ORF8'}->{'start'} = 'ATGAAATTTCTTGTTTTCTTAGGAA';
$gene_to_termini{'ORF8'}->{'end'} = 'TTGTTTTAGATTTCATCTAA';

# ORF8a and ORF8b?

$gene_to_termini{'N'}->{'start'} = 'ATGTCTGATAATGGACCCCAAAATCAGCG';
$gene_to_termini{'N'}->{'end'} = 'GCTGACTCAACTCAGGCCTAA';

# N: NOL1 (Wuhan-Hu-1 28274-28282)
$gene_to_termini{'N_NOL1'}->{'start'} = 'ATGTCTGAT'; # whole thing
$gene_to_termini{'N_NOL1'}->{'end'} = 'ATGTCTGAT';

# N: OL ORF9b (Wuhan-Hu-1 28283-28579)
$gene_to_termini{'N_OL_ORF9b'}->{'start'} = 'AATGGACCCCAAAATCAG';
$gene_to_termini{'N_OL_ORF9b'}->{'end'} = 'GTGACGGTAAAATGAAA';

# N: NOL2 (Wuhan-Hu-1 28580-28732)
$gene_to_termini{'N_NOL2'}->{'start'} = 'GATCTCAGTCCAAGATGGTATTT';
$gene_to_termini{'N_NOL2'}->{'end'} = 'ACCCGCAATCCTGCTAAC';

# N: OL ORF9c (Wuhan-Hu-1 28733-28957)
$gene_to_termini{'N_OL_ORF9c'}->{'start'} = 'AATGCTGCAATCGTGCTACAACTTCCTC';
$gene_to_termini{'N_OL_ORF9c'}->{'end'} = 'CTGCTGCTTGACAGATTGAAC';

# N: NOL3 (Wuhan-Hu-1 28958-29533)
$gene_to_termini{'N_NOL3'}->{'start'} = 'CAGCTTGAGAGCAAAATGTCTGG';
$gene_to_termini{'N_NOL3'}->{'end'} = 'GCTGACTCAACTCAGGCCTAA';

$gene_to_termini{'ORF9b'}->{'start'} = 'ATGGACCCCAA'; #AATCAGCGAAATG';
$gene_to_termini{'ORF9b'}->{'end'} = 'ACCAGACGAATTCGTGGTGGTGACGGTAAAATGA';

$gene_to_termini{'ORF9c'}->{'start'} = 'ATGCTGCAATCGTGCTACAACTT';
$gene_to_termini{'ORF9c'}->{'end'} = 'CTGCTGCTTGACAGATTGA';

$gene_to_termini{'ORF10'}->{'start'} = 'ATGGGCTATATAAACGTTTTCGC';
$gene_to_termini{'ORF10'}->{'end'} = 'TAACTTTAATCTCACATAG';


#print "\n################################################################################";
#print "\nANALYSIS MODE=$mode\...\n";
my %gene_to_sites;

my @gene_names_sorted = qw/ORF1ab_1 ORF1ab_2 ORF1ab_1_NOL ORF1ab_1_OL_ORF1ab_2 ORF1ab_2_NOL 
						ORF1a nsp1 nsp2 nsp3 nsp4 nsp5 nsp6 nsp7 nsp8 nsp9 nsp10 nsp11 
						nsp12_1 nsp12_2 nsp13 nsp14 nsp15 nsp16 S ORF3a ORF3a_NOL1 
						ORF3a_OL_ORF3d ORF3a_NOL2 ORF3d E M ORF6 ORF7a ORF7a_NOL ORF7b ORF7b_NOL 
						ORF8p N N_NOL1 N_OL_ORF9b N_NOL2 N_OL_ORF9c N_NOL3 ORF9b ORF9c 
						ORF10/;
						
#my @gene_names_sorted = qw/ORF1ab_1 ORF1ab_2 ORF1a nsp1 nsp2 nsp3 nsp4 nsp5 nsp6 nsp7 nsp8 nsp9 nsp10 nsp11 nsp12_1 nsp12_2 
#						nsp13 nsp14 nsp15 nsp16 S S_confirm ORF3a ORF3d E M ORF6 ORF7a ORF7b ORF8p N ORF9b ORF9c ORF10/;


GENE_LOOP: foreach my $gene_name (@gene_names_sorted) {
#GENE_LOOP: foreach my $gene_name (sort keys %gene_to_termini) {
	my $start_seq = $gene_to_termini{$gene_name}->{'start'};
	my $end_seq = $gene_to_termini{$gene_name}->{'end'};
	
	SEQ_LOOP: foreach my $sequence_name (sort keys %header2sequence) {
		my $sequence = $header2sequence{$sequence_name};
		
		# Find index of start
		my $start_index = index($sequence, $start_seq);
		
		# Find index of end
		my $end_index = index($sequence, $end_seq);
		
		if($start_index > -1 && $end_index > -1) {
			# Compute positions
			my $start_site = $start_index + 1;
			my $end_site = $end_index + length($end_seq);
			
			# Print GTF line, go to next gene
			print "SARSCoV2_exhaustive\tSNPGenie\tCDS\t$start_site\t$end_site\t.\t+\t0\tgene_id \"$gene_name\";\n";
			
			# If S gene, print the 9 preceding codons
			if($gene_name eq 'S') {
				print "SARSCoV2_exhaustive\tSNPGenie\tCDS\t" . 
					($start_site - (9 * 3)) . "\t" .
					($start_site - 1) . "\t" .
					".\t+\t0\tgene_id \"$gene_name\_altStart\";\n";
			}
			
			next GENE_LOOP;
		}
		
	}
	
	print "\n### WARNING: didn't find $gene_name\n";
}

exit;


