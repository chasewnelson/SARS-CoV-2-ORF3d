#! /usr/bin/perl

# PROGRAM: takes in one aligned FAST MSA; outputs haplotypes

# Copyright (C) 2020 Chase W. Nelson

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

# DATE CREATED: May 2020
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# CONTACT2: cwnelson88@gmail.com
# AFFILIATION1: Sackler Institute for Comparative Genomics, American Museum of Natural
#     History, New York, NY 10024, USA
# AFFILIATION2: Special Volunteer, Division of Cancer Epidemiology & Genetics, National
#     Cancer Institute, National Institutes of Health, Rockville, MD 20850, USA
# AFFILIATION3: BigPlant Consortium, Center for Genomics and Systems Biology, New York 
#     University, New York, NY 10003, USA

use strict;
#use warnings;
use Data::Dumper;


#########################################################################################
# INITIALIZE VARIABLES (all optional and/or defaulted)
my $fasta_file;

if(@ARGV == 1) {
	$fasta_file = $ARGV[0];
} else {
	die "\n### WARNING: must be only one unnamed argument: a multiple sequence alignment in FASTA format. TERMINATED\n\n";
}

# Get the time
my $time1 = time;
my $local_time1 = localtime;

#print "\n##########################################################################################";
#print "\nHaplotype analysis initiated at local time $local_time1";
#print "\n##########################################################################################\n\n";

unless(-f "$fasta_file") {
	die "\n### WARNING: FASTA input file does not exist.\n".
		"### Script terminated.\n\n";
}


#########################################################################################
# STORE FASTA SEQUENCE & FULL-SEQUENCE HAPLOTYPES
my %haplotype_data;

my $header = '';
my $seq = '';
my $seq_num = 0;
my $last_seq_length;

open(IN_FASTA, "$fasta_file") or die "Could not open FASTA file $fasta_file\n";

#print "\n\nReading in sequences from FASTA...\n\n";

while(<IN_FASTA>) {
	chomp;
	if(/>/) {
		if($seq_num == 0) {
			$header = $_;
			$seq_num++;
		} else { # finished a new sequence
			$seq = uc($seq);
			$haplotype_data{$seq}->{count}++;
			
			if($haplotype_data{$seq}->{count} == 1) {
				$haplotype_data{$seq}->{first_ID} = $header;
			}
			
			# New header
			$header = $_;
			$seq_num ++;
			
			my $this_seq_length = length($seq);
			
			#print "\nseq $seq_num is of length $this_seq_length\n";
			
			if($last_seq_length && ($last_seq_length != $this_seq_length)) {
				die "\n\nDIE: The sequences must be aligned, i.e., must be the same length. TERMINATED.\n\n";
			} else {
				$last_seq_length = $this_seq_length;
				#print "\nseq: $seq\n";
				$seq = '';
			}
		}
	} else {
		$seq .= $_;
	}
}

# Finished last sequence
$seq = uc($seq);
$haplotype_data{$seq}->{count}++;
			
if($haplotype_data{$seq}->{count} == 1) {
	$haplotype_data{$seq}->{first_ID} = $header;
}

close IN_FASTA;

#$seq_num = scalar(@seqs_arr);
my $seq_length = length($seq);


#########################################################################################
# PRINT haplotype data
foreach my $haplotype (sort keys %haplotype_data) {
	print ">" . $haplotype_data{$haplotype}->{first_ID} . "_n" . $haplotype_data{$haplotype}->{count} . "\n";
	print "$haplotype\n";
}

#print "\n\n";
#print "\nHERE WE END WITHOUT ADO.\n\n";

exit;


#########################################################################################
#########################################################################################
####################################                 ####################################
####################################   SUBROUTINES   ####################################
####################################                 ####################################
#########################################################################################
#########################################################################################


#########################################################################################
sub get_product_names_from_gtf {
	my ($cds_file) = @_;
	#print "\n\n$cds_file\n\n";
	my %products_hash;
	open (CURRINFILE, $cds_file) or die "\n## Cannot open $cds_file. TERMINATED.\n\n";
	while (<CURRINFILE>) {
		if($_ =~ /CDS\t\d+\t\d+\t[\.\d+]\t\+/) { # Must be on the + strand
			if($_ =~/gene_id \"gene\:([\w\s\.\-\:']+)\"/) { # transcript_id not a problem
				$products_hash{$1} = 1;
			} elsif($_ =~ /gene_id \"([\w\s\.\-\:']+ [\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
			} elsif($_ =~/gene_id \"([\w\s\.\-\:']+)\"/) {
				$products_hash{$1} = 1;
			} else {
				die "\n\n## WARNING: CDS annotation(s) in $cds_file does not have a ".
					"gene_id. SNPGenie terminated.\n\n";
			}
		}
	}
	close CURRINFILE;
	
	my @product_names = keys %products_hash;
	#print "\n@product_names\n\n";
	return @product_names;
}

#########################################################################################
sub fisher_yates_shuffle {
	my $array = shift;
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}

#########################################################################################
sub translate_nt_seq { # also requires &get_amino_acid
	my ($nt_seq) = @_;
	
	my $input_length = length($nt_seq);
	
	if ($input_length < 3) {
		die "\n\n## WARNING: All coding sequences be at least 3 nucleotides.\n\n";
	}
	
	if(($input_length % 3) > 0) {
		die "\n\n## WARNING: All coding sequences must be a multiple of 3 (complete codons).\n\n";
	}
	
	#print "len is ".length($input_seq)."\n$input_seq";
	
	my $amino_acid_seq;
	
	for(my $index = 0; $index < $input_length; $index+=3) {
		my $codon = substr($nt_seq, $index, 3);
		$amino_acid_seq .= &get_amino_acid($codon);
	}
	
#	print "\nAmino Acid Sequence:\n$amino_acid_seq\n\n";
	
	return $amino_acid_seq;
	
}

#########################################################################################
# Get the amino acid (single-letter code) encoded by a given DNA or RNA codon
# Returns an array with:
#	returned[0] = number of nonsynonymous sites
#	returned[1] = number of synonymous sites
sub get_amino_acid {
	#my ($codon) = @_;
	my ($codon) = uc($_[0]); # uc returns uppercase
	$codon =~ tr/U/T/;
	
	my $amino_acid;
	
	# Establish genetic code for use with synonymous sites; DNA or RNA
	my %code = (
		"AAA"=>"K","AAC"=>"N","AAG"=>"K","AAT"=>"N","ACA"=>"T","ACC"=>"T","ACG"=>"T",
		"ACT"=>"T","AGA"=>"R","AGC"=>"S","AGG"=>"R","AGT"=>"S","ATA"=>"I","ATC"=>"I",
		"ATG"=>"M","ATT"=>"I","CAA"=>"Q","CAC"=>"H","CAG"=>"Q","CAT"=>"H","CCA"=>"P",
		"CCC"=>"P","CCG"=>"P","CCT"=>"P","CGA"=>"R","CGC"=>"R","CGG"=>"R","CGT"=>"R",
		"CTA"=>"L","CTC"=>"L","CTG"=>"L","CTT"=>"L","GAA"=>"E","GAC"=>"D","GAG"=>"E",
		"GAT"=>"D","GCA"=>"A","GCC"=>"A","GCG"=>"A","GCT"=>"A","GGA"=>"G","GGC"=>"G",
		"GGG"=>"G","GGT"=>"G","GTA"=>"V","GTC"=>"V","GTG"=>"V","GTT"=>"V","TAA"=>"*",
		"TAC"=>"Y","TAG"=>"*","TAT"=>"Y","TCA"=>"S","TCC"=>"S","TCG"=>"S","TCT"=>"S",
		"TGA"=>"*","TGC"=>"C","TGG"=>"W","TGT"=>"C","TTA"=>"L","TTC"=>"F","TTG"=>"L",
		"TTT"=>"F"
	);
	
	$amino_acid = $code{$codon};
	
	if($amino_acid eq '') {
		$amino_acid = 'X'; # previously '?'
	}
	
	return $amino_acid;
}

#########################################################################################
sub get_header_names {
	# Originally assumed that we've received a tempfile ending in "_snpg9temp.txt"
	# However, we're now calling it at least once before creating the tempfile to 
	# see what kind of processing (e.g., Geneious to CLC) is needed prior to tempfile
	# creation. Must include capability to get headers for .CSV file
	my ($curr_snp_report_filename,$filename) = @_;
	#print "\n$curr_snp_report_filename\n";
	
	#my $newline_char = &detect_newline_char($curr_snp_report_filename);
	#my $old_newline = $/;
	#$/ = $newline_char;
	
	my $seen_tab_delimited = 0;
	my $seen_comma_delimited = 0;
	my $seen_vcf_tab_delimited = 0;
	my @line_arr;
	
	my $line = 0;
	open (CURRINFILE, $curr_snp_report_filename);
	#seek(CURRINFILE,0,0);
	while (<CURRINFILE>) {
		#print "$_";
		if($line == 0) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t\w+\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,\w+,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
				#print "COMMA!!!!!";
				last;
			}

			$line++;
		} elsif($line > 0 && $_ =~ /^##/) {
			$line++;
		} elsif($line > 0 && ($_ =~ /^#CHROM/)) {
			chomp;
			# CHOMP for 3 operating systems
			if($_ =~ /\r\n$/) {
				$_ =~ s/\r\n//;
			} elsif($_ =~ /\r$/) {
				$_ =~ s/\r//;
			} elsif($_ =~ /\n$/) {
				$_ =~ s/\n//;
			}
			
			if($_ =~/\t/) { # it's TAB-delimited
				@line_arr = split(/\t/,$_);
				#print "TAB!!!!!";
				last;
			} elsif($_ =~/,/) { # it's COMMA-delimited
				@line_arr = split(/,/,$_);
				#print "COMMA!!!!!";
				last;
			} else {
				chdir('SNPGenie_Results');
				open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
				# FILE | PRODUCT | SITE | CODON | WARNING
				
				# No change OR error should occur if the file does not, in fact, end
				# with this SUFFIX
				my $file_nm = $curr_snp_report_filename;
				#$file_nm =~ s/_snpg9temp.txt/.txt/;
				$file_nm =~ s/_\w\w\w\w.txt/.txt/;
				
				print ERROR_FILE "$filename\tN/A\tN/A\t".
					"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
				close ERROR_FILE;
				chdir('..');
				
				#unlink $curr_snp_report_filename;
				
				die "\n\n## WARNING: The SNP Report $filename is ".
					"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
					"terminated\n\n";
			}
		} else {
			chdir('SNPGenie_Results');
			open(ERROR_FILE,">>SNPGenie\_LOG\.txt");
			# FILE | PRODUCT | SITE | CODON | WARNING
			
			# No change OR error should occur if the file does not, in fact, end
			# with this SUFFIX
			my $file_nm = $curr_snp_report_filename;
			#$file_nm =~ s/_snpg9temp.txt/.txt/;
			$file_nm =~ s/_\w\w\w\w.txt/.txt/;
			
			print ERROR_FILE "$filename\tN/A\tN/A\t".
				"File not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPGenie terminated.\n";
			close ERROR_FILE;
			chdir('..');
			
			#unlink $curr_snp_report_filename;
			
			die "\n\n## WARNING: The SNP Report $filename is ".
				"not TAB(\\t)- or COMMA-delimited, or there is only one column. SNPgenie ".
				"terminated\n\n";
		}
	}
	seek(CURRINFILE,0,0);
	close CURRINFILE;
	#$/ = $old_newline;
	return @line_arr;
}
