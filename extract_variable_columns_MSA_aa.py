#! /usr/bin/env python

###############################################################################
## LICENSE
##
## Copyright (C) 2020
##
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
###############################################################################

# DATE: 2020
# AUTHOR: Chase W. Nelson
# CONTACT: cnelson@gate.sinica.edu.tw
# CITATION: https://github.com/chasewnelson/SARS-CoV-2-ORF3d

from Bio import SeqIO, AlignIO
import re
import sys
import os

Usage = "Script for examining FASTA (arg 1 file) to identify variant sites with minor alleles >= a given frequency (arg 2)"

# Check two arguments
if not len(sys.argv) == 3:
    raise SystemExit("\n### TERMINATED: need 2 unnamed arguments: (1) infile MSA and (2) minor allele freq cutoff\n")
infile = sys.argv[1]
minor_allele_freq = float(sys.argv[2])

# Create outfile name
outfile = str(infile)
outfile = outfile.replace(".fasta", "")
outfile = outfile + "_variants_MAF" + str(minor_allele_freq) + ".fasta"
outfile_name = str(outfile)

outfile_tsv = str(infile)
outfile_tsv = outfile_tsv.replace(".fasta", "")
outfile_tsv = outfile_tsv + "_variants_MAF" + str(minor_allele_freq) + ".tsv"

print("Will write output to: " + outfile + " and " + outfile_tsv)

# Check infile exists
if not os.path.isfile(infile):
    raise SystemExit("\n### TERMINATED: infile does not exist\n")

# Open infile as multiple sequence alignment (MSA)
alignment = AlignIO.read(infile, 'fasta')
print("Number of sequences: %i" % len(alignment))
aln_num_seqs = len(alignment)
aln_num_cols = alignment.get_alignment_length()

alignment_names = []

for rec in alignment:
    alignment_names.append(str(rec.id))


############################################################################
### EXAMINE EACH COLUMN; IF VARIABLE, KEEP AND LABEL
# Loop through all columns and count the number of gaps
variable_cols = []
invariable_cols = []
for i in range(aln_num_cols):  # number of columns
    # print(i)
    this_col_letters = alignment[:, i].upper()
    this_col_dict = {}

    for this_letter in list(this_col_letters):
        # print(this_letter)
        if this_col_dict.get(this_letter) is None:
            this_col_dict[this_letter] = 1
        else:
            this_col_dict[this_letter] += 1

    num_alleles = 0
    minor_allele = 'X'
    minor_allele_count = sum(this_col_dict.values())
    defined_allele_count = 0
    for key in this_col_dict.keys():
        # print(key + "=" + str(this_col_dict.get(key)))
        if key != "-" and key != "X" and key != "*":
            num_alleles += 1
            defined_allele_count += this_col_dict[key]

            if this_col_dict[key] < minor_allele_count:
                minor_allele = key
                minor_allele_count = this_col_dict[key]

    # print(str(num_alleles))

    if num_alleles > 1 and ((minor_allele_count / aln_num_seqs) >= minor_allele_freq):  # variable column
        variable_cols.append(i)
    else:
        invariable_cols.append(i)

print("Variable columns: " + str(variable_cols))
# print("Invariable columns: " + str(invariable_cols))


# INCLUDE ONLY columns that are variable
for i in reversed(range(aln_num_cols)):
    if (i in invariable_cols):
        alignment = alignment[:, :i] + alignment[:, i + 1:]

# Write FASTA output
with open(outfile, 'w') as outfile:
    outfile.writelines(alignment.format('fasta'))  # new_align

# Open and write header to TSV output
OUTFILE_TSV = open(outfile_tsv, "w")
OUTFILE_TSV.write("ID")
for this_site_index in variable_cols:
    this_site = this_site_index + 1
    OUTFILE_TSV.write("\t" + str(this_site))
OUTFILE_TSV.write("\n")

# Open FASTA to write TSV columns
recs = SeqIO.parse(outfile_name, "fasta")
for rec in recs:
    this_line = str(rec.id)
    this_seq = list(str(rec.seq).upper())

    for this_col_nt in this_seq:
        this_line = this_line + "\t" + this_col_nt

    OUTFILE_TSV.write(this_line + "\n")

OUTFILE_TSV.close()


