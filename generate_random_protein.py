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

import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Usage = "Script to generate random protein sequences of a given length from an input proteome in FASTA format"

# Check two arguments
if not len(sys.argv) == 4:
    raise SystemExit("\n### TERMINATED: need 3 unnamed arguments:\n" + \
                     "(1) infile peptide MSA (FASTA)\n" + \
                     "(2) length of peptide to generate\n" + \
                     "(3) number of peptides to generate\n")

# Input
proteome_file = sys.argv[1]
peptide_length = sys.argv[2]
peptide_length = int(peptide_length)
num_peptides = sys.argv[3]
num_peptides = int(num_peptides)

# Check infile exists
if not os.path.isfile(proteome_file):
    raise SystemExit("\n### TERMINATED: infile does not exist\n")

# Create outfile name
outfile_name = str(proteome_file)
outfile_name = outfile_name.replace(".fasta", "")
outfile_name = outfile_name + "_random_l" + str(peptide_length) + "n" + str(num_peptides) + ".fasta"
print("Will write output to: " + outfile_name)

# Open infile, e.g., import SARS-CoV-2 (Wuhan-Hu-1) ORF3d amino acid sequence
recs = SeqIO.parse(proteome_file, "fasta")

proteome = ""
for rec in recs:
    proteome = proteome + str(rec.seq)

# Produce a FASTA of random peptides using the proteome
# Following: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc284
with open(outfile_name, "w") as outfile_hdl:
    for i in range(num_peptides):
        this_random_peptide = random.choices(proteome, k = peptide_length)
        shuffled_rec = SeqRecord(
            Seq("".join(this_random_peptide), rec.seq.alphabet),
            id = "s%i" % (i + 1),
            description = "random peptide drawn from %s" % proteome_file,
        )
        outfile_hdl.write(shuffled_rec.format("fasta"))

