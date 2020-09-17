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

from Bio import SeqIO
import os
import sys

Usage = "Script for extracting a subset of sequences (arg 1 file, single column) from a FASTA (arg 2 file)"

# Check two argument
if not len(sys.argv) == 3:
    raise SystemExit("\n### TERMINATED: need 2 unnamed arguments: (1) infile with 1 column of IDs; and (2) FASTA\n")
infile_IDs = sys.argv[1]
infile_FASTA = sys.argv[2]

# Create outfile name
outfile = str(infile_IDs)
outfile = outfile.replace(".txt", "")
outfile = outfile.replace(".tsv", "")
outfile = outfile + "_seqs.fasta"
print("Will write output to: " + outfile + "\n")

# Check infiles exist
if not os.path.isfile(infile_IDs):
    raise SystemExit("\n### TERMINATED: infile 1 (IDs) does not exist\n")

if not os.path.isfile(infile_FASTA):
    raise SystemExit("\n### TERMINATED: infile 2 (FASTA) does not exist\n")

# Store IDs
desired_IDs = []
with open(infile_IDs, "r") as f:
    lines = (line.rstrip() for line in f)

    for line in lines:
        desired_IDs.append(line)

# Open FASTA and get desired records
recs = SeqIO.parse(infile_FASTA, "fasta")
outfile = open(outfile, "w")

for rec in recs:
    ID = rec.id

    if ID in desired_IDs:
        seq = rec.seq
        seq = str(seq)
        seq = seq.upper()

        # Write to FASTA file
        outfile.write(">" + ID + "\n")
        outfile.write(str(seq) + "\n")

outfile.close()


