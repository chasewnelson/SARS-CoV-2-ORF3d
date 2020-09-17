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
from datetime import date, datetime, time
import sys
import os

Usage = "Script to extract SARS-CoV-2 GISAID sequences by location"

# Check argument number
if not len(sys.argv) == 5:
    raise SystemExit("\n### TERMINATED: need 3 unnamed arguments:\n" + \
                     "    # (1) .TSV with date metadata\n" + \
                     "    # (2) .FASTA alignment\n" + \
                     "    # (3) name of desired output file (.fasta)\n" + \
                     "    # (4) string to detect in Location field\n\n")

infile_metadata = sys.argv[1]
infile_FASTA = sys.argv[2]
outfile_name = sys.argv[3]
location_keyword = sys.argv[4]

print("Will extraction sequences with Location keyword: " + location_keyword + "\n")

# Check infiles exist
if not os.path.isfile(infile_metadata):
    raise SystemExit("\n### TERMINATED: infile 1 (IDs) does not exist\n")

if not os.path.isfile(infile_FASTA):
    raise SystemExit("\n### TERMINATED: infile 2 (FASTA) does not exist\n")

print("Will write output to: " + outfile_name + "\n")

# Build a dictionary of locations
ID_to_location = {}

# for now, ASSUME col 1 is ID, and col 3 (idx 2) is location
# location examples:
# North America / USA / Washington / King County
# North America / USA / New Hampshire
# Asia / China / Guandong / Shenzhen
# Asia / China / Guangdong

line_num = 0
with open(infile_metadata, "r") as f:
    lines = (line.rstrip() for line in f)

    for line in lines:
        line_num += 1
        line_list = line.split("\t")
        #this_date = datetime(2019, 12, 20)

        if line_num == 1:
            # print("Accession ID")
            # print(line_list[0])

            if not line_list[0] == "Accession ID" and line_list[2] == "Location":
                raise SystemExit(
                    "\n### TERMINATED: header must be col1=\"Accession ID\" and col3=\"Location\"\n")
        else:
            location_worked = False
            seq_ID = line_list[0]
            seq_location = line_list[2]

            if seq_location != "":
                ID_to_location[seq_ID] = seq_location
            else:
                print("Empty location: seq_ID=" + str(seq_ID))

print("\n")


# Open FASTA for reading
recs = SeqIO.parse(infile_FASTA, "fasta")

# Open FASTA for writing
outfile_hdl = open(outfile_name, "w")

seqs_from_location_count = 0

for rec in recs:
    this_ID = str(rec.id)

    if ID_to_location[this_ID] is not None and location_keyword in ID_to_location[this_ID]:
        seqs_from_location_count += 1
        this_seq = str(rec.seq).upper()

        # Write to FASTA file
        outfile_hdl.write(">" + this_ID + "\n")
        outfile_hdl.write(str(this_seq) + "\n")

# Close output file
outfile_hdl.close()

print("Sequences extracted from " + location_keyword + ": " + str(seqs_from_location_count) + "\n")


