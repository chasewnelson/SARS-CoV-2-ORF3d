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
import sys

Usage = "Script for tallying epitope coverage for one product by 9- or 15-mers"

# Check argument number
if not len(sys.argv) == 3:
    raise SystemExit("\n### TERMINATED: need 2 unnamed arguments:\n" + \
                     "    # (1) .TSV file with 5 columns: ID/NB/product/codon_start/codon_end\n" + \
                     "    # (2) LENGTH of linear peptides: 9 (MHC I) or 15 (MHC II)\n\n")

# Input file
random_ORF_file = sys.argv[1]
epitope_len = int(sys.argv[2])

# Gather data
epitope_data = {}

# Check file exists
if os.path.isfile(random_ORF_file):
    with open(random_ORF_file, "r") as random_ORF_hdl:
        lines = (line.rstrip() for line in random_ORF_hdl)
        this_file_prefix = str(random_ORF_file).rstrip(".tsv")

        for line in lines:
            if not line.startswith("ID"):  # skip header
                line_list = line.split("\t")

                # Extract information for this site
                this_ID = line_list[0]
                this_NB = int(line_list[1])
                this_start = int(line_list[3])

                if this_ID not in epitope_data.keys():
                    epitope_data[this_ID] = {}
                    epitope_data[this_ID][this_start] = this_NB
                else:
                    epitope_data[this_ID][this_start] = this_NB

# Loop and compute
for this_ID in epitope_data.keys():
    starts_sorted = list(epitope_data[this_ID].keys())
    starts_sorted.sort(key=int)

    # Initialize
    site_NB_coverage = {}

    for i in range(min(starts_sorted), max(starts_sorted) + epitope_len):
        site_NB_coverage[i] = 0

    for this_start in starts_sorted:
        this_end = this_start + epitope_len - 1
        this_NB = epitope_data[this_ID][this_start]

        for i in range(this_start, this_end + 1):
            site_NB_coverage[i] += this_NB

    for this_start in starts_sorted:
        this_end = this_start + epitope_len - 1
        this_NB = epitope_data[this_ID][this_start]
        this_NB_coverage = site_NB_coverage[this_start]

        # Print to STDOUT
        print(this_ID + "\t" + str(this_start) + "\t" + str(this_end) + "\t" + str(this_NB) + "\t" + str(this_NB_coverage))


