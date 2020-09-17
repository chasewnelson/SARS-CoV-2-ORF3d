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
# CITATION: https://github.com/chasewnelson/

from Bio import SeqIO
from datetime import date, datetime, time
import sys
import os

Usage = "Track the frequency of alleles at specific sites in a sliding time window"

# Check argument number
if not len(sys.argv) == 7:
    raise SystemExit("\n### TERMINATED: need 6 unnamed arguments:\n" + \
                     "    # (1) .TSV with date metadata\n" + \
                     "    # (2) .FASTA alignment\n" + \
                     "    # (3) .TSV with a single column of sites to track\n" + \
                     "    # (4) name of desired results file (.TSV)\n" + \
                     "    # (5) window size\n" + \
                     "    # (6) step size\n\n")

infile_metadata = sys.argv[1]
infile_FASTA = sys.argv[2]
infile_sites = sys.argv[3]
outfile_name = sys.argv[4]
window_size = int(sys.argv[5])
step_size = int(sys.argv[6])

print("Date metadata from " + str(infile_metadata))
print("FASTA alignment from " + str(infile_FASTA))
print("Sites to track from " + str(infile_sites))
print("RESULTS will be placed in " + str(outfile_name))
print("Will perform a time-based sliding window: window_size=" + str(window_size) + ", step_size=" + str(step_size) + " (days)\n")

# Check infiles exist
if not os.path.isfile(infile_metadata):
    raise SystemExit("\n### TERMINATED: infile 1 (IDs) does not exist\n")

if not os.path.isfile(infile_FASTA):
    raise SystemExit("\n### TERMINATED: infile 2 (FASTA) does not exist\n")

if not os.path.isfile(infile_sites):
    raise SystemExit("\n### TERMINATED: infile 3 (one column of sites) does not exist\n")

# Create outfile_prefix name
print("Will write output to: " + outfile_name + "\n")

# Build a dictionary of dates/times
ID_to_date = {}
min_date = datetime(2019, 12, 20)

# for now, ASSUME col 1 is ID, and col
line_num = 0
with open(infile_metadata, "r") as f:
    lines = (line.rstrip() for line in f)

    for line in lines:
        line_num += 1
        line_list = line.split("\t")
        this_date = datetime(2019, 12, 20)

        if line_num == 1:
            # print("Accession ID")
            # print(line_list[0])

            if not line_list[0] == "Accession ID" and line_list[3] == "Collection date":
                raise SystemExit(
                    "\n### TERMINATED: header must be col1=\"Accession ID\" and col4=\"Collection date\"\n")
        else:
            date_worked = False
            seq_ID = line_list[0]
            seq_date = line_list[3]
            try:
                this_date = datetime.strptime(seq_date, "%Y-%m-%d")  # shortcut %F didn't work
                date_worked = True
            except:
                try:
                    # ID_to_date[seq_ID] =
                    this_date = datetime.strptime(seq_date, "%m/%d/%y")  # shortcut %D didn't work for
                    date_worked = True
                except:
                    print("Improper date: seq_ID=" + str(seq_ID) + ",seq_date=" + str(seq_date))

            if this_date >= min_date and date_worked:
                ID_to_date[seq_ID] = this_date

print("\n")

# find day 0
day0_date = min(ID_to_date.values())
day0_ID = list(ID_to_date.keys())[list(ID_to_date.values()).index(day0_date)]

# Get number of times from time 0
ID_to_days_elapsed = {}

for this_seq_ID, this_date in zip(ID_to_date.keys(), ID_to_date.values()):
    this_delta = this_date - day0_date
    ID_to_days_elapsed[this_seq_ID] = this_delta.days

### Gather sites to track
tracked_sites_list = []
with open(infile_sites) as infile_sites_hdl:
    lines = [x.rstrip() for x in infile_sites_hdl]

    for line in lines:
        tracked_sites_list.append(int(line))

print("Sites to track=" + ",".join(map(str, tracked_sites_list)))

print("Sequence counts per window:\n")

time_site_nts_dict = {}

# Gather frequency data
this_window_start = 0
while (this_window_start + window_size) <= max(ID_to_days_elapsed.values()):
    this_window_end = this_window_start + window_size
    this_window_midpoint = (this_window_start + this_window_end) / 2

    # Open FASTA for reading
    recs = SeqIO.parse(infile_FASTA, "fasta")

    # Initialize entry for this window if needed
    time_site_nts_dict[this_window_midpoint] = {}

    this_window_seq_count = 0
    for rec in recs:
        this_ID = str(rec.id)

        if this_ID in ID_to_days_elapsed.keys():
            this_ID_days_elapsed = ID_to_days_elapsed[this_ID]

            # If this sequence belongs to this window, add its allele(s)
            if this_ID_days_elapsed >= this_window_start and this_ID_days_elapsed <= this_window_end:
                this_window_seq_count += 1
                this_seq = str(rec.seq).upper()

                # Extract column for alignment subset including sequence of this timepoint
                for this_site in tracked_sites_list:
                    if this_site not in time_site_nts_dict[this_window_midpoint]:
                        time_site_nts_dict[this_window_midpoint][this_site] = str(this_seq[this_site - 1])
                    else:
                        time_site_nts_dict[this_window_midpoint][this_site] += str(this_seq[this_site - 1])

    print("this_window_start=" + str(this_window_start) + "," + "this_window_end=" + str(this_window_end) + "," + \
          "this_window_midpoint=" + str(this_window_midpoint) + "," + "this_window_seq_count=" + str(this_window_seq_count))

    # Update window
    this_window_start += step_size

# Open TSV for writing
outfile_hdl = open(outfile_name, "w")

# Write header
outfile_hdl.write("\t".join(map(str, ['this_midpoint', 'this_site', 'defined_allele_count', 'defined_nt_count', \
                                     'major_nt', 'major_nt_count', 'major_nt_freq', 'minor_nt', 'minor_nt_count', 'minor_nt_freq', \
                                     "A_count", "C_count", "G_count", "T_count", \
                                     "A_freq", "C_freq", "G_freq", "T_freq"])) + \
                          "\n")

for this_midpoint in time_site_nts_dict.keys():
    for this_site in time_site_nts_dict[this_midpoint].keys():
        this_col_letters = str(time_site_nts_dict[this_midpoint][this_site])

        this_col_nt_counts = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}

        for this_letter in list(this_col_letters):
            # print(this_letter)
            if this_col_nt_counts.get(this_letter) is None:
                this_col_nt_counts[this_letter] = 1
            else:
                this_col_nt_counts[this_letter] += 1

        minor_nt = '-'
        minor_nt_count = sum(this_col_nt_counts.values()) # max it out first
        major_nt = '-'
        major_nt_count = 0
        defined_nt_count = 0
        defined_allele_count = 0

        for key in this_col_nt_counts.keys():
            if key == "A" or key == "C" or key == "G" or key == "T":
                defined_nt_count += this_col_nt_counts[key]
                defined_allele_count += 1

                if this_col_nt_counts[key] < minor_nt_count and this_col_nt_counts[key] > 0:
                    minor_nt = key
                    minor_nt_count = this_col_nt_counts[key]

                if this_col_nt_counts[key] > major_nt_count:
                    major_nt = key
                    major_nt_count = this_col_nt_counts[key]

        major_nt_freq = major_nt_count / defined_nt_count
        minor_nt_freq = minor_nt_count / defined_nt_count

        this_col_nt_freqs = {"A": 0, "C": 0, "G": 0, "T": 0}
        for key in this_col_nt_counts.keys():
            if key == "A" or key == "C" or key == "G" or key == "T":
                this_col_nt_freqs[key] = this_col_nt_counts[key] / defined_nt_count

        # Write to TSV file
        outfile_hdl.write("\t".join(map(str, [this_midpoint, this_site, defined_allele_count, defined_nt_count, \
                                     major_nt, major_nt_count, major_nt_freq, minor_nt, minor_nt_count, minor_nt_freq, \
                                     this_col_nt_counts["A"], this_col_nt_counts["C"], this_col_nt_counts["G"], this_col_nt_counts["T"], \
                                     this_col_nt_freqs["A"], this_col_nt_freqs["C"], this_col_nt_freqs["G"], this_col_nt_freqs["T"]])) + \
                          "\n")

# Close output file
outfile_hdl.close()

raise SystemExit("\n### All done; we stopped here, dear.\n")


