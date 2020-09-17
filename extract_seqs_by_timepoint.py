#! /usr/bin/env python

from Bio import SeqIO
from datetime import date, datetime, time
import os
import sys

#Usage = "Script for extracting a subset of sequences (arg 1 file, single column) from a FASTA (arg 2 file)"

# Check two argument
if not len(sys.argv) == 3:
    raise SystemExit("\n### TERMINATED: need 2 unnamed arguments: (1) .TSV with date metadata; and (2) .FASTA alignment\n")
infile_metadata = sys.argv[1]
infile_FASTA = sys.argv[2]
window_size = 14
step_size = 7

print("Will perform a time-based sliding window: window_size=" + str(window_size) + ", step_size=" + str(step_size) + " (days)\n")

# Check infiles exist
if not os.path.isfile(infile_metadata):
    raise SystemExit("\n### TERMINATED: infile 1 (IDs) does not exist\n")

if not os.path.isfile(infile_FASTA):
    raise SystemExit("\n### TERMINATED: infile 2 (FASTA) does not exist\n")

# Create outfile_prefix name
outfile_prefix = str(infile_FASTA)
outfile_prefix = outfile_prefix.replace(".fasta", "")
outfile_prefix = outfile_prefix.replace(".txt", "")
outfile_prefix = outfile_prefix.replace(".tsv", "")
print("Will write output to: " + outfile_prefix + "_[DATE].fasta\n")

# Build a dictionary of dates/times
ID_to_date = {}
min_date = datetime(2019, 12, 20)

# for now, ASSUME col 1 (idx0) is ID, and col 4 (idx 3) is date
line_num = 0
with open(infile_metadata, "r") as f:
    lines = (line.rstrip() for line in f)

    for line in lines:
        line_num += 1
        line_list = line.split("\t")
        this_date = datetime(2019, 12, 20)

        if line_num == 1:

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
                    # pass
                    print("Improper date: seq_ID=" + str(seq_ID) + ",seq_date=" + str(seq_date)) # "\n"

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

print("Sequence counts per window:\n")

this_window_start = 0
while (this_window_start + window_size) <= max(ID_to_days_elapsed.values()):
    this_window_end = this_window_start + window_size

    # Open FASTA for reading
    recs = SeqIO.parse(infile_FASTA, "fasta")

    # Create directory for FASTA file
    new_dir = outfile_prefix + "_" + str(this_window_start) + "to" + str(this_window_end)
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    # Open FASTA for writing
    this_outfile_name = new_dir + "/" + outfile_prefix + "_" + str(this_window_start) + "to" + str(
        this_window_end) + ".fasta"
    this_outfile_hdl = open(this_outfile_name, "w")

    this_window_seq_count = 0
    for rec in recs:
        this_ID = str(rec.id)

        if this_ID in ID_to_days_elapsed.keys():
            this_ID_days_elapsed = ID_to_days_elapsed[this_ID]

            if this_ID_days_elapsed >= this_window_start and this_ID_days_elapsed <= this_window_end:
                this_window_seq_count += 1
                this_seq = str(rec.seq).upper()

                # Write to FASTA file
                this_outfile_hdl.write(">" + this_ID + "\n")
                this_outfile_hdl.write(str(this_seq) + "\n")

    # Close output file
    this_outfile_hdl.close()

    print("this_window_start=" + str(this_window_start) + "," + "this_window_end=" + str(this_window_end) + "," + \
          "this_window_seq_count=" + str(this_window_seq_count))

    # Update window
    this_window_start += step_size


