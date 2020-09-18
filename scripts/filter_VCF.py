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
from scipy.stats import binom
import numpy as np
import os
import random
import re
import sys

Usage = "Script to dynamically filter within-host variants to control for a user-defined false-discovery rate"

#n = total reads
#p = error rate / 3
#x = number of reads of minor allele
#FDR = (1 - binomcdf(n, p, x -1) ) * length_of_the_genome

# Check argument number
if not len(sys.argv) == 6:
    raise SystemExit("\n### TERMINATED: need 5 unnamed arguments:\n" + \
                     #"    # (1) .VCF SNP report\n" + \
                     "    # (1) analysis-wide FDR cutoff\n" + \
                     "    # (2) minor allele frequency cutoff (min allowed)\n" + \
                     "    # (3) genome length\n" + \
                     "    # (4) sequencing error rate per site (assumes all nts equally probable)\n" + \
                     "    # (5) number of samples in analysis\n\n")

#infile_VCF = sys.argv[1]
FDR_cutoff = float(sys.argv[1])
MAF_cutoff = float(sys.argv[2])
genome_len = int(sys.argv[3])
error_per_site = float(sys.argv[4])
num_samples = int(sys.argv[5])


###############################################################################
### Gather VCF files in working directory
VCF_filenames = [name for name in os.listdir(".") if os.path.isfile(name) and name.endswith(".vcf")]

print("Analysis-wide FDR cutoff: " + str(FDR_cutoff))
print("Minor allele frequency cutoff (min allowed): " + str(MAF_cutoff))
print("Genome length: " + str(genome_len))
print("Sequencing error rate per site (assumes all nts equally probable): " + str(error_per_site))
print("Number of samples in analysis: " + str(num_samples))
print("VCF files to examine, gathered from working directory: " + str(VCF_filenames))
print("Will write output to files with names of the form: <VCF_name>_filtered.vcf\n")


###############################################################################
### Extract SNP data from VCF
# Prepare regex
VCF_DP_pattern = r"DP=\d+" # DP=266
VCF_DP_regex = re.compile(VCF_DP_pattern)
VCF_AF_pattern = r"AF=[\d\.]+" # AF=0.015038
VCF_AF_regex = re.compile(VCF_AF_pattern)
sample_count = 0
total_variant_count = 0
total_variant_pass_count = 0;
variants_pass_freqs = []
variants_pass_ACs = []
variants_pass_DPs = []
variants_fail_freqs = []
variants_fail_ACs = []
variants_fail_DPs = []

for this_VCF in VCF_filenames:

    # Check file exists
    if os.path.isfile(this_VCF):  # this_VCF

        # Create outfile name
        this_outfile_prefix = str(this_VCF)
        this_outfile_prefix = this_outfile_prefix.replace(".vcf", "")
        this_outfile_name = this_outfile_prefix + "_filtered.vcf"
        print("\nFiltered VCF output to: " + this_outfile_name + "\nEliminated records:")

        # Open output file for writing
        this_outfile_hdl = open(this_outfile_name, "w")
        new_metadata_lines = "##FILTER=<ID=MAF_cutoff,Description=\"Minor allele frequency cutoff (min allowed): " + \
                             str(round(MAF_cutoff, 5)) + "\">\n" + \
                             "##FILTER=<ID=FDR_cutoff,Description=\"Analysis-wide FDR cutoff (max): " + \
                             str(round(FDR_cutoff, 5)) + "\">\n" + \
                             "##FILTER=<ID=error_per_site,Description=\"Sequencing error rate per site: " + \
                             str(round(error_per_site, 5)) + "\">"  # no newline

        # Keep track of number pass
        variant_count = 0
        variant_pass_count = 0;

        with open(this_VCF, "r") as f:
            lines = (line.rstrip() for line in f)
            this_VCF_ID = str(this_VCF).rstrip(".vcf")
            sample_count += 1

            for line in lines:
                if line.startswith("##"):  # simply re-print existing metadata headers
                    # Write to FASTA file
                    this_outfile_hdl.write(line + "\n")
                elif line.startswith("#"): # add new metadata, then print header line
                    this_outfile_hdl.write(new_metadata_lines + "\n")
                    this_outfile_hdl.write(line + "\n")
                else:
                    variant_count += 1
                    total_variant_count += 1
                    line_list = line.split("\t")

                    # Extract information for this SNP
                    this_CHROM = line_list[0]
                    this_POS = int(line_list[1])
                    this_INFO = line_list[7]
                    this_DP = VCF_DP_regex.search(this_INFO)
                    this_DP = this_INFO[this_DP.start():this_DP.end()]
                    this_DP = int(this_DP.replace("DP=", ""))
                    this_AF = VCF_AF_regex.search(this_INFO)
                    this_AF = this_INFO[this_AF.start():this_AF.end()]  # AF=0.367003
                    this_AF = float(this_AF.replace("AF=", ""))

                    # Make sure it's the MINOR allele frequency
                    this_AF_minor = float(this_AF)
                    if this_AF_minor > 0.5:
                        this_AF_minor = 1 - this_AF_minor

                    # Allele count to be the nearest integer
                    this_AC_minor = round(this_DP * this_AF_minor)

                    # Calculate analysis-wide FDR for this allele count and coverage
                    FDR_per_genome_dataset = (1 - binom.cdf(this_AC_minor - 1, this_DP, error_per_site / 3)) * genome_len * num_samples

                    if FDR_per_genome_dataset > FDR_cutoff or this_AF_minor < MAF_cutoff: # exclude
                        variants_fail_freqs.append(this_AF_minor)
                        variants_fail_ACs.append(this_AC_minor)
                        variants_fail_DPs.append(this_DP)
                        print("Excluded the following record: " + str(this_VCF) + ", " + str(this_CHROM) + ":" + \
                            str(this_POS) + ",DP=" + str(this_DP) + ",AC=" + str(this_AC_minor) + \
                              ",FDR=" + str(round(FDR_per_genome_dataset, 5)) + ",MAF=" + \
                            str(round(this_AF_minor, 5)))
                    else: # include (print)
                        variants_pass_freqs.append(this_AF_minor)
                        variants_pass_ACs.append(this_AC_minor)
                        variants_pass_DPs.append(this_DP)
                        variant_pass_count += 1
                        total_variant_pass_count += 1
                        this_outfile_hdl.write(line + "\n")

        this_outfile_hdl.close()
        print("Variants pass: " + str(variant_pass_count) + "/" + str(variant_count))
    else:
        raise SystemExit("\n### TERMINATED: file doesn't exist: " + str(this_VCF) + "\n")


###############################################################################
### FINISH
mean_freq_pass = np.mean(np.array(variants_pass_freqs))
median_freq_pass = np.median(np.array(variants_pass_freqs))
min_freq_pass = np.min(np.array(variants_pass_freqs))
max_freq_pass = np.max(np.array(variants_pass_freqs))

mean_AC_pass = np.mean(np.array(variants_pass_ACs))
median_AC_pass = np.median(np.array(variants_pass_ACs))
min_AC_pass = np.min(np.array(variants_pass_ACs))
max_AC_pass = np.max(np.array(variants_pass_ACs))

mean_DP_pass = np.mean(np.array(variants_pass_DPs))
median_DP_pass = np.median(np.array(variants_pass_DPs))
min_DP_pass = np.min(np.array(variants_pass_DPs))
max_DP_pass = np.max(np.array(variants_pass_DPs))

mean_freq_fail = np.mean(np.array(variants_fail_freqs))
median_freq_fail = np.median(np.array(variants_fail_freqs))
min_freq_fail = np.min(np.array(variants_fail_freqs))
max_freq_fail = np.max(np.array(variants_fail_freqs))

mean_AC_fail = np.mean(np.array(variants_fail_ACs))
median_AC_fail = np.median(np.array(variants_fail_ACs))
min_AC_fail = np.min(np.array(variants_fail_ACs))
max_AC_fail = np.max(np.array(variants_fail_ACs))

mean_DP_fail = np.mean(np.array(variants_fail_DPs))
median_DP_fail = np.median(np.array(variants_fail_DPs))
min_DP_fail = np.min(np.array(variants_fail_DPs))
max_DP_fail = np.max(np.array(variants_fail_DPs))


print("\n")
print("###############################################################################")
print("Total samples examined: " + str(sample_count))
print("Total variants examined: " + str(total_variant_count))
print("Total variants pass: " + str(total_variant_pass_count))
print("Fraction variants pass: " + str(total_variant_pass_count / total_variant_count))

print("min / mean / median / max FREQUENCY of passing alleles: " + str(min_freq_pass) + " / " + \
      str(mean_freq_pass) + " / " + str(median_freq_pass) + " / " + str(max_freq_pass))

print("min / mean / median / max AC of passing alleles: " + str(min_AC_pass) + " / " + \
      str(mean_AC_pass) + " / " + str(median_AC_pass) + " / " + str(max_AC_pass))

print("min / mean / median / max DP of passing alleles: " + str(min_DP_pass) + " / " + \
      str(mean_DP_pass) + " / " + str(median_DP_pass) + " / " + str(max_DP_pass))

print("min / mean / median / max FREQUENCY of failing alleles: " + str(min_freq_fail) + " / " + \
      str(mean_freq_fail) + " / " + str(median_freq_fail) + " / " + str(max_freq_fail))

print("min / mean / median / max AC of failing alleles: " + str(min_AC_fail) + " / " + \
      str(mean_AC_fail) + " / " + str(median_AC_fail) + " / " + str(max_AC_fail))

print("min / mean / median / max DP of failing alleles: " + str(min_DP_fail) + " / " + \
      str(mean_DP_fail) + " / " + str(median_DP_fail) + " / " + str(max_DP_fail))

print("###############################################################################\n")
raise SystemExit("\n### All done; we stopped here, dear.\n")


