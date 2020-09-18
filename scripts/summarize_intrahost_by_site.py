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
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import os
import re
import sys

Usage = "Script to create a genome database cataloguing within-host variants"

# Check argument number
if not len(sys.argv) == 3:
    raise SystemExit("\n### TERMINATED: need 2 unnamed arguments:\n" + \
                     "    # (1) FASTA reference (1 sequence)\n" + \
                     "    # (2) GTF file with comprehensive gene annotations\n\n")

infile_FASTA = sys.argv[1]
infile_GTF = sys.argv[2]

# Define outfile name
outfile_prefix = str(infile_FASTA)
outfile_prefix = outfile_prefix.replace(".vcf", "")
outfile_prefix = outfile_prefix.replace(".fasta", "")
outfile_prefix = outfile_prefix.replace(".txt", "")
outfile_prefix = outfile_prefix.replace(".tsv", "")
outfile_name = outfile_prefix + "_site_database" + ".tsv"
print("Will write output to: " + outfile_name + "\n")


###############################################################################
### Gather VCF files in working directory
VCF_filenames = [name for name in os.listdir(".") if os.path.isfile(name) and name.endswith(".vcf")]
print("VCF files to examine, gathered from working directory: " + str(VCF_filenames))
print("Will write output to files with names of the form: <VCF_name>_filtered.vcf\n")


###############################################################################
### Record reference sequence
# Open FASTA for reading
rec = SeqIO.read(infile_FASTA, "fasta") # will throw error if more than one record!
ref_seq = str(rec.seq)
ref_seq_len = len(ref_seq)


###############################################################################
### Build dictionary of data by site
site_data = {}

for this_nt_idx in range(ref_seq_len):
    this_nt = ref_seq[this_nt_idx]
    this_site = this_nt_idx + 1
    site_data[this_site] = {}

    # Add reference and an entry for each nt
    site_data[this_site]['REF'] = this_nt
    site_data[this_site]['A'] = "."
    site_data[this_site]['C'] = "."
    site_data[this_site]['G'] = "."
    site_data[this_site]['T'] = "."
    site_data[this_site]['gene1'] = "."
    site_data[this_site]['gene2'] = "."
    site_data[this_site]['gene1_CDS'] = "."
    site_data[this_site]['gene2_CDS'] = "."
    site_data[this_site]['gene1_codon'] = "."
    site_data[this_site]['gene2_codon'] = "."
    site_data[this_site]['gene1_codon_position'] = "."
    site_data[this_site]['gene2_codon_position'] = "."
    site_data[this_site]['gene1_aa_REF'] = "."
    site_data[this_site]['gene2_aa_REF'] = "."
    site_data[this_site]['gene1_codon_ALT_A'] = "."
    site_data[this_site]['gene2_codon_ALT_A'] = "."
    site_data[this_site]['gene1_codon_ALT_C'] = "."
    site_data[this_site]['gene2_codon_ALT_C'] = "."
    site_data[this_site]['gene1_codon_ALT_G'] = "."
    site_data[this_site]['gene2_codon_ALT_G'] = "."
    site_data[this_site]['gene1_codon_ALT_T'] = "."
    site_data[this_site]['gene2_codon_ALT_T'] = "."
    site_data[this_site]['gene1_aa_ALT_A'] = "."
    site_data[this_site]['gene2_aa_ALT_A'] = "."
    site_data[this_site]['gene1_aa_ALT_C'] = "."
    site_data[this_site]['gene2_aa_ALT_C'] = "."
    site_data[this_site]['gene1_aa_ALT_G'] = "."
    site_data[this_site]['gene2_aa_ALT_G'] = "."
    site_data[this_site]['gene1_aa_ALT_T'] = "."
    site_data[this_site]['gene2_aa_ALT_T'] = "."

# site_data[14409]['REF']
# site_data[14409]['A']


###############################################################################
### Loop GTF file to add genes, including up to one OLG
VCF_gene_id_pattern = r'gene_id "([\w\d\-\_\.\'\(\)\"]+)"'
VCF_gene_id_regex = re.compile(VCF_gene_id_pattern)

gene_coordinates = {}

with open(infile_GTF, "r") as infile_GTF_hdl:
    lines = (line.rstrip() for line in infile_GTF_hdl)

    for line in lines:

        line_list = line.split("\t")

        if len(line_list) == 9:
            this_start = int(line_list[3])
            this_end = int(line_list[4])
            this_name_field = str(line_list[8])

            # Extract name of gene
            this_name = VCF_gene_id_regex.search(this_name_field)
            this_name = this_name_field[this_name.start():this_name.end()]
            this_name = this_name.replace("gene_id \"", "")
            this_name = this_name[:-1]
            # print(this_name)

            if this_name not in gene_coordinates.keys():
                gene_coordinates[this_name] = {}
                gene_coordinates[this_name][this_start] = this_end
            elif this_start not in gene_coordinates[this_name].keys():
                gene_coordinates[this_name][this_start] = this_end
            else:
                raise SystemExit("\n### TERMINATED: multiple segments of " + str(this_name) + \
                                 " begin at the same site.\n")

# print(gene_coordinates)
# {'ORF1ab': {266: 13468, 13468: 21555}, 'S': {21563: 25384}, 'ORF3a': {25393: 26220}, 'ORF3bp': {25524: 25697}, 'E': {26245: 26472}, 'M': {26523: 27191}, 'ORF6': {27202: 27387}, 'ORF7a': {27394: 27759}, 'ORF7b': {27756: 27887}, 'ORF8p': {27894: 28259}, 'N': {28274: 29533}, 'ORF9b': {28284: 28577}, 'ORF9c': {28734: 28955}, 'ORF10': {29558: 29674}}
infile_GTF_hdl.close()

###############################################################################
### Loop each gene to ensure multiple of 3, and build CDS sequence
gene_CDS = {}
for this_gene in list(gene_coordinates.keys()):
    # print(this_gene)
    total_length = 0

    this_gene_starts = list(gene_coordinates[this_gene].keys())
    this_gene_starts.sort(key=int)

    for this_start in this_gene_starts:
        # print(this_start)
        this_end = gene_coordinates[this_gene][this_start]
        segment_length = (this_end - this_start + 1)
        total_length = total_length + segment_length
        segment_seq = ref_seq[(this_start - 1):this_end]

        if this_gene not in gene_CDS.keys():
            gene_CDS[this_gene] = str(segment_seq)
        else:
            gene_CDS[this_gene] = str(gene_CDS[this_gene]) + str(segment_seq)
    # print("")

    if total_length % 3 == 0:
        print("The length of " + this_gene + " is " + str(total_length) + ", a multiple of 3.")
    else:
        raise SystemExit("\n### TERMINATED: The length of " + this_gene + " is a NOT multiple of 3.\n")
print("\n")


###############################################################################
### Loop each gene and add to site data
for this_gene in gene_coordinates.keys():
    # print(this_gene)
    this_gene_starts = list(gene_coordinates[this_gene].keys())
    this_gene_starts.sort(key=int)

    # set up variables
    curr_codon_position = 1

    # global curr_codon
    curr_codon = "NNN"
    curr_aa_REF = "X"
    curr_CDS_position = 1
    curr_CDS = gene_CDS[this_gene]

    # Loop!
    for this_start in this_gene_starts:
        this_end = gene_coordinates[this_gene][this_start]

        for this_site in range(this_start, this_end + 1):
            curr_gene_key = "gene1"
            if site_data[this_site]['gene1'] != ".":
                curr_gene_key = "gene2"

            # Store gene at site
            site_data[this_site][curr_gene_key] = str(this_gene)

            # prepare ALT codons and amino acids
            curr_codon_ALT_A = "."
            curr_codon_ALT_C = "."
            curr_codon_ALT_G = "."
            curr_codon_ALT_T = "."
            curr_aa_ALT_A = "."
            curr_aa_ALT_C = "."
            curr_aa_ALT_G = "."
            curr_aa_ALT_T = "."

            if curr_codon_position == 1:

                curr_codon = str(curr_CDS[(curr_CDS_position - 1):(curr_CDS_position + 2)])
                curr_codon_Seq = Seq(curr_codon, generic_dna)
                curr_aa_REF = curr_codon_Seq.translate()

                # generate ALT codons
                curr_codon_ALT_A = str('A' + curr_codon[1:3])
                curr_codon_ALT_A_Seq = Seq(curr_codon_ALT_A, generic_dna)
                curr_aa_ALT_A = curr_codon_ALT_A_Seq.translate()
                curr_codon_ALT_C = str('C' + curr_codon[1:3])
                curr_codon_ALT_C_Seq = Seq(curr_codon_ALT_C, generic_dna)
                curr_aa_ALT_C = curr_codon_ALT_C_Seq.translate()
                curr_codon_ALT_G = str('G' + curr_codon[1:3])
                curr_codon_ALT_G_Seq = Seq(curr_codon_ALT_G, generic_dna)
                curr_aa_ALT_G = curr_codon_ALT_G_Seq.translate()
                curr_codon_ALT_T = str('T' + curr_codon[1:3])
                curr_codon_ALT_T_Seq = Seq(curr_codon_ALT_T, generic_dna)
                curr_aa_ALT_T = curr_codon_ALT_T_Seq.translate()
            elif curr_codon_position == 2:
                curr_codon_ALT_A = str(curr_codon[0:1] + 'A' + curr_codon[2:3])
                curr_codon_ALT_A_Seq = Seq(curr_codon_ALT_A, generic_dna)
                curr_aa_ALT_A = curr_codon_ALT_A_Seq.translate()
                curr_codon_ALT_C = str(curr_codon[0:1] + 'C' + curr_codon[2:3])
                curr_codon_ALT_C_Seq = Seq(curr_codon_ALT_C, generic_dna)
                curr_aa_ALT_C = curr_codon_ALT_C_Seq.translate()
                curr_codon_ALT_G = str(curr_codon[0:1] + 'G' + curr_codon[2:3])
                curr_codon_ALT_G_Seq = Seq(curr_codon_ALT_G, generic_dna)
                curr_aa_ALT_G = curr_codon_ALT_G_Seq.translate()
                curr_codon_ALT_T = str(curr_codon[0:1] + 'T' + curr_codon[2:3])
                curr_codon_ALT_T_Seq = Seq(curr_codon_ALT_T, generic_dna)
                curr_aa_ALT_T = curr_codon_ALT_T_Seq.translate()
            elif curr_codon_position == 3:
                curr_codon_ALT_A = str(curr_codon[0:2] + 'A')
                curr_codon_ALT_A_Seq = Seq(curr_codon_ALT_A, generic_dna)
                curr_aa_ALT_A = curr_codon_ALT_A_Seq.translate()
                curr_codon_ALT_C = str(curr_codon[0:2] + 'C')
                curr_codon_ALT_C_Seq = Seq(curr_codon_ALT_C, generic_dna)
                curr_aa_ALT_C = curr_codon_ALT_C_Seq.translate()
                curr_codon_ALT_G = str(curr_codon[0:2] + 'G')
                curr_codon_ALT_G_Seq = Seq(curr_codon_ALT_G, generic_dna)
                curr_aa_ALT_G = curr_codon_ALT_G_Seq.translate()
                curr_codon_ALT_T = str(curr_codon[0:2] + 'T')
                curr_codon_ALT_T_Seq = Seq(curr_codon_ALT_T, generic_dna)
                curr_aa_ALT_T = curr_codon_ALT_T_Seq.translate()

            # Store data
            site_data[this_site][curr_gene_key + '_CDS'] = int(curr_CDS_position)
            site_data[this_site][curr_gene_key + '_codon'] = str(curr_codon)
            site_data[this_site][curr_gene_key + '_codon_position'] = int(curr_codon_position)
            site_data[this_site][curr_gene_key + '_aa_REF'] = str(curr_aa_REF)

            site_data[this_site][str(curr_gene_key + '_codon_ALT_A')] = str(curr_codon_ALT_A)
            site_data[this_site][curr_gene_key + '_codon_ALT_C'] = str(curr_codon_ALT_C)
            site_data[this_site][curr_gene_key + '_codon_ALT_G'] = str(curr_codon_ALT_G)
            site_data[this_site][curr_gene_key + '_codon_ALT_T'] = str(curr_codon_ALT_T)
            site_data[this_site][curr_gene_key + '_aa_ALT_A'] = str(curr_aa_ALT_A)
            site_data[this_site][curr_gene_key + '_aa_ALT_C'] = str(curr_aa_ALT_C)
            site_data[this_site][curr_gene_key + '_aa_ALT_G'] = str(curr_aa_ALT_G)
            site_data[this_site][curr_gene_key + '_aa_ALT_T'] = str(curr_aa_ALT_T)

            # Update CDS and codon positions
            curr_CDS_position = curr_CDS_position + 1
            if curr_codon_position == 3:
                curr_codon_position = 1
            else:
                curr_codon_position = curr_codon_position + 1

#print(site_data[14409])


###############################################################################
### Independently, gather all intrahost data fro VCF files in working directory
VCF_data = {}
VCF_ID_list = []

# Prepare regex
VCF_DP_pattern = r"DP=\d+"  # DP=266
VCF_DP_regex = re.compile(VCF_DP_pattern)
VCF_AF_pattern = r"AF=[\d\.]+"  # AF=0.015038
VCF_AF_regex = re.compile(VCF_AF_pattern)
sample_count = 0
total_variant_count = 0
total_variant_pass_count = 0
variants_pass_freqs = []
variants_pass_ACs = []
variants_pass_DPs = []
variants_fail_freqs = []
variants_fail_ACs = []
variants_fail_DPs = []

for this_VCF in VCF_filenames:

    # Check file exists
    if os.path.isfile(this_VCF):  # this_VCF
        # print(this_VCF + " exists!")
        with open(this_VCF, "r") as this_VCF_hdl:
            lines = (line.rstrip() for line in this_VCF_hdl)
            this_VCF_ID = str(this_VCF).rstrip(".vcf")
            VCF_ID_list.append(this_VCF_ID)

            for line in lines:
                if not line.startswith("#"):  # skip metadata and headers
                    line_list = line.split("\t")

                    # Extract information for this SNP
                    this_POS = int(line_list[1])
                    this_REF = line_list[3]
                    this_ALT = str(line_list[4])
                    this_INFO = line_list[7]
                    this_DP = VCF_DP_regex.search(this_INFO)
                    this_DP = this_INFO[this_DP.start():this_DP.end()]  # AF=0.367003
                    this_DP = int(this_DP.replace("DP=", ""))
                    this_AF = VCF_AF_regex.search(this_INFO)
                    this_AF = this_INFO[this_AF.start():this_AF.end()]  # AF=0.367003
                    this_AF = float(this_AF.replace("AF=", ""))

                    this_ALT_count = int(round(this_AF * this_DP))
                    this_REF_count = int(this_DP - this_ALT_count)

                    if this_POS not in VCF_data.keys():
                        VCF_data[this_POS] = {}
                        VCF_data[this_POS][this_VCF_ID] = {}
                        VCF_data[this_POS][this_VCF_ID]['REF'] = this_REF
                        VCF_data[this_POS][this_VCF_ID]['ALT'] = this_ALT
                        VCF_data[this_POS][this_VCF_ID]['REF_count'] = this_REF_count
                        VCF_data[this_POS][this_VCF_ID]['ALT_count'] = this_ALT_count
                    elif this_VCF_ID not in VCF_data[this_POS].keys():
                        VCF_data[this_POS][this_VCF_ID] = {}
                        VCF_data[this_POS][this_VCF_ID]['REF'] = this_REF
                        VCF_data[this_POS][this_VCF_ID]['ALT'] = this_ALT
                        VCF_data[this_POS][this_VCF_ID]['REF_count'] = this_REF_count
                        VCF_data[this_POS][this_VCF_ID]['ALT_count'] = this_ALT_count
                    elif 'ALT2' not in VCF_data[this_POS][this_VCF_ID].keys():
                        VCF_data[this_POS][this_VCF_ID]['ALT2'] = this_ALT
                        VCF_data[this_POS][this_VCF_ID]['ALT2_count'] = this_ALT_count
                        VCF_data[this_POS][this_VCF_ID]['REF_count'] = VCF_data[this_POS][this_VCF_ID]['REF_count'] - this_ALT_count
                    elif 'ALT3' not in VCF_data[this_POS][this_VCF_ID].keys():
                        VCF_data[this_POS][this_VCF_ID]['ALT3'] = this_ALT
                        VCF_data[this_POS][this_VCF_ID]['ALT3_count'] = this_ALT_count
                        VCF_data[this_POS][this_VCF_ID]['REF_count'] = VCF_data[this_POS][this_VCF_ID]['REF_count'] - this_ALT_count
                    else:
                        raise SystemExit("\n### TERMINATED: sample " + str(this_VCF_ID) + \
                                         " contains two records for site " + str(this_POS) + \
                                         ". Further coding needed.\n")


    else:
        raise SystemExit("\n### TERMINATED: file doesn't exist: " + str(this_VCF) + "\n")


###############################################################################
### Print output file

# Sort sites
sites_list = list(site_data.keys())
sites_list.sort(key=int)

# Sort VCF IDs
VCF_ID_list.sort()

# Open TSV for writing
outfile_hdl = open(outfile_name, "w")

# Print header
header_list = ['this_site',
               'REF',
               'gene1',
               'gene2',
               'gene1_CDS',
               'gene2_CDS',
               'gene1_codon',
               'gene2_codon',
               'gene1_codon_position',
               'gene2_codon_position',
               'gene1_REF_aa',
               'gene2_REF_aa',
               'ALT',
               'gene1_ALT_codon',
               'gene2_ALT_codon',
               'gene1_ALT_aa',
               'gene2_ALT_aa'
               ]
header_list.extend(VCF_ID_list)
header_line = "\t".join(map(str, header_list))
outfile_hdl.write(header_line + "\n")

# Write data
for this_site in sites_list:
    this_line_list_start = [this_site,
                            site_data[this_site]['REF'],
                            site_data[this_site]['gene1'],
                            site_data[this_site]['gene2'],
                            site_data[this_site]['gene1_CDS'],
                            site_data[this_site]['gene2_CDS'],
                            site_data[this_site]['gene1_codon'],
                            site_data[this_site]['gene2_codon'],
                            site_data[this_site]['gene1_codon_position'],
                            site_data[this_site]['gene2_codon_position'],
                            site_data[this_site]['gene1_aa_REF'],
                            site_data[this_site]['gene2_aa_REF']]

    this_line_list_A = [site_data[this_site]['gene1_codon_ALT_A'],
                        site_data[this_site]['gene2_codon_ALT_A'],
                        site_data[this_site]['gene1_aa_ALT_A'],
                        site_data[this_site]['gene2_aa_ALT_A']]

    this_line_list_C = [site_data[this_site]['gene1_codon_ALT_C'],
                        site_data[this_site]['gene2_codon_ALT_C'],
                        site_data[this_site]['gene1_aa_ALT_C'],
                        site_data[this_site]['gene2_aa_ALT_C']]

    this_line_list_G = [site_data[this_site]['gene1_codon_ALT_G'],
                        site_data[this_site]['gene2_codon_ALT_G'],
                        site_data[this_site]['gene1_aa_ALT_G'],
                        site_data[this_site]['gene2_aa_ALT_G']]

    this_line_list_T = [site_data[this_site]['gene1_codon_ALT_T'],
                        site_data[this_site]['gene2_codon_ALT_T'],
                        site_data[this_site]['gene1_aa_ALT_T'],
                        site_data[this_site]['gene2_aa_ALT_T']]

    if this_site in VCF_data.keys():  # there was at least one intrahost variant at this site
        for this_VCF_ID in VCF_ID_list:
            if this_VCF_ID in VCF_data[this_site].keys():  # this sample had a variant here
                # same REF allele and a defined ALT nucleotide
                if VCF_data[this_site][this_VCF_ID]['REF'] == site_data[this_site]['REF']:
                    if VCF_data[this_site][this_VCF_ID]['ALT'] in ['A', 'C', 'G', 'T']:
                        if VCF_data[this_site][this_VCF_ID]['ALT'] == 'A':
                            this_line_list_A.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                    "," + str(VCF_data[this_site][this_VCF_ID]['ALT_count']))
                        elif VCF_data[this_site][this_VCF_ID]['ALT'] == 'C':
                            this_line_list_C.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                    "," + str(VCF_data[this_site][this_VCF_ID]['ALT_count']))
                        elif VCF_data[this_site][this_VCF_ID]['ALT'] == 'G':
                            this_line_list_G.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                    "," + str(VCF_data[this_site][this_VCF_ID]['ALT_count']))
                        elif VCF_data[this_site][this_VCF_ID]['ALT'] == 'T':
                            this_line_list_T.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                    "," + str(VCF_data[this_site][this_VCF_ID]['ALT_count']))
                        else:
                            print("Something went wrong at site=" + this_site + ", VCF=" + this_VCF_ID + "\n")
                            this_line_list_A.append(".")
                            this_line_list_C.append(".")
                            this_line_list_G.append(".")
                            this_line_list_T.append(".")

                        # Multiallelic sites
                        if 'ALT2' in VCF_data[this_site][this_VCF_ID].keys() and VCF_data[this_site][this_VCF_ID]['ALT2'] in ['A', 'C', 'G', 'T']:
                            if VCF_data[this_site][this_VCF_ID]['ALT2'] == 'A':
                                this_line_list_A.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT2_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT2'] == 'C':
                                this_line_list_C.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT2_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT2'] == 'G':
                                this_line_list_G.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT2_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT2'] == 'T':
                                this_line_list_T.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT2_count']))
                            else:
                                print("Something went wrong at site=" + this_site + ", VCF=" + this_VCF_ID + "\n")
                                this_line_list_A.append(".")
                                this_line_list_C.append(".")
                                this_line_list_G.append(".")
                                this_line_list_T.append(".")
                                
                        if 'ALT3' in VCF_data[this_site][this_VCF_ID].keys() and VCF_data[this_site][this_VCF_ID]['ALT3'] in ['A', 'C', 'G', 'T']:
                            if VCF_data[this_site][this_VCF_ID]['ALT3'] == 'A':
                                this_line_list_A.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT3_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT3'] == 'C':
                                this_line_list_C.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT3_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT3'] == 'G':
                                this_line_list_G.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT3_count']))
                            elif VCF_data[this_site][this_VCF_ID]['ALT3'] == 'T':
                                this_line_list_T.append(str(VCF_data[this_site][this_VCF_ID]['REF_count']) + \
                                                        "," + str(VCF_data[this_site][this_VCF_ID]['ALT3_count']))
                            else:
                                print("Something went wrong at site=" + this_site + ", VCF=" + this_VCF_ID + "\n")
                                this_line_list_A.append(".")
                                this_line_list_C.append(".")
                                this_line_list_G.append(".")
                                this_line_list_T.append(".")
                else:  # something went wrong
                    print("Something went wrong at site=" + this_site + ", VCF=" + this_VCF_ID + "\n")
                    this_line_list_A.append(".")
                    this_line_list_C.append(".")
                    this_line_list_G.append(".")
                    this_line_list_T.append(".")
            else:  # no variant here
                this_line_list_A.append(".")
                this_line_list_C.append(".")
                this_line_list_G.append(".")
                this_line_list_T.append(".")
    else:  # no variants, all blank
        for this_VCF_ID in VCF_ID_list:  # TODO: find a repeat function for this
            this_line_list_A.append(".")
            this_line_list_C.append(".")
            this_line_list_G.append(".")
            this_line_list_T.append(".")

    # print(this_line_list)
    this_line_start = "\t".join(map(str, this_line_list_start))
    this_line_A = this_line_start + "\tA\t" + "\t".join(map(str, this_line_list_A))
    this_line_C = this_line_start + "\tC\t" + "\t".join(map(str, this_line_list_C))
    this_line_G = this_line_start + "\tG\t" + "\t".join(map(str, this_line_list_G))
    this_line_T = this_line_start + "\tT\t" + "\t".join(map(str, this_line_list_T))

    # Write to output file
    outfile_hdl.write(this_line_A + "\n")
    outfile_hdl.write(this_line_C + "\n")
    outfile_hdl.write(this_line_G + "\n")
    outfile_hdl.write(this_line_T + "\n")

# Close TSV output file
outfile_hdl.write("\n")
outfile_hdl.close()


