<img src="https://github.com/chasewnelson/SARS-CoV-2-ORF3d/blob/master/images/cover_image.png?raw=true" title="Overlapping gene products" alt="Overlapping gene products" align="left" size="small">

# SARS-CoV-2 *ORF3d*

Supplementary scripts for Nelson et al. (2020) paper on SARS-CoV-2 *ORF3d*:


Nelson CW, Ardern Z, Goldberg TL, Meng C, Kuo C-H, Ludwig C, Kolokotronis S-O, Wei X. 2020. 
<a target="_blank" href="https://elifesciences.org/articles/59633">Dynamically evolving novel overlapping gene as a factor in the SARS-CoV-2 pandemic.</a> *eLife* **9**: e59633. DOI: 10.7554/eLife.59633


## <a name="contents"></a>Contents

* [Supplementary data](#supplementary-data)
* [Supplementary scripts](#supplementary-scripts)
	* [**Figure 1**. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members](#figure-1).
		* `fig1B.bash`
		* `ORF_length.R`
	* [**Figure 2**. Re-analysis of SARS-CoV-2 gene expression in publicly available ribosome profiling and mass spectrometry datasets](#figure-2).
		* `aligned_fasta2haplotypes.pl`
		* `riboseq_sliding_window.R`
		* `riboseq.R`
		* `mass-spec_vs_riboseq.R`
	* [**Figure 3**. SARS-CoV-2 protein sequence properties](#figure-3).
		* `generate_random_protein.py`
		* `tally_epitope_coverage.py`
		* `epitope_MHCI.R`
		* `epitope_MHCII.R`
		* `hydrophobicity_profiles_ORF3a.R`
	* [**Figure 5**. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels](#figure-5).
		* `SARS-CoV-2_locate_genes.pl`
		* `generate_seqs_from_VCF.py`
		* `three_levels_diversity.R`
		* `selection_vs_expression.R`
		* `extract_seqs_by_timepoint.py`
		* `extract_seqs_by_location.py`
		* `extract_positions_by_timepoint.py`
		* `temporal_pi.R`
	* [**Figure 6**. Between-taxa sliding window analysis of natural selection on overlapping frames of *ORF3a*](#figure-6).
		* `selection_sliding_windows.R`
	* [**Figure 7**. Pandemic spread of the EP+1 haplotype and the hitchhiking of *ORF3d*-LOF](#figure-7).
		* `extract_variable_columns_MSA.py`
		* `extract_variable_columns_MSA_aa.py`
		* `Fig7.m`
		* `Fig7b.m`
	* [**Figure 8**. High-frequency within-host mutations](#figure-8).
		* `filter_vcf.py`
		* `summarize_intrahost_by_site.py`
		* `Fig8.m`
		* `Fig8_supp1.m`
	* [**Additional scripts**](#additional-scripts).
		* `extract_fasta_by_sites.pl`
		* `extract_seq_subset.py`
		* `translate_nt_seq_file.pl`
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)


## <a name="supplementary-data"></a>Supplementary data

All supplementary are available at <a target="_blank" href="https://zenodo.org/record/4052729">Zenodo</a> under record ID 4052729.

For easy access, the four most important supplementary data files are available in the `/data/` directory of this repository:

1. `SARS-related-CoV_ALN.fasta`: whole-genome multiple sequence alignment of *n*=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus* (between-taxa analysis). See [manuscript](#citation) for details. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.

2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within **SARS-related-CoV_ALN.fasta**.

3. `SARS-CoV-2.gtf`: Gene Transfer Format (GTF) file giving gene positions within the SARS-CoV-2 <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1 genome</a> genome. Records are ordered by start site in the genome. Full gene regions are given, as well as just those segments of the gene overlapping (OL) or not overlapping (NOL) other genes. For example, the positions of *ORF3d* are 25524-25697; the positions of ORF3a codons overlapping *ORF3d* (ORF3a_OL_ORF3d) are 25522-25698; and the positions of *ORF3a* codons involved in the *ORF3a*/*ORF3c*/*ORF3d* triple overlap (ORF3a_OL_ORF3c_OL_ORF3d) are 25522-25584. *ORF3b* refers to the first short (23 codon) ORF in SARS-CoV-2 and is considered OL; the remainder of the region homologous to *ORF3b* in SARS-CoV is not considered OL.

4. `Supplementary_Tables.xlsx`: Supplementary Tables referred to in the [manuscript](#citation).


## <a name="supplementary-scripts"></a>Supplementary scripts

Scripts are arranged by Figure, and therefore by analysis. Although we are not able to provide all input data files because of the <a target="_blank" href="https://www.gisaid.org/">GISAID</a> privacy agreement, users can apply for their own access. Every effort has been made to describe all steps and input, here or in each script's comments. Where applicable, scripts are arranged in the order they should be executed. The scripts are of two types: 

1. ***Command-line***. These scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual analysis***. These R scripts document the bulk of our data analyses and visualizations. They are intended to be executed manually line-by-line in R/RStudio. The user should replace path names and arguments with the appropriate values for the user's analysis and directories. Attention has been drawn to lines or variables that should be modified by the user with the flag: `CHANGE THIS`.


### <a name="figure-1"></a>Figure 1. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members

* `fig1b_genome_annotation.bash` (*command-line script*)
	* **Description**. Produces Figure 1B using PyGenomeTracks.
	* **Requirements**. PyGenomeTracks, Seqkit.
	* **Input**. The following files need to be in the working directory: 
		1. `SARS-related-CoV_ALN.fasta`: multiple alignment `.fasta` file of n=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus*. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.
		2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within `SARS-related-CoV_ALN.fasta`.
		3. `sbc_rename2.nw`: Newick tree for sarbecovirus alignment
		4. `parameters_input.txt`
		5. `parameters_input2.txt`
	* **Output**. 
		1. `Fig1b.png`
	* **Example**:

			fig1b_genome_annotation.bash

* `ORF_length.R` (*manual analysis script*)
	* **Description**. Analyze the genome-wide ORF length results of the Schlub et al. (2018) codon permutation method and produce Figure 1—figure supplement 1.
	* **Requirements**. tidyverse, RColorBrewer, ggrepel, patchwork.
	* **Input**. 
		1. `frameshift_results.txt`, produced by `frameshift_analysis.bash`.
	* **Output**. 
		1. `ORF_length_heatmap_data.tsv`, p-values for each internal overlapping gene within a previously annotated gene, redundatly annotated on site-by-site basis.


### <a name="figure-2"></a>Figure 2. Re-analysis of SARS-CoV-2 gene expression in publicly available ribosome profiling and mass spectrometry datasets

* `aligned_fasta2haplotypes.pl` (*command-line script*)
	* **Description**. Script to determine all existing haplotypes in a set of sequences. Output used as a list of peptide search queries in mass spectrometry.
	* **Input**. Unnamed arguments in the following order: 
		1. A `.fasta` file containing a multiple sequence amino acid alignment of one protein product of interest
	* **Output**. 
		1. To STDOUT, prints a two-column `.tsv` table: column 1 contain contains the number of occurrences (*n*) of the haplotype; column 2 contains the haplotype sequence itself
	* **Example**:

			aligned_fasta2haplotypes.pl SARSCOV2_ORF3d_aa.fasta

* `riboseq_sliding_window.R` (*command-line script*)
	* **Description**. Sliding window script to calculate proportion of reads in each frame for a specified read length and window size, separately for each condition (treatment/time).
	* **Input**. Unnamed arguments in the following order: 
		1. `INFILE_FRAMES`: file with condition metadata for samples (here, **riboseq_sample_metadata.tsv**). Columns should be in the following order: sample, condition, time, and treatment, where condition is the combination of time and treatment (see script and file)
		2. `INFILE_READS`: file with read data; columns described within script; produced by concatenating the output of ribosome profiling read mapping at the command line as follows:

				grep -E "\s" SRR*.txt | gsed -E s/_all/\\t/ | gsed -E s/ntd-reads_framed_zerobased.txt:/\\t/ | gsed -E s/\\s+/\\t/ > mapped_reads_by_readlength.tsv


		3. `WIN_SIZE`: size of the sliding window (integer)
		4. `READ_LEN`: length of reads to consider (integer)
	* **Output**. 
		1. Table in `.tsv` format giving the sum of reads in each frame for each condition and position (start of window).
	* **Example**:

			Rscript riboseq_sliding_window.R riboseq_sample_metadata.tsv mapped_reads_by_readlength.tsv 30 30

* `riboseq.R` (*manual analysis script*)
	* **Description**. Analyze all ribosome profiling data to generate for analyses underlying Figure 2 and its copious supplement.
	* **Requirements**. R libraries patchwork, RColorBrewer, scales, tidyverse.
	* **Input**. 
		1. `mapped_reads_by_readlength.tsv`, ribosome profiling read depth by sample, read length, position, and frame, after mapping of data from Finkel et al. (2020)
		2. `riboseq_sample_metadata.tsv`, file with condition metadata for samples with columns in the following order: sample, condition, time, and treatment, where condition is the combination of time and treatment (see file)
		3. `mapped_reads_by_readlength_sumSamples_win*read*.tsv`, results of the sliding window script **riboseq\_sliding\_window.R** for various window sizes (win\*) and read length (read\*) summarizing the proportion of reads mapping to each frame for each condition (treatment/time)
	* **Output**. 
		1. Figures and statistics

* `mass-spec_vs_riboseq.R` (*manual analysis script*)
	* **Description**. Analyze and compare expression estimates from mass spectrometry to those from ribosome profiling read depth for Figure 2 and supplement.
	* **Requirements**. R libraries RColorBrewer, scales, tidyverse.
	* **Input**. 
		1. `riboseq_upstream_peaks.tsv`, ribosome profile profiling read depth near gene start sites, shown in Figure 2A, derived from analysis of data from Finkel et al. (2020) and wrangled in the script **riboseq.R**
		2. `mapped_reads_by_readlength_wCDS.tsv`, ribosome profiling read depth by sample, read length, position, and frame. Source data underlying Figure 2C, derived from analysis of data from Finkel et al. (2020), produced in the script **riboseq.R**.
		3. `expression_data_by_gene.tsv`, summary of expression estimates from both ribosome profiling and mass spectrometry for all gene detectable by mass spectrometry
	* **Output**. 
		1. Figures and statistics


### <a name="figure-3"></a>Figure 3. SARS-CoV-2 protein sequence properties

* `generate_random_protein.py` (*command-line script*)
	* **Description**. Script to generate random protein sequences of a given length from an input proteome in `.fasta` format.
	* **Requirements**. Python packages Bio, Bio.Seq, Bio.SeqRecord, os, random, sys
	* **Input**. Unnamed arguments in the following order: 
		1. input file with peptide sequence in `.fasta` format
		2. length of peptide to generate
		3. number of peptides to generate
	* **Output**. 
		1. Multiple sequence alignment in `.fasta` format containing the randomized protein sequence(s).
	* **Example**:

			generate_random_protein.py ORF3d_aa.fasta 57 1000

* `tally_epitope_coverage.py` (*command-line script*)
	* **Description**. Script for tallying epitope coverage for one protein product in a sliding window.
	* **Requirements**. Python packages os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. input file in `.tsv` format; NetMHCpan output file with 5 columns in this order: ID, NB (number MHC alleles bound based on NetMHCpan output), product, codon\_start, codon\_end 
		2. length (integer) of linear peptides used in the epitope analysis (9 for MHC I with NetMHCIIpan; 15 for MHC II with NetMHCIIpan)
	* **Output**. 
		1. Table in `.tsv` format giving the sum of bound epitopes overlapping each site.
	* **Example**:

			tally_epitope_coverage.py ORF3d_random.tsv 9

* `epitope_MHCI.R` (*manual analysis script*)
	* **Description**. Analyze MHC class I epitopes for Figure 3A.
	* **Requirements**. R libraries, ggrepel, patchwork, RColorBrewer, tidyverse.
	* **Input**. 
		1. `NetMHCpan_Wuhan-Hu-1.tsv`, NetMHCpan results for the basic protein set of the <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1 genome</a>, plus *ORF3d*
		2. `NetMHCpan_Wuhan-Hu-1_additional.tsv`, NetMHCpan results for *ORF3c*, *ORF3d-2*, and *ORF3b* (short)
		3. `MHCI_epitope_summary_tally.tsv`, result of **tally_epitope_coverage.py** applied to the proteins encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
		4. `MHCI_short_unannot_ORFs.tsv`, result of **tally_epitope_coverage.py** applied to the short unannotated ORFs encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
		5. `NetMHCpan_*_random_tally.tsv`, result of **tally_epitope_coverage.py** applied to 1000 randomized peptides based on each of the proteins encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
	* **Output**. 
		1. `MHCI_epitope_summary_test.tsv`, compiled results for MHC class I analysis, used for Figure 3A
		2. Figures and statistics

* `epitope_MHCII.R` (*manual analysis script*)
	* **Description**. Analyze MHC class II epitopes for Figure 3A.
	* **Requirements**. R libraries, ggrepel, patchwork, RColorBrewer, tidyverse.
	* **Input**. 
		1. `NetMHCIIpan_Wuhan-Hu-1.tsv`, NetMHCIIpan results for the basic protein set of the <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1 genome</a>, plus *ORF3d*
		2. `NetMHCIIpan_Wuhan-Hu-1_additional.tsv`, NetMHCIIpan results for *ORF3c*, *ORF3d-2*, and *ORF3b* (short)
		3. `MHCII_epitope_summary_tally.tsv`, result of **tally_epitope_coverage.py** applied to the proteins encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
		4. `MHCII_short_unannot_ORFs.tsv`, result of **tally_epitope_coverage.py** applied to the short unannotated ORFs encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
		5. `NetMHCIIpan_*_random_tally.tsv`, result of **tally_epitope_coverage.py** applied to 1000 randomized peptides based on each of the proteins encoded by <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>
	* **Output**. 
		1. `MHCII_epitope_summary_test.tsv`, compiled results for MHC class II analysis, used for Figure 3A
		2. Figures and statistics

* `hydrophobicity_profiles_ORF3a.R` (*manual analysis script*)
	* **Description**. Analyze hydrophobicity profiles of the peptides encoded by the three forward-strand frames of *ORF3a* for Figure 3B, Figure 3—figure supplement 2, and Figure 3—figure supplement 3.
	* **Requirements**. R libraries patchwork, RColorBrewer, scales, tidyverse.
	* **Input**. 
		1. `hydrophobicity_profiles_ORF3a.csv`, output of the <a target="_blank" href="http://volpes.univie.ac.at/">VOLPES server</a> using the unitless hydrophobicity scale “FAC1” (Factor 1) with a sliding window of 25 amino acids
	* **Output**. 
		1. `hydrophobicity_profiles_ORF3a_corr.tsv`, correlations between hydrophobicity profiles in a sliding window of 25 residues
		2. Figures and statistics


### <a name="figure-5"></a>Figure 5. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels

* `SARS-CoV-2_locate_genes.pl` (*command-line script*)
	* **Description**. Script to locate gene start and stop sites by finding the first sequences beginning and ending with (hardcoded) nucleotide sequences taken from the reference sequence <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>. The code itself is a useful resource, as it contains the beginning and ending of most genes.
	* **Requirements**. Perl
	* **Input**. Three unnamed arguments in the following order: 
		1. A `.fasta` file containing a multiple sequence alignment of SARS-CoV-2 whole-genome nucleotide sequences. Note that the script may fail for alignments with gaps or highly diverged from the <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a> genotype.	
	* **Output**. 
		1. To STDOUT, prints a `.gtf` file containing gene coordinates. Note that some genes may be missed if the alignment contains sequences highly diverged from SARS-CoV-2.
	* **Example**:

			SARS-CoV-2_locate_genes.pl SARS-CoV-2_ALN.fasta

* `generate_seqs_from_VCF.py` (*command-line script*)
	* **Description**. Script to generate a `.fasta` with randomly interspersed variants (from VCF), for use with <a target="_blank" href="https://github.com/chasewnelson/OLGenie">OLGenie</a> (Nelson et al. 2020) or any software that requires a multiple sequence alignment input and does not depend upon linkage.
	* **Requirements**. Python packages Bio, os, random, re, sys
	* **Input**. Three unnamed arguments in the following order: 
		1. A `.fasta` file containing exactly one (1) reference sequence (*e.g.*, SARS-CoV-2 <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>)
		2. a `.vcf` file containing the variants of interest
		3. number of sequences to generate (integer)	* **Output**. 
		1. A new `*.fasta` multiple sequence alignment file with variants randomy interpersed at the frequencies defined in the VCF file.
	* **Example**:

			generate_seqs_from_VCF.py reference.fasta variants.vcf 1000

* `three_levels_diversity.R` (*manual analysis script*)
	* **Description**. Input, wrangle, and bootstrap SNPGenie and OLGenie results for each gene and evolutionary level (between-taxa, between-host, and within-host) to produce the analysis underlying Figure 5.
	* **Requirements**. R libraries boot, RColorBrewer, scales, tidyverse.
	* **Input**.
		1. `between-taxa_codon_OL_category.tsv`, metadata file containing the OLG category for each codon in each product in the *Severe acute respiratory syndrome-related coronavirus* alignment, **SARS-related-CoV_ALN.fasta**
		2. `between-taxa_SNPGenie.tsv`, SNPGenie (non-OLG *d*<sub>N</sub>/*d*<sub>S</sub>, , Jukes-Cantor corrected) results for the alignment, **SARS-related-CoV_ALN.fasta**
		3. `between-taxa_SNPGenie_taxaRestricted.tsv`, SNPGenie results for our alignment, **SARS-related-CoV_ALN.fasta**, but limiting to taxa that have intact *ORF3d* (SARS-CoV-2 and pangolin-CoV GX/P5L) or lacking full-length *ORF3b* (9 sequences; see Figure 4)
		4. `between-taxa_OLGenie.tsv`, OLGenie (OLG *d*<sub>N</sub>/*d*<sub>S</sub>, Jukes-Cantor corrected) results for most products and frames of annotated genes in **SARS-related-CoV_ALN.fasta**
		5. `between-taxa_OLGenie2.tsv`, OLGenie results for *ORF3a* in the ss13 frames to analyze *ORF3c* for **SARS-related-CoV_ALN.fasta**
		6. `between-taxa_OLGenie_addNN.tsv`, curated OLGenie results for **SARS-related-CoV_ALN.fasta**, incorporating the Wei-Zhang method for codon 71 of *ORF3a* (NN change in *ORF3d*)
		7. `between-host_SNPGenie.tsv`, SNPGenie (non-OLG *π*<sub>N</sub>/*π*<sub>S</sub>) results for our SARS-CoV-2 (GISAID) alignment
		8. `between-host_OLGenie.tsv`, OLGenie (OLG *π*<sub>N</sub>/*π*<sub>S</sub>) results for our SARS-CoV-2 (GISAID) alignment
		9. `within-host_SNPGenie.tsv`, SNPGenie results for our deeply sequenced within-host samples, before taking means by codon
		10. `within-host_OLGenie.tsv`, OLGenie results for our deeply sequenced within-host samples, based on pseudo-alignments produced with **generate\_seqs\_from\_VCF.py**, before taking means by codon
	* **Output**. 
		1. `between-taxa_FINAL.tsv`, final, bootstrapped, OLG-aware Jukes-Cantor corrected *d*<sub>N</sub>/*d*<sub>S</sub> results for between-taxa (*Severe acute respiratory syndrome-related coronavirus*) data
		2. `between-host_FINAL.tsv`, final, bootstrapped, OLG-aware *π*<sub>N</sub>/*π*<sub>S</sub> results for between-host (GISAID) data
		3. `within-host_FINAL.tsv`: final, bootstrapped, OLG-aware *π*<sub>N</sub>/*π*<sub>S</sub> results for within-host (SRR) data

* `selection_vs_expression.R` (*manual analysis script*)
	* **Description**. Compare selection (*π*<sub>N</sub>/*π*<sub>S</sub> and *d*<sub>N</sub>/*d*<sub>S</sub>) to expression level for Figure 5—figure supplement 1.
	* **Requirements**. R libraries RColorBrewer, scales, tidyverse
	* **Input**. 
		1. `expression_data_by_gene_frame.tsv`, summary of expression estimates from both ribosome profiling and mass spectrometry for all gene detectable by mass spectrometry, by frame (source data of Figure 2B)
		2. `selection_three_levels.txt`, final *π*<sub>N</sub>/*π*<sub>S</sub> or *d*<sub>N</sub>/*d*<sub>S</sub> estimates underlying Figure 5 for each gene be evolutionary level (between-taxa, between-host, within-host) and region type (non-overlapping, overlapping).
	* **Output**. 
		1. Figures and statistics

* `extract_seqs_by_timepoint.py` (*command-line script*)
	* **Description**. Script to extract SARS-CoV-2 GISAID sequences by timepoint in a sliding window for analysis in Figure 5—figure supplement 2. Sliding window size (14 days) and step size (7 days) may be changed in the code. Day 0 is taken to be the earliest date on or following 2019/12/20; sequences sampled at earlier dates will be excluded.
	* **Requirements**. Python packages Bio, datetime, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. the GISAID acknowledgements table, saved as a `.tsv` file with a one-row header; this contains the sampling dates of each sequence. The script expects the sequence ID in column 1 (index 0) and the date in column 4 (index 3)
		2. a `.fasta` multiple sequence alignment file, where the headers are the GISAID IDs, *i.e.*, they match the accession IDs in the GISAID table
	* **Output**. 
		1. One `*.fasta` multiple sequence alignment file for each 14 day window, containing just those sequences sampled during that period. For example, the first file will have the name `*_0to14.fasta` (window 1), the second file will have the name `*_7to21.fasta` (window 2), and so on.
	* **Example**:

			extract_seqs_by_timepoint.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN.fasta

* `extract_seqs_by_location.py` (*command-line script*)
	* **Description**. Script to extract SARS-CoV-2 GISAID sequences by location for analysis in Figure 5—figure supplement 2.
	* **Requirements**. Python packages Bio, datetime, os, sys
	* **Input**. Four unnamed arguments in the following order: 
		1. the GISAID acknowledgements table, saved as a `.tsv` file with a one-row header; this contains the sampling dates of each sequence. The script expects the sequence ID in column 1 (index 0), location in column 3 (index 2), and the date in column 4 (index 3)
		2. a `.fasta` multiple sequence alignment file, where the headers are the GISAID IDs, *i.e.*, they match the accession IDs in the GISAID table
		3. name of desired output file
		4. location, *i.e.*, string to detect in `Location` columns
	* **Output**. 
		1. One `*.fasta` multiple sequence alignment file with the name given by input argument 3 above, containing only those sequences matching the location given in argument 4 above.
	* **Example**:

			extract_seqs_by_location.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN.fasta SARS-CoV-2_ALN_Asia.fasta Asia

* `extract_positions_by_timepoint.py` (*command-line script*)
	* **Description**. Track the frequency of alleles at specific sites in a sliding time window for allele trajectory analysis in Figure 5—figure supplement 2.
	* **Requirements**. Python packages Bio, datetime, os, sys
	* **Input**. Six unnamed arguments in the following order: 
		1. the GISAID acknowledgements table, saved as a `.tsv` file with a one-row header; this contains the sampling dates of each sequence. The script expects the sequence ID in column 1 (index 0), and the date in column 4 (index 3)
		2. a `.fasta` multiple sequence alignment file, where the headers are the GISAID IDs, *i.e.*, they match the accession IDs in the GISAID table
		3. a `.tsv` file with a single column: a list of sites (1-based genome positions) to track
		4. name of desired output file (`.tsv`)
		5. window size in days (int), the size of the windows in which to calculate allele frequencies
		6. step size in days (int), how many days to move forward for each successive window
	* **Output**. 
		1. One `.tsv` table with the name given by input argument 4 above, where each row is a unique combination of timepoint and genome position, and columns report the numbers and frequencies of the alleles for that timepoint and position.
	* **Example**:

			extract_positions_by_timepoint.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN_Asia.fasta sites_to_track.tsv SARS-CoV-2_ALN_Asia_tracked.fasta 14 1

* `temporal_pi.R` (*manual analysis script*)
	* **Description**. Calculate and plot nucleotide diversity (*π*) and location entropy as a function of time for Figure 5—figure supplement 2. The user must first use <a target="_blank" href="https://github.com/chasewnelson/SNPGenie">SNPGenie</a> (`snpgenie_within_group.pl`) to analyze the results of `extract_seqs_by_timepoint.py`, separately for each window.
	* **Requirements**. R libraries boot, patchwork, RColorBrewer, scales, tidyverse
	* **Input**.
		1. `within_group_codon_results.tsv`, combined across timepoints
		2. `timepoint_seq_IDs.txt`, a table with sequence IDs (column `ID`,  *i.e.*, EPI_*) for each time window (column `time_period`, *e.g.*, "0to14")
		3. `gisaid_cov2020_acknowledgement_table.tsv`
	* **Output**. 
		1. Figures and statistics




### <a name="figure-6"></a>Figure 6. Between-taxa sliding window analysis of natural selection on overlapping frames of *ORF3a*

* `selection_sliding_windows.R` (*manual analysis script*) 
	* **Description**. Visualize overlapping gene *d*<sub>N</sub>/*d*<sub>S</sub>, calculated with OLGenie, in sliding windows for pairs of taxa shown in Figure 6 and supplement
	* **Requirements**. R libraries boot, RColorBrewer, tidyverse.
	* **Input**.
		1. `SARS-CoV-2-ref_ORF3a_ss12_windows_prepangolin.txt`, OLGenie (*d*<sub>NN</sub>/*d*<sub>NS</sub>) results in a sliding window (window size=50 codons, step size=1 codon) for SARS-CoV-2 against each other taxon (every pair) for *ORF3a* reading frame ss12 (i.e., the frame of *ORF3d*)
		2. `SARS-CoV-2-ref_ORF3a_ss12_windows_pangolin.txt`, updated OLGenie results in a sliding window for SARS-CoV-2 against pangolin GX/P5L, *i.e.*, the only other taxon in which *ORF3d* is intact, meant to replace the data in **SARS-CoV-2-ref_ORF3a_ss12_windows.txt**
		3. `SARS-CoV-2-ref_N_ss13_windows.tsv`, OLGenie results in a sliding window for SARS-CoV-2 against each other taxon (every pair) for *N* reading frame ss13 (i.e., the frame of *ORF9b* and *ORF9c*)
		4. `SARS-CoV-ref_ORF3a_ss13_windows.txt`, OLGenie results in a sliding window for SARS-CoV against each other taxon (every pair) for *ORF3a* reading frame ss13 (i.e., the frame of *ORF3c* and *ORF3b*)
		5. `SARS-CoV-ref_N_ss13_windows.txt`, OLGenie results in a sliding window for SARS-CoV against each other taxon (every pair) for *N* reading frame ss13 (i.e., the frame of *ORF9b* and *ORF9c*)
	* **Output**. 
		1. Figures and statistics


### <a name="figure-7"></a>Figure 7. Pandemic spread of the EP+1 haplotype and the hitchhiking of *ORF3d*-LOF

* `extract_variable_columns_MSA.py` (*command-line script*)
	* **Description**. Script for examining a multiple sequence alignment to identify variable (segregating) sites
	* **Requirements**. Python packages Bio, os, re, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. a `.fasta` file containing aligned nucleotide sequences 
		2. a minimum minor allele frequency to allow (numeric)
	* **Output**. 
		1. To a file named `*_variants_MAF[\d\.].fasta`, prints a `.fasta` multiple sequence alignment of just the variable sites
		2. To a file named `*_variants_MAF[\d\.].tsv`, prints a `.tsv` table giving the nucleotide present in each sequence (row) at each variable site (column)
		2. To STDOUT, reports the positions of the variable sites
	* **Example**:

			extract_variable_columns_MSA.py SARS-CoV-2_ALN.fasta 0.02

* `extract_variable_columns_MSA_aa.py` (*command-line script*)
	* **Description**. Same as `extract_variable_columns_MSA.py`, but for amino acid sequences, necessary to account for STOP codons (potential non-word characters).
	* **Example**:

			extract_variable_columns_MSA_aa.py SARS-CoV-2_ORF3d_aa_ALN.fasta 0.02

* `Fig7.m`
	* **Input**
		1. `SARSCOV2_MAFFT_processed_variants_MAF0.02_ManualCheck_Fig7.xlsx`, file from output, with manual correction of the entry with different ordering format of the collection date

* `Fig7b.m`
	* **Description**. Calculate the haplotype frequency trajectory and output figure.
	* **Output** 
		1. `Fig7b.jpg` --> input in Fig7.m
		2. `Fig7a.jpg` (artwork, input in Fig7.m)
		3. `Fig7.m` (figure collage): output Fig7.jpg


### <a name="figure-8"></a>Figure 8. High-frequency within-host mutations

* `filter_VCF.py` (*command-line script*)
	* **Description**. Script to dynamically filter within-host variants using a binomial cutoff to control for a user-defined false-discovery rate (FDR). Automatically detects and analyzes all `.vcf` (variant call format) files in the working directory. For use with <a target="_blank" href="https://github.com/chasewnelson/SNPGenie">SNPGenie</a> (Nelson et al. 2015) or any software requiring FDR-filtered variants in `.vcf` files.
	* **Requirements**. Python packages Bio, numpy, os, random, re, scipy.stats, sys
	* **Input**. One or more `.vcf` files in the working directory, and five unnamed arguments in the following order: 
		1. analysis-wide FDR cutoff (integer): the maximum absolute number of false-positive variants called across all VCF files examined
		2. the minimum allowed minor allele frequency (numeric), independent of FDR
		3. genome length (integer)
		4. sequencing error rate per site (numeric): assumes all nucleotide changes are equally probable	
		5. number of samples in analysis (integer), i.e., total number of deeply sequenced within-host samples (here, the number of VCF files)
	* **Output**. 
		1. Creates a new (`*_filtered.vcf`) VCF file for each VCF file in the working directory, including only those variants that pass the FDR threshold.
		2. Prints summary statistics to STDOUT regarding the FREQUENCY, AC (allele count), and DP (total read depth at site) of passing and failing variants.
	* **Example**:

			filter_VCF.py 1 0 29903 0.002 401

* `summarize_intrahost_by_site.py` (*command-line script*)
	* **Description**. Script to create a genome database cataloguing within-host variants. Automatically detects and analyzes all `.vcf` (variant call format) files in the working directory.
	* **Requirements**. Python packages Bio, Bio.alphabet, Bio.seq, os, re, sys
	* **Input**. One or more `.vcf` files in the working directory, and two unnamed arguments in the following order: 
		1. A `.fasta` file containing exactly one (1) reference sequence (*e.g.*, SARS-CoV-2 <a target="_blank" href="https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3">Wuhan-Hu-1</a>)
		2. A `.gtf` file containing containing genes to be annotated in the output table, with up to two genes overlapping each site (others will be ignored)
	* **Output**. 
		1. An output file by the name `*_site_database.tsv`. There are four rows for each position in the genome (defined by input 1), corresponding to each of the four possible nucleotide changes (including self-nucleotide). For example, a position with A in the reference (REF), there will four possible single nucleotide changes (ALT): A (self), C, G, and T. Each row is also labelled with up to two genes overlapping the site, and the codon, codon position, and amino acid encoded by each gene. Finally, each column following the metadata is a sample, giving the number of `REF,ALT` reads in that same at that position, if its VCF file contains a record. Note that the beginning of the file is largely unpopulated, as the first rows correspond to the 5'-UTR region lacking genes and coverage.
	* **Example**:

			summarize_intrahost_by_site.py NC_045512.fasta NC_045512.gtf

* `Fig8.m`
	* **Description**. Calculate the frequency of alternative allele and partition by the type of mutation and output figure.
	* **Input**
		1. `NC_045512_site_database_altA.tsv`
		2. `NC_045512_site_database_altC.tsv`
		3. `NC_045512_site_database_altG.tsv`
		4. `NC_045512_site_database_altT.tsv`


### <a name="additional-scripts"></a>Additional scripts

* `extract_fasta_by_sites.pl` (*command-line script*)
	* **Description**. Script to extract (excise; cut out) regions of a multiple sequence alignment, e.g., pull out specific genes.
	* **Requirements**. Perl
	* **Input**. Two unnamed arguments in the following order: 
		1. a `.fasta` file containing containing one or more aligned sequences 
		2. a`.gtf` file containing CDS products to extract. They need not really be CDS, but should be labelled as such in the file.	Only works for forward-strand (+) products.
	* **Output**. 
		1. one `.fasta` file for each CDS record in the GTF file. Just the DNA segment corresponding to the CDS coordinates of each record will be present in the resulting `.fasta` files.
	* **Example**:

			extract_fasta_by_sites.pl Wuhan_Hu_1.fasta Wuhan_Hu_1.gtf

* `extract_seq_subset.py` (*command-line script*)
	* **Description**.  Script for extracting a subset of sequences from a `.fasta` based on header ID.
	* **Requirements**. Python packages Bio, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. a text file containing exactly one column of sequence IDs (`.fasta` headers) 
		2. a multiple sequence alignment of nucleotides in `.fasta` format	* **Output**. 
		1. A multiple sequence alignment in `.fasta` format based on the file given by argument 2, but including only those sequences with headers provided in the file given by argument 1.
	* **Example**:

			extract_seq_subset.py seq_ID_list.txt SARS-COV-2_ALN.fasta

* `translate_nt_seq_file.pl` (*command-line script*)
	* **Description**. Script to translate a file of (un-aligned) nucleotide sequences and print the proteins.
	* **Requirements**. Perl
	* **Input**. One unnamed argument: 
		1. one nucleotide `.fasta` file containing the original (un-aligned) coding nucleotide sequences (complete codon sets, *i.e.*, multiples of 3)
	* **Output**. 
		1. To STDOUT, print the translated sequences.
	* **Example**:

			translate_nt_seq_file.pl coding_nt_seqs.fasta


## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Postdoctoral Research Fellowship from Academia Sinica (to C.W.N. under P.I. Wen-Hsiung Li); funding from the Bavarian State Government and National Philanthropic Trust (to Z.A. under P.I. Siegfried Scherer); NSF IOS grants #1755370 and #1758800 (to S.-O.K.); and the University of Wisconsin-Madison John D. MacArthur Professorship Chair (to T.L.G). Copyright-free images were obtained from Pixabay. The authors thank the GISAID platform and the originating and submitting laboratories who kindly uploaded SARS-CoV-2 sequences to the GISAID EpiCov™ Database for public access (Supplement). The authors thank Maciej F. Boni, Reed A. Cartwright, John Flynn, Kyle Friend, Dan Graur, Robert S. Harbert, Cheryl Hayashi, David G. Karlin, Niloufar Kavian, Kin-Hang (Raven) Kok, Wen-Hsiung Li, Meiyeh Lu, David A. Matthews, Lisa Mirabello, Apurva Narechania, Felix Li Jin, and attendees of the UC Berkeley popgen journal club for useful information and discussion; Andrew E. Firth, Alexander Gorbalenya, Irwin Jungreis, Manolis Kellis, Raven Kok, Angelo Pavesi, Kei Sato, Manuela Sironi, and Noam Stern-Ginossar for an invaluable discussion regarding standardizing nomenclature; Helen Piontkivska, Patricia Wittkopp, Antonis Rokas, and one anonymous reviewer for critical suggestions; Ming-Hsueh Lin for immense feedback on figures; and special thanks to Priya Moorjani, Jacob Tennessen, Montgomery Slatkin, Yun S. Song, Jianzhi George Zhang, Xueying Li, Hongxiang Zheng, Qinqin Yu, Meredith Yeager, and Michael Dean for commenting on earlier drafts of the manuscript.


## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Ardern Z, Goldberg TL, Meng C, Kuo C-H, Ludwig C, Kolokotronis S-O, Wei X. 2020. 
<a target="_blank" href="https://elifesciences.org/articles/59633">Dynamically evolving novel overlapping gene as a factor in the SARS-CoV-2 pandemic.</a> *eLife* **9**: e59633. DOI: 10.7554/eLife.59633

and this page:

>https://github.com/chasewnelson/SARS-CoV-2-ORF3d


## <a name="contact"></a>Contact and troubleshooting

If you have questions about our scripts or study, please first thoroughly read the documentation and in-line comments relevant to the script of interest. If these do not answer your question, please click on the <a target="_blank" href="https://github.com/chasewnelson/SARS-CoV-2-ORF3d/issues">Issues</a> tab at the top of this page and search to see if your question has already been answered; if not, begin a new issue, so that others might benefit from the discussion.

Other queries should be addressed to the corresponding authors: 

*  Chase W. Nelson, cnelson <**AT**> gate <**DOT**> sinica <**DOT**> edu <**DOT**> tw
*  Zachary Ardern, zachary <**DOT**> ardern <**AT**> tum <**DOT**> de
*  Xinzhu (April) Wei, aprilwei <**AT**> berkeley <**DOT**> edu


## <a name="references"></a>References

* Finkel Y, Mizrahi O, Nachshon A, Weingarten-Gabbay S, Morgenstern D, Yahalom-Ronen Y, Tamir H, Achdout H, Stein D, Israeli O, Adi B-D, Melamed S, Weiss S, Israely T, Paran N, Schwartz M, Stern-Ginossar N. 2020. <a target="_blank" href="https://www.nature.com/articles/s41586-020-2739-1">The coding capacity of SARS-CoV-2</a>. *Nature*, **in press**, doi: https://doi.org/10.1038/s41586-020-2739-1
* Nelson CW, Ardern Z, Wei X. 2020. <a target="_blank" href="https://academic.oup.com/mbe/article/37/8/2440/5815567">OLGenie: estimating natural selection to predict functional overlapping genes</a>. *Molecular Biology and Evolution* **37**(8):2440–2449. doi: https://doi.org/10.1093/molbev/msaa087
* Nelson CW, Ardern Z, Goldberg TL, Meng C, Kuo C-H, Ludwig C, Kolokotronis S-O, Wei X. 2020. 
<a target="_blank" href="https://elifesciences.org/articles/59633">Dynamically evolving novel overlapping gene as a factor in the SARS-CoV-2 pandemic.</a> *eLife* **9**: e59633. doi: 10.7554/eLife.59633
* Nelson CW, Moncla LH, Hughes AL. 2015. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/31/22/3709/241742">SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data</a>. *Bioinformatics* **31**(22):3709–3711. doi: https://doi.org/10.1093/bioinformatics/btv449
* Schlub TE, Buchmann JP, Holmes EC. 2018. <a target="_blank" href="https://academic.oup.com/mbe/article/35/10/2572/5067730">A simple method to detect candidate overlapping genes in viruses using single genome sequences</a>. *Molecular Biology and Evolution* **35**(10):2572–2581. doi: https://doi.org/10.1093/molbev/msy155
