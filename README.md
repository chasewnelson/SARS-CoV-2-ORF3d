<img src="https://github.com/chasewnelson/SARS-CoV-2-ORF3d/blob/master/images/cover_image.png?raw=true" title="Overlapping gene products" alt="Overlapping gene products" align="left" size="small">

# SARS-CoV-2 *ORF3d*

Supplementary data and scripts for Nelson et al. (2020) paper on SARS-CoV-2 *ORF3d*.


## <a name="contents"></a>Contents

* [Supplementary data](#supplementary-data)
	* `SARS-related-CoV_ALN.fasta`: whole-genome multiple sequence alignment of *n*=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus* (between-taxa analysis)
	* `SARS-related-CoV_ALN.gtf`: a `.gtf` file giving gene positions within `SARS-related-CoV_ALN.fasta`
	* `Supplementary_Tables.xlsx`: Supplementary Tables referred to in the [manuscript](#citation)
* [Supplementary scripts](#supplementary-scripts)
	* [**Figure 1**. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members](#figure-1).
		* `fig1B.bash`
		* `ORF_length.R`: analyze ORF lengths for analysis in Figure 1—figure supplement 1
	* [**Figure 2**. Re-analysis of SARS-CoV-2 gene expression in publicly available ribosome profiling and mass spectrometry datasets](#figure-2).
		* `aligned_fasta2haplotypes.pl`: determine all existing haplotypes in a set of sequences to compile a list of peptide search queries
		* `riboseq_sliding_window.R`: calculate proportion of ribosome profiling reads in each frame
	* [**Figure 3**. SARS-CoV-2 protein sequence properties](#figure-3).
		* `generate_random_protein.py`: generate random proteins given amino acid content
		* `tally_epitope_coverage.py`: tally epitope coverage in a sliding window
		* `epitope_MHCI.R`: analyze MHC class I epitopes for Figure 3A
		* `epitope_MHCII.R`: analyze MHC class II epitopes for Figure 3A
	* [**Figure 4**. Amino acid variation in proteins encoded by genes overlapping *ORF3a* in viruses of the species *Severe acute respiratory syndrome-related coronavirus*](#figure-4).
	* [**Figure 5**. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels](#figure-5).
		* `SARS-CoV-2_locate_genes.pl`: automatically find the coordinates of each SARS-CoV-2 gene in a nucleotide multiple sequence alignment
		* `generate_seqs_from_VCF.py`: generate <a target="_blank" href="https://github.com/chasewnelson/OLGenie">OLGenie</a> input for within-host analysis
		* `extract_seqs_by_timepoint.py`: extract SARS-CoV-2 GISAID sequences by timepoint in a sliding window for analysis in Figure 5—figure supplement 2
		* `extract_seqs_by_location.py`: extract sequences by location for location-specific timepoint analyses in Figure 5—figure supplement 2
		* `temporal_pi.R`: calculate and plot nucleotide diversity (*π*) and location entropy as a function of time for Figure 5—figure supplement 2
		
	* [**Figure 6**. Between-taxa sliding window analysis of natural selection on overlapping frames of *ORF3a*](#figure-6).
	* [**Figure 7**. Pandemic spread of the EP+1 haplotype and the hitchhiking of *ORF3d*-LOF](#figure-7).
		* `extract_variable_columns_MSA.py`: identify variable sites in a nucleotide multiple sequence alignment
		* `extract_variable_columns_MSA_aa.py`: identify variable sites in an amino acid multiple sequence alignment
	* [**Figure 8**. High-frequency within-host mutations](#figure-8).
		* `filter_vcf.py`: apply a binomial false-discovery rate correction to within-host variants
		* `summarize_intrahost_by_site.py`: create a genome database cataloguing within-host variants
	* [**Additional scripts**](#additional-scripts).
		* `extract_fasta_by_sites.pl`: extracts gene regions of a multiple sequence alignment
		* `extract_seq_subset.py`: extract a subset of sequences from a `.fasta` file
		* `translate_nt_seq_file.pl`: translate a file of protein-coding nucleotide sequences
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)


## <a name="supplementary-data"></a>Supplementary data

The following supplementary data are provided in this GitHub repository:

1. `SARS-related-CoV_ALN.fasta`: whole-genome multiple sequence alignment of *n*=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus* (between-taxa analysis). See [manuscript](#citation) for details. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.

2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within `SARS-related-CoV_ALN.fasta`.

3. `Supplementary_Tables.xlsx`: Supplementary Tables referred to in the [manuscript](#citation).


## <a name="supplementary-scripts"></a>Supplementary scripts

Scripts are arranged by Figure, and therefore by analysis. Although we are not able to provide all input data files because of the <a target="_blank" href="https://www.gisaid.org/">GISAID</a> privacy agreement, users can apply for their own access. Every effort has been made to describe all steps and input, here or in each script's comments. Where applicable, scripts are arranged in the order they should be executed. The scripts are of two types: 

1. ***Command-line*** scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual analysis*** scripts are R scripts or Python Jupyter Notebooks documenting the bulk of our data analyses and visualizations. These are intended to be executed manually line-by-line in R/RStudio or Jupyter. The use should replace path names and arguments with the appropriate values for the user's analysis and directories. Attention has been drawn to lines or variables that should be modified by the user with the flag: `<-- CHANGE THIS`


### TEMPLATE
* `name.py` (*command-line script*)
	* **Description**. XXX.
	* **Requirements**. XXX.
	* **Input**.
		1. `frameshift_results.txt`
	* **Intermediate files**.
		1. `ORF_length_heatmap_data.tsv`
	* **Output**. XXX.
	* **Example**:

```Shell
XXX
```
		

### <a name="figure-1"></a>Figure 1. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members

* `fig1B.bash` (*command-line script*)
	* **Description**. Produces Figure1B using PyGenomeTracks.
	* **Requirements**. PyGenomeTracks, Seqkit.
	* **Input**. The following files need to be in the working directory: 
		1. `SARS-related-CoV_ALN.fasta`: multiple alignment `.fasta` file of n=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus*. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.
		2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within `SARS-related-CoV_ALN.fasta`.
		3. `sbc_rename2.nw`: Newick tree for sarbecovirus alignment **!!TODO: Zac, I need this <--**
		4. `parameters_input.txt` **!!TODO: Zac, I need this <--**
		5. `parameters_input2.txt` **!!TODO: Zac, I need this <--**
	* **Output**. 
		1. `Fig1b.png`
	* **Example**:

```Shell
fig1B.bash
```

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

```Shell
aligned_fasta2haplotypes.pl SARSCOV2_ORF3d_aa.fasta
```

* `riboseq_sliding_window.R` (*command-line script*)
	* **Description**. Sliding window script to calculate proportion of reads in each frame for a specified read length and window size, separately for each treatment combination.
	* **Input**. Unnamed arguments in the following order: 
		1. `INFILE_FRAMES`: file with condition metadata for samples. Columns described within script
		2. `INFILE_READS`: file with read data. Columns described within script
		3. `WIN_SIZE`: size of the sliding window (integer)
		4. `READ_LEN`: length of reads to consider (integer)
	* **Output**. 
		1. Table in `.tsv` format giving the sum of reads in each frame for each condition and position (start of window).
	* **Example**:

```Shell
Rscript riboseq_sliding_window.R frames_table.txt mapped_reads_by_readlength_ALL.tsv 30 30
```


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

```Shell
generate_random_protein.py ORF3d_aa.fasta 57 1000
```

* `tally_epitope_coverage.py` (*command-line script*)
	* **Description**. Script for tallying epitope coverage for one protein product in a sliding window.
	* **Requirements**. Python packages os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. input file in `.tsv` format; NetMHCpan output file with 5 columns in this order: ID, NB (number MHC alleles bound based on NetMHCpan output), product, codon_start, codon_end 
		2. length (integer) of linear peptides used in the epitope analysis (9 for MHC I with NetMHCIIpan; 15 for MHC II with NetMHCIIpan)
	* **Output**. 
		1. Table in `.tsv` format giving the sum of bound epitopes overlapping each site.
	* **Example**:

```Shell
tally_epitope_coverage.py ORF3d_random.tsv 9
```

* `epitope_MHCI.R` (*manual analysis script*)
	* **Description**. Analyze MHC class I epitopes for Figure 3A.
	* **Requirements**. R libraries, ggrepel, patchwork, RColorBrewer, tidyverse.
	* **Input**. 
		1. `frameshift_results.txt`, produced by `frameshift_analysis.bash`.
	* **Output**. 
		1. `ORF_length_heatmap_data.tsv`, p-values for each internal overlapping gene within a previously annotated gene, redundatly annotated on site-by-site basis.






### <a name="figure-5"></a>Figure 5. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels

* `SARS-CoV-2_locate_genes.pl` (*command-line script*)
	* **Description**. Script to locate gene start and stop sites by finding the first sequences beginning and ending with (hardcoded) nucleotide sequences taken from the reference sequence: https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3. The code itself is a useful resource, as it contains the beginning and ending of most genes.
	* **Requirements**. Perl
	* **Input**. Three unnamed arguments in the following order: 
		1. A `.fasta` file containing a multiple sequence alignment of SARS-CoV-2 whole-genome nucleotide sequences. Note that the script may fail for alignments with gaps or highly diverged from the Wuhan-Hu-1 genotype.	
	* **Output**. 
		1. To STDOUT, prints a `.gtf` file containing gene coordinates. Note that some genes may be missed if the alignment contains sequences highly diverged from SARS-CoV-2.
	* **Example**:

```Shell
SARS-CoV-2_locate_genes.pl SARS-CoV-2_ALN.fasta
```

* `generate_seqs_from_VCF.py` (*command-line script*)
	* **Description**. Script to generate a `.fasta` with randomly interspersed variants (from VCF), for use with <a target="_blank" href="https://github.com/chasewnelson/OLGenie">OLGenie</a> (Nelson et al. 2020) or any software that requires a multiple sequence alignment input and does not depend upon linkage.
	* **Requirements**. Python packages Bio, os, random, re, sys
	* **Input**. Three unnamed arguments in the following order: 
		1. A `.fasta` file containing exactly one (1) reference sequence (*e.g.*, SARS-CoV-2 Wuhan-Hu-1)
		2. a `.vcf` file containing the variants of interest
		3. number of sequences to generate (integer)	* **Output**. 
		1. A new `*.fasta` multiple sequence alignment file with variants randomy interpersed at the frequencies defined in the VCF file.
	* **Example**:

```Shell
generate_seqs_from_VCF.py reference.fasta variants.vcf 1000
```

* `extract_seqs_by_timepoint.py` (*command-line script*)
	* **Description**. Script to extract SARS-CoV-2 GISAID sequences by timepoint in a sliding window for analysis in Figure 5—figure supplement 2. Sliding window size (14 days) and step size (7 days) may be changed in the code. Day 0 is taken to be the earliest date on or following 2019/12/20; sequences sampled at earlier dates will be excluded.
	* **Requirements**. Python packages Bio, datetime, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. the GISAID acknowledgements table, saved as a `.tsv` file with a one-row header; this contains the sampling dates of each sequence. The script expects the sequence ID in column 1 (index 0) and the date in column 4 (index 3)
		2. a `.fasta` multiple sequence alignment file, where the headers are the GISAID IDs, *i.e.*, they match the accession IDs in the GISAID table
	* **Output**. 
		1. One `*.fasta` multiple sequence alignment file for each 14 day window, containing just those sequences sampled during that period. For example, the first file will have the name `*_0to14.fasta` (window 1), the second file will have the name `*_7to21.fasta` (window 2), and so on.
	* **Example**:

```Shell
extract_seqs_by_timepoint.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN.fasta
```

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

```Shell
extract_seqs_by_location.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN.fasta SARS-CoV-2_ALN_Asia.fasta Asia
```

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

```Shell
extract_positions_by_timepoint.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN_Asia.fasta sites_to_track.tsv SARS-CoV-2_ALN_Asia_tracked.fasta 14 1
```

* `temporal_pi.R` (*manual analysis script*)
	* **Description**. Calculate and plot nucleotide diversity (*π*) and location entropy as a function of time for Figure 5—figure supplement 2. The user must first use <a target="_blank" href="https://github.com/chasewnelson/SNPGenie">SNPGenie</a> (`snpgenie_within_group.pl`) to analyze the results of `extract_seqs_by_timepoint.py`, separately for each window.
	* **Requirements**. R libraries boot, patchwork, RColorBrewer, scales, tidyverse
	* **Input**. Three unnamed arguments in the following order: 
		1. `within_group_codon_results.tsv`, combined across timepoints
		2. `timepoint_seq_IDs.txt`, a table with sequence IDs (column `ID`,  *i.e.*, EPI_*) for each time window (column `time_period`, *e.g.*, "0to14")
		3. `gisaid_cov2020_acknowledgement_table.tsv`
	* **Output**. 
		1. Figures and statistics related to Figure 5—figure supplement 2


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

```Shell
extract_variable_columns_MSA.py SARS-CoV-2_ALN.fasta 0.02
```

* `extract_variable_columns_MSA_aa.py` (*command-line script*)
	* **Description**. Same as `extract_variable_columns_MSA.py`, but for amino acid sequences, necessary to account for STOP codons (potential non-word characters).
	* **Example**:

```Shell
extract_variable_columns_MSA_aa.py SARS-CoV-2_ORF3d_aa_ALN.fasta 0.02
```



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

```Shell
filter_VCF.py 1 0 29903 0.002 401
```

* `summarize_intrahost_by_site.py` (*command-line script*)
	* **Description**. Script to create a genome database cataloguing within-host variants. Automatically detects and analyzes all `.vcf` (variant call format) files in the working directory.
	* **Requirements**. Python packages Bio, Bio.alphabet, Bio.seq, os, re, sys
	* **Input**. One or more `.vcf` files in the working directory, and two unnamed arguments in the following order: 
		1. A `.fasta` file containing exactly one (1) reference sequence (*e.g.*, SARS-CoV-2 Wuhan-Hu-1)
		2. A `.gtf` file containing containing genes to be annotated in the output table, with up to two genes overlapping each site (others will be ignored)
	* **Output**. 
		1. An output file by the name `*_site_database.tsv`. There are four rows for each position in the genome (defined by input 1), corresponding to each of the four possible nucleotide changes (including self-nucleotide). For example, a position with A in the reference (REF), there will four possible single nucleotide changes (ALT): A (self), C, G, and T. Each row is also labelled with up to two genes overlapping the site, and the codon, codon position, and amino acid encoded by each gene. Finally, each column following the metadata is a sample, giving the number of `REF,ALT` reads in that same at that position, if its VCF file contains a record. Note that the beginning of the file is largely unpopulated, as the first rows correspond to the 5'-UTR region lacking genes and coverage.
	* **Example**:

```Shell
summarize_intrahost_by_site.py NC_045512.fasta NC_045512.gtf
```


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

```Shell
extract_fasta_by_sites.pl Wuhan_Hu_1.fasta Wuhan_Hu_1.gtf
```

* `extract_seq_subset.py` (*command-line script*)
	* **Description**.  Script for extracting a subset of sequences from a `.fasta` based on header ID.
	* **Requirements**. Python packages Bio, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. a text file containing exactly one column of sequence IDs (`.fasta` headers) 
		2. a multiple sequence alignment of nucleotides in `.fasta` format	* **Output**. 
		1. A multiple sequence alignment in `.fasta` format based on the file given by argument 2, but including only those sequences with headers provided in the file given by argument 1.
	* **Example**:

```Shell
extract_seq_subset.py seq_ID_list.txt SARS-COV-2_ALN.fasta
```

* `translate_nt_seq_file.pl` (*command-line script*)
	* **Description**. Script to translate a file of (un-aligned) nucleotide sequences and print the proteins.
	* **Requirements**. Perl
	* **Input**. One unnamed argument: 
		1. one nucleotide `.fasta` file containing the original (un-aligned) coding nucleotide sequences (complete codon sets, *i.e.*, multiples of 3)
	* **Output**. 
		1. To STDOUT, print the translated sequences.
	* **Example**:

```Shell
translate_nt_seq_file.pl coding_nt_seqs.fasta
```


## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Postdoctoral Research Fellowship from Academia Sinica (to C.W.N. under P.I. Wen-Hsiung Li); funding from the Bavarian State Government and National Philanthropic Trust (to Z.A. under P.I. Siegfried Scherer); NSF IOS grants #1755370 and #1758800 (to S.-O.K.); and the University of Wisconsin-Madison John D. MacArthur Professorship Chair (to T.L.G). Copyright-free images were obtained from Pixabay. The authors thank the GISAID platform and the originating and submitting laboratories who kindly uploaded SARS-CoV-2 sequences to the GISAID EpiCov™ Database for public access (Supplement). The authors thank Maciej F. Boni, Reed A. Cartwright, John Flynn, Kyle Friend, Dan Graur, Robert S. Harbert, Cheryl Hayashi, David G. Karlin, Niloufar Kavian, Kin-Hang (Raven) Kok, Wen-Hsiung Li, Meiyeh Lu, David A. Matthews, Lisa Mirabello, Apurva Narechania, Felix Li Jin, and attendees of the UC Berkeley popgen journal club for useful information and discussion; Andrew E. Firth, Alexander Gorbalenya, Irwin Jungreis, Manolis Kellis, Raven Kok, Angelo Pavesi, Kei Sato, Manuela Sironi, and Noam Stern-Ginossar for an invaluable discussion regarding standardizing nomenclature; Helen Piontkivska, Patricia Wittkopp, Antonis Rokas, and one anonymous reviewer for critical suggestions; Ming-Hsueh Lin for immense feedback on figures; and special thanks to Priya Moorjani, Jacob Tennessen, Montgomery Slatkin, Yun S. Song, Jianzhi George Zhang, Xueying Li, Hongxiang Zheng, Qinqin Yu, Meredith Yeager, and Michael Dean for commenting on earlier drafts of the manuscript.


## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Ardern Z, Goldberg TL, Meng C, Kuo C-H, Ludwig C, Kolokotronis S-O, Wei X. 2020. Dynamically evolving novel overlapping gene as a factor in the SARS-CoV-2 pandemic. *eLife*, in review. bioRxiv doi: <a target="_blank" rel="noopener noreferrer" href="https://doi.org/10.1101/2020.05.21.109280">https://doi.org/10.1101/2020.05.21.109280</a>

and this page:

>https://github.com/chasewnelson/SARS-CoV-2-ORF3d


## <a name="contact"></a>Contact and troubleshooting

If you have questions about our scripts or study, please first thoroughly read the documentation and in-line comments relevant to the script of interest. If these do not answer your question, please click on the <a target="_blank" href="https://github.com/chasewnelson/SARS-CoV-2-ORF3d/issues">Issues</a> tab at the top of this page and search to see if your question has already been answered; if not, begin a new issue, so that others might benefit from the discussion.

Other queries should be addressed to the corresponding authors: 

*  Chase W. Nelson, cnelson <**AT**> gate <**DOT**> sinica <**DOT**> edu <**DOT**> tw
*  Zachary Ardern, zachary <**DOT**> ardern <**AT**> tum <**DOT**> de
*  Xinzhu (April) Wei, aprilwei <**AT**> berkeley <**DOT**> edu


## <a name="references"></a>References

* Nelson CW, Ardern Z, Wei X. 2020. <a target="_blank" href="https://academic.oup.com/mbe/article/37/8/2440/5815567">OLGenie: estimating natural selection to predict functional overlapping genes</a>. *Molecular Biology and Evolution* **37**(8):2440–2449. doi: https://doi.org/10.1093/molbev/msaa087
* Nelson CW, Moncla LH, Hughes AL. 2015. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/31/22/3709/241742">SNPGenie: estimating evolutionary parameters to detect natural selection using pooled next-generation sequencing data</a>. *Bioinformatics* **31**(22):3709–3711. doi: https://doi.org/10.1093/bioinformatics/btv449
* Schlub TE, Buchmann JP, Holmes EC. 2018. <a target="_blank" href="https://academic.oup.com/mbe/article/35/10/2572/5067730">A simple method to detect candidate overlapping genes in viruses using single genome sequences</a>. *Molecular Biology and Evolution* **35**(10):2572–2581. doi: https://doi.org/10.1093/molbev/msy155
