# SARS-CoV-2 *ORF3d*

Supplementary scripts for Nelson et al. (2020) paper on SARS-CoV-2 *ORF3d*.


## <a name="contents"></a>Contents

* [Supplementary data](#supplementary-data)
* [Supplementary scripts](#supplementary-scripts)
	* [**Figure 1**. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members](#figure-1).
		* `fig1B.bash`
		* `ORF_length.R`
	* [**Figure 2**. Re-analysis of SARS-CoV-2 gene expression in publicly available ribosome profiling and mass spectrometry datasets](#figure-2).
		* `riboseq_sliding_window.R`
	* [**Figure 3**. SARS-CoV-2 protein sequence properties](#figure-3).
		* `generate_random_protein.py`
		* `tally_epitope_coverage.py`
	* [**Figure 4**. Amino acid variation in proteins encoded by genes overlapping *ORF3a* in viruses of the species *Severe acute respiratory syndrome-related coronavirus*](#figure-4).
	* [**Figure 5**. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels](#figure-5).
		* `generate_seqs_from_VCF.py`
	* [**Figure 6**. Between-taxa sliding window analysis of natural selection on overlapping frames of *ORF3a*](#figure-6).
	* [**Figure 7**. Pandemic spread of the EP+1 haplotype and the hitchhiking of *ORF3d*-LOF](#figure-7).
	* [**Figure 8**. High-frequency within-host mutations](#figure-8).
		* `filter_vcf.py`
		* `summarize_intrahost_by_site.py`
	* [**Additional scripts**](#additional-scripts).
		* `extract_seq_subset.py`
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)


## <a name="supplementary-data"></a>Supplementary data

The following supplementary data are provided in this GitHub repository:

1. `SARS-related-CoV_ALN.fasta`: the whole-genome multiple sequence alignment (between-taxa analysis) of *n*=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus*. See [manuscript](#citation) for details. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.
2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within `SARS-related-CoV_ALN.fasta`.
3. `Supplementary_Tables.xlsx`: Supplementary Tables referred to in the [manuscript](#citation).


## <a name="supplementary-scripts"></a>Supplementary scripts

Scripts are arranged by Figure, and therefore by analysis. The scripts are of two types: 

1. ***Command-line*** scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual analysis*** scripts are R scripts or Python Jupyter Notebooks documenting the bulk of our data analyses and visualizations. These are intended to be executed manually line-by-line in R/RStudio or Jupyter. The use should replace path names and arguments with the appropriate values for the user's analysis and directories.


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

		`XXX`
		

### <a name="figure-1"></a>Figure 1. Gene repertoire and evolutionary relationships of *Severe acute respiratory syndrome-related coronavirus* species members

* `fig1B.bash` (*command-line script*)
	* **Description**. Produces Figure1B using PyGenomeTracks.
	* **Requirements**. PyGenomeTracks, Seqkit.
	* **Input**. The following files need to be in the working directory: 
		1. `SARS-related-CoV_ALN.fasta`: multiple alignment FASTA file of n=21 genomes of the species *Severe acute respiratory syndrome-related coronavirus*. Note that the pangolin-CoV GD/1 sequence has been masked as `N`, because GISAID permission is required for data access.
		2. `SARS-related-CoV_ALN.gtf`: Gene Transfer Format (GTF) file giving gene positions within `SARS-related-CoV_ALN.fasta`.
		3. `sbc_rename2.nw`: Newick tree for sarbecovirus alignment **!!TODO: Zac, I need this <--**
		4. `parameters_input.txt` **!!TODO: Zac, I need this <--**
		5. `parameters_input2.txt` **!!TODO: Zac, I need this <--**
	* **Output**. 
		1. `Fig1b.png`
	* **Example**:

		`fig1B.bash`


* `ORF_length.R` (*manual analysis script*)
	* **Description**. Analyze the genome-wide ORF length results of the Schlub et al. (2018) codon permutation method and produce Figure 1—figure supplement 1.
	* **Requirements**. tidyverse, RColorBrewer, ggrepel, patchwork.
	* **Input**. 
		1. `frameshift_results.txt`, produced by `frameshift_analysis.bash`.
	* **Output**. 
		1. `ORF_length_heatmap_data.tsv`, p-values for each internal overlapping gene within a previously annotated gene, redundatly annotated on site-by-site basis.


### <a name="figure-2"></a>Figure 2. Re-analysis of SARS-CoV-2 gene expression in publicly available ribosome profiling and mass spectrometry datasets

* `riboseq_sliding_window.R` (*command-line script*)
	* **Description**. Sliding window script to calculate proportion of reads in each frame for a specified read length and window size, separately for each treatment combination.
	* **Input**. Unnamed arguments in the following order: 
		1. `INFILE_FRAMES`: file with condition metadata for samples
		2. `INFILE_READS`: file with read data
		3. `WIN_SIZE`: size of the sliding window
		4. `READ_LEN`: length of reads to consider.
	* **Output**. 
		1. Table in `.tsv` format giving the sum of reads in each frame for each condition and position (start of window).
	* **Example**:

		`Rscript riboseq_sliding_window.R frames_table.txt mapped_reads_by_readlength_ALL.tsv 30 30`


### <a name="figure-3"></a>Figure 3. SARS-CoV-2 protein sequence properties

* `generate_random_protein.py` (*command-line script*)
	* **Description**. Script to generate random protein sequences of a given length from an input proteome in FASTA format.
	* **Requirements**. Python packages Bio, Bio.Seq, Bio.SeqRecord, os, random, sys
	* **Input**. Unnamed arguments in the following order: 
		1. input file with peptide sequence in FASTA format
		2. length of peptide to generate
		3. number of peptides to generate
	* **Output**. 
		1. Multiple sequence alignment in FASTA format containing the randomized protein sequence(s).
	* **Example**:

		`generate_random_protein.py ORF3d_aa.fasta 57 1000`


* `tally_epitope_coverage.py` (*command-line script*)
	* **Description**. Script for tallying epitope coverage for one protein product in a sliding window.
	* **Requirements**. Python packages os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. input file in `.tsv` format; NetMHCPan output file with 5 columns in this order: ID, NB (number MHC alleles bound based on NetMHCPan output), product, codon_start, codon_end 
		2. length (integer) of linear peptides used in the epitope analysis (9 for MHC I with NetMHCPan; 15 for MHC II with NetMHCIIPan)
	* **Output**. 
		1. Table in `.tsv` format giving the sum of bound epitopes overlapping each site.
	* **Example**:

		`tally_epitope_coverage.py ORF3d_random.tsv 9`


### <a name="figure-5"></a>Figure 5. Natural selection analysis of viral nucleotide differences at three hierarchical evolutionary levels

* `generate_seqs_from_VCF.py` (*command-line script*)
	* **Description**. Script to generate a FASTA with randomly interspersed variants (from VCF), for use with OLGenie or any software that requires a multiple sequence alignment input and does not depend upon linkage.
	* **Requirements**. Python packages Bio, os, random, re, sys
	* **Input**. Three unnamed arguments in the following order: 
		1. A FASTA file containing exactly one (1) reference sequence (*e.g.*, SARS-CoV-2 Wuhan-Hu-1)
		2. a `.vcf` file containing the variants of interest
		3. number of sequences to generate (integer)	* **Output**. 
		1. A new `*.fasta` multiple sequence alignment file with variants randomy interpersed at the frequencies defined in the VCF file.
	* **Example**:

		`generate_seqs_from_VCF.py reference.fasta variants.vcf 1000`


### <a name="figure-8"></a>Figure 8. High-frequency within-host mutations

* `filter_VCF.py` (*command-line script*)
	* **Description**. Script to dynamically filter within-host variants using a binomial cutoff to control for a user-defined false-discovery rate (FDR). Automatically detects and analyzes all `.vcf` (variant call format) files in the working directory.
	* **Requirements**. Python packages Bio, numpy, os, random, re, scipy.stats, sys
	* **Input**. Five unnamed arguments in the following order: 
		1. analysis-wide FDR cutoff (integer): the maximum absolute number of false-positive variants called across all VCF files examined
		2. the minimum allowed minor allele frequency (numeric), independent of FDR
		3. genome length (integer)
		4. sequencing error rate per site (numeric): assumes all nucleotide changes are equally probable	
		5. number of samples in analysis (integer), i.e., total number of deeply sequenced within-host samples (here, the number of VCF files)
	* **Output**. 
		1. Creates a new (`*_filtered.vcf`) VCF file for each VCF file in the working directory, including only those variants that pass the FDR threshold.
		2. Prints summary statistics to STDOUT regarding the FREQUENCY, AC (allele count), and DP (total read depth at site) of passing and failing variants.
	* **Example**:

		`filter_VCF.py 1 0 29903 0.002 401`


```
* `summarize_intrahost_by_site.py` (*command-line script*)
	* **Description**.  Script for extracting a subset of sequences from a FASTA based on header ID.
	* **Requirements**. Python packages Bio, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. a text file containing exactly one column of sequence IDs (FASTA headers) 
		2. a multiple sequence alignment of nucleotides in FASTA format	* **Output**. 
		1. A multiple sequence alignment in FASTA format based on the file given by argument 2, but including only those sequences with headers provided in the file given by argument 1.
	* **Example**:

		`extract_seq_subset.py seq_ID_list.txt SARS-COV-2_ALN.fasta`
```


### <a name="additional-scripts"></a>Additional scripts

* `extract_seq_subset.py` (*command-line script*)
	* **Description**.  Script for extracting a subset of sequences from a FASTA based on header ID.
	* **Requirements**. Python packages Bio, os, sys
	* **Input**. Two unnamed arguments in the following order: 
		1. a text file containing exactly one column of sequence IDs (FASTA headers) 
		2. a multiple sequence alignment of nucleotides in FASTA format	* **Output**. 
		1. A multiple sequence alignment in FASTA format based on the file given by argument 2, but including only those sequences with headers provided in the file given by argument 1.
	* **Example**:

		`extract_seq_subset.py seq_ID_list.txt SARS-COV-2_ALN.fasta`


## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Postdoctoral Research Fellowship from Academia Sinica (to C.W.N. under P.I. Wen-Hsiung Li); funding from the Bavarian State Government and National Philanthropic Trust (to Z.A. under P.I. Siegfried Scherer); NSF IOS grants #1755370 and #1758800 (to S.-O.K.); and the University of Wisconsin-Madison John D. MacArthur Professorship Chair (to T.L.G). The authors thank the GISAID platform and the originating and submitting laboratories who kindly uploaded SARS-CoV-2 sequences to the GISAID EpiCov™ Database for public access (Supplement). The authors thank Maciej F. Boni, Reed A. Cartwright, John Flynn, Kyle Friend, Dan Graur, Robert S. Harbert, Cheryl Hayashi, David G. Karlin, Niloufar Kavian, Kin-Hang (Raven) Kok, Wen-Hsiung Li, Meiyeh Lu, David A. Matthews, Lisa Mirabello, Apurva Narechania, Felix Li Jin, and attendees of the UC Berkeley popgen journal club for useful information and discussion; Andrew E. Firth, Alexander Gorbalenya, Irwin Jungreis, Manolis Kellis, Raven Kok, Angelo Pavesi, Kei Sato, Manuela Sironi, and Noam Stern-Ginossar for an invaluable discussion regarding standardizing nomenclature; Helen Piontkivska, Patricia Wittkopp, Antonis Rokas, and one anonymous reviewer for critical suggestions; Ming-Hsueh Lin for immense feedback on figures; and special thanks to Priya Moorjani, Jacob Tennessen, Montgomery Slatkin, Yun S. Song, Jianzhi George Zhang, Xueying Li, Hongxiang Zheng, Qinqin Yu, Meredith Yeager, and Michael Dean for commenting on earlier drafts of the manuscript.


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
* Schlub TE, Buchmann JP, Holmes EC. 2018. <a target="_blank" href="https://academic.oup.com/mbe/article/35/10/2572/5067730">A simple method to detect candidate overlapping genes in viruses using single genome sequences</a>. *Molecular Biology and Evolution* **35**(10):2572–2581. doi: https://doi.org/10.1093/molbev/msy155
