# SARS-CoV-2 *ORF3d*

Supplementary scripts for Nelson et al. (2020) paper on SARS-CoV-2 *ORF3d*.


## <a name="contents"></a>Contents

* [Supplementary scripts](#supplementary-scripts)
	* [Figure 1](#figure-1)
	* [Figure 2](#figure-2)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)


## <a name="supplementary-scripts"></a>Supplementary scripts

Scripts are arranged by Figure, and therefore by analysis. The scripts are of two types: 

1. ***Command-line*** scripts are intended to be executed from the bash command line with the specified arguments. 

2. ***Manual execution*** scripts are R scripts or Python Jupyter Notebooks documenting the bulk of our data analyses and visualizations. These are intended to be executed manually line-by-line in R/RStudio or Jupyter. The use should replace path names and arguments with the appropriate values for the user's analysis and directories.


### <a name="figure-1"></a>Figure 1

1. `fig1B.bash` (*command-line script*)
	* **Decription**. Produces Figure1B using “PyGenomeTracks”.
	* **Input**. The following files need to be in the working directory: (1) multi fasta file Sarbecovirus_n21_ALN_FINAL_v2.fasta; (2) GTF file with gene positions, Sarbecovirus_n21_ALN_FINAL_v3.gtf; (3) Newick tree for sarbecovirus alignment, sbc_rename2.nw; (4) parameters file parameters_input.txt; and (5) parameter file parameters_input2.txt.
	* **Output**. Fig1b.png.
	* **Requirements**. PyGenomeTracks, Seqkit.
	* **Example**:

		`fig1B.bash`


### <a name="figure-2"></a>Figure 2

1. `SARSCOV2_ribo_seq_batch.R`. (*command-line script*)
	* **Description**. Sliding window script to calculate proportion of reads in each frame for a specified read length and window size, separately for each treatment combination.
	* **Input**. Unnamed arguments in this order: (1) `INFILE_FRAMES`, file with condition metadata for samples; (2) `INFILE_READS`, file with read data; (3) `WIN_SIZE`, size of the sliding window; and (4) `READ_LEN`, length of reads to consider.
	* **Output**. A `.TSV` table giving the sum of reads in each frame for each condition and position (start of window).
	* **Example**:

		`Rscript SARSCOV2_ribo_seq_batch.R frames_table.txt mapped_reads_by_readlength_ALL.tsv 29 28`


## <a name="acknowledgments"></a>Acknowledgments

This work was supported by a Postdoctoral Research Fellowship from Academia Sinica (to C.W.N. under P.I. Wen-Hsiung Li); funding from the Bavarian State Government and National Philanthropic Trust (to Z.A. under P.I. Siegfried Scherer); NSF IOS grants #1755370 and #1758800 (to S.-O.K.); and the University of Wisconsin-Madison John D. MacArthur Professorship Chair (to T.L.G). The authors thank the GISAID platform and the originating and submitting laboratories who kindly uploaded SARS-CoV-2 sequences to the GISAID EpiCov™ Database for public access (Supplement). The authors thank Maciej F. Boni, Reed A. Cartwright, John Flynn, Kyle Friend, Dan Graur, Robert S. Harbert, Cheryl Hayashi, David G. Karlin, Niloufar Kavian, Kin-Hang (Raven) Kok, Wen-Hsiung Li, Meiyeh Lu, David A. Matthews, Lisa Mirabello, Apurva Narechania, Felix Li Jin, and attendees of the UC Berkeley popgen journal club for useful information and discussion; Andrew E. Firth, Alexander Gorbalenya, Irwin Jungreis, Manolis Kellis, Raven Kok, Angelo Pavesi, Kei Sato, Manuela Sironi, and Noam Stern-Ginossar for an invaluable discussion regarding standardizing nomenclature; Helen Piontkivska, Patricia Wittkopp, Antonis Rokas, and one anonymous reviewer for critical suggestions; Ming-Hsueh Lin for immense feedback on figures; and special thanks to Priya Moorjani, Jacob Tennessen, Montgomery Slatkin, Yun S. Song, Jianzhi George Zhang, Xueying Li, Hongxiang Zheng, Qinqin Yu, Meredith Yeager, and Michael Dean for commenting on earlier drafts of the manuscript.


## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Ardern Z, Goldberg TL, Meng C, Kuo C-H, Ludwig C, Kolokotronis S-O, Wei X. 2020. Dynamically evolving novel overlapping gene as a factor in the SARS-CoV-2 pandemic. *eLife*, in review. bioRxiv doi: <a target="_blank" rel="noopener noreferrer" href="https://doi.org/10.1101/2020.05.21.109280">https://doi.org/10.1101/2020.05.21.109280</a>

and this page:

>https://github.com/chasewnelson/SARS-CoV-2-ORF3d


## <a name="contact"></a>Contact

If you have questions about our scripts or study, please click on the <a target="_blank" href="https://github.com/chasewnelson/SARS-CoV-2-ORF3d/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to the corresponding authors: 

*  Chase W. Nelson, cnelson <**AT**> gate <**DOT**> sinica <**DOT**> edu <**DOT**> tw
*  Zachary Ardern, zachary <**DOT**> ardern <**AT**> tum <**DOT**> de
*  Xinzhu (April) Wei, aprilwei <**AT**> berkeley <**DOT**> edu


## <a name="references"></a>References

* Nelson CW, Ardern Z, Wei X. 2020. <a target="_blank" href="https://academic.oup.com/mbe/article/37/8/2440/5815567">OLGenie: Estimating Natural Selection to Predict Functional Overlapping Genes</a>. *Molecular Biology and Evolution* **37**(8):2440–2449. doi: https://doi.org/10.1093/molbev/msaa087