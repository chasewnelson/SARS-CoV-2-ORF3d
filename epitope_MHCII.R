### EPITOPE ANALYSIS of NetMHCIIpan results

library(tidyverse)
library(boot)
library(scales)
library(RColorBrewer)


# set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

##############################
### MHC class II epitopes
(epitope_summary_MHCII <- read_tsv("data/NetMHCIIpan_Wuhan-Hu-1.txt"))
#248,820 x 13

## strong binder = SB = % Rank 0.5
## weak binder = WB = % Rank 2

# inventory genes
unique(epitope_summary_MHCII$ID)

# RENAME ORF3d the new name
epitope_summary_MHCII[epitope_summary_MHCII$ID == 'ORF3b', ]$ID <- 'ORF3d'

# CHANGE TO 0-based: MHCII, unlike MHCI, is 1-based!
epitope_summary_MHCII$Pos <- epitope_summary_MHCII$Pos - 1

# FILTER DOWN by randomly choosing one allele (redundant info)
(epitope_summary_MHCII <- filter(epitope_summary_MHCII, Target == "HLA-DPA10201-DPB10401"))

# BOOKKEEPING
(epitope_summary_MHCII_codon_counts <- epitope_summary_MHCII %>%
    group_by(ID) %>%
    summarise(
      residues = n()
    ))


###############################
### ADD IN THE NEW CANDIDATE PEPTIDES
(epitope_data_newPeptides <- read_tsv("data/NetMHCIIpan_Wuhan-Hu-1_additional.txt"))

# CHANGE TO 0-based: MHCII, unlike MHCI, is 1-based!
epitope_data_newPeptides$Pos <- epitope_data_newPeptides$Pos - 1

# Summarize 
(epitope_summary_MHCII_newPeptides <- epitope_data_newPeptides %>%
    group_by(Pos, Peptide, ID, Peptide) %>% # Peptide should be redundant, just making sure
    summarise(
      Ave = mean(Score), # for MHC I, called `1-log50k`
      NB = mean(NB),
      SB = mean(SB),
      WB = mean(WB)
    ))

# RENAME genes
unique(epitope_summary_MHCII_newPeptides$ID)
# "3d_2_25596-25697" "S_iO1_21744-21863" "3c_25457-25582" "3b_25814-25882" "S_iO2_21768-21863"

epitope_summary_MHCII_newPeptides[epitope_summary_MHCII_newPeptides$ID == "3d_2_25596-25697", ]$ID <- 'ORF3d2'
epitope_summary_MHCII_newPeptides[epitope_summary_MHCII_newPeptides$ID == "S_iO1_21744-21863", ]$ID <- 'S-iORF1'
epitope_summary_MHCII_newPeptides[epitope_summary_MHCII_newPeptides$ID == "3c_25457-25582", ]$ID <- 'ORF3c'
epitope_summary_MHCII_newPeptides[epitope_summary_MHCII_newPeptides$ID == "3b_25814-25882", ]$ID <- 'ORF3b' # this is the first short ORF3b in SARS-CoV-2
epitope_summary_MHCII_newPeptides[epitope_summary_MHCII_newPeptides$ID == "S_iO2_21768-21863", ]$ID <- 'S-iORF2'

# Examine
(epitope_summary_MHCII_newPeptides_codon_counts <- epitope_summary_MHCII_newPeptides %>%
    group_by(ID) %>%
    summarise(
      residues = n()
    ))


##############################
### COMBINE THE OLD AND NEW ANALYSES

# Reorder or remove superfluous columns
epitope_summary_MHCII <- dplyr::select(epitope_summary_MHCII, Pos, Peptide, ID, NB, SB, WB) 
epitope_summary_MHCII_newPeptides <- dplyr::select(epitope_summary_MHCII_newPeptides, Pos, Peptide, ID, NB, SB, WB) 

# add together
(epitope_summary_MHCII <- rbind(as.tibble(epitope_summary_MHCII), as.tibble(epitope_summary_MHCII_newPeptides)))

# number codons
epitope_summary_MHCII$codon_start = epitope_summary_MHCII$Pos + 1
epitope_summary_MHCII$codon_end = epitope_summary_MHCII$codon_start + 15 - 1

# reorder columns
epitope_summary_MHCII <- dplyr::select(epitope_summary_MHCII, ID, Pos, codon_start, codon_end, Peptide, Ave, NB, SB, WB)
epitope_summary_MHCII

# SAVE a version with 5 columns for epitope tallying: ID/NB/product/codon_start/codon_end
epitope_summary_MHCII_5col <- epitope_summary_MHCII
epitope_summary_MHCII_5col$product <- epitope_summary_MHCII_5col$ID
epitope_summary_MHCII_5col <- dplyr::select(epitope_summary_MHCII_5col, ID, NB, product, codon_start, codon_end)

### SAVE <-- CHANGE THIS
#write_tsv(epitope_summary_MHCII, "data/MHCII_epitope_summary.tsv") # INCLUDING NEW PEPTIDES


##################################################
### SLIDING WINDOW in Python to calculate epitope COVERAGE of each site
##################################################

# DO: tally_epitope_coverage.py MHCII_epitope_summary.tsv 15 > MHCII_epitope_summary_tally.tsv


####################
### RELOAD AFTER RUNNING PYTHON
(epitope_results <- read_tsv("data/MHCII_epitope_summary_tally.tsv", 
                             col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))) 
# 9,671 x 5

# rename column
names(epitope_results)[names(epitope_results) == "ID"] <- 'product'
unique(epitope_results$product)

# FILTER and RENAME desired products
epitope_gene_IDs <- c("S", "S-iORF1", "S-iORF2", 
                      "ORF3a", "ORF3c", "ORF3d", "ORF3d2", "ORF3b",
                      "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF9c", "ORF10") # not renamed ORF8p as ORF8 
epitope_gene_names <- c("S", "S-i1", "S-i2", 
                        "3a", "3c", "3d", "3d-2", "3b",
                        "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10")
(epitope_results <- filter(epitope_results, product %in% c(epitope_gene_IDs)))
# 2,794 x 5

# FACTOR gene
epitope_results$product <- factor(epitope_results$product,
                                  levels = epitope_gene_IDs,
                                  labels = epitope_gene_names)

### SUMMARIZE COVERAGE FOR EACH GENE
(epitope_summary_MHCII <- epitope_results %>%
    group_by(product) %>%
    summarise(
      nonamer_count = n(),
      mean_NB = mean(NB_coverage),
      sd_NB = sd(NB_coverage),
      SE_NB = sd_NB/ sqrt(nonamer_count),
      sum_NB = sum(NB_coverage)
    ))

(epitope_summary_MHCII_hits <- filter(epitope_results, NB_coverage > 0) %>%
    group_by(product) %>%
    summarise(
      nonamer_epitope_count = n()
    )) # for each product, the number of residues hitting at least one epitope

# join number of nonamers for each gene
(epitope_summary_MHCII <- left_join(epitope_summary_MHCII, epitope_summary_MHCII_hits, by = "product"))
epitope_summary_MHCII[is.na(epitope_summary_MHCII$nonamer_epitope_count), ]$nonamer_epitope_count <- 0

# calculate proportion of the product hit
epitope_summary_MHCII$prop_product_hit <- epitope_summary_MHCII$nonamer_epitope_count / epitope_summary_MHCII$nonamer_count
epitope_summary_MHCII$NB_per_codon <- epitope_summary_MHCII$sum_NB / epitope_summary_MHCII$nonamer_count 

### PLOT PROPORTION BY GENE
ggplot(epitope_summary_MHCII, mapping = aes(x = product, y = prop_product_hit)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  ylab("Predicted epitope coverage") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10),
        panel.border = element_rect(),
        strip.background = element_blank()) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1))
# 3d-2 is 0!


##################################################
### ADD TWO CLASSES OF NEGATIVE CONTROLS!
##################################################

####################
# (1) NONFUNCTIONAL ORFS. result expected for real ORFs that have been evolving in the genome without functional constraint
(nonfunctional_ORFs <- read_tsv("data/MHCII_short_unannot_ORFs.tsv"))
#84,890 x 13
nonfunctional_ORFs <- dplyr::select(nonfunctional_ORFs, Pos, ID, Target, NB)
names(nonfunctional_ORFs)[names(nonfunctional_ORFs) == "Target"] <- "Allele"
nonfunctional_ORFs$Pos <- nonfunctional_ORFs$Pos - 1 # convert to 0-based for compatibility with MHC I approach
(nonfunctional_ORFs <- filter(nonfunctional_ORFs, Allele == "HLA-DPA10201-DPB10401"))
nonfunctional_ORFs$codon_start <- nonfunctional_ORFs$Pos + 1
nonfunctional_ORFs$codon_end <- nonfunctional_ORFs$codon_start + 15 - 1

# SAVE a version with 5 columns for epitope tallying: ID/NB/product/codon_start/codon_end
nonfunctional_ORFs_5col <- nonfunctional_ORFs
nonfunctional_ORFs_5col$product <- nonfunctional_ORFs_5col$ID
nonfunctional_ORFs_5col <- dplyr::select(nonfunctional_ORFs_5col, ID, NB, product, codon_start, codon_end)
nonfunctional_ORFs_5col
#write_tsv(nonfunctional_ORFs_5col, "data/MHCII_short_unannot_ORFs_stripped.tsv") # INCLUDING NEW PEPTIDES


##################################################
### SLIDING WINDOW in Python to calculate epitope COVERAGE of each site
##################################################

# IN: /Users/cwnelson88/Desktop/SARS-CoV-2/epitopes_MHCII/result_mhc_ii/
# DO: tally_epitope_coverage.py MHCII_short_unannot_ORFs_stripped.tsv 15 > MHCII_short_unannot_ORFs_tally.tsv


##################################################
## RELOAD TALLY after Python
(nonfunctional_ORFs <- read_tsv("data/MHCII_short_unannot_ORFs_tally.tsv",
                                col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage")))
# 3,265 x 5
names(nonfunctional_ORFs)[names(nonfunctional_ORFs) == "ID"] <- 'product'
length(unique(nonfunctional_ORFs$product)) # n=103
unique(nonfunctional_ORFs$product)
#[1] "s26342-26181" "s21038-20844" "s19856-19713" "s18185-17979" "s17543-17445" "s28081-28191" "s16295-16200" "s14027-13905" "s13550-13359" "s10856-10755"
#...

### NB
(nonfunctional_ORFs_summary <- nonfunctional_ORFs %>%
    group_by(product) %>%
    summarise(
      nonamer_count = n(),
      mean_NB = mean(NB_coverage),
      sd_NB = sd(NB_coverage),
      SE_NB = sd_NB/ sqrt(nonamer_count),
      sum_NB = sum(NB_coverage)
    ))

(nonfunctional_ORFs_summary_hits <- filter(nonfunctional_ORFs, NB_coverage > 0) %>%
    group_by(product) %>%
    summarise(
      nonamer_nonfunctional_ORFs_count = n()
    ))

# join
(nonfunctional_ORFs_summary <- left_join(nonfunctional_ORFs_summary, nonfunctional_ORFs_summary_hits, by = "product"))
nonfunctional_ORFs_summary[is.na(nonfunctional_ORFs_summary$nonamer_nonfunctional_ORFs_count), ]$nonamer_nonfunctional_ORFs_count <- 0

# prop hit by epitope
nonfunctional_ORFs_summary$prop_product_hit <- nonfunctional_ORFs_summary$nonamer_nonfunctional_ORFs_count / nonfunctional_ORFs_summary$nonamer_count
nonfunctional_ORFs_summary$NB_per_codon <- nonfunctional_ORFs_summary$sum_NB / nonfunctional_ORFs_summary$nonamer_count
nonfunctional_ORFs_summary # 103 x 9

# Define confidence bounds
(nonfunctional_ORFs_epitopes_per9mer_mean <- mean(nonfunctional_ORFs_summary$NB_per_codon))
(nonfunctional_ORFs_epitopes_per9mer_SD <- sd(nonfunctional_ORFs_summary$NB_per_codon))
(nonfunctional_ORFs_epitopes_per9mer_SE <- sd(nonfunctional_ORFs_summary$NB_per_codon) / sqrt(nrow(nonfunctional_ORFs_summary)))

# 95% CI
(nonfunctional_ORFs_epitopes_per9mer_p025 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.025))
(nonfunctional_ORFs_epitopes_per9mer_p975 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.975))


####################
## (2) RANDOM ORFS. result expected for entirely random ORFs that have the same amino acid content of each protein

# Because they are very large, we cut them down to size, e.g.,
#DO: cat NC_045512.2_E_aa_random_l75n1000_26_allele_netMHCIIpan_parsed.txt | gsed -E 's/_rand[a-zA-Z0-9_]+//' | cut -f1,3,11 > NC_045512.2_E_aa_random_l75n1000_26_allele_netMHCIIpan_parsed_cut.txt
#DO: cut down to only one allele (the 26 are redundant)
#DO: save as, e.g., E_random.tsv


################################################################################
### TALLY each site's overlap using the Python script AGAIN, because R can't handle the damned loop. For example:

# DO: tally_epitope_coverage.py E_random.tsv 15 > NetMHCIIpan_E_random_tally.tsv

S_random_tally <- read_tsv("data/NetMHCIIpan_S_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
S_random_tally$product <- "S"
ORF3a_random_tally <- read_tsv("data/NetMHCIIpan_ORF3a_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3a_random_tally$product <- "ORF3a"
ORF3c_random_tally <- read_tsv("data/NetMHCIIpan_ORF3c_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3c_random_tally$product <- "ORF3c"
ORF3d_random_tally <- read_tsv("data/NetMHCIIpan_ORF3d_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage")) # used to be called ORF3c
ORF3d_random_tally$product <- "ORF3d"
ORF3d2_random_tally <- read_tsv("data/NetMHCIIpan_ORF3d2_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3d2_random_tally$product <- "ORF3d-2"
ORF3b_random_tally <- read_tsv("data/NetMHCIIpan_ORF3b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3b_random_tally$product <- "ORF3b"
E_random_tally <- read_tsv("data/NetMHCIIpan_E_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
E_random_tally$product <- "E"
M_random_tally <- read_tsv("data/NetMHCIIpan_M_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
M_random_tally$product <- "M"
ORF6_random_tally <- read_tsv("data/NetMHCIIpan_ORF6_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF6_random_tally$product <- "ORF6"
ORF7a_random_tally <- read_tsv("data/NetMHCIIpan_ORF7a_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF7a_random_tally$product <- "ORF7a"
ORF7b_random_tally <- read_tsv("data/NetMHCIIpan_ORF7b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF7b_random_tally$product <- "ORF7b"
ORF8_random_tally <- read_tsv("data/NetMHCIIpan_ORF8_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF8_random_tally$product <- "ORF8"
N_random_tally <- read_tsv("data/NetMHCIIpan_N_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
N_random_tally$product <- "N"
ORF9b_random_tally <- read_tsv("data/NetMHCIIpan_ORF9b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF9b_random_tally$product <- "ORF9b"
ORF9c_random_tally <- read_tsv("data/NetMHCIIpan_ORF9c_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF9c_random_tally$product <- "ORF9c"
ORF10_random_tally <- read_tsv("data/NetMHCIIpan_ORF10_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF10_random_tally$product <- "ORF10"


### Combine
ALL_ORFs_random <- rbind(S_random_tally, ORF3a_random_tally, ORF3c_random_tally, ORF3d_random_tally, ORF3d2_random_tally, ORF3b_random_tally,
                         E_random_tally, M_random_tally, ORF6_random_tally, ORF7a_random_tally, ORF7b_random_tally, ORF8_random_tally,
                         N_random_tally, ORF9b_random_tally, ORF9c_random_tally, ORF10_random_tally)
ALL_ORFs_random # 2,747,000 x 6
unique(ALL_ORFs_random$ID) #"s1"    "s2"    "s3"    "s4"    "s5"    "s6"    "s7"    "s8"    "s9" ...

### SUMMARIZE BY REPLICATE (ID) and GENE
ALL_ORFs_random

(ALL_ORFs_random_summary <- ALL_ORFs_random %>%
    group_by(product, ID) %>%
    summarise(
      nonamer_count = n(),
      mean_NB = mean(NB_coverage),
      sd_NB = sd(NB_coverage),
      SE_NB = sd_NB/ sqrt(nonamer_count),
      sum_NB = sum(NB_coverage)
    ))

(ALL_ORFs_random_summary_hits <- filter(ALL_ORFs_random, NB_coverage > 0) %>%
    group_by(product, ID) %>%
    summarise(
      nonamer_ALL_ORFs_random_count = n()
    ))

# join
(ALL_ORFs_random_summary <- left_join(ALL_ORFs_random_summary, ALL_ORFs_random_summary_hits, by = c("product", "ID")))
ALL_ORFs_random_summary[is.na(ALL_ORFs_random_summary$nonamer_ALL_ORFs_random_count), ]$nonamer_ALL_ORFs_random_count <- 0
# 16,000 x 8

ALL_ORFs_random_summary$prop_product_hit <- ALL_ORFs_random_summary$nonamer_ALL_ORFs_random_count / ALL_ORFs_random_summary$nonamer_count
ALL_ORFs_random_summary$NB_per_codon <- ALL_ORFs_random_summary$sum_NB / ALL_ORFs_random_summary$nonamer_count

# Summarize by gene
(ALL_ORFs_random_summary_byGene <- ALL_ORFs_random_summary %>%
    group_by(product) %>%
    summarise(
      random_replicate_count = n(),
      random_NB_per_codon_mean = mean(NB_per_codon),
      random_NB_per_codon_sd = sd(NB_per_codon),
      random_NB_per_codon_SE = random_NB_per_codon_sd / sqrt(random_replicate_count),
      random_NB_per_codon_Q1 = quantile(NB_per_codon, 0.025),
      random_NB_per_codon_Q3 = quantile(NB_per_codon, 0.975)
    ))


# FACTOR genes to short form to match names
random_gene_IDs <- c("ORF1ab","S", "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3b",
                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF9c", "ORF10")
random_gene_names <- c("1ab", "S", "3a", "3c", "3d", "3d-2", "3b",
                       "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10")

# Recode product names
ALL_ORFs_random_summary_byGene$product <- factor(ALL_ORFs_random_summary_byGene$product,
                                                 levels = random_gene_IDs,
                                                 labels = random_gene_names)
ALL_ORFs_random_summary_byGene
#View(ALL_ORFs_random_summary_byGene)



#######################################################
### JOIN negative control to test data! (ORF1ab may be lost)
epitope_summary_MHCII <- left_join(epitope_summary_MHCII, ALL_ORFs_random_summary_byGene, by = "product")
epitope_summary_MHCII # 18 x 15


#######################################################
###  Pivot to LONG and transfer values to columns
(epitope_summary_MHCII_LONG <- epitope_summary_MHCII %>%
   pivot_longer(cols = c('NB_per_codon', 'random_NB_per_codon_mean'), names_to = "sample_type", values_to = "NB_per_codon")) 

# Define 95% CI for observed ORFs
epitope_summary_MHCII_LONG$NB_per_codon_LOW <- epitope_summary_MHCII_LONG$NB_per_codon - 1.96 * epitope_summary_MHCII_LONG$SE_NB
epitope_summary_MHCII_LONG$NB_per_codon_HIGH <- epitope_summary_MHCII_LONG$NB_per_codon + 1.96 * epitope_summary_MHCII_LONG$SE_NB

# Define 95% CI for randomized ORFs
epitope_summary_MHCII_LONG[epitope_summary_MHCII_LONG$sample_type == "random_NB_per_codon_mean", ]$NB_per_codon_LOW <- 
  epitope_summary_MHCII_LONG[epitope_summary_MHCII_LONG$sample_type == "random_NB_per_codon_mean", ]$random_NB_per_codon_Q1 # not Q1, but 2.5%ile
epitope_summary_MHCII_LONG[epitope_summary_MHCII_LONG$sample_type == "random_NB_per_codon_mean", ]$NB_per_codon_HIGH <- 
  epitope_summary_MHCII_LONG[epitope_summary_MHCII_LONG$sample_type == "random_NB_per_codon_mean", ]$random_NB_per_codon_Q3 # not Q3, but 97.5%ile

# FACTOR sample type
epitope_summary_MHCII_LONG$sample_type <- factor(epitope_summary_MHCII_LONG$sample_type,
                                                 levels = c("NB_per_codon", "random_NB_per_codon_mean"),
                                                 labels = c("Observed", "Expected (randomized peptide)"))

# FACTOR gene
length(unique(epitope_summary_MHCII_LONG$product))
unique(epitope_summary_MHCII_LONG$product) 
epitope_summary_MHCII_LONG$product <- factor(epitope_summary_MHCII_LONG$product,
                                             levels = random_gene_names,
                                             labels = random_gene_names)

# subset for use
epitope_gene_names_subset <- setdiff(epitope_gene_names, c("3d-2"))

# 95% CI for short unannotated ORFs
(nonfunctional_ORFs_epitopes_per9mer_p025 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.025))
(nonfunctional_ORFs_epitopes_per9mer_p975 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.975))

### ADD TO PLOT
(epitopes_per9mer_byProduct_MHCII_PLOT <- ggplot(filter(epitope_summary_MHCII_LONG, product %in% epitope_gene_names_subset),
                                                 mapping = aes(x = product, y = NB_per_codon, fill = sample_type)) +
    
    # nonfunctional ORF control
    geom_rect(xmin = -Inf, xmax = Inf, 
              ymin = nonfunctional_ORFs_epitopes_per9mer_p025, ymax = nonfunctional_ORFs_epitopes_per9mer_p975, 
              fill = "#E9E9E9", linetype = 0, inherit.aes = FALSE) + 
    geom_hline(yintercept = nonfunctional_ORFs_epitopes_per9mer_mean, color = "grey", linetype = "dashed", size = 0.25) + 
    
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8) +
    geom_errorbar(mapping = aes(ymin = NB_per_codon_LOW, ymax = NB_per_codon_HIGH), position = position_dodge(0.8), width = 0.175, size = 0.175) +
    xlab("Protein") +
    ylab("Epitopes per residue") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          axis.text.x = element_text(size = rel(0.85)), 
          axis.text.y = element_text(size = rel(0.85)),
          axis.title.x = element_text(size = rel(0.85)),
          axis.title.y = element_text(size = rel(0.85)),
          panel.border = element_rect(),
          strip.background = element_blank()) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("#0072B2", "gray"))
)


## Quantify which are significant
NUM_REPLICATES <- 1000
epitope_summary_MHCII_test <- epitope_summary_MHCII
epitope_summary_MHCII_test <- dplyr::select(epitope_summary_MHCII_test, product, NB_per_codon)

# Recode randomized product names
ALL_ORFs_random_summary$product <- factor(ALL_ORFs_random_summary$product,
                                          levels = random_gene_IDs,
                                          labels = random_gene_names)

epitope_summary_MHCII_test$random_NB_gt <- 0
epitope_summary_MHCII_test$random_NB_eq <- 0
epitope_summary_MHCII_test$random_NB_lt <- 0

for (this_product in unique(epitope_summary_MHCII_test$product)) {
  NB_per_codon_observed <- epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$NB_per_codon
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$random_NB_gt <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon > NB_per_codon_observed)
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$random_NB_eq <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon == NB_per_codon_observed)
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$random_NB_lt <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon < NB_per_codon_observed)
}

# one-sided ASL P-values
epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P <- 1
epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P <- 1
epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P <- epitope_summary_MHCII_test$random_NB_lt / (epitope_summary_MHCII_test$random_NB_lt + epitope_summary_MHCII_test$random_NB_eq + epitope_summary_MHCII_test$random_NB_gt)
epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P <- epitope_summary_MHCII_test$random_NB_gt / (epitope_summary_MHCII_test$random_NB_lt + epitope_summary_MHCII_test$random_NB_eq + epitope_summary_MHCII_test$random_NB_gt)


# two-sided ASL P-values
epitope_summary_MHCII_test$ASL_Obs_ne_Exp_P <- 1
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P < epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_ne_Exp_P <- 
  2 * epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P < epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_lt_Exp_P
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P > epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_ne_Exp_P <- 
  2 * epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_Obs_lt_Exp_P > epitope_summary_MHCII_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_gt_Exp_P
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_Obs_ne_Exp_P == 0, ]$ASL_Obs_ne_Exp_P <- 1 / NUM_REPLICATES

### Finally, test versus nonfunctionals
epitope_summary_MHCII_test$nonfunctional_NB_gt <- 0
epitope_summary_MHCII_test$nonfunctional_NB_eq <- 0
epitope_summary_MHCII_test$nonfunctional_NB_lt <- 0

for (this_product in unique(epitope_summary_MHCII_test$product)) {
  NB_per_codon_observed <- epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$NB_per_codon
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$nonfunctional_NB_gt <- sum(nonfunctional_ORFs_summary$NB_per_codon > NB_per_codon_observed)
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$nonfunctional_NB_eq <- sum(nonfunctional_ORFs_summary$NB_per_codon == NB_per_codon_observed)
  epitope_summary_MHCII_test[epitope_summary_MHCII_test$product == this_product, ]$nonfunctional_NB_lt <- sum(nonfunctional_ORFs_summary$NB_per_codon < NB_per_codon_observed)
}


# one-sided ASL_nonfunctional P-values
epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P <- 1
epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P <- 1
epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P <- epitope_summary_MHCII_test$nonfunctional_NB_lt / (epitope_summary_MHCII_test$nonfunctional_NB_lt + epitope_summary_MHCII_test$nonfunctional_NB_eq + epitope_summary_MHCII_test$nonfunctional_NB_gt)
epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P <- epitope_summary_MHCII_test$nonfunctional_NB_gt / (epitope_summary_MHCII_test$nonfunctional_NB_lt + epitope_summary_MHCII_test$nonfunctional_NB_eq + epitope_summary_MHCII_test$nonfunctional_NB_gt)


# two-sided ASL_nonfunctional P-values
epitope_summary_MHCII_test$ASL_nonfunctional_Obs_ne_Exp_P <- 1
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P < epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 
  2 * epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P < epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_lt_Exp_P
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P > epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 
  2 * epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_nonfunctional_Obs_lt_Exp_P > epitope_summary_MHCII_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_gt_Exp_P
epitope_summary_MHCII_test[epitope_summary_MHCII_test$ASL_nonfunctional_Obs_ne_Exp_P == 0, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 1 / NUM_REPLICATES


### SAVE
write_tsv(epitope_summary_MHCII_test, "data/MHCII_epitope_summary_test.tsv")



############################################################################################################
# NOW SIMPLY COMBINE MHCI and MHCII


