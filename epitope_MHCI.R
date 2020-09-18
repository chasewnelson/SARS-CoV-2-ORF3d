### EPITOPE ANALYSIS of NetMHCpan results

library(tidyverse)
library(boot)
library(scales)
library(RColorBrewer)

# set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")


##############################
### MHC class I epitopes -- summarized all alleles by hand and added codon position in 20200413_NetMHCpan_WuhanHu1_ORFs_nt_aa.xlsx
(epitope_summary <- read_tsv("data/NetMHCpan_Wuhan-Hu-1.tsv"))
#29,075 x 9

##"Ave" is mean `1-log50k`. Looks.
## strong binder = SB = % Rank 0.5
## weak binder = WB = % Rank 2

# inventory genes
unique(epitope_summary$ID)

# RENAME after figuring out which ID is which gene, because names were truncated by NetMHCpan
epitope_summary[epitope_summary$ID == 'N_NOLa', ]$ID <- 'N_NOL2'
epitope_summary[epitope_summary$ID == 'N_NOLb', ]$ID <- 'N_NOL3'
epitope_summary[epitope_summary$ID == 'N_OL_a', ]$ID <- 'N_OL_ORF9b'
epitope_summary[epitope_summary$ID == 'N_OL_b', ]$ID <- 'N_OL_ORF9c'
epitope_summary[epitope_summary$ID == 'ORF1a_a', ]$ID <- 'ORF1a'
epitope_summary[epitope_summary$ID == 'ORF1a_b', ]$ID <- 'ORF1ab_1'
epitope_summary[epitope_summary$ID == 'ORF1a_c', ]$ID <- 'ORF1ab_1_NOL'
epitope_summary[epitope_summary$ID == 'ORF1a_d', ]$ID <- 'ORF1ab_2'
epitope_summary[epitope_summary$ID == 'ORF1a_e', ]$ID <- 'ORF1ab_2_NOL'
epitope_summary[epitope_summary$ID == 'ORF3a_a', ]$ID <- 'ORF3a'
epitope_summary[epitope_summary$ID == 'ORF3a_b', ]$ID <- 'ORF3a_NOL1'
epitope_summary[epitope_summary$ID == 'ORF3a_c', ]$ID <- 'ORF3a_NOL2'
epitope_summary[epitope_summary$ID == 'ORF3a_d', ]$ID <- 'ORF3a_OL_ORF3bp'
epitope_summary[epitope_summary$ID == 'ORF7a_a', ]$ID <- 'ORF7a'
epitope_summary[epitope_summary$ID == 'ORF7a_b', ]$ID <- 'ORF7a_NOL'
epitope_summary[epitope_summary$ID == 'ORF7b_a', ]$ID <- 'ORF7b'
epitope_summary[epitope_summary$ID == 'ORF7b_b', ]$ID <- 'ORF7b_NOL'
epitope_summary[epitope_summary$ID == 'ORF3b', ]$ID <- 'ORF3d' # previously used ORF3b as the name for ORF3d

# BOOKKEEPING: Because names were truncated, ambiguous
(epitope_summary_codon_counts <- epitope_summary %>%
    group_by(ID) %>%
    summarise(
      residues = n()
    ))


###############################
### ADD IN THE NEW CANDIDATE PEPTIDES
(epitope_data_newPeptides <- read_tsv("data/NetMHCpan_Wuhan-Hu-1_additional.tsv"))
# 1572 x 12

# Summarize 
(epitope_summary_newPeptides <- epitope_data_newPeptides %>%
    group_by(Pos, Peptide, ID, core) %>%
    summarise(
      Ave = mean(`1-log50k`),
      NB = mean(NB),
      SB = mean(SB),
      WB = mean(WB)
    )) # 131 x 10

# RENAME after figuring out which ID is which gene, because names were truncated by NetMHCpan
unique(epitope_summary_newPeptides$ID)

epitope_summary_newPeptides[epitope_summary_newPeptides$ID == "3d_2_25596-2569", ]$ID <- 'ORF3d2'
epitope_summary_newPeptides[epitope_summary_newPeptides$ID == "S_iO1_21744-218", ]$ID <- 'S-iORF1'
epitope_summary_newPeptides[epitope_summary_newPeptides$ID == "3c_25457-25582", ]$ID <- 'ORF3c'
epitope_summary_newPeptides[epitope_summary_newPeptides$ID == "3b_25814-25882", ]$ID <- 'ORF3b' # this is the first short ORF3b in SARS-CoV-2
epitope_summary_newPeptides[epitope_summary_newPeptides$ID == "S_iO2_21768-218", ]$ID <- 'S-iORF2'

# Examine
(epitope_summary_newPeptides_codon_counts <- epitope_summary_newPeptides %>%
    group_by(ID) %>%
    summarise(
      residues = n()
    ))


##############################
### JOIN THE OLD AND NEW PEPTIDE ANALYSES
epitope_summary <- dplyr::select(epitope_summary, Pos, Peptide, ID, Ave, NB, SB, WB)

# Reorder or remove superfluous columns
epitope_summary_newPeptides <- dplyr::select(epitope_summary_newPeptides, Pos, Peptide, ID, Ave, NB, SB, WB)

# add together
(epitope_summary <- rbind(as_tibble(epitope_summary), as.tibble(epitope_summary_newPeptides)))
#29,206 x 7

# number codons, convert to 1-based
epitope_summary$codon_start = epitope_summary$Pos + 1
epitope_summary$codon_end = epitope_summary$codon_start + 9 - 1

# reorder columns
epitope_summary <- dplyr::select(epitope_summary, ID, Pos, codon_start, codon_end, Peptide, Ave, NB, SB, WB)
epitope_summary # 29,206 x 9

### SAVE AND RELOAD
# SAVE for epitope tallying: ID/NB/product/codon_start/codon_end
epitope_summary_5col <- epitope_summary
epitope_summary_5col$product <- epitope_summary_5col$ID
epitope_summary_5col <- dplyr::select(epitope_summary_5col, ID, NB, product, codon_start, codon_end)
write_tsv(epitope_summary_5col, "data/MHCI_epitope_summary.tsv") # INCLUDING NEW PEPTIDES


##################################################
### SLIDING WINDOW in Python to calculate epitope COVERAGE of each site
##################################################

# DO: tally_epitope_coverage.py MHCI_epitope_summary.tsv 9 > MHCI_epitope_summary_tally.tsv 


####################
### RELOAD AFTER RUNNING PYTHON
(epitope_results <- read_tsv("data/MHCI_epitope_summary_tally.tsv", 
                             col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))) 
# 29,206 x 5

# rename column
names(epitope_results)[names(epitope_results) == "ID"] <- 'product'
unique(epitope_results$product)

# FILTER and RENAME desired products
epitope_gene_IDs <- c("S", "S-iORF1", "S-iORF2", 
                      "ORF3a", "ORF3c", "ORF3d", "ORF3d2", "ORF3b",
                      "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8p", "N", "ORF9b", "ORF9c", "ORF10") # ORF8p at first
epitope_gene_names <- c("S", "S-i1", "S-i2", 
                        "3a", "3c", "3d", "3d-2", "3b",
                        "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10")
(epitope_results <- filter(epitope_results, product %in% c(epitope_gene_IDs)))
# 2902 x 5

# FACTOR gene
epitope_results$product <- factor(epitope_results$product,
                                  levels = epitope_gene_IDs,
                                  labels = epitope_gene_names)

### SUMMARIZE COVERAGE FOR EACH GENE
(epitope_summary <- epitope_results %>%
    group_by(product) %>%
    summarise(
      nonamer_count = n(),
      mean_NB = mean(NB_coverage),
      sd_NB = sd(NB_coverage),
      SE_NB = sd_NB/ sqrt(nonamer_count),
      sum_NB = sum(NB_coverage)
    ))

(epitope_summary_hits <- filter(epitope_results, NB_coverage > 0) %>%
    group_by(product) %>%
    summarise(
      nonamer_epitope_count = n()
    ))

# join number of nonamers for each gene
(epitope_summary <- left_join(epitope_summary, epitope_summary_hits, by = "product"))
epitope_summary[is.na(epitope_summary$nonamer_epitope_count), ]$nonamer_epitope_count <- 0 # superlfuous, THIS time

# calculate proportion of the product hit
epitope_summary$prop_product_hit <- epitope_summary$nonamer_epitope_count / epitope_summary$nonamer_count
epitope_summary$NB_per_codon <- epitope_summary$sum_NB / epitope_summary$nonamer_count
epitope_summary

# PLOT PROPORTION BY GENE
ggplot(epitope_summary, mapping = aes(x = product, y = prop_product_hit)) +
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
# 3d-2 and N have the lowest proportion of residues covered by an eptiope


##################################################
### ADD TWO CLASSES OF NEGATIVE CONTROLS! 
##################################################

####################
# (1) NONFUNCTIONAL ORFS. result expected for real ORFs that have been evolving in the genome without functional constraint
#(nonfunctional_ORFs <- read_tsv("/Users/cwnelson88/Desktop/SARS-CoV-2/NetMHCpan/NetMHCpan_Epitope_Negative_Controls/SARSCOV2-ORFs-90-ATG-longest_unannotated_noStop_netMHCpan_parsed_cut.txt"))
#nonfunctional_ORFs <- filter(nonfunctional_ORFs, Allele == "HLA-A01:01") # because they contain redundant information
#nonfunctional_ORFs$codon_start <- nonfunctional_ORFs$Pos + 1
#nonfunctional_ORFs$codon_end <- nonfunctional_ORFs$Pos + 9
#nonfunctional_ORFs$NB_coverage <- 0

# did the below with Python instead later (tally_*.py)
# uncomment to repeat here, but LENGTHY <-- CHANGE THIS
#for (this_ID in unique(nonfunctional_ORFs$ID)) { 
#  #this_ID <- "s26342-26181"
#  cat(this_ID, "\n")
#  
#  for (this_codon_start in sort(unique(nonfunctional_ORFs[nonfunctional_ORFs$ID == this_ID, ]$codon_start))) {
#    cat('.')
#    this_codon_end <- nonfunctional_ORFs[nonfunctional_ORFs$ID == this_ID & nonfunctional_ORFs$codon_start == this_codon_start, ]$codon_end
#    this_NB <- nonfunctional_ORFs[nonfunctional_ORFs$ID == this_ID & nonfunctional_ORFs$codon_start == this_codon_start, ]$NB
#    
#    for (i in seq(this_codon_start, this_codon_end, 1)) {
#      nonfunctional_ORFs[nonfunctional_ORFs$ID == this_ID & nonfunctional_ORFs$codon_start == i, ]$NB_coverage <- 
#        nonfunctional_ORFs[nonfunctional_ORFs$ID == this_ID & nonfunctional_ORFs$codon_start == i, ]$NB_coverage + this_NB
#    }
#  }
#  cat("\n")
#}

## SAVE AND RELOAD
#write_tsv(nonfunctional_ORFs, "data/MHCI_short_unannot_ORFs.tsv")

## RELOAD IMPORTANT
nonfunctional_ORFs <- read_tsv("data/MHCI_short_unannot_ORFs.tsv")
nonfunctional_ORFs
names(nonfunctional_ORFs)[names(nonfunctional_ORFs) == "ID"] <- 'product'

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

# join nonamers hit and nonamer counts
(nonfunctional_ORFs_summary <- left_join(nonfunctional_ORFs_summary, nonfunctional_ORFs_summary_hits, by = "product"))
nonfunctional_ORFs_summary[is.na(nonfunctional_ORFs_summary$nonamer_nonfunctional_ORFs_count), ]$nonamer_nonfunctional_ORFs_count <- 0 # superfluous THIS time

# proportions
nonfunctional_ORFs_summary$prop_product_hit <- nonfunctional_ORFs_summary$nonamer_nonfunctional_ORFs_count / nonfunctional_ORFs_summary$nonamer_count
nonfunctional_ORFs_summary$NB_per_codon <- nonfunctional_ORFs_summary$sum_NB / nonfunctional_ORFs_summary$nonamer_count
nonfunctional_ORFs_summary

# Define confidence bounds -- 95% confidence interval
(nonfunctional_ORFs_epitopes_per9mer_mean <- mean(nonfunctional_ORFs_summary$NB_per_codon))
(nonfunctional_ORFs_epitopes_per9mer_SD <- sd(nonfunctional_ORFs_summary$NB_per_codon))
(nonfunctional_ORFs_epitopes_per9mer_SE <- sd(nonfunctional_ORFs_summary$NB_per_codon) / sqrt(nrow(nonfunctional_ORFs_summary)))
(nonfunctional_ORFs_epitopes_per9mer_p025 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.025))
(nonfunctional_ORFs_epitopes_per9mer_p975 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.975))


####################
## (2) RANDOM ORFs. result expected for entirely random ORFs that have the same amino acid content of each protein

# Because they are very large, we cut them down to size, e.g.,
# DO: cat NC_045512.2_S_aa_random_netMHCpan_parsed.txt | gsed -E 's/_rand[a-zA-Z0-9_]+//' | cut -f1,3,10 > S_random.tsv
# DO: convert to 1-based positions by adding 1
# DO: cut down to only one allele (the 12 are redundant)
# DO: save as, e.g., S_random.tsv


################################################################################
### TALLY each site's overlap using the Python script because R can't handle the damned loop and I'm too stupid to vectorize. For example:
# DO, for each protein:
#tally_epitope_coverage.py S_random.tsv 9 > S_random_tally.tsv

# LOAD AND PROCESS
S_random_tally <- read_tsv("data/NetMHCpan_S_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
S_random_tally$product <- "S"
ORF3a_random_tally <- read_tsv("data/NetMHCpan_ORF3a_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3a_random_tally$product <- "ORF3a"
ORF3c_random_tally <- read_tsv("data/NetMHCpan_ORF3c_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3c_random_tally$product <- "ORF3c"
ORF3d_random_tally <- read_tsv("data/NetMHCpan_ORF3d_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage")) # used to be called ORF3c
ORF3d_random_tally$product <- "ORF3d"
ORF3d2_random_tally <- read_tsv("data/NetMHCpan_ORF3d2_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3d2_random_tally$product <- "ORF3d-2"
ORF3b_random_tally <- read_tsv("data/NetMHCpan_ORF3b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF3b_random_tally$product <- "ORF3b"
E_random_tally <- read_tsv("data/NetMHCpan_E_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
E_random_tally$product <- "E"
M_random_tally <- read_tsv("data/NetMHCpan_M_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
M_random_tally$product <- "M"
ORF6_random_tally <- read_tsv("data/NetMHCpan_ORF6_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF6_random_tally$product <- "ORF6"
ORF7a_random_tally <- read_tsv("data/NetMHCpan_ORF7a_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF7a_random_tally$product <- "ORF7a"
ORF7b_random_tally <- read_tsv("data/NetMHCpan_ORF7b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF7b_random_tally$product <- "ORF7b"
ORF8_random_tally <- read_tsv("data/NetMHCpan_ORF8_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF8_random_tally$product <- "ORF8"
N_random_tally <- read_tsv("data/NetMHCpan_N_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
N_random_tally$product <- "N"
ORF9b_random_tally <- read_tsv("data/NetMHCpan_ORF9b_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF9b_random_tally$product <- "ORF9b"
ORF9c_random_tally <- read_tsv("data/NetMHCpan_ORF9c_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF9c_random_tally$product <- "ORF9c"
ORF10_random_tally <- read_tsv("data/NetMHCpan_ORF10_random_tally.tsv", col_names = c("ID", "codon_start", "codon_end", "NB", "NB_coverage"))
ORF10_random_tally$product <- "ORF10"


### Combine
ALL_ORFs_random <- rbind(S_random_tally, ORF3a_random_tally, ORF3c_random_tally, ORF3d_random_tally, ORF3d2_random_tally, ORF3b_random_tally,
                         E_random_tally, M_random_tally, ORF6_random_tally, ORF7a_random_tally, ORF7b_random_tally, ORF8_random_tally,
                         N_random_tally, ORF9b_random_tally, ORF9c_random_tally, ORF10_random_tally)
ALL_ORFs_random # 2,625,000 x 6
unique(ALL_ORFs_random$ID) #"s1"    "s2"    "s3"    "s4"    "s5"    "s6"    "s7"    "s8"    "s9" ...
unique(ALL_ORFs_random$product)


### SUMMARIZE BY REPLICATE (ID) and GENE
(ALL_ORFs_random_summary <- ALL_ORFs_random %>%
    group_by(product, ID) %>%
    summarise(
      nonamer_count = n(),
      mean_NB = mean(NB_coverage),
      sd_NB = sd(NB_coverage),
      SE_NB = sd_NB/ sqrt(nonamer_count),
      sum_NB = sum(NB_coverage)
    ))
#16,000 x 7

(ALL_ORFs_random_summary_hits <- filter(ALL_ORFs_random, NB_coverage > 0) %>%
    group_by(product, ID) %>%
    summarise(
      nonamer_ALL_ORFs_random_count = n()
    ))
# 15,993 x 3

# join nonamer hits and count
(ALL_ORFs_random_summary <- left_join(ALL_ORFs_random_summary, ALL_ORFs_random_summary_hits, by = c("product", "ID")))
ALL_ORFs_random_summary[is.na(ALL_ORFs_random_summary$nonamer_ALL_ORFs_random_count), ]$nonamer_ALL_ORFs_random_count <- 0 

ALL_ORFs_random_summary$prop_product_hit <- ALL_ORFs_random_summary$nonamer_ALL_ORFs_random_count / ALL_ORFs_random_summary$nonamer_count
ALL_ORFs_random_summary$NB_per_codon <- ALL_ORFs_random_summary$sum_NB / ALL_ORFs_random_summary$nonamer_count 
ALL_ORFs_random_summary 

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
random_gene_IDs <- c("ORF1ab",
                     "S", "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3b",
                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF9c", "ORF10")
random_gene_names <- c("1ab", 
                       "S", "3a", "3c", "3d", "3d-2", "3b",
                       "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10")

# Recode product names
ALL_ORFs_random_summary_byGene$product <- factor(ALL_ORFs_random_summary_byGene$product,
                                                 levels = random_gene_IDs,
                                                 labels = random_gene_names)
ALL_ORFs_random_summary_byGene



#######################################################
### JOIN negative control to test data!
epitope_summary <- left_join(epitope_summary, ALL_ORFs_random_summary_byGene, by = "product")
epitope_summary # 18 x 15


#######################################################
###  Pivot to LONG and transfer values to columns
(epitope_summary_LONG <- epitope_summary %>%
   pivot_longer(cols = c('NB_per_codon', 'random_NB_per_codon_mean'), names_to = "sample_type", values_to = "NB_per_codon")) 

# Define 95% CI for observed ORFs
epitope_summary_LONG$NB_per_codon_LOW <- epitope_summary_LONG$NB_per_codon - 1.96 * epitope_summary_LONG$SE_NB
epitope_summary_LONG$NB_per_codon_HIGH <- epitope_summary_LONG$NB_per_codon + 1.96 * epitope_summary_LONG$SE_NB

# Define 95% CI for randomized ORFs
epitope_summary_LONG[epitope_summary_LONG$sample_type == "random_NB_per_codon_mean", ]$NB_per_codon_LOW <- 
  epitope_summary_LONG[epitope_summary_LONG$sample_type == "random_NB_per_codon_mean", ]$random_NB_per_codon_Q1 # not Q1, actually 2.5%ile
epitope_summary_LONG[epitope_summary_LONG$sample_type == "random_NB_per_codon_mean", ]$NB_per_codon_HIGH <- 
  epitope_summary_LONG[epitope_summary_LONG$sample_type == "random_NB_per_codon_mean", ]$random_NB_per_codon_Q3 # not Q3, actually 97.5%ile

# FACTOR sample type
epitope_summary_LONG$sample_type <- factor(epitope_summary_LONG$sample_type,
                                           levels = c("NB_per_codon", "random_NB_per_codon_mean"),
                                           labels = c("Observed", "Expected (randomized peptide)"))
#View(epitope_summary_LONG)

# FACTOR gene
length(unique(epitope_summary_LONG$product)) # 18
epitope_summary_LONG$product <- factor(epitope_summary_LONG$product,
                                       levels = random_gene_names,
                                       labels = random_gene_names)

# subset for use
epitope_gene_names_subset <- setdiff(epitope_gene_names, c("3d-2"))

# short unannotated ORFs
(nonfunctional_ORFs_epitopes_per9mer_p025 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.025))
(nonfunctional_ORFs_epitopes_per9mer_p975 <- quantile(nonfunctional_ORFs_summary$NB_per_codon, 0.975))


### ADD TO PLOT
(epitopes_per9mer_byProduct_MHCI_PLOT <- ggplot(filter(epitope_summary_LONG, product %in% epitope_gene_names_subset),
                                                mapping = aes(x = product, y = NB_per_codon, fill = sample_type)) + 
    
    # nonfunctional ORF control
    geom_rect(xmin = -Inf, xmax = Inf,
              ymin = nonfunctional_ORFs_epitopes_per9mer_p025, ymax = nonfunctional_ORFs_epitopes_per9mer_p975,
              fill = "#E9E9E9", linetype = 0, inherit.aes = FALSE) +
    geom_hline(yintercept = nonfunctional_ORFs_epitopes_per9mer_mean, color = "grey", linetype = "dashed", size = 0.25) + # grey
    
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.8) +
    
    geom_errorbar(mapping = aes(ymin = NB_per_codon_LOW, ymax = NB_per_codon_HIGH), position = position_dodge(0.8), width = 0.2, size = 0.175) +
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
epitope_summary_test <- epitope_summary
epitope_summary_test <- dplyr::select(epitope_summary_test, product, NB_per_codon)

# Recode randomized product names
ALL_ORFs_random_summary$product <- factor(ALL_ORFs_random_summary$product,
                                          levels = random_gene_IDs,
                                          labels = random_gene_names)
ALL_ORFs_random_summary

epitope_summary_test$random_NB_gt <- 0
epitope_summary_test$random_NB_eq <- 0
epitope_summary_test$random_NB_lt <- 0

for (this_product in unique(epitope_summary_test$product)) {
  NB_per_codon_observed <- epitope_summary_test[epitope_summary_test$product == this_product, ]$NB_per_codon
  epitope_summary_test[epitope_summary_test$product == this_product, ]$random_NB_gt <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon > NB_per_codon_observed)
  epitope_summary_test[epitope_summary_test$product == this_product, ]$random_NB_eq <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon == NB_per_codon_observed)
  epitope_summary_test[epitope_summary_test$product == this_product, ]$random_NB_lt <- sum(ALL_ORFs_random_summary[ALL_ORFs_random_summary$product == this_product, ]$NB_per_codon < NB_per_codon_observed)
}

# one-sided ASL P-values
epitope_summary_test$ASL_Obs_lt_Exp_P <- 1
epitope_summary_test$ASL_Obs_gt_Exp_P <- 1
epitope_summary_test$ASL_Obs_lt_Exp_P <- epitope_summary_test$random_NB_lt / (epitope_summary_test$random_NB_lt + epitope_summary_test$random_NB_eq + epitope_summary_test$random_NB_gt)
epitope_summary_test$ASL_Obs_gt_Exp_P <- epitope_summary_test$random_NB_gt / (epitope_summary_test$random_NB_lt + epitope_summary_test$random_NB_eq + epitope_summary_test$random_NB_gt)


# two-sided ASL P-values
epitope_summary_test$ASL_Obs_ne_Exp_P <- 1
epitope_summary_test[epitope_summary_test$ASL_Obs_lt_Exp_P < epitope_summary_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_ne_Exp_P <- 
  2 * epitope_summary_test[epitope_summary_test$ASL_Obs_lt_Exp_P < epitope_summary_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_lt_Exp_P
epitope_summary_test[epitope_summary_test$ASL_Obs_lt_Exp_P > epitope_summary_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_ne_Exp_P <- 
  2 * epitope_summary_test[epitope_summary_test$ASL_Obs_lt_Exp_P > epitope_summary_test$ASL_Obs_gt_Exp_P, ]$ASL_Obs_gt_Exp_P
epitope_summary_test[epitope_summary_test$ASL_Obs_ne_Exp_P == 0, ]$ASL_Obs_ne_Exp_P <- 1 / NUM_REPLICATES


### Finally, test versus nonfunctionals
epitope_summary_test$nonfunctional_NB_gt <- 0
epitope_summary_test$nonfunctional_NB_eq <- 0
epitope_summary_test$nonfunctional_NB_lt <- 0

for (this_product in unique(epitope_summary_test$product)) {
  NB_per_codon_observed <- epitope_summary_test[epitope_summary_test$product == this_product, ]$NB_per_codon
  epitope_summary_test[epitope_summary_test$product == this_product, ]$nonfunctional_NB_gt <- sum(nonfunctional_ORFs_summary$NB_per_codon > NB_per_codon_observed)
  epitope_summary_test[epitope_summary_test$product == this_product, ]$nonfunctional_NB_eq <- sum(nonfunctional_ORFs_summary$NB_per_codon == NB_per_codon_observed)
  epitope_summary_test[epitope_summary_test$product == this_product, ]$nonfunctional_NB_lt <- sum(nonfunctional_ORFs_summary$NB_per_codon < NB_per_codon_observed)
}


# one-sided ASL_nonfunctional P-values
epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P <- 1
epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P <- 1
epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P <- epitope_summary_test$nonfunctional_NB_lt / (epitope_summary_test$nonfunctional_NB_lt + epitope_summary_test$nonfunctional_NB_eq + epitope_summary_test$nonfunctional_NB_gt)
epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P <- epitope_summary_test$nonfunctional_NB_gt / (epitope_summary_test$nonfunctional_NB_lt + epitope_summary_test$nonfunctional_NB_eq + epitope_summary_test$nonfunctional_NB_gt)


# two-sided ASL_nonfunctional P-values
epitope_summary_test$ASL_nonfunctional_Obs_ne_Exp_P <- 1
epitope_summary_test[epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P < epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 
  2 * epitope_summary_test[epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P < epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_lt_Exp_P
epitope_summary_test[epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P > epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 
  2 * epitope_summary_test[epitope_summary_test$ASL_nonfunctional_Obs_lt_Exp_P > epitope_summary_test$ASL_nonfunctional_Obs_gt_Exp_P, ]$ASL_nonfunctional_Obs_gt_Exp_P
epitope_summary_test[epitope_summary_test$ASL_nonfunctional_Obs_ne_Exp_P == 0, ]$ASL_nonfunctional_Obs_ne_Exp_P <- 1 / NUM_REPLICATES


### SAVE
write_tsv(epitope_summary_test, "data/MHCI_epitope_summary_test.tsv")


###  MHC class II in ANOTHER SCRIPT: epitope_MHCII.R




