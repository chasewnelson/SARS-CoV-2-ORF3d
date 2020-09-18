###############################################################################
### SARSCOV2: between-taxa, between-host, and within-host diversity (THREE LEVELS)

library(tidyverse)
library(boot)
library(RColorBrewer)
library(scales)


# Set working directory 
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")


####################################################################################################
####################################################################################################
### BETWEEN-HOST (GISAID), called interhost herein
####################################################################################################
####################################################################################################

####################################################################################################
# IMPORT **SNPGenie** after manually classifying every codon as NOL, OL, or dOL (double overlapping; ORF3c/ORF3d region)
(interhost_results <- read_tsv("data/between-host_SNPGenie.tsv"))
## 9,985 x 10

# Rename OLG categories
interhost_results[interhost_results$OL_category == 'NOL', ]$OL_category <- "Non-OLG"
interhost_results[interhost_results$OL_category == 'OL', ]$OL_category <- "OLG"
interhost_results[interhost_results$OL_category == 'dOL', ]$OL_category <- "Double OLG"

# Another column for double checking
interhost_results$OL_category2 <- "Non-OLG"

# Save names to remember which columns to keep later
interhost_results_names <- names(interhost_results)


####################################################################################################
# IMPORT OLGenie (all genes and frames) & COMBINE
(interhost_results_OLGenie <- read_tsv("data/between-host_OLGenie.tsv"))

# Process
interhost_results_OLGenie$product <- str_replace(string = interhost_results_OLGenie$product, pattern = "../SARSCOV2_MAFFT_processed_", replacement = "") # for GISAID
interhost_results_OLGenie$product <- str_replace(string = interhost_results_OLGenie$product, pattern = ".fasta", replacement = "")

# Limit to analyzed gene regions, using reference frame perspective for OLGs
products_to_analyze <- c("S", "ORF3a", "N") # we won't consider the S-iORF to be legit, so won't exclude from S or show in figure
(interhost_results_OLGenie <- filter(interhost_results_OLGenie, product %in% products_to_analyze))

### Simply JOIN by product/codon, add coluns at the end; replace SNPGenie values with OLGenie values at the proper positions; delete superfluous columns
### Genes with ss13 overlap ###
(interhost_results <- left_join(interhost_results, filter(interhost_results_OLGenie, frame == 'ss13'), by = c('product', 'codon')))
# 9,985 x 24

# Reference gene with OLG in ss13: ORF3a, codons 22-43 (ORF3c before double overlap)
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$N_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$NN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$S_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$SN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$N_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$NN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$S_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$SN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(22, 43), ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(22, 43), ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(22, 43), ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(22, 43), ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 22:43, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: ORF3a/ORF3c (ORF3a codons 22-43, trim to 23-42)
interhost_results_ORF3c <- filter(interhost_results, product == "ORF3a", codon %in% 23:42)
interhost_results_ORF3c$product <- "ORF3c"
interhost_results_ORF3c$N_sites <- interhost_results_ORF3c$NN_sites
interhost_results_ORF3c$S_sites <- interhost_results_ORF3c$NS_sites # NS; ORF3a perspective
interhost_results_ORF3c$N_diffs <- interhost_results_ORF3c$NN_diffs
interhost_results_ORF3c$S_diffs <- interhost_results_ORF3c$NS_diffs # NS; ORF3a perspective
interhost_results_ORF3c$codon <- interhost_results_ORF3c$codon - min(interhost_results_ORF3c$codon) + 2 # starts at 2
interhost_results_ORF3c$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF3c"), interhost_results_ORF3c)
rm(interhost_results_ORF3c)

# Reference gene with OLG in ss13: ORF3a, codons 141-164 (ORF3b, first ORF)
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$N_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$NN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$S_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$SN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$N_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$NN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$S_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$SN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(141, 164), ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(141, 164), ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(141, 164), ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(141, 164), ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 141:164, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: ORF3a/ORF3b (ORF3a codons 141-164, trim to 142-163)
interhost_results_ORF3b <- filter(interhost_results, product == "ORF3a", codon %in% 142:163)
interhost_results_ORF3b$product <- "ORF3b"
interhost_results_ORF3b$N_sites <- interhost_results_ORF3b$NN_sites
interhost_results_ORF3b$S_sites <- interhost_results_ORF3b$NS_sites # NS; ORF3a perspective
interhost_results_ORF3b$N_diffs <- interhost_results_ORF3b$NN_diffs
interhost_results_ORF3b$S_diffs <- interhost_results_ORF3b$NS_diffs # NS; ORF3a perspective
interhost_results_ORF3b$codon <- interhost_results_ORF3b$codon - min(interhost_results_ORF3b$codon) + 2 # starts at 2
interhost_results_ORF3b$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF3b"), interhost_results_ORF3b)
rm(interhost_results_ORF3b)

# Alternate gene in ss13: S/S.iORF1 (S codons 61-101, trim to 62-100)
interhost_results_SiORF1 <- filter(interhost_results, product == "S", codon %in% 62:100)
interhost_results_SiORF1$product <- "SiORF1"
interhost_results_SiORF1$N_sites <- interhost_results_SiORF1$NN_sites
interhost_results_SiORF1$S_sites <- interhost_results_SiORF1$NS_sites # NS; S perspective
interhost_results_SiORF1$N_diffs <- interhost_results_SiORF1$NN_diffs
interhost_results_SiORF1$S_diffs <- interhost_results_SiORF1$NS_diffs # NS; S perspective
interhost_results_SiORF1$codon <- interhost_results_SiORF1$codon - min(interhost_results_SiORF1$codon) + 2 # starts at 2
interhost_results_SiORF1$OL_category <- "OLG" # ADDITIONAL
interhost_results_SiORF1$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "SiORF1"), interhost_results_SiORF1)
rm(interhost_results_SiORF1)

# Reference gene with OLG in ss13: N, codons 4-102 (ORF9b)
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$N_sites <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$NN_sites
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$S_sites <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$SN_sites
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$N_diffs <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$NN_diffs
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$S_diffs <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$SN_diffs
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(4, 102), ]$N_sites <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(4, 102), ]$S_sites <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(4, 102), ]$N_diffs <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(4, 102), ]$S_diffs <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 4:102, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9b (N codons 4-102, trim to 5-101)
interhost_results_ORF9b <- filter(interhost_results, product == "N", codon %in% 5:101)
interhost_results_ORF9b$product <- "ORF9b"
interhost_results_ORF9b$N_sites <- interhost_results_ORF9b$NN_sites
interhost_results_ORF9b$S_sites <- interhost_results_ORF9b$NS_sites # NS; N perspective
interhost_results_ORF9b$N_diffs <- interhost_results_ORF9b$NN_diffs
interhost_results_ORF9b$S_diffs <- interhost_results_ORF9b$NS_diffs # NS; N perspective
interhost_results_ORF9b$codon <- interhost_results_ORF9b$codon - min(interhost_results_ORF9b$codon) + 2 # starts at 2
interhost_results_ORF9b$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF9b"), interhost_results_ORF9b)
rm(interhost_results_ORF9b)

# Reference gene with OLG in ss13: N, codons 154-228 (ORF9c)
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$N_sites <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$NN_sites
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$S_sites <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$SN_sites
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$N_diffs <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$NN_diffs
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$S_diffs <- 
  interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$SN_diffs
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(154, 228), ]$N_sites <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(154, 228), ]$S_sites <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(154, 228), ]$N_diffs <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% c(154, 228), ]$S_diffs <- NA
interhost_results[interhost_results$product == "N" & interhost_results$codon %in% 154:228, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9c (N codons 154-228, trim to 155-227)
interhost_results_ORF9c <- filter(interhost_results, product == "N", codon %in% 155:227)
interhost_results_ORF9c$product <- "ORF9c"
interhost_results_ORF9c$N_sites <- interhost_results_ORF9c$NN_sites
interhost_results_ORF9c$S_sites <- interhost_results_ORF9c$NS_sites # NS; N perspective
interhost_results_ORF9c$N_diffs <- interhost_results_ORF9c$NN_diffs
interhost_results_ORF9c$S_diffs <- interhost_results_ORF9c$NS_diffs # NS; N perspective
interhost_results_ORF9c$codon <- interhost_results_ORF9c$codon - min(interhost_results_ORF9c$codon) + 2 # starts at 2
interhost_results_ORF9c$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF9c"), interhost_results_ORF9c)
rm(interhost_results_ORF9c)

# strip columns
(interhost_results <- dplyr::select(interhost_results, interhost_results_names))
# 10,064 x 11


### Genes with ss12 overlap ###
(interhost_results <- left_join(interhost_results, filter(interhost_results_OLGenie, frame == 'ss12'), by = c('product', 'codon')))
# 10,064 x 29

# Reference gene with OLG in ss12: ORF3a, codons 65-102 (ORF3d after double overlap)
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$N_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$NN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$S_sites <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$SN_sites
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$N_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$NN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$S_diffs <- 
  interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$SN_diffs
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(65, 102), ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(65, 102), ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(65, 102), ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% c(65, 102), ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 65:102, ]$OL_category2 <- "OLG"

# Alternate gene in ss12: ORF3a/ORF3d (ORF3a codons 65-102, trim to 66-101)
interhost_results_ORF3d <- filter(interhost_results, product == "ORF3a", codon %in% 66:101)
interhost_results_ORF3d$product <- "ORF3d"
interhost_results_ORF3d$N_sites <- interhost_results_ORF3d$NN_sites
interhost_results_ORF3d$S_sites <- interhost_results_ORF3d$NS_sites # NS; ORF3a perspective
interhost_results_ORF3d$N_diffs <- interhost_results_ORF3d$NN_diffs
interhost_results_ORF3d$S_diffs <- interhost_results_ORF3d$NS_diffs # NS; ORF3a perspective
interhost_results_ORF3d$codon <- interhost_results_ORF3d$codon - min(interhost_results_ORF3d$codon) + 23 # starts at 65-44+2=23
interhost_results_ORF3d$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF3d"), interhost_results_ORF3d)
rm(interhost_results_ORF3d)

# Alternate gene in ss12: ORF3a/ORF3d2 (ORF3a codons 68-102, trim to 69-101)
interhost_results_ORF3d2 <- filter(interhost_results, product == "ORF3a", codon %in% 69:101)
interhost_results_ORF3d2$product <- "ORF3d2"
interhost_results_ORF3d2$N_sites <- interhost_results_ORF3d2$NN_sites
interhost_results_ORF3d2$S_sites <- interhost_results_ORF3d2$NS_sites # NS; ORF3a perspective
interhost_results_ORF3d2$N_diffs <- interhost_results_ORF3d2$NN_diffs
interhost_results_ORF3d2$S_diffs <- interhost_results_ORF3d2$NS_diffs # NS; ORF3a perspective
interhost_results_ORF3d2$codon <- interhost_results_ORF3d2$codon - min(interhost_results_ORF3d2$codon) + 2 # starts at 2
interhost_results_ORF3d2$OL_category2 <- "OLG"
interhost_results <- rbind(filter(interhost_results, product != "ORF3d2"), interhost_results_ORF3d2)
rm(interhost_results_ORF3d2)

# strip columns
(interhost_results <- dplyr::select(interhost_results, interhost_results_names))


### DOUBLE OVERLAP REGION
# Reference gene with Double OLG: ORF3a, codons 44-64 (ORF3a/ORF3c/ORF3d)
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 44:64, ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 44:64, ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 44:64, ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 44:64, ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF3a" & interhost_results$codon %in% 44:64, ]$OL_category2 <- "Double OLG"


### OTHER OLG REGIONS
# ORF1ab/ORF1a: ORF1ab codons 4401-4407
interhost_results[interhost_results$product == "ORF1ab" & interhost_results$codon %in% 4401:4407, ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF1ab" & interhost_results$codon %in% 4401:4407, ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF1ab" & interhost_results$codon %in% 4401:4407, ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF1ab" & interhost_results$codon %in% 4401:4407, ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF1ab" & interhost_results$codon %in% 4401:4407, ]$OL_category2 <- "OLG"

# ORF7a: codons 121-122 (overlap with ORF7b)
interhost_results[interhost_results$product == "ORF7a" & interhost_results$codon %in% 121:122, ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF7a" & interhost_results$codon %in% 121:122, ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF7a" & interhost_results$codon %in% 121:122, ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF7a" & interhost_results$codon %in% 121:122, ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF7a" & interhost_results$codon %in% 121:122, ]$OL_category2 <- "OLG"

# ORF7b: codons 1-2 (overlap with ORF7a)
interhost_results[interhost_results$product == "ORF7b" & interhost_results$codon %in% 1:2, ]$N_sites <- NA
interhost_results[interhost_results$product == "ORF7b" & interhost_results$codon %in% 1:2, ]$S_sites <- NA
interhost_results[interhost_results$product == "ORF7b" & interhost_results$codon %in% 1:2, ]$N_diffs <- NA
interhost_results[interhost_results$product == "ORF7b" & interhost_results$codon %in% 1:2, ]$S_diffs <- NA
interhost_results[interhost_results$product == "ORF7b" & interhost_results$codon %in% 1:2, ]$OL_category2 <- "OLG"


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
interhost_results$num_defined_seqs <- 6


####################################################################################################
### BOOTSTRAP EACH STATISTIC

############################################################################################################
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT
dNdS_diff_boot_fun <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  
  (dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE))
  (dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE))
  (dNdS <- dN / dS)
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0) # 345
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0) # 0
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0) # 655
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}

############################################################################################################
### ANALYSIS VARIABLES
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6


############################################################################################################
### INITIALIZE DATA FRAME: INTERHOST
interhost_results_bootstrap <- data.frame(gene_name = character(),
                                          OL_category = character(), # ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric())


### LOOP EACH DATA SUBSET (4 minutes with 6 CPUs)
for (this_gene in sort(unique(interhost_results$product))) {
  #this_gene <- 'ORF10'
  
  # ELIMINATE FOR SNPGENIE-ONLY <-- CHANGE
  for (this_OL_cat in sort(unique(interhost_results$OL_category))) {
    
    # Filter; CHANGE: ELIMINATE OL_category for SNPGENIE ONLY <-- CHANGE
    this_data <- filter(interhost_results, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS, 
                        OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG")) # OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG")
    
    if(nrow(this_data) >= 9) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      #if(PREPEND_TO_OUTPUT != '') {
      #  summary_data <- paste(PREPEND_TO_OUTPUT, summary_data, sep = "\t")
      #}
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      interhost_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                         OL_category = this_OL_cat, # ELIMINATE FOR SNPGENIE ONLY  <-- CHANGE
                                                         num_bootstraps = NBOOTSTRAPS,
                                                         min_defined_codons = MIN_DEFINED_CODONS,
                                                         num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                         N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                         S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                         N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                         S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                         num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                         dN = as.numeric(boot_dNdS_vector['dN']),
                                                         dS = as.numeric(boot_dNdS_vector['dS']),
                                                         dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                         dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                         boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                         boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                         boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                         boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                         boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                         P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                         boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                         boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                         boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                         ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                         ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
      
      # Add the 4 new rows to results
      interhost_results_bootstrap <- rbind(interhost_results_bootstrap, interhost_results_bootstrap_ADDITION)
    }
  } # ELIMINATE FOR SNPGENIE-ONLY <-- CHANGE
}
# 13 warnings: NAs introduced by coercion (OKAY)

names(interhost_results_bootstrap) <- c('gene_name', 
                                        'OL_category', # ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

### FILTER UNUSED GENES
interhost_results_bootstrap <- filter(interhost_results_bootstrap, ! gene_name %in% c("ORF3d2", "SiORF1")) # <-- CHANGE


### Manual 2-sided ASL P-value
interhost_results_bootstrap$P_ALS <- NA
interhost_results_bootstrap[interhost_results_bootstrap$ASL_dN_gt_dS_P < interhost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * interhost_results_bootstrap[interhost_results_bootstrap$ASL_dN_gt_dS_P < interhost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
interhost_results_bootstrap[interhost_results_bootstrap$ASL_dN_gt_dS_P > interhost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * interhost_results_bootstrap[interhost_results_bootstrap$ASL_dN_gt_dS_P > interhost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
interhost_results_bootstrap[interhost_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS


### SAVE
#write_tsv(interhost_results_bootstrap, "data/between-host_FINAL.tsv")



####################################################################################################
####################################################################################################
### WITHIN-HOST (SRR*), called 'intrahost' herein
####################################################################################################
####################################################################################################


####################################################################################################
# LOAD with numbered codons
(intrahost_results <- read_tsv("data/within-host_SNPGenie.tsv"))

# check codon counts
intrahost_results %>%
  group_by(product) %>%
  summarise(
    highest_codon_num = max(codon)
  ) # confirmed

# rename products
intrahost_results[intrahost_results$product == 'ORF3bp', ]$product <- "ORF3d" # for continuity; will not use anyway (replaced with reference frame perspective)
intrahost_results[intrahost_results$product == 'ORF8p', ]$product <- "ORF8"


### SUMMARIZE RESULTS BY CODON MEANS, ACROSS ALL SAMPLES: SPECIFIC TO INTRAHOST
(intrahost_results_codonMeans <- intrahost_results %>%
    group_by(product, codon) %>% ## diff
    summarise(
      N_sites = mean(N_sites, na.rm = TRUE),
      S_sites = mean(S_sites, na.rm = TRUE),
      N_diffs = mean(N_diffs, na.rm = TRUE),
      S_diffs = mean(S_diffs, na.rm = TRUE)
    )) # 9,985 rows

# clear workspace
intrahost_results <- intrahost_results_codonMeans
rm(intrahost_results_codonMeans)

# JOIN OL CATEGORY # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
codon_OL_category <- dplyr::select(interhost_results, product, product_order, codon, OL_category) # <-- CHANGE?
intrahost_results <- left_join(intrahost_results, codon_OL_category, by = c('product', 'codon')) # <-- CHANGE?
intrahost_results # 9,985 x 8

# Another column for double checking
intrahost_results$OL_category2 <- "Non-OLG"

# Bookkeeping
(intrahost_results %>% 
       group_by(product, OL_category) %>%
       summarise(
         count = n()
       ))

# Save names to remember which columns to keep later
intrahost_results_names <- names(intrahost_results)

# Save the SNPGenie-only results
#write_tsv(intrahost_results, "data/within-host_SNPGenie_byCodon.tsv")


####################################################################################################
# IMPORT OLGenie & COMBINE <-- SKIP IF SNPGENIE ONLY
(intrahost_results_OLGenie <- read_tsv("data/within-host_OLGenie.tsv"))
# 644,808 x 15

# Process
intrahost_results_OLGenie$product <- str_replace(string = intrahost_results_OLGenie$product, pattern = "../SARSCOV2_MAFFT_processed_", replacement = "")
intrahost_results_OLGenie$product <- str_replace(string = intrahost_results_OLGenie$product, pattern = "../SRR\\d+_bbduk_removed_lofreq3_filtered_nSeqs1000_", replacement = "") # for INTRAHOST
intrahost_results_OLGenie$product <- str_replace(string = intrahost_results_OLGenie$product, pattern = ".fasta", replacement = "")
intrahost_results_OLGenie[intrahost_results_OLGenie$product == "ORF3bp", ]$product <- "ORF3d" # for continuity; won't use

### SUMMARIZE RESULTS BY CODON MEANS, SPECIFIC TO INTRAHOST
intrahost_results_OLGenie_codonMeans <- intrahost_results_OLGenie %>%
  group_by(product, frame, codon) %>% ## diff
  summarise(
    NN_sites = mean(NN_sites, na.rm = TRUE),
    SN_sites = mean(SN_sites, na.rm = TRUE),
    NS_sites = mean(NS_sites, na.rm = TRUE),
    SS_sites = mean(SS_sites, na.rm = TRUE),
    NN_diffs = mean(NN_diffs, na.rm = TRUE),
    SN_diffs = mean(SN_diffs, na.rm = TRUE),
    NS_diffs = mean(NS_diffs, na.rm = TRUE),
    SS_diffs = mean(SS_diffs, na.rm = TRUE),
    count = n() # 401 samples
  ) 


# clear workspace
intrahost_results_OLGenie <- intrahost_results_OLGenie_codonMeans
rm(intrahost_results_OLGenie_codonMeans)

# Limit to analyzed gene regions, using reference frame perspective for OLGs
products_to_analyze <- c("S", "ORF3a", "N") # we won't consider the S-iORF to be legit, so won't exclude from S or show in figure
(intrahost_results_OLGenie <- filter(intrahost_results_OLGenie, product %in% products_to_analyze))

# Check the needed frames exist (ss12 and ss13 for ORF3a; ss13 for N)
intrahost_results_OLGenie %>% 
  group_by(product, frame) %>%
  summarise(
    count = n()
  ) # YES


####################################################################################################
### Simply JOIN by product/codon, add coluns at the end; replace SNPGenie values with OLGenie values at the proper positions; delete superfluous columns

### Genes with ss13 overlap ###
intrahost_results <- left_join(intrahost_results, filter(intrahost_results_OLGenie, frame == 'ss13'), by = c('product', 'codon'))

# Reference gene with OLG in ss13: ORF3a, codons 22-43 (ORF3c before double overlap)
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$N_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$NN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$S_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$SN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$N_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$NN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$S_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$SN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(22, 43), ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(22, 43), ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(22, 43), ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(22, 43), ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 22:43, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: ORF3a/ORF3c (ORF3a codons 22-43, trim to 23-42)
intrahost_results_ORF3c <- filter(intrahost_results, product == "ORF3a", codon %in% 23:42)
intrahost_results_ORF3c$product <- "ORF3c"
intrahost_results_ORF3c$N_sites <- intrahost_results_ORF3c$NN_sites
intrahost_results_ORF3c$S_sites <- intrahost_results_ORF3c$NS_sites # NS; ORF3a perspective
intrahost_results_ORF3c$N_diffs <- intrahost_results_ORF3c$NN_diffs
intrahost_results_ORF3c$S_diffs <- intrahost_results_ORF3c$NS_diffs # NS; ORF3a perspective
intrahost_results_ORF3c$codon <- intrahost_results_ORF3c$codon - min(intrahost_results_ORF3c$codon) + 2 # starts at 2
intrahost_results_ORF3c$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF3c"), intrahost_results_ORF3c)
rm(intrahost_results_ORF3c)

# Reference gene with OLG in ss13: ORF3a, codons 141-164 (ORF3b, first ORF)
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$N_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$NN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$S_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$SN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$N_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$NN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$S_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$SN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(141, 164), ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(141, 164), ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(141, 164), ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(141, 164), ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 141:164, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: ORF3a/ORF3b (ORF3a codons 141-164, trim to 142-163)
intrahost_results_ORF3b <- filter(intrahost_results, product == "ORF3a", codon %in% 142:163)
intrahost_results_ORF3b$product <- "ORF3b"
intrahost_results_ORF3b$N_sites <- intrahost_results_ORF3b$NN_sites
intrahost_results_ORF3b$S_sites <- intrahost_results_ORF3b$NS_sites # NS; ORF3a perspective
intrahost_results_ORF3b$N_diffs <- intrahost_results_ORF3b$NN_diffs
intrahost_results_ORF3b$S_diffs <- intrahost_results_ORF3b$NS_diffs # NS; ORF3a perspective
intrahost_results_ORF3b$codon <- intrahost_results_ORF3b$codon - min(intrahost_results_ORF3b$codon) + 2 # starts at 2
intrahost_results_ORF3b$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF3b"), intrahost_results_ORF3b)
rm(intrahost_results_ORF3b)

# Reference gene with OLG in ss13: N, codons 4-102 (ORF9b)
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$N_sites <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$NN_sites
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$S_sites <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$SN_sites
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$N_diffs <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$NN_diffs
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$S_diffs <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$SN_diffs
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(4, 102), ]$N_sites <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(4, 102), ]$S_sites <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(4, 102), ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(4, 102), ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 4:102, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9b (N codons 4-102, trim to 5-101)
intrahost_results_ORF9b <- filter(intrahost_results, product == "N", codon %in% 5:101)
intrahost_results_ORF9b$product <- "ORF9b"
intrahost_results_ORF9b$N_sites <- intrahost_results_ORF9b$NN_sites
intrahost_results_ORF9b$S_sites <- intrahost_results_ORF9b$NS_sites # NS; N perspective
intrahost_results_ORF9b$N_diffs <- intrahost_results_ORF9b$NN_diffs
intrahost_results_ORF9b$S_diffs <- intrahost_results_ORF9b$NS_diffs # NS; N perspective
intrahost_results_ORF9b$codon <- intrahost_results_ORF9b$codon - min(intrahost_results_ORF9b$codon) + 2 # starts at 2
intrahost_results_ORF9b$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF9b"), intrahost_results_ORF9b)
rm(intrahost_results_ORF9b)

# Reference gene with OLG in ss13: N, codons 154-228 (ORF9c)
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$N_sites <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$NN_sites
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$S_sites <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$SN_sites
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$N_diffs <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$NN_diffs
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$S_diffs <- 
  intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$SN_diffs
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(154, 228), ]$N_sites <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(154, 228), ]$S_sites <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(154, 228), ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% c(154, 228), ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "N" & intrahost_results$codon %in% 154:228, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9c (N codons 154-228, trim to 155-227)
intrahost_results_ORF9c <- filter(intrahost_results, product == "N", codon %in% 155:227)
intrahost_results_ORF9c$product <- "ORF9c"
intrahost_results_ORF9c$N_sites <- intrahost_results_ORF9c$NN_sites
intrahost_results_ORF9c$S_sites <- intrahost_results_ORF9c$NS_sites # NS; N perspective
intrahost_results_ORF9c$N_diffs <- intrahost_results_ORF9c$NN_diffs
intrahost_results_ORF9c$S_diffs <- intrahost_results_ORF9c$NS_diffs # NS; N perspective
intrahost_results_ORF9c$codon <- intrahost_results_ORF9c$codon - min(intrahost_results_ORF9c$codon) + 2 # starts at 2
intrahost_results_ORF9c$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF9c"), intrahost_results_ORF9c)
rm(intrahost_results_ORF9c)

# strip columns
intrahost_results <- dplyr::select(intrahost_results, intrahost_results_names)

## Genes with ss12 overlap
(intrahost_results <- left_join(intrahost_results, filter(intrahost_results_OLGenie, frame == 'ss12'), by = c('product', 'codon')))

# Reference gene with OLG in ss12: ORF3a, codons 65-102 (ORF3d after double overlap)
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$N_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$NN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$S_sites <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$SN_sites
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$N_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$NN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$S_diffs <- 
  intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$SN_diffs
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(65, 102), ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(65, 102), ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(65, 102), ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% c(65, 102), ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 65:102, ]$OL_category2 <- "OLG"

# Alternate gene in ss12: ORF3a/ORF3d (ORF3a codons 65-102, trim to 66-101)
intrahost_results_ORF3d <- filter(intrahost_results, product == "ORF3a", codon %in% 66:101)
intrahost_results_ORF3d$product <- "ORF3d"
intrahost_results_ORF3d$N_sites <- intrahost_results_ORF3d$NN_sites
intrahost_results_ORF3d$S_sites <- intrahost_results_ORF3d$NS_sites # NS; ORF3a perspective
intrahost_results_ORF3d$N_diffs <- intrahost_results_ORF3d$NN_diffs
intrahost_results_ORF3d$S_diffs <- intrahost_results_ORF3d$NS_diffs # NS; ORF3a perspective
intrahost_results_ORF3d$codon <- intrahost_results_ORF3d$codon - min(intrahost_results_ORF3d$codon) + 23 # starts at 65-44+2=23
intrahost_results_ORF3d$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF3d"), intrahost_results_ORF3d)
rm(intrahost_results_ORF3d)

# Alternate gene in ss12: ORF3a/ORF3d2 (ORF3a codons 68-102, trim to 69-101)
intrahost_results_ORF3d2 <- filter(intrahost_results, product == "ORF3a", codon %in% 69:101)
intrahost_results_ORF3d2$product <- "ORF3d2"
intrahost_results_ORF3d2$N_sites <- intrahost_results_ORF3d2$NN_sites
intrahost_results_ORF3d2$S_sites <- intrahost_results_ORF3d2$NS_sites # NS; ORF3a perspective
intrahost_results_ORF3d2$N_diffs <- intrahost_results_ORF3d2$NN_diffs
intrahost_results_ORF3d2$S_diffs <- intrahost_results_ORF3d2$NS_diffs # NS; ORF3a perspective
intrahost_results_ORF3d2$codon <- intrahost_results_ORF3d2$codon - min(intrahost_results_ORF3d2$codon) + 2 # starts at 2
intrahost_results_ORF3d2$OL_category2 <- "OLG"
intrahost_results <- rbind(filter(intrahost_results, product != "ORF3d2"), intrahost_results_ORF3d2)
rm(intrahost_results_ORF3d2)

# strip columns
intrahost_results <- dplyr::select(intrahost_results, intrahost_results_names)


### DOUBLE OVERLAP REGION
# Reference gene with Double OLG: ORF3a, codons 44-64 (ORF3a/ORF3c/ORF3d)
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 44:64, ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 44:64, ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 44:64, ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 44:64, ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF3a" & intrahost_results$codon %in% 44:64, ]$OL_category2 <- "Double OLG"


### OTHER OLG REGIONS
# ORF1ab/ORF1a: ORF1ab codons 4401-4407
intrahost_results[intrahost_results$product == "ORF1ab" & intrahost_results$codon %in% 4401:4407, ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF1ab" & intrahost_results$codon %in% 4401:4407, ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF1ab" & intrahost_results$codon %in% 4401:4407, ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF1ab" & intrahost_results$codon %in% 4401:4407, ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF1ab" & intrahost_results$codon %in% 4401:4407, ]$OL_category2 <- "OLG"

# ORF7a: codons 121-122 (overlap with ORF7b)
intrahost_results[intrahost_results$product == "ORF7a" & intrahost_results$codon %in% 121:122, ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF7a" & intrahost_results$codon %in% 121:122, ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF7a" & intrahost_results$codon %in% 121:122, ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF7a" & intrahost_results$codon %in% 121:122, ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF7a" & intrahost_results$codon %in% 121:122, ]$OL_category2 <- "OLG"

# ORF7b: codons 1-2 (overlap with ORF7a)
intrahost_results[intrahost_results$product == "ORF7b" & intrahost_results$codon %in% 1:2, ]$N_sites <- NA
intrahost_results[intrahost_results$product == "ORF7b" & intrahost_results$codon %in% 1:2, ]$S_sites <- NA
intrahost_results[intrahost_results$product == "ORF7b" & intrahost_results$codon %in% 1:2, ]$N_diffs <- NA
intrahost_results[intrahost_results$product == "ORF7b" & intrahost_results$codon %in% 1:2, ]$S_diffs <- NA
intrahost_results[intrahost_results$product == "ORF7b" & intrahost_results$codon %in% 1:2, ]$OL_category2 <- "OLG"


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
intrahost_results$num_defined_seqs <- 6


############################################################################################################
# BOOTSTRAP FUNCTION: same as before

############################################################################################################
### ANALYSIS VARIABLES: same as before
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6

############################################################################################################
### INITIALIZE DATA FRAME: intrahost
intrahost_results_bootstrap <- data.frame(gene_name = character(),
                                          OL_category = character(), # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
                                          num_bootstraps = integer(),
                                          min_defined_codons = integer(),
                                          num_codons = integer(),
                                          N_sites = numeric(),
                                          S_sites = numeric(),
                                          N_diffs = numeric(),
                                          S_diffs = numeric(),
                                          num_replicates = integer(),
                                          dN = numeric(),
                                          dS = numeric(),
                                          dNdS = numeric(),
                                          dN_m_dS = numeric(),
                                          boot_dN_SE = numeric(),
                                          boot_dS_SE = numeric(),
                                          boot_dN_over_dS_SE = numeric(),
                                          boot_dN_over_dS_P = numeric(),
                                          boot_dN_m_dS_SE = numeric(),
                                          P_value = numeric(),
                                          boot_dN_gt_dS_count = integer(), 
                                          boot_dN_eq_dS_count = integer(), 
                                          boot_dN_lt_dS_count = integer(), 
                                          ASL_dN_gt_dS_P = numeric(), 
                                          ASL_dN_lt_dS_P = numeric())


### LOOP EACH DATA SUBSET (4 minutes with 6 CPUs)
for (this_gene in sort(unique(intrahost_results$product))) {
  #this_gene <- 'ORF10'
  
  # ELIMINATE IF SNPGENIE ONLY
  for (this_OL_cat in sort(unique(intrahost_results$OL_category))) {
    
    # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
    # Filter by gene, frame, and minimum number of defined codons
    this_data <- filter(intrahost_results, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS,
                        OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG")) # OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG") <-- CHANGE
    #this_data <- filter(intrahost_results, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS) # <-- CHANGE
    
    if(nrow(this_data) >= 9) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      #if(PREPEND_TO_OUTPUT != '') {
      #  summary_data <- paste(PREPEND_TO_OUTPUT, summary_data, sep = "\t")
      #}
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      intrahost_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                         OL_category = this_OL_cat, # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
                                                         num_bootstraps = NBOOTSTRAPS,
                                                         min_defined_codons = MIN_DEFINED_CODONS,
                                                         num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                         N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                         S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                         N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                         S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                         num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                         dN = as.numeric(boot_dNdS_vector['dN']),
                                                         dS = as.numeric(boot_dNdS_vector['dS']),
                                                         dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                         dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                         boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                         boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                         boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                         boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                         boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                         P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                         boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                         boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                         boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                         ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                         ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
      
      # Add the 4 new rows to results
      intrahost_results_bootstrap <- rbind(intrahost_results_bootstrap, intrahost_results_bootstrap_ADDITION)
    }
  } # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
}
# warnings: NAs introduced by coercion (OKAY)

names(intrahost_results_bootstrap) <- c('gene_name', 
                                        'OL_category', # ELIMINATE IF SNPGENIE ONLY <-- CHANGE
                                        'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                        'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                        'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                        'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                        'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

### FILTER UNUSED GENES
unique(intrahost_results_bootstrap$gene_name)
intrahost_results_bootstrap <- filter(intrahost_results_bootstrap, ! gene_name %in% c("ORF3d2", "SiORF1")) # <-- CHANGE


### Manual 2-sided ASL P-value
intrahost_results_bootstrap$P_ALS <- NA
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P < intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$ASL_dN_gt_dS_P) & intrahost_results_bootstrap$ASL_dN_gt_dS_P > intrahost_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
intrahost_results_bootstrap[! is.na(intrahost_results_bootstrap$P_ALS) & intrahost_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS


### SAVE
#write_tsv(intrahost_results_bootstrap, "data/within-host_FINAL.tsv")



####################################################################################################
####################################################################################################
### BETWEEN-TAXA (Severe acute respiratory syndrome-related coronaviruses, n=21), herein called 'interspecies'
####################################################################################################
####################################################################################################

####################################################################################################
# IMPORT **SNPGenie** 
(interspecies_results <- read_tsv("data/between-taxa_SNPGenie.tsv"))

# rename ORF3bp
interspecies_results[interspecies_results$product == "ORF3bp", ]$product <- "ORF3d" # for continuity; won't use anyway

### JOIN OLG CATEGORY (add ORF8b ourselves); manually classified every codon as NOL, OL, or dOL (double overlapping; ORF3c/ORF3d region)
(codon_OL_category <- read_tsv("data/between-taxa_codon_OL_category.tsv"))

#join
interspecies_results <- left_join(interspecies_results, codon_OL_category, by = c('product', 'codon'))

# classify ORF8b as NOL
interspecies_results[interspecies_results$product == "ORF8b", ]$OL_category <- "NOL"

# Rename OLG categories
interspecies_results[interspecies_results$OL_category == 'NOL', ]$OL_category <- "Non-OLG"
interspecies_results[interspecies_results$OL_category == 'OL', ]$OL_category <- "OLG"
interspecies_results[interspecies_results$OL_category == 'dOL', ]$OL_category <- "Double OLG"


# Reclassify the ORF3a/ORF3b region, codons 165-276, as "Non-OLG" because we only included non-SARS, non-ORF3b; will replace first ORF3b later
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$OL_category <- "Non-OLG"

# Another column for double checking
interspecies_results$OL_category2 <- "Non-OLG"

# Save names to remember which columns to keep later
interspecies_results_names <- names(interspecies_results)


####################################################################################################
# IMPORT TAXA-RESTRICTED SNPGenie & COMBINE: SNPGenie_ORF3a_allMinus12SARS

## ORF3c: all sequences
## ORF3d: SARS-CoV-2 and pangolin-CoV GX/P5L
## ORF3b: BELOW, we include the SNPGenie result for only those 21-12=9 sequences LACKING full-length ORF3b (for ORF3a Non-OLG regions)
## ORF9b: all
## ORF9c: all

# input results restricted to taxa with functional ORFs
interspecies_results_Rest<- read_tsv("data/between-taxa_SNPGenie_taxaRestricted.tsv")

### Simply JOIN by product/codon; replace SNPGenie values with restricted values at the proper positions; delete superfluous columns
(interspecies_results <- left_join(interspecies_results, interspecies_results_Rest, by = c('product', 'codon')))
#10,419 x 16

# Replace values in the ORF3b region with those without full-length ORF3b
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$N_sites_Rest
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$S_sites_Rest
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$N_diffs_Rest
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon >= 165, ]$S_diffs_Rest

# strip superfluous columns
interspecies_results <- dplyr::select(interspecies_results, interspecies_results_names)



####################################################################################################
# IMPORT OLGenie & COMBINE

# Import full analysis
(interspecies_results_OLGenie <- read_tsv("data/between-taxa_OLGenie.tsv"))
#920 x 15

# Import ORF3a ss13
(interspecies_results_OLGenie2 <- read_tsv("data/between-taxa_OLGenie2.tsv"))

# combine
interspecies_results_OLGenie <- rbind(interspecies_results_OLGenie, interspecies_results_OLGenie2)
rm(interspecies_results_OLGenie2)

# Process
interspecies_results_OLGenie$product <- str_replace(string = interspecies_results_OLGenie$product, pattern = ".*_([a-zA-Z0-9]+).fasta", replacement = "\\1") # for Sarbecovirus

# rename ORF3bp, one of many old names for ORF3d; but won't use anyway
interspecies_results_OLGenie[interspecies_results_OLGenie$product == "ORF3bp", ]$product <- "ORF3d"

# Limit to analyzed gene regions, using reference frame perspective for OLGs
products_to_analyze <- c("S", "ORF3a", "N") # we won't consider the S-iORF to be legit, so won't exclude from S or show in figure
(interspecies_results_OLGenie <- filter(interspecies_results_OLGenie, product %in% products_to_analyze))
unique(paste0(interspecies_results_OLGenie$product, "-", interspecies_results_OLGenie$frame))


### Simply JOIN by PAIR/product/codon (pair being particular analysis, e.g., "all" sequences); 
### add columns at the end; replace SNPGenie values with OLGenie values at the proper positions; delete superfluous columns

### Genes with ss13 overlap ###
(interspecies_results <- left_join(interspecies_results, filter(interspecies_results_OLGenie, frame == 'ss13'), by = c('product', 'codon')))
# 10,419 x 25

# Reference gene with OLG in ss13: ORF3a, codons 22-43 (ORF3c before double overlap)
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$NN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$SN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$NN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$SN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(22, 43), ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(22, 43), ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(22, 43), ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(22, 43), ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 22:43, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: ORF3a/ORF3c (ORF3a codons 22-43, trim to 23-42)
interspecies_results_ORF3c <- filter(interspecies_results, product == "ORF3a", codon %in% 23:42)
interspecies_results_ORF3c$product <- "ORF3c"
interspecies_results_ORF3c$N_sites <- interspecies_results_ORF3c$NN_sites
interspecies_results_ORF3c$S_sites <- interspecies_results_ORF3c$NS_sites # NS; ORF3a perspective
interspecies_results_ORF3c$N_diffs <- interspecies_results_ORF3c$NN_diffs
interspecies_results_ORF3c$S_diffs <- interspecies_results_ORF3c$NS_diffs # NS; ORF3a perspective
interspecies_results_ORF3c$codon <- interspecies_results_ORF3c$codon - min(interspecies_results_ORF3c$codon) + 2 # starts at 2
interspecies_results_ORF3c$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF3c"), interspecies_results_ORF3c)
rm(interspecies_results_ORF3c)

# Reference gene with OLG in ss13: ORF3a, codons 141-164 (ORF3b, first ORF)
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$NN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$SN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$NN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$SN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(141, 164), ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(141, 164), ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(141, 164), ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(141, 164), ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 141:164, ]$OL_category2 <- "OLG"
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon > 164, ]$OL_category2 <- "Non-OLG" 
# ^last is additional because only included subset; already joined non-SARS (ORF3b) for non-OLG estimate

# Alternate gene in ss13: ORF3a/ORF3b (ORF3a codons 141-164, trim to 142-163)
interspecies_results_ORF3b <- filter(interspecies_results, product == "ORF3a", codon %in% 142:163)
interspecies_results_ORF3b$product <- "ORF3b"
interspecies_results_ORF3b$N_sites <- interspecies_results_ORF3b$NN_sites
interspecies_results_ORF3b$S_sites <- interspecies_results_ORF3b$NS_sites # NS; ORF3a perspective
interspecies_results_ORF3b$N_diffs <- interspecies_results_ORF3b$NN_diffs
interspecies_results_ORF3b$S_diffs <- interspecies_results_ORF3b$NS_diffs # NS; ORF3a perspective
interspecies_results_ORF3b$codon <- interspecies_results_ORF3b$codon - min(interspecies_results_ORF3b$codon) + 2 # starts at 2
interspecies_results_ORF3b$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF3b"), interspecies_results_ORF3b)
rm(interspecies_results_ORF3b)

# Reference gene with OLG in ss13: N, codons 4-103 (ORF9b)
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$NN_sites
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$SN_sites
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$NN_diffs
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$SN_diffs
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(4, 103), ]$N_sites <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(4, 103), ]$S_sites <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(4, 103), ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(4, 103), ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 4:103, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9b (N codons 4-103, trim to 5-102)
interspecies_results_ORF9b <- filter(interspecies_results, product == "N", codon %in% 5:102)
interspecies_results_ORF9b$product <- "ORF9b"
interspecies_results_ORF9b$N_sites <- interspecies_results_ORF9b$NN_sites
interspecies_results_ORF9b$S_sites <- interspecies_results_ORF9b$NS_sites # NS; N perspective
interspecies_results_ORF9b$N_diffs <- interspecies_results_ORF9b$NN_diffs
interspecies_results_ORF9b$S_diffs <- interspecies_results_ORF9b$NS_diffs # NS; N perspective
interspecies_results_ORF9b$codon <- interspecies_results_ORF9b$codon - min(interspecies_results_ORF9b$codon) + 2 # starts at 2
interspecies_results_ORF9b$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF9b"), interspecies_results_ORF9b)
rm(interspecies_results_ORF9b)

# Reference gene with OLG in ss13: N, codons 155-229 (ORF9c)
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$NN_sites
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$SN_sites
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$NN_diffs
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$SN_diffs
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(155, 229), ]$N_sites <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(155, 229), ]$S_sites <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(155, 229), ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% c(155, 229), ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "N" & interspecies_results$codon %in% 155:229, ]$OL_category2 <- "OLG"

# Alternate gene in ss13: N/ORF9c (N codons 155-229, trim to 156-228)
interspecies_results_ORF9c <- filter(interspecies_results, product == "N", codon %in% 156:228)
interspecies_results_ORF9c$product <- "ORF9c"
interspecies_results_ORF9c$N_sites <- interspecies_results_ORF9c$NN_sites
interspecies_results_ORF9c$S_sites <- interspecies_results_ORF9c$NS_sites # NS; N perspective
interspecies_results_ORF9c$N_diffs <- interspecies_results_ORF9c$NN_diffs
interspecies_results_ORF9c$S_diffs <- interspecies_results_ORF9c$NS_diffs # NS; N perspective
interspecies_results_ORF9c$codon <- interspecies_results_ORF9c$codon - min(interspecies_results_ORF9c$codon) + 2 # starts at 2
interspecies_results_ORF9c$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF9c"), interspecies_results_ORF9c)
rm(interspecies_results_ORF9c)

# strip columns
interspecies_results <- dplyr::select(interspecies_results, interspecies_results_names)


### NO genes with ss12 overlap shared by the full set of genomes; see immediately below



####################################################################################################
# IMPORT OLGenie for ORF3d analysis of human/P5L with Wei-Zhang discovered NN & COMBINE (P4L later) <-- CHANGE
(interspecies_results_OLGenie <- read_tsv("data/between-taxa_OLGenie_addNN.tsv"))
# 274 x 15

# Process
interspecies_results_OLGenie$product <- "ORF3a"
unique(paste0(interspecies_results_OLGenie$product, "-", interspecies_results_OLGenie$frame))
# "ORF3a-ss12"


### Simply JOIN by pair/product/codon; replace SNPGenie values with OLGenie values at the proper positions; delete superfluous columns

### Genes with ss12 overlap ###
(interspecies_results <- left_join(interspecies_results, filter(interspecies_results_OLGenie, frame == 'ss12'), by = c('product', 'codon')))

# Reference gene with OLG in ss12: ORF3a, codons 65-102 (ORF3d after double overlap)
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$N_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$NN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$S_sites <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$SN_sites
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$N_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$NN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$S_diffs <- 
  interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$SN_diffs
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(65, 102), ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(65, 102), ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(65, 102), ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% c(65, 102), ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 65:102, ]$OL_category2 <- "OLG"

# Alternate gene in ss12: ORF3a/ORF3d (ORF3a codons 65-102, trim to 66-101)
interspecies_results_ORF3d <- filter(interspecies_results, product == "ORF3a", codon %in% 66:101)
interspecies_results_ORF3d$product <- "ORF3d"
interspecies_results_ORF3d$N_sites <- interspecies_results_ORF3d$NN_sites
interspecies_results_ORF3d$S_sites <- interspecies_results_ORF3d$NS_sites # NS; ORF3a perspective
interspecies_results_ORF3d$N_diffs <- interspecies_results_ORF3d$NN_diffs
interspecies_results_ORF3d$S_diffs <- interspecies_results_ORF3d$NS_diffs # NS; ORF3a perspective
interspecies_results_ORF3d$codon <- interspecies_results_ORF3d$codon - min(interspecies_results_ORF3d$codon) + 23 # starts at 65-44+2=23
interspecies_results_ORF3d$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF3d"), interspecies_results_ORF3d)
rm(interspecies_results_ORF3d)

# Alternate gene in ss12: ORF3a/ORF3d2 (ORF3a codons 68-102, trim to 69-101)
interspecies_results_ORF3d2 <- filter(interspecies_results, product == "ORF3a", codon %in% 69:101)
interspecies_results_ORF3d2$product <- "ORF3d2"
interspecies_results_ORF3d2$N_sites <- interspecies_results_ORF3d2$NN_sites
interspecies_results_ORF3d2$S_sites <- interspecies_results_ORF3d2$NS_sites # NS; ORF3a perspective
interspecies_results_ORF3d2$N_diffs <- interspecies_results_ORF3d2$NN_diffs
interspecies_results_ORF3d2$S_diffs <- interspecies_results_ORF3d2$NS_diffs # NS; ORF3a perspective
interspecies_results_ORF3d2$codon <- interspecies_results_ORF3d2$codon - min(interspecies_results_ORF3d2$codon) + 2 # starts at 2
interspecies_results_ORF3d2$OL_category2 <- "OLG"
interspecies_results <- rbind(filter(interspecies_results, product != "ORF3d2"), interspecies_results_ORF3d2)
rm(interspecies_results_ORF3d2)

# strip columns
interspecies_results <- dplyr::select(interspecies_results, interspecies_results_names)


### DOUBLE OVERLAP REGION
# Reference gene with Double OLG: ORF3a, codons 44-64 (ORF3a/ORF3c/ORF3d)
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 44:64, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 44:64, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 44:64, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 44:64, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF3a" & interspecies_results$codon %in% 44:64, ]$OL_category2 <- "Double OLG"


### OTHER OLG REGIONS
# ORF1ab/ORF1a: ORF1ab codons 4460-4466 <-- also differs from SARS-CoV-2
interspecies_results[interspecies_results$product == "ORF1ab" & interspecies_results$codon %in% 4460:4466, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF1ab" & interspecies_results$codon %in% 4460:4466, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF1ab" & interspecies_results$codon %in% 4460:4466, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF1ab" & interspecies_results$codon %in% 4460:4466, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF1ab" & interspecies_results$codon %in% 4460:4466, ]$OL_category2 <- "OLG"

# ORF7a: codons 122-124 (overlap with ORF7b)  <-- also differs from SARS-CoV-2
interspecies_results[interspecies_results$product == "ORF7a" & interspecies_results$codon %in% 122:124, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF7a" & interspecies_results$codon %in% 122:124, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF7a" & interspecies_results$codon %in% 122:124, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF7a" & interspecies_results$codon %in% 122:124, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF7a" & interspecies_results$codon %in% 122:124, ]$OL_category2 <- "OLG"

# ORF7b: codons 1-3 (overlap with ORF7a)  <-- also differs from SARS-CoV-2
interspecies_results[interspecies_results$product == "ORF7b" & interspecies_results$codon %in% 1:3, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF7b" & interspecies_results$codon %in% 1:3, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF7b" & interspecies_results$codon %in% 1:3, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF7b" & interspecies_results$codon %in% 1:3, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF7b" & interspecies_results$codon %in% 1:3, ]$OL_category2 <- "OLG"


### ELIMINATE DATA SUBSETS
#E: trim beginning because OLG with ORF3b SARSr taxa? Eliminate codons 1-13
interspecies_results[interspecies_results$product == "E" & interspecies_results$codon %in% 1:13, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "E" & interspecies_results$codon %in% 1:13, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "E" & interspecies_results$codon %in% 1:13, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "E" & interspecies_results$codon %in% 1:13, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "E" & interspecies_results$codon %in% 1:13, ]$OL_category2 <- "OLG"

#ORF6: trim last 3, codon 62-64
interspecies_results[interspecies_results$product == "ORF6" & interspecies_results$codon %in% 62:64, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF6" & interspecies_results$codon %in% 62:64, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF6" & interspecies_results$codon %in% 62:64, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF6" & interspecies_results$codon %in% 62:64, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF6" & interspecies_results$codon %in% 62:64, ]$OL_category2 <- "Non-OLG"

#9c: trim last 3 codons, codons 72-74
interspecies_results[interspecies_results$product == "ORF9c" & interspecies_results$codon %in% 72:74, ]$N_sites <- NA
interspecies_results[interspecies_results$product == "ORF9c" & interspecies_results$codon %in% 72:74, ]$S_sites <- NA
interspecies_results[interspecies_results$product == "ORF9c" & interspecies_results$codon %in% 72:74, ]$N_diffs <- NA
interspecies_results[interspecies_results$product == "ORF9c" & interspecies_results$codon %in% 72:74, ]$S_diffs <- NA
interspecies_results[interspecies_results$product == "ORF9c" & interspecies_results$codon %in% 72:74, ]$OL_category2 <- "OLG"

#ORF10: no filtering


###############################################################################
# BOOTSTRAP PROCESS

# MANUALLY ENSURE THIS IS THE CASE BEFOREHAND
interspecies_results$num_defined_seqs <- 6


####################################################################################################
### BOOTSTRAP EACH STATISTIC with ***JUKES-CANTOR CORRECTION***

############################################################################################################
# *BASIC* BOOTSTRAP FUNCTION (dN - dS) for CODON UNIT
dNdS_diff_boot_fun_JC <- function(codon_results, numerator, denominator, num_replicates, num_cpus) {
  
  # Function for dN
  dN_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    return(dN)
  }
  
  # Function for dN
  dS_function <- function(D, indices) {
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    return(dS)
  }
  
  # Function for dN - dS
  dN_m_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_m_dS <- dN - dS
    return(dN_m_dS)
  }
  
  # Function for dN/dS
  dN_over_dS_function <- function(D, indices) {
    dN <- sum(D[indices, paste0(numerator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(numerator, "_sites")], na.rm = TRUE)
    dN <- -3/4 * log(1 - (4/3) * dN)
    dS <- sum(D[indices, paste0(denominator, "_diffs")], na.rm = TRUE) / sum(D[indices, paste0(denominator, "_sites")], na.rm = TRUE)
    dS <- -3/4 * log(1 - (4/3) * dS)
    dN_over_dS <- dN / dS
    return(dN_over_dS)
  }
  
  # CREATE FUNCTION FOR dN/dS TO CALCULATE ITS SE
  
  dN <- sum(as.vector(codon_results[ , paste0(numerator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(numerator, "_sites")]), na.rm = TRUE)
  dN <- -3/4 * log(1 - (4/3) * dN)
  dS <- sum(as.vector(codon_results[ , paste0(denominator, "_diffs")]), na.rm = TRUE) / sum(as.vector(codon_results[ , paste0(denominator, "_sites")]), na.rm = TRUE)
  dS <- -3/4 * log(1 - (4/3) * dS)
  dNdS <- dN / dS
  
  # Run the BOOTSTRAPS
  # boot dN
  (boot_dN <- boot(data = codon_results, R = num_replicates, statistic = dN_function, parallel = 'multicore', ncpus = num_cpus))
  (dN <- boot_dN$t0)
  (boot_dN_SE <- sd(boot_dN$t))
  
  # boot dS
  (boot_dS <- boot(data = codon_results, R = num_replicates, statistic = dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dS <- boot_dS$t0)
  (boot_dS_SE <- sd(boot_dS$t))
  
  # boot dN - dS
  (boot_dN_m_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_m_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_m_dS <- boot_dN_m_dS$t0)
  (boot_dN_m_dS_SE <- sd(boot_dN_m_dS$t))
  (boot_dN_m_dS_Z <- dN_m_dS / boot_dN_m_dS_SE)
  (boot_dN_m_dS_P <- 2 * pnorm(-abs(boot_dN_m_dS_Z)))
  
  # boot dN/dS
  (boot_dN_over_dS <- boot(data = codon_results, R = num_replicates, statistic = dN_over_dS_function, parallel = 'multicore', ncpus = num_cpus))
  (dN_over_dS <- boot_dN_over_dS$t0)
  (boot_dN_over_dS_SE <- sd(boot_dN_over_dS$t))
  (boot_dN_over_dS_Z <- dN_over_dS / boot_dN_over_dS_SE)
  (boot_dN_over_dS_P <- 2 * pnorm(-abs(boot_dN_over_dS_Z)))
  
  ### NEW: ASL (acheived significance level)
  boot_dN_gt_dS_count <- sum(boot_dN_m_dS$t > 0) # 345
  boot_dN_eq_dS_count <- sum(boot_dN_m_dS$t == 0) # 0
  boot_dN_lt_dS_count <- sum(boot_dN_m_dS$t < 0) # 655
  ASL_dN_gt_dS_P <- boot_dN_lt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  ASL_dN_lt_dS_P <- boot_dN_gt_dS_count / (boot_dN_gt_dS_count + boot_dN_eq_dS_count + boot_dN_lt_dS_count)
  
  return(paste(num_replicates, dN, dS, dNdS, dN_m_dS, boot_dN_SE, boot_dS_SE, boot_dN_over_dS_SE, boot_dN_over_dS_P, 
               boot_dN_m_dS_SE, boot_dN_m_dS_P, 
               boot_dN_gt_dS_count, boot_dN_eq_dS_count, boot_dN_lt_dS_count, ASL_dN_gt_dS_P, ASL_dN_lt_dS_P,
               sep = "\t"))
}


############################################################################################################
### ANALYSIS VARIABLES
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6


############################################################################################################
### INITIALIZE DATA FRAME: interspecies
interspecies_results_bootstrap <- data.frame(gene_name = character(),
                                             OL_category = character(), # ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
                                             num_bootstraps = integer(),
                                             min_defined_codons = integer(),
                                             num_codons = integer(),
                                             N_sites = numeric(),
                                             S_sites = numeric(),
                                             N_diffs = numeric(),
                                             S_diffs = numeric(),
                                             num_replicates = integer(),
                                             dN = numeric(),
                                             dS = numeric(),
                                             dNdS = numeric(),
                                             dN_m_dS = numeric(),
                                             boot_dN_SE = numeric(),
                                             boot_dS_SE = numeric(),
                                             boot_dN_over_dS_SE = numeric(),
                                             boot_dN_over_dS_P = numeric(),
                                             boot_dN_m_dS_SE = numeric(),
                                             P_value = numeric(),
                                             boot_dN_gt_dS_count = integer(), 
                                             boot_dN_eq_dS_count = integer(), 
                                             boot_dN_lt_dS_count = integer(), 
                                             ASL_dN_gt_dS_P = numeric(), 
                                             ASL_dN_lt_dS_P = numeric())


### LOOP EACH DATA SUBSET
for (this_gene in sort(unique(interspecies_results$product))) {
  #this_gene <- 'ORF10'
  
  # ELIMINATE FOR SNPGENIE-ONLY <-- CHANGE
  for (this_OL_cat in sort(unique(interspecies_results$OL_category))) { # <-- CHANGE
    
    # Filter; ELIMINATE OL_category for SNPGENIE ONLY <-- CHANGE
    this_data <- filter(interspecies_results, product == this_gene, num_defined_seqs >= MIN_DEFINED_CODONS, 
                        OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG")) # OL_category == this_OL_cat, OL_category %in% c("OLG", "Non-OLG")
    
    if(nrow(this_data) >= 9) {
      # LEADING SUMMARY COLUMNS:
      N_sites <- sum(this_data$N_sites, na.rm = T)
      S_sites <- sum(this_data$S_sites, na.rm = T)
      N_diffs <- sum(this_data$N_diffs, na.rm = T)
      S_diffs <- sum(this_data$S_diffs, na.rm = T)
      
      summary_data <- paste(nrow(this_data),
                            N_sites, S_sites, 
                            N_diffs, S_diffs, 
                            sep = "\t")
      
      #if(PREPEND_TO_OUTPUT != '') {
      #  summary_data <- paste(PREPEND_TO_OUTPUT, summary_data, sep = "\t")
      #}
      
      # BOOTSTRAP THE ONE RATIO
      boot_dNdS <- dNdS_diff_boot_fun_JC(this_data, 'N', 'S', NBOOTSTRAPS, NCPUS)
      
      # RECORD HEADER
      boot_vector_names <- c('num_codons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs',
                             'num_replicates', 
                             'dN', 'dS', 'dNdS', 'dN_m_dS', 'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                             'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')
      
      boot_dNdS_vector <- unlist(str_split(string = paste(summary_data, boot_dNdS, sep = "\t"), pattern = "\t"))
      
      # Add names
      names(boot_dNdS_vector) <- boot_vector_names
      
      # Prepare additional rows
      interspecies_results_bootstrap_ADDITION <- data.frame(gene_name = this_gene,
                                                            OL_category = this_OL_cat, # ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
                                                            num_bootstraps = NBOOTSTRAPS,
                                                            min_defined_codons = MIN_DEFINED_CODONS,
                                                            num_codons = as.integer(boot_dNdS_vector['num_codons']),
                                                            N_sites = as.numeric(boot_dNdS_vector['N_sites']),
                                                            S_sites = as.numeric(boot_dNdS_vector['S_sites']),
                                                            N_diffs = as.numeric(boot_dNdS_vector['N_diffs']),
                                                            S_diffs = as.numeric(boot_dNdS_vector['S_diffs']),
                                                            num_replicates = as.integer(boot_dNdS_vector['num_replicates']),
                                                            dN = as.numeric(boot_dNdS_vector['dN']),
                                                            dS = as.numeric(boot_dNdS_vector['dS']),
                                                            dNdS = as.numeric(boot_dNdS_vector['dNdS']),
                                                            dN_m_dS = as.numeric(boot_dNdS_vector['dN_m_dS']),
                                                            boot_dN_SE = as.numeric(boot_dNdS_vector['boot_dN_SE']),
                                                            boot_dS_SE = as.numeric(boot_dNdS_vector['boot_dS_SE']),
                                                            boot_dN_over_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_over_dS_SE']),
                                                            boot_dN_over_dS_P = as.numeric(boot_dNdS_vector['boot_dN_over_dS_P']),
                                                            boot_dN_m_dS_SE = as.numeric(boot_dNdS_vector['boot_dN_m_dS_SE']),
                                                            P_value = as.numeric(boot_dNdS_vector['P_value']),
                                                            boot_dN_gt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_gt_dS_count']),
                                                            boot_dN_eq_dS_count = as.numeric(boot_dNdS_vector['boot_dN_eq_dS_count']),
                                                            boot_dN_lt_dS_count = as.numeric(boot_dNdS_vector['boot_dN_lt_dS_count']),
                                                            ASL_dN_gt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_gt_dS_P']),
                                                            ASL_dN_lt_dS_P = as.numeric(boot_dNdS_vector['ASL_dN_lt_dS_P']))
      
      # Add the 4 new rows to results
      interspecies_results_bootstrap <- rbind(interspecies_results_bootstrap, interspecies_results_bootstrap_ADDITION)
    }
  } # ELIMINATE FOR SNPGENIE-ONLY <-- CHANGE
}
# warnings: NAs introduced by coercion (OKAY)
# warnings()

names(interspecies_results_bootstrap) <- c('gene_name', 
                                           'OL_category', # ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
                                           'num_bootstraps', 'min_defined_codons', 'num_codons', 
                                           'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                           'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                           'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                           'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

### FILTER UNUSED GENES
interspecies_results_bootstrap <- filter(interspecies_results_bootstrap, ! gene_name %in% c("ORF3d2", "SiORF1", "ORF8b"),
                                         ! (gene_name == "E" & OL_category == "OLG")) # <-- CHANGE

### Manual 2-sided ASL P-value
interspecies_results_bootstrap$P_ALS <- NA
interspecies_results_bootstrap[! is.na(interspecies_results_bootstrap$ASL_dN_gt_dS_P) & interspecies_results_bootstrap$ASL_dN_gt_dS_P < interspecies_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * interspecies_results_bootstrap[! is.na(interspecies_results_bootstrap$ASL_dN_gt_dS_P) & interspecies_results_bootstrap$ASL_dN_gt_dS_P < interspecies_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_gt_dS_P
interspecies_results_bootstrap[! is.na(interspecies_results_bootstrap$ASL_dN_gt_dS_P) & interspecies_results_bootstrap$ASL_dN_gt_dS_P > interspecies_results_bootstrap$ASL_dN_lt_dS_P, ]$P_ALS <- 
  2 * interspecies_results_bootstrap[! is.na(interspecies_results_bootstrap$ASL_dN_gt_dS_P) & interspecies_results_bootstrap$ASL_dN_gt_dS_P > interspecies_results_bootstrap$ASL_dN_lt_dS_P, ]$ASL_dN_lt_dS_P
interspecies_results_bootstrap[! is.na(interspecies_results_bootstrap$ASL_dN_gt_dS_P) & interspecies_results_bootstrap$P_ALS == 0, ]$P_ALS <- 1 / NBOOTSTRAPS

### SAVE
#write_tsv(interspecies_results_bootstrap, "data/between-taxa_FINAL.tsv")




####################################################################################################
####################################################################################################
### COMBINE INTRAHOST/INTERHOST/INTERSPECIES
####################################################################################################
####################################################################################################

### RELOAD ALL DATASETS AND ADD LEVELS

####################
### INTRAHOST
intrahost_results_bootstrap <- read_tsv("data/within-host_FINAL.tsv")
intrahost_results_bootstrap_LONG <- intrahost_results_bootstrap %>%
  pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value")

intrahost_results_bootstrap_LONG$d_measure <- factor(intrahost_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
intrahost_results_bootstrap_LONG$d_SE_min <- NA
intrahost_results_bootstrap_LONG$d_SE_max <- NA
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_min <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value - intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_max <- 
  intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value + intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
intrahost_results_bootstrap_LONG[intrahost_results_bootstrap_LONG$d_SE_min < 0 & ! is.na(intrahost_results_bootstrap_LONG$d_SE_min), ]$d_SE_min <- 0 # OK if fail because none negative

# Error bars -- ELIMINATE OL_category FOR SNPGENIE ONLY <-- CHANGE
intrahost_error_bar_colors <- dplyr::select(intrahost_results_bootstrap_LONG, gene_name, OL_category, d_measure, d_value, d_SE_min, d_SE_max)
intrahost_error_bar_colors$this_color <- NA
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dN', ]$this_color <- 'pink'
intrahost_error_bar_colors[intrahost_error_bar_colors$d_measure == 'dS', ]$this_color <- 'lightblue'

# Analysis level
intrahost_results_bootstrap_LONG$level <- 'Intrahost'
intrahost_error_bar_colors$level <- "Intrahost"


####################
### INTERHOST
interhost_results_bootstrap <- read_tsv("data/between-host_FINAL.tsv")
interhost_results_bootstrap_LONG <- interhost_results_bootstrap %>%
  pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value")

interhost_results_bootstrap_LONG$d_measure <- factor(interhost_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
interhost_results_bootstrap_LONG$d_SE_min <- NA
interhost_results_bootstrap_LONG$d_SE_max <- NA
interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_min <- 
  interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value - interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_max <- 
  interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$d_value + interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_min <- 
  interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value - interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_max <- 
  interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$d_value + interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
interhost_results_bootstrap_LONG[interhost_results_bootstrap_LONG$d_SE_min < 0 & ! is.na(interhost_results_bootstrap_LONG$d_SE_min), ]$d_SE_min <- 0 # OK if fail because none negative

# Error bars -- ELIMINATE OL_category FOR SNPGENIE ONLY <-- CHANGE
interhost_error_bar_colors <- dplyr::select(interhost_results_bootstrap_LONG, gene_name, OL_category, d_measure, d_value, d_SE_min, d_SE_max)
interhost_error_bar_colors$this_color <- NA
interhost_error_bar_colors[interhost_error_bar_colors$d_measure == 'dN', ]$this_color <- 'lightblue'
interhost_error_bar_colors[interhost_error_bar_colors$d_measure == 'dS', ]$this_color <- 'pink'

# Analysis level
interhost_results_bootstrap_LONG$level <- 'Interhost'
interhost_error_bar_colors$level <- 'Interhost'


####################
### INTERSPECIES
interspecies_results_bootstrap <- read_tsv("data/between-taxa_FINAL.tsv")
(interspecies_results_bootstrap_LONG <- interspecies_results_bootstrap %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))

interspecies_results_bootstrap_LONG$d_measure <- factor(interspecies_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
interspecies_results_bootstrap_LONG$d_SE_min <- NA
interspecies_results_bootstrap_LONG$d_SE_max <- NA
interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_min <- 
  interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$d_value - interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$d_SE_max <- 
  interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$d_value + interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dN', ]$boot_dN_SE
interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_min <- 
  interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$d_value - interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$d_SE_max <- 
  interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$d_value + interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_measure == 'dS', ]$boot_dS_SE
interspecies_results_bootstrap_LONG[interspecies_results_bootstrap_LONG$d_SE_min < 0 & ! is.na(interspecies_results_bootstrap_LONG$d_SE_min), ]$d_SE_min <- 0 # OK if fail because none negative

# OLGENIE -- ELIMINATE OL_category FOR SNPGENIE ONLY <-- CHANGE
interspecies_error_bar_colors <- dplyr::select(interspecies_results_bootstrap_LONG, gene_name, OL_category, d_measure, d_value, d_SE_min, d_SE_max)
interspecies_error_bar_colors$this_color <- NA
interspecies_error_bar_colors[interspecies_error_bar_colors$d_measure == 'dN', ]$this_color <- 'lightblue'
interspecies_error_bar_colors[interspecies_error_bar_colors$d_measure == 'dS', ]$this_color <- 'pink'

# Analysis level
interspecies_results_bootstrap_LONG$level <- 'Interspecies'
interspecies_error_bar_colors$level <- "Interspecies"



####################################################################################################
# Combine data
allLevels_results_bootstrap_LONG <- rbind(intrahost_results_bootstrap_LONG, interhost_results_bootstrap_LONG, interspecies_results_bootstrap_LONG)
allLevels_error_bar_colors <- rbind(intrahost_error_bar_colors, interhost_error_bar_colors, interspecies_error_bar_colors)

# Add factor levels and ordering
allLevels_results_bootstrap_LONG$level <- factor(allLevels_results_bootstrap_LONG$level, levels = rev(c("Intrahost", "Interhost", "Interspecies")))
allLevels_error_bar_colors$level <- factor(allLevels_error_bar_colors$level, levels = rev(c("Intrahost", "Interhost", "Interspecies")))

# dN and dS
allLevels_results_bootstrap_LONG$d_measure <- factor(allLevels_results_bootstrap_LONG$d_measure, levels = c('dN', 'dS'))
allLevels_error_bar_colors$d_measure <- factor(allLevels_error_bar_colors$d_measure, levels = c('dN', 'dS'))

# ORDER AND NAME THE GENES
unique(allLevels_results_bootstrap_LONG$gene_name)
gene_ids_sorted <- c('ORF1ab', 'S', 'ORF3a', 'ORF3c', 'ORF3d', 'ORF3b', 'E', 'M', 'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'N', 'ORF9b', 'ORF9c', 'ORF10')
gene_names_sorted <- c('1ab', 'S', '3a', '3c', '3d', '3b', 'E', 'M', '6', '7a', '7b', '8', 'N', '9b', '9c', '10')
allLevels_results_bootstrap_LONG$gene_name <- factor(allLevels_results_bootstrap_LONG$gene_name, 
                                                     levels = gene_ids_sorted, 
                                                     labels = gene_names_sorted)
allLevels_error_bar_colors$gene_name <- factor(allLevels_error_bar_colors$gene_name, 
                                               levels = gene_ids_sorted, 
                                               labels = gene_names_sorted)

# ORDER OL CATEGORY -- ELIMINATE FOR SNPGENIE ONLY <-- CHANGE
allLevels_results_bootstrap_LONG$OL_category <- factor(allLevels_results_bootstrap_LONG$OL_category, 
                                                       levels = c('Non-OLG', 'OLG'),
                                                       labels = c("Non-OLG regions", "OLG regions"))
allLevels_error_bar_colors$OL_category <- factor(allLevels_error_bar_colors$OL_category, 
                                                 levels = c('Non-OLG', 'OLG'),
                                                 labels = c("Non-OLG regions", "OLG regions"))


### PREPARE PLOT
# Replace ORF3c value so it doesn't wash out everything else <-- CHANGE
sorted_dNdS <- sort(unique(allLevels_results_bootstrap_LONG$dNdS), decreasing = TRUE)
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$gene_name == "3c" & allLevels_results_bootstrap_LONG$level == "Interhost", ]$dNdS <- 
  1 / (1.25 * sorted_dNdS[1]) #1 / (3 * sorted_dNdS[1])

### normalized dNdS based on ratio (negative reciprocal for negative selection)
allLevels_results_bootstrap_LONG$dNdS_norm <- NA
# Intrahost
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Intrahost", ]$dNdS_norm <- 
  allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Intrahost", ]$dNdS
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Intrahost", ]$dNdS_norm <- 
  (-(1 / allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Intrahost", ]$dNdS))
# Interhost
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Interhost", ]$dNdS_norm <- 
  allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Interhost", ]$dNdS
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Interhost", ]$dNdS_norm <- 
  (-(1 / allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Interhost", ]$dNdS))
# Interspecies
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Interspecies", ]$dNdS_norm <- 
  allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS > 1 & allLevels_results_bootstrap_LONG$level == "Interspecies", ]$dNdS
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Interspecies", ]$dNdS_norm <- 
  (-(1 / allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS < 1 & allLevels_results_bootstrap_LONG$level == "Interspecies", ]$dNdS))
allLevels_results_bootstrap_LONG[allLevels_results_bootstrap_LONG$dNdS_norm == -Inf & allLevels_results_bootstrap_LONG$level == "Interspecies", ]$dNdS_norm <- 
  -(max(abs(allLevels_results_bootstrap_LONG[! is.infinite(allLevels_results_bootstrap_LONG$dNdS_norm), ]$dNdS_norm)))

### PLOT diversity chart ###
gene_names_allLevels <- gene_names_sorted # old: c('1ab', 'S', '3a', "3c", 'E', 'M', '6', '7a', '7b', '8', 'N', '9b', '9c', '10') 
pi_corr_factor <- 100
max_dNdS <- max(allLevels_results_bootstrap_LONG$dNdS)
max_d_SE_max <- max(allLevels_results_bootstrap_LONG$d_SE_max)


### PLOT
(allLevels_diversity_PLOT <- ggplot(data = filter(allLevels_results_bootstrap_LONG, gene_name %in% gene_names_allLevels),
                                    mapping = aes(x = gene_name, y = d_value * pi_corr_factor, group = d_measure)) +
    # Backdrop boxes based on which pi is higher
    geom_bar(data = filter(allLevels_results_bootstrap_LONG, d_measure == 'dN', gene_name %in% gene_names_allLevels), mapping = aes(y = Inf, fill = dNdS_norm), stat = 'identity') + # dN_m_dS_norm '#F7F7FF' BLUES: F7F7FF/F4F4FF/F2F2FF
    
    geom_errorbar(data = filter(allLevels_error_bar_colors, gene_name %in% gene_names_allLevels), 
                  mapping =  aes(ymin = d_SE_min * pi_corr_factor, ymax = d_SE_max * pi_corr_factor),
                  color = rev(rep(c('pink', 'lightblue'), 54)),
                  position = position_dodge(width = 0.5), width = 0, size = 1) +
    geom_point(data = filter(allLevels_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_min * pi_corr_factor),
               color = rev(rep(c('pink', 'lightblue'), 54)),
               position = position_dodge(width = 0.5), size = 0.29) + # size = 1
    geom_point(data = filter(allLevels_error_bar_colors, gene_name %in% gene_names_allLevels), 
               mapping =  aes(y = d_SE_max * pi_corr_factor),
               color = rev(rep(c('pink', 'lightblue'), 54)),
               position = position_dodge(width = 0.5), size = 0.29) +
    geom_point(stat = 'identity', position = position_dodge(width = 0.5), pch = 21, size = 2, stroke = 0, 
               fill = rev(rep(c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2]), 54))) +
    #ggtitle("Non-OLG method only") + # USE FOR SNPGENIE ONLY
    # ELIMINATE IF SNPGENIE ONLY
    facet_grid(level ~ OL_category, scales = 'free', space = 'free_x',) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.x = element_text(size = 9, face = "italic"),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9.5),
          panel.border = element_rect(),
          strip.text = element_text(size = 9.5),
          strip.background = element_blank()
          ) +
    xlab("") + 
    ylab(bquote('Differences per site ('*'x 10'^'-2'*')')) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), expand = expand_scale(mult = c(0, 0.1))) + 
    scale_fill_gradient2(low = brewer.pal(9, "Blues")[5], mid = 'white', high = brewer.pal(9, "Reds")[5], midpoint = 0)) 


