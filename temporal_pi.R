###############################################################################
### Plot nucleotide diversity (π) as a function of time
###############################################################################

library(tidyverse)
library(boot)
library(patchwork)
library(scales)
library(RColorBrewer)

# Set the working directory <-- CHANGE THIS
setwd("SARS-CoV-2-ORF3d")

# FIRST, run the python script and put all the sequences from the given time periods in separate folders:
#extract_seqs_by_timepoint.py gisaid_cov2020_acknowledgement_table.tsv SARS-CoV-2_ALN.fasta

# NEXT, run SNPGenie in each subdirectory:
#snpgenie_within_group.pl --fasta_file_name=SARS-CoV-2_ALN_0to14.fasta --gtf_file_name=SARS-CoV-2_ALN.gtf

# COMBINE all timepoint results, skipping headers:
# DO: cat SARS-CoV-2_ALN_*to*/within_group_codon_results.txt | grep -v variability > within_group_codon_results_TIMEPOINTS.txt

# Input the combined SNPGenie results <-- CHANGE THIS
codon_results <- read_tsv(file = "within_group_codon_results_TIMEPOINTS.txt",
                          col_names = c('file', 'product', 'codon', 'variability', 'comparisons', 'N_sites', 'S_sites', 'N_diffs', 'S_diffs'))

# Extract timepoint
codon_results$start_time <- codon_results$file
codon_results$start_time <- str_replace(string = codon_results$start_time,  pattern = ".+_(\\d+)to(\\d+).+", replacement = "\\1")
codon_results$start_time <- as.integer(codon_results$start_time)
codon_results$end_time <- codon_results$file
codon_results$end_time <- str_replace(string = codon_results$end_time,  pattern = ".+_(\\d+)to(\\d+).\\w+", replacement = "\\2")
codon_results$end_time <- as.integer(codon_results$end_time)
codon_results$midpoint_time <- (codon_results$start_time + codon_results$end_time) / 2
sort(unique(codon_results$start_time))
sort(unique(codon_results$midpoint_time))
sort(unique(codon_results$end_time))

# MANUALLY ADD REAL DATES
times_to_dates <- data.frame(timepoint = unique(sort(codon_results$midpoint_time)))
time0 = as.Date("2019-12-24")
times_to_dates$date <- time0 # time 0
times_to_dates[times_to_dates$timepoint == 7, ]$date <- time0 + 7
times_to_dates[times_to_dates$timepoint == 14, ]$date <- time0 + 14
times_to_dates[times_to_dates$timepoint == 21, ]$date <- time0 + 21
times_to_dates[times_to_dates$timepoint == 28, ]$date <- time0 + 28
times_to_dates[times_to_dates$timepoint == 35, ]$date <- time0 + 35
times_to_dates[times_to_dates$timepoint == 42, ]$date <- time0 + 42
times_to_dates[times_to_dates$timepoint == 49, ]$date <- time0 + 49
times_to_dates[times_to_dates$timepoint == 56, ]$date <- time0 + 56
times_to_dates[times_to_dates$timepoint == 63, ]$date <- time0 + 63
times_to_dates[times_to_dates$timepoint == 70, ]$date <- time0 + 70
times_to_dates[times_to_dates$timepoint == 77, ]$date <- time0 + 77
times_to_dates[times_to_dates$timepoint == 84, ]$date <- time0 + 84
times_to_dates[times_to_dates$timepoint == 91, ]$date <- time0 + 91


# MANUALLY ENSURE DEFINED SEQUENCES SATISFIED
codon_results$num_defined_seqs <- 6

# Nonredunant gene names
unique(codon_results$product) 
gene_set_NOL <- c("ORF1ab_1_NOL", "ORF1ab_2_NOL", "S", "ORF3a_NOL1", "ORF3a_NOL2", "E", "M","ORF6", "ORF7a_NOL", "ORF7b_NOL", "ORF8p",
                  "N_NOL1", "N_NOL2", "N_NOL3", "ORF10")
gene_set_OL <- c("ORF1ab_1_OL_ORF1ab_2", "ORF3a_OL_ORF3bp", "N_OL_ORF9b", "N_OL_ORF9c")

# Which genes to consider? CHANGE THIS <--
codon_results$region_type <- NA
codon_results[codon_results$product %in% gene_set_NOL, ]$region_type <- 'NOL'
codon_results[codon_results$product %in% gene_set_OL, ]$region_type <- 'OL'

### SUMMARIZE RESULTS BY GENE AND FRAME
(codon_results_summary <- codon_results %>% 
    group_by(region_type, midpoint_time) %>% # product
    summarise(
      N_sites = sum(N_sites, na.rm = TRUE),
      S_sites = sum(S_sites, na.rm = TRUE),
      N_diffs = sum(N_diffs, na.rm = TRUE),
      S_diffs = sum(S_diffs, na.rm = TRUE)
    ))

codon_results_summary$dN <- codon_results_summary$N_diffs / codon_results_summary$N_sites
codon_results_summary$dS <- codon_results_summary$S_diffs / codon_results_summary$S_sites
codon_results_summary$dNdS <- codon_results_summary$dN / codon_results_summary$dS

(codon_results_summary_LONG <- codon_results_summary %>%
    pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value"))


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
  
  # count results for ASL 
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


### ANALYSIS VARIABLES
MIN_DEFINED_CODONS <- 6
NBOOTSTRAPS <- 10000
NCPUS <- 6

### INITIALIZE DATA FRAME
bootstrap_gene_results <- data.frame(timepoint = integer(), 
                                     region_type = character(),
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


### LOOP EACH DATA SUBSET; takes a LONG TIME because it's a loop and this IS R
for(this_timepoint in sort(unique(codon_results$midpoint_time))) {
  #this_timepoint <- 7
  
  for(this_region_type in sort(unique(codon_results[! is.na(codon_results$region_type), ]$region_type))) { 
    
    # Filter by gene, frame, and minimum number of defined codons
    this_data <- filter(codon_results, midpoint_time == this_timepoint, region_type == this_region_type) 
    
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
      bootstrap_gene_results_ADDITION <- data.frame(timepoint = this_timepoint, # gene_name = this_gene,
                                                    region_type = this_region_type,
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
      bootstrap_gene_results <- rbind(bootstrap_gene_results, bootstrap_gene_results_ADDITION)
    }
  }
}
# 13 warnings: NAs introduced by coercion (OKAY)

names(bootstrap_gene_results) <- c('timepoint', 'region_type', 'num_bootstraps', 'min_defined_codons', 'num_codons', # 'gene_name'
                                   'N_sites', 'S_sites', 'N_diffs', 'S_diffs', 
                                   'num_replicates', 'dN', 'dS', 'dNdS', 'dN_m_dS', 
                                   'boot_dN_SE', 'boot_dS_SE', 'boot_dN_over_dS_SE', 'boot_dN_over_dS_P', 'boot_dN_m_dS_SE', 'P_value',
                                   'boot_dN_gt_dS_count', 'boot_dN_eq_dS_count', 'boot_dN_lt_dS_count', 'ASL_dN_gt_dS_P', 'ASL_dN_lt_dS_P')

# ADD DATE
bootstrap_gene_results <- left_join(x = bootstrap_gene_results, y = times_to_dates, by = "timepoint")

# SAVE <-- CHANGE THIS
#write_tsv(bootstrap_gene_results, "within_group_codon_results_TIMEPOINTS_BOOTSTRAP.tsv")


####################################################################################################
# RELOAD <-- CHANGE THIS
bootstrap_gene_results <- read_tsv("within_group_codon_results_TIMEPOINTS_BOOTSTRAP.tsv")

### SUMMARIZE RESULTS BY TIMEPOINT AND FRAME, NOW WITH BOOTSTRAPS
bootstrap_gene_results_LONG <- bootstrap_gene_results %>%
  pivot_longer(cols = c('dN', 'dS'), names_to = "d_measure", values_to = "d_value")
# so the 2 measures of d (dN and dS) are on different lines

# Get SE values for dN and dS on their new appropriate rows
bootstrap_gene_results_LONG$d_SE <- NA
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$d_SE <- bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$boot_dN_SE
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$d_SE <- bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$boot_dS_SE


####################################################################################################
### PLOT dN and dS by timepoint
bootstrap_gene_results_LONG$d_measure <- factor(bootstrap_gene_results_LONG$d_measure, levels = c('dN', 'dS'))

# Add error bars
bootstrap_gene_results_LONG$d_SE_min <- NA
bootstrap_gene_results_LONG$d_SE_max <- NA
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$d_SE_min <- 
  bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$d_value - bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$boot_dN_SE
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$d_SE_max <- 
  bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$d_value + bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dN', ]$boot_dN_SE
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$d_SE_min <- 
  bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$d_value - bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$boot_dS_SE
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$d_SE_max <- 
  bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$d_value + bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_measure == 'dS', ]$boot_dS_SE
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$d_SE_min < 0 & ! is.na(bootstrap_gene_results_LONG$d_SE_min), ]$d_SE_min <- 0 # GOOD IF FAIL because none negative

# significance for P-value
bootstrap_gene_results_LONG$significance <- NA
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$P_value < 0.05, ]$significance <- '*'
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$P_value < 0.01, ]$significance <- '**'
bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$P_value < 0.001, ]$significance <- '***'

### Now try adding the number of countries/locations for each window
# Need:
#(1) all the EPIs for that time;
#(2) extract their metadata
#(3) that's it

# CREATE AND LOAD a table with sequence IDs (col ID,  i.e., EPI_*) for each time window (col time_period) <-- CHANGE THIS
timepoint_seq_IDs <- read_tsv("timepoint_seq_IDs.txt")
timepoint_seq_IDs$surrogate_key <- 1:nrow(timepoint_seq_IDs)

# regex time out
timepoint_seq_IDs$start_time <- timepoint_seq_IDs$time_period
timepoint_seq_IDs$start_time <- str_replace(string = timepoint_seq_IDs$start_time,  pattern = "(\\d+)to(\\d+)", replacement = "\\1")
timepoint_seq_IDs$start_time <- as.integer(timepoint_seq_IDs$start_time)
timepoint_seq_IDs$end_time <- timepoint_seq_IDs$time_period
timepoint_seq_IDs$end_time <- str_replace(string = timepoint_seq_IDs$end_time,  pattern = "(\\d+)to(\\d+)", replacement = "\\2")
timepoint_seq_IDs$end_time <- as.integer(timepoint_seq_IDs$end_time)
timepoint_seq_IDs$timepoint <- (timepoint_seq_IDs$start_time + timepoint_seq_IDs$end_time) / 2

# regex ID out
timepoint_seq_IDs$ID <- str_replace(string = timepoint_seq_IDs$ID,  pattern = ".*EPI_", replacement = "EPI_")
nrow(timepoint_seq_IDs) # 7885

# add date
timepoint_seq_IDs <- left_join(x = timepoint_seq_IDs, y = times_to_dates, by = "timepoint")

# LOAD and JOIN GISAID acknowledgments table for metadata, primarily location <-- CHANGE THIS
GISAID_metadata <- read_tsv("gisaid_cov2020_acknowledgement_table.tsv")
names(GISAID_metadata)[names(GISAID_metadata) == "Accession ID"] <- "ID"
nrow(GISAID_metadata) # 5818

timepoint_seq_IDs <- left_join(x = timepoint_seq_IDs, y = GISAID_metadata, by = "ID")
nrow(timepoint_seq_IDs) # 7885

(timepoint_location_data <- timepoint_seq_IDs %>%
    group_by(timepoint, date, Location) %>%
    summarise(
      count = n()
    ))

# Examine unique locations, looks for duplicates
unique(timepoint_location_data$Location)

# Summarize COUNTS and ENTROPY
(timepoint_location_counts <- timepoint_seq_IDs %>%
    group_by(timepoint, Location) %>%
    summarise(
      location_n = n()
    ) %>%
    group_by(timepoint) %>%
    summarise(
      num_locations = n(),
      location_entropy = -sum(((location_n) / sum(location_n)) * log((location_n) / sum(location_n)))
    ))

# Extract COUNTRY
timepoint_seq_IDs$country <- timepoint_seq_IDs$Location
timepoint_seq_IDs$country <- str_replace(string = timepoint_seq_IDs$country, pattern = "([-\\w ]+) / ([-\\w ]+).*", replacement = "\\1 / \\2")
timepoint_seq_IDs$country <- str_replace(string = timepoint_seq_IDs$country, pattern = "\\s+$", replacement = "")

# Examine and clean
#unique(timepoint_seq_IDs$country)
timepoint_seq_IDs[timepoint_seq_IDs$country == "North America/USA/Virginia", ]$country <- "North America / USA"
timepoint_seq_IDs[timepoint_seq_IDs$country == "Asia/ Malaysia/Selangor", ]$country <- "Asia / Malaysia"
timepoint_seq_IDs[timepoint_seq_IDs$country == "North America /USA / Virginia", ]$country <- "North America / USA"
timepoint_seq_IDs[timepoint_seq_IDs$country == "Asia/China/Zhejiang/Hangzhou", ]$country <- "Asia / China"

# Entropy by country
(timepoint_country_counts <- timepoint_seq_IDs %>%
    group_by(timepoint, country) %>%
    summarise(
      country_n = n()
    ) %>%
    group_by(timepoint) %>%
    summarise(
      num_countries = n(),
      country_entropy = -sum(((country_n) / sum(country_n)) * log((country_n) / sum(country_n)))
    ))

# JOIN
timepoint_location_counts <- left_join(timepoint_location_counts, timepoint_country_counts, by = "timepoint")
timepoint_country_counts

# add date
(timepoint_location_counts <- left_join(x = timepoint_location_counts, y = times_to_dates, by = "timepoint"))
as.Date("2020-03-24") + 7

# Save final <-- CHANGE THIS
#write_tsv(timepoint_location_counts, "timepoint_location_country_counts.tsv")

# Reload final <-- CHANGE THIS
timepoint_location_counts <- read_tsv("timepoint_location_country_counts.tsv")

### PLOT again with unique locations
(max_location_count <- max(timepoint_location_counts$num_locations))
(max_location_entropy <- max(timepoint_location_counts$location_entropy))
(max_country_count <- max(timepoint_location_counts$num_countries))
(max_country_entropy <- max(timepoint_location_counts$country_entropy))
(max_piS <- max(bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == 'NOL', ]$d_SE_max))
timepoint_location_counts$d_measure <- "dN" # dummy
timepoint_location_counts$region_type <- "NOL" # dummy

# PLOT
(SARSCOV2_pi_TIMELINE_entropy <- ggplot(data = filter(bootstrap_gene_results_LONG, region_type == 'NOL'), mapping = aes(x = timepoint, y = d_value, color = d_measure)) +
    geom_line() +
    geom_ribbon(mapping = aes(ymin = d_value - d_SE, ymax = d_value + d_SE, fill = d_measure), alpha = 0.075, linetype = 0) +
    geom_text(mapping = aes(x = 91, y = 3.893860e-04), color = 'black', label = expression(italic('π')['S']), hjust = -1) +
    geom_text(mapping = aes(x = 91, y = 1.914569e-04), color = 'black', label = expression(italic('π')['N']),  hjust = -1) +
    
    # Location entropy
    geom_line(mapping = aes(x = timepoint, y = location_entropy * max_piS / (max_location_entropy)), inherit.aes = FALSE, linetype = 'dotted', data = timepoint_location_counts) +
    geom_text(mapping = aes(x = 63, y = 4.25  * max_piS / (max_location_entropy)), hjust = 1.05, vjust = 0, color = 'black', label = 'Location entropy', inherit.aes = FALSE) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          strip.background = element_blank()) +
    xlab("Days since first sample") + ylab("Nucleotide diversity") + 
    scale_color_manual(values = c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2]),
                       labels = c(expression(italic('π')['N']), expression(italic('π')['S']))))


### THREE PLOTS: 
## (1) piN and piS
## (2) piN/piS
## (3) LOCATIONS

## (1) piN and piS
(max_SE_value <- max(bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == "NOL", ]$d_value + bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == "NOL", ]$d_SE))
pi_normalization <- 10000
(SARSCOV2_pi_TIMELINE_plot1 <- ggplot(data = filter(bootstrap_gene_results_LONG, region_type == 'NOL'), mapping = aes(x = date, y = pi_normalization*d_value, color = d_measure)) +
    geom_line() +
    geom_ribbon(mapping = aes(ymin = pi_normalization*(d_value - d_SE), ymax = pi_normalization*(d_value + d_SE), fill = d_measure), alpha = 0.075, linetype = 0) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 9),
          strip.background = element_blank()) +
    scale_x_date(labels = date_format("%b %d"),
                 expand = expand_scale(mult = c(0, 0)),
                 limits = c(as.Date(time0 + 7), as.Date(time0 + 91 + 7)),
                 breaks = seq(as.Date(time0 + 7), as.Date(time0 + 91 + 7), by = "14 day")) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, pi_normalization*max_SE_value), expand = expand_scale(mult = c(0, 0.01))) +
    xlab("") + ylab("") + # ylab("Nucleotide diversity") +
    scale_color_manual(values = c(brewer.pal(9, 'Set1')[1], brewer.pal(9, 'Set1')[2])))

# (2) piN/piS
bootstrap_gene_results_LONG[is.na(bootstrap_gene_results_LONG$boot_dN_over_dS_SE), ]$boot_dN_over_dS_SE <- Inf
(SARSCOV2_pi_TIMELINE_plot2 <- ggplot(data = filter(bootstrap_gene_results_LONG, region_type == 'NOL', d_measure == "dN"), mapping = aes(x = date, y = dNdS)) +
    geom_line(mapping = aes(x = date, y = dNdS), color = 'black') +
    geom_ribbon(mapping = aes(ymin = dNdS - boot_dN_over_dS_SE, ymax = dNdS + boot_dN_over_dS_SE), alpha = 0.075, linetype = 0) +
    geom_abline(slope = 0, intercept = 1, linetype = "dashed", color = "grey") +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 9),
          strip.background = element_blank()) +
    xlab("") + ylab("") + #ylab('πN / πS') +
    scale_x_date(labels = date_format("%b %d"),
                 expand = expand_scale(mult = c(0, 0)),
                 limits = c(as.Date(time0 + 7), as.Date(time0 + 91 + 7)),
                 breaks = seq(as.Date(time0 + 7), as.Date(time0 + 91 + 7), by = "14 day")) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), limits = c(0, 1.1), expand = expand_scale(mult = c(0, 0.01))))

mean(bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == "NOL" & bootstrap_gene_results_LONG$d_measure == "dN", ]$dNdS)
sd(bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == "NOL" & bootstrap_gene_results_LONG$d_measure == "dN", ]$dNdS) / 
  sqrt(length(bootstrap_gene_results_LONG[bootstrap_gene_results_LONG$region_type == "NOL" & bootstrap_gene_results_LONG$d_measure == "dN", ]$dNdS) - 1)
# 0.4610758 +- 0.02985601 SE


# (3) LOCATION AND COUNTRY ENTROPIES
(SARSCOV2_pi_TIMELINE_plot3 <- ggplot(data = timepoint_location_counts, mapping = aes(x = date, y = location_entropy)) +
    geom_line(color = brewer.pal(8, 'Accent')[5]) + 
    geom_line(mapping = aes(x = date, y = country_entropy), data = timepoint_location_counts, inherit.aes = FALSE, color = brewer.pal(8, 'Accent')[1]) +
    theme_classic() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
          axis.text.y = element_text(size = 9),
          strip.background = element_blank()) +
    xlab("") + ylab("") +
    scale_x_date(labels = date_format("%b %d"),
                 expand = expand_scale(mult = c(0, 0)),
                 limits = c(as.Date(time0 + 7), as.Date(time0 + 91 + 7)),
                 breaks = seq(as.Date(time0 + 7), as.Date(time0 + 91 + 7), by = "14 day")) +
    scale_y_continuous(breaks = scales::pretty_breaks(3), expand = expand_scale(mult = c(0, 0.01))))

### HERE, run AF_trajectories to get whatever trajectory you want; can add at the end

# PLOT
SARSCOV2_pi_TIMELINE_plot1 / SARSCOV2_pi_TIMELINE_plot2 / SARSCOV2_pi_TIMELINE_plot3 #/ trajectory_plot_14404 / trajectory_plot_23399 / trajectory_plot_25559


