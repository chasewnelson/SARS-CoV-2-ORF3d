### SARS-CoV-2 ribosome profiling vs. proteomics comparison
library(tidyverse)
library(RColorBrewer)
library(scales)

### set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

####################################################################################################
### Ribo-seq UPSTREAM PEAKS vs. mass spectrometry iBAQ%

### load & wrangle Ribo-seq depth MAX UPSTREAM PEAKS
riboseq_peak_coverage_maxima <- read_tsv("data/riboseq_upstream_peaks.tsv")
(riboseq_peak_coverage_maxima <- filter(riboseq_peak_coverage_maxima, peak == "Peak"))

### load & wrangle iBAQ%
(proteome_data2 <- read_tsv("data/mass-spec_data.tsv"))
proteome_data2[proteome_data2$`Protein Names` == "ORF9b and ORF9b-short_form", ]$`Protein Names` <- "ORF9b"
proteomic_gene_IDs2 <- c("ORF1ab polyprotein", "surface glycoprotein", "ORF3a protein", "candidate-from-permutation-analysis",
                        "envelope protein", "membrane glycoprotein", 
                        "ORF6 protein", "ORF7a protein", "ORF7b", "ORF8 protein", 
                        "nucleocapsid phosphoprotein", "ORF9b", "ORF9b-short_form")
proteomic_gene_names2 <- c("ORF1ab", "S", "ORF3a", "ORF1ab-sas",
                          "E", "M", 
                          "ORF6", "ORF7a", "ORF7b", "ORF8", 
                          "N", "ORF9b", "ORF9b-short")

proteome_data2$gene <- factor(proteome_data2$`Protein Names`,
                                  levels = proteomic_gene_IDs2,
                                  labels = proteomic_gene_names2)

# Finally, mean of the two datasets
(proteome_data2_means <- proteome_data2 %>%
  group_by(gene) %>%
  summarise(
    mean_iBAQ_pct = mean(`iBAQ%`)
  )) 

### COMBINE the two
(expression_data <- left_join(proteome_data2_means, riboseq_peak_coverage_maxima, by = "gene"))

# Some e.g., ORF1ab-sas aren't shared, so eliminate
(expression_data <- filter(expression_data, ! is.na(reads_per_million)))

### FACTOR GENE BY EXPRESSION (proteomics)
gene_order_proteomics <- arrange(expression_data, desc(mean_iBAQ_pct))$gene
expression_data$gene <- factor(expression_data$gene, levels = gene_order_proteomics)

### CORRELATE
cor.test(expression_data$mean_iBAQ_pct, expression_data$reads_per_million, method = 'spearman') # 0.8909091 p-value = 0.0004284

### PLOT exploration
(expression_data_PLOT <- ggplot(data = expression_data, mapping = aes(x = log10(mean_iBAQ_pct), y = reads_per_million, color = gene, shape = gene)) +
    geom_point(size = 3) +
    xlab("Mass spec expression (log10[iBAQ%])") +
    ylab("Ribo-seq peak (reads per million mapped)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.75, "line"),
          legend.text = element_text(size = 9, face = "italic"),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(size = 9), 
          axis.title.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.border = element_rect(),
          strip.background = element_blank()) +
    scale_shape_manual(values = c(1:11)) +
    scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[0:5], 'yellow', brewer.pal(11, 'RdYlBu')[7:11]))) 



################################################################################
### Now instead compare to Ribo-seq MEAN DEPTH PER GENE, BY FRAME
mapped_reads_by_readlength_wCDS <- read_tsv("data/mapped_reads_by_readlength_wCDS.tsv")

# limit to reads of length 30, 5 hpi
mapped_reads_by_readlength_wCDS <- filter(mapped_reads_by_readlength_wCDS, read_length == 30, time == "5 hpi")

### vvv CHANGE THIS vvv ### <-- CHANGE
# in-frame only! #
mapped_reads_by_readlength_wCDS <- filter(mapped_reads_by_readlength_wCDS, frame == "Codon pos 1")
###

# add condition_combn
mapped_reads_by_readlength_wCDS$condition_combn <- paste0(mapped_reads_by_readlength_wCDS$treatment, " / ", mapped_reads_by_readlength_wCDS$time)

# SUMMARIZE read_count by time, treatment, sample, gene, region_type
(mapped_reads_by_readlength_wCDS_SUMMARY <- mapped_reads_by_readlength_wCDS %>%
  group_by(condition_combn, sample, gene, region_type) %>%
  summarize(
    mean_depth = mean(read_count)
  ))

# SUMMARIZE FURTHER by taking the means of the samples
(mapped_reads_by_readlength_wCDS_SUMMARY <- mapped_reads_by_readlength_wCDS_SUMMARY %>%
  group_by(condition_combn, gene, region_type) %>%
  summarize(
    mean_depth = mean(mean_depth)
  )) 

### COMBINE WITH EARLIER EXPRESSION DATA

# set aside proteomics for joining
expression_data_PROTEOMICS <- dplyr::select(expression_data, gene, mean_iBAQ_pct)

# extract upstream peaks
expression_data_PEAKS <- dplyr::select(expression_data, gene, reads_per_million)
unique(expression_data_PEAKS$gene)
names(expression_data_PEAKS) <- c("gene", "riboseq")
expression_data_PEAKS$condition_combn <- "Upstream peak"

# excise treatment means
riboseq_summary_allFrames_PROCESSED <- filter(mapped_reads_by_readlength_wCDS_SUMMARY, region_type == "gene")
riboseq_summary_allFrames_PROCESSED <- rbind(riboseq_summary_allFrames_PROCESSED, filter(mapped_reads_by_readlength_wCDS_SUMMARY, gene == "ORF9b")) # add ORF9b
riboseq_summary_allFrames_PROCESSED[riboseq_summary_allFrames_PROCESSED$gene == "ORF1ab_1", ]$gene <- "ORF1ab" # rename first half of ORF1ab
unique(riboseq_summary_allFrames_PROCESSED$gene)
riboseq_summary_allFrames_PROCESSED <- dplyr::select(riboseq_summary_allFrames_PROCESSED, gene, mean_depth, condition_combn)
names(riboseq_summary_allFrames_PROCESSED) <- c("gene", "riboseq", "condition_combn")

### COMBINE ###
expression_data_COMBINED <- rbind(as.data.frame(expression_data_PEAKS), as.data.frame(riboseq_summary_allFrames_PROCESSED))

# join
expression_data_COMBINED <- left_join(expression_data_COMBINED, expression_data_PROTEOMICS, by = "gene")

# remove NA
expression_data_COMBINED <- filter(expression_data_COMBINED, ! is.na(mean_iBAQ_pct))


### FACTOR GENE BY EXPRESSION (proteomics)
gene_order_proteomics <- arrange(expression_data, desc(mean_iBAQ_pct))$gene
expression_data_COMBINED$gene <- factor(expression_data_COMBINED$gene, levels = gene_order_proteomics)

### FACTOR CONDITION
expression_data_COMBINED$condition_combn <- factor(expression_data_COMBINED$condition_combn,
                                                   levels = c("Upstream peak", "CHX / 5 hpi", "Harr / 5 hpi", "LTM / 5 hpi", "mRNA / 5 hpi"))

##############################
### CORRELATE after limiting to Codon pos 1***
cor.test(filter(expression_data_COMBINED, condition_combn == "Upstream peak")$mean_iBAQ_pct, 
         filter(expression_data_COMBINED, condition_combn == "Upstream peak")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(expression_data_COMBINED, condition_combn == "Upstream peak")$mean_iBAQ_pct and filter(expression_data_COMBINED, condition_combn == "Upstream peak")$riboseq
#S = 24, p-value = 0.0004284
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8909091 

cor.test(filter(expression_data_COMBINED, condition_combn == "CHX / 5 hpi")$mean_iBAQ_pct, 
         filter(expression_data_COMBINED, condition_combn == "CHX / 5 hpi")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(expression_data_COMBINED, condition_combn == "CHX / 5 hpi")$mean_iBAQ_pct and filter(expression_data_COMBINED, condition_combn == "CHX / 5 hpi")$riboseq
#S = 82, p-value = 0.04399
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.6272727 

cor.test(filter(expression_data_COMBINED, condition_combn == "Harr / 5 hpi")$mean_iBAQ_pct, 
         filter(expression_data_COMBINED, condition_combn == "Harr / 5 hpi")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(expression_data_COMBINED, condition_combn == "Harr / 5 hpi")$mean_iBAQ_pct and filter(expression_data_COMBINED, condition_combn == "Harr / 5 hpi")$riboseq
#S = 88, p-value = 0.05618
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.6 

cor.test(filter(expression_data_COMBINED, condition_combn == "LTM / 5 hpi")$mean_iBAQ_pct, 
         filter(expression_data_COMBINED, condition_combn == "LTM / 5 hpi")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(expression_data_COMBINED, condition_combn == "LTM / 5 hpi")$mean_iBAQ_pct and filter(expression_data_COMBINED, condition_combn == "LTM / 5 hpi")$riboseq
#S = 86, p-value = 0.05188
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.6090909 

cor.test(filter(expression_data_COMBINED, condition_combn == "mRNA / 5 hpi")$mean_iBAQ_pct, 
         filter(expression_data_COMBINED, condition_combn == "mRNA / 5 hpi")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(expression_data_COMBINED, condition_combn == "mRNA / 5 hpi")$mean_iBAQ_pct and filter(expression_data_COMBINED, condition_combn == "mRNA / 5 hpi")$riboseq
#S = 156, p-value = 0.3864
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.2909091 


### PLOT all conditions
(expression_data_COMBINED_PLOT <- ggplot(data = expression_data_COMBINED, mapping = aes(x = log10(mean_iBAQ_pct), y = riboseq, color = gene, shape = gene)) +
    geom_point(size = 2) + 
    geom_smooth(method = 'lm', se = FALSE, size = 0.4, color = "darkgrey") + 
    xlab("Mass spec expression (log10[iBAQ%])") +
    ylab("Ribo-seq coverage (reads per million mapped)") +
    facet_grid(condition_combn ~ ., scales = "free_y") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.75, "line"),
          legend.text = element_text(size = 9, face = "italic"),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(size = 9), 
          axis.title.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.border = element_rect(),
          strip.background = element_blank()) +
    scale_shape_manual(values = c(1:11)) +
    scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[0:5], 'yellow', brewer.pal(11, 'RdYlBu')[7:11]))) 



### FIGURES AND SUPPLEMENT, ALL SITES, ALL CODON POSITIONS

### RELOAD load & wrangle Ribo-seq MEAN DEPTH PER GENE
mapped_reads_by_readlength_wCDS <- read_tsv("data/mapped_reads_by_readlength_wCDS.tsv")

# limit to reads of length 30, 5 hpi
mapped_reads_by_readlength_wCDS <- filter(mapped_reads_by_readlength_wCDS, read_length == 30, time == "5 hpi") # 61,974 x 34

# add condition_combn
mapped_reads_by_readlength_wCDS$condition_combn <- paste0(mapped_reads_by_readlength_wCDS$treatment, " / ", mapped_reads_by_readlength_wCDS$time)

### SUMMARY ALL FRAME
# SUMMARIZE read_count by time, treatment, sample, gene, region_type
(riboseq_summary_allFrame <- mapped_reads_by_readlength_wCDS %>%
    group_by(condition_combn, sample, gene, frame, region_type) %>%
    summarize(
      mean_depth = mean(read_count)
    )) 

# SUMMARIZE FURTHER by taking the means of the samples
(riboseq_summary_allFrame <- riboseq_summary_allFrame %>%
    group_by(condition_combn, gene, region_type) %>%
    summarize(
      mean_depth = mean(mean_depth)
    )) # 84 x 4

riboseq_summary_allFrame$frame <- "All frames"
(riboseq_summary_allFrame <- dplyr::select(riboseq_summary_allFrame, condition_combn, gene, frame, region_type, mean_depth))
# 84 x 5

### SUMMARY BY FRAME
# SUMMARIZE read_count by time, treatment, sample, gene, frame, region_type
(riboseq_summary_byFrame <- mapped_reads_by_readlength_wCDS %>%
    group_by(condition_combn, sample, gene, frame, region_type) %>%
    summarize(
      mean_depth = mean(read_count)
    )) 

# SUMMARIZE FURTHER by taking the means of the samples
(riboseq_summary_byFrame <- riboseq_summary_byFrame %>%
    group_by(condition_combn, gene, frame, region_type) %>%
    summarize(
      mean_depth = mean(mean_depth)
    )) # 236 x 5

# combine
riboseq_summary_allFrames <- rbind(riboseq_summary_allFrame, riboseq_summary_byFrame)


### COMBINE WITH EARLIER EXPRESSION DATA

# set aside proteomics for joining
expression_data_PROTEOMICS <- dplyr::select(expression_data, gene, mean_iBAQ_pct)

# excise treatment means
riboseq_summary_allFrames_PROCESSED <- filter(riboseq_summary_allFrames, region_type == "gene")
riboseq_summary_allFrames_PROCESSED <- rbind(riboseq_summary_allFrames_PROCESSED, filter(riboseq_summary_allFrames, gene == "ORF9b")) # add ORF9b
riboseq_summary_allFrames_PROCESSED[riboseq_summary_allFrames_PROCESSED$gene == "ORF1ab_1", ]$gene <- "ORF1ab" # rename first half of ORF1ab
unique(riboseq_summary_allFrames_PROCESSED$gene)

# join
riboseq_summary_allFrames_PROCESSED <- left_join(riboseq_summary_allFrames_PROCESSED, expression_data_PROTEOMICS, by = "gene")

# remove NA
riboseq_summary_allFrames_PROCESSED <- filter(riboseq_summary_allFrames_PROCESSED, ! is.na(mean_iBAQ_pct))

### FACTOR GENE BY EXPRESSION (proteomics)
gene_order_proteomics <- arrange(expression_data, desc(mean_iBAQ_pct))$gene
riboseq_summary_allFrames_PROCESSED$gene <- factor(riboseq_summary_allFrames_PROCESSED$gene, levels = gene_order_proteomics)

### FACTOR CONDITION
riboseq_summary_allFrames_PROCESSED$condition_combn <- factor(riboseq_summary_allFrames_PROCESSED$condition_combn,
                                                   levels = c("CHX / 5 hpi", "Harr / 5 hpi", "LTM / 5 hpi", "mRNA / 5 hpi"))

### FACTOR FRAME
riboseq_summary_allFrames_PROCESSED$frame <- factor(riboseq_summary_allFrames_PROCESSED$frame,
                                                              levels = c("All frames", "Codon pos 1", "Codon pos 2", "Codon pos 3"))

# Rename Gene
names(riboseq_summary_allFrames_PROCESSED)[names(riboseq_summary_allFrames_PROCESSED) == "gene"] <- "Gene"

# PLOT Figure 2—figure supplement 2
(riboseq_summary_allFrames_PROCESSED_PLOT <- ggplot(data = riboseq_summary_allFrames_PROCESSED, mapping = aes(x = log10(mean_iBAQ_pct), y = mean_depth, color = Gene, shape = Gene)) +
    geom_point(size = 2, stroke = 0.8) + 
    geom_smooth(method = 'lm', se = FALSE, size = 0.4, color = "darkgrey") + 
    xlab("log10(iBAQ%) (mass spectrometry)") +
    ylab("Reads per million mapped (ribosomal profiling))") +
    facet_grid(condition_combn ~ frame, scales = "free_y") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.75, "line"),
          legend.text = element_text(size = 9, face = "italic"),
          strip.text = element_text(size = 9),
          axis.text.x = element_text(size = 9), 
          axis.title.x = element_text(size = 9),
          axis.text.y = element_text(size = 9),
          axis.title.y = element_text(size = 9),
          panel.border = element_rect(),
          strip.background = element_blank()) +
    scale_x_continuous(expand = expand_scale(mult = c(0.075, 0.075))) +
    scale_y_continuous(expand = expand_scale(mult = c(0.075, 0.075))) +
    scale_shape_manual(values = c(1:10, 13)) +
    scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[0:5], brewer.pal(11, 'BrBG')[5], brewer.pal(11, 'RdYlBu')[7:11])))


### CORRELATIONS

## ALL FRAMES ###
cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "All frames")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "All frames")$mean_depth, 
         method = "spearman")
#p-value = 0.08745, r=0.5454545 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "All frames")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "All frames")$mean_depth, 
         method = "spearman")
#p-value = 0.2139, r=0.4090909

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "All frames")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "All frames")$mean_depth, 
         method = "spearman")
#p-value = 0.1069, r=0.5181818

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "All frames")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "All frames")$mean_depth, 
         method = "spearman")
#p-value = 0.4021, r=0.2818182 


## CODON POS 1 ##
cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 1")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 1")$mean_depth, 
         method = "spearman")
#p-value = 0.04399, r=0.6272727

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 1")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 1")$mean_depth, 
         method = "spearman")
#p-value = 0.05618, r=0.6 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 1")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 1")$mean_depth, 
         method = "spearman")
#p-value = 0.05188, r=0.6090909 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 1")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 1")$mean_depth, 
         method = "spearman")
#p-value = 0.3864, r=0.2909091 


## CODON POS 2 ##
cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 2")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 2")$mean_depth, 
         method = "spearman")
#p-value = 0.08155, r=0.5545455 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 2")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 2")$mean_depth, 
         method = "spearman")
#p-value = 0.2484, r=0.3818182 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 2")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 2")$mean_depth, 
         method = "spearman")
#p-value = 0.2365, r=0.3909091 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 2")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 2")$mean_depth, 
         method = "spearman")
#p-value = 0.4021, r=0.2818182 


## CODON POS 3 ##
cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 3")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "CHX / 5 hpi", frame == "Codon pos 3")$mean_depth, 
         method = "spearman")
#p-value = 0.1456, r=0.4727273 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 3")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "Harr / 5 hpi", frame == "Codon pos 3")$mean_depth, 
         method = "spearman")
#p-value = 0.1456, r=0.4727273 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 3")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "LTM / 5 hpi", frame == "Codon pos 3")$mean_depth, 
         method = "spearman")
#p-value = 0.2994, r=0.3454545 

cor.test(filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 3")$mean_iBAQ_pct, 
         filter(riboseq_summary_allFrames_PROCESSED, condition_combn == "mRNA / 5 hpi", frame == "Codon pos 3")$mean_depth, 
         method = "spearman")
#p-value = 0.313, r=0.3363636 





#########################################################################################################
### FINAL FIGURE 2B: (1) upstream peak, (2) In frame, and (3) out-of-frame, average of Harr and LTM

### RELOAD & wrangle Ribo-seq MEAN DEPTH PER GENE
mapped_reads_by_readlength_wCDS <- read_tsv("data/mapped_reads_by_readlength_wCDS.tsv")

# limit to reads of length 30, 5 hpi
mapped_reads_by_readlength_wCDS <- filter(mapped_reads_by_readlength_wCDS, read_length == 30, time == "5 hpi")

# add condition_combn
mapped_reads_by_readlength_wCDS$condition_combn <- paste0(mapped_reads_by_readlength_wCDS$treatment, " / ", mapped_reads_by_readlength_wCDS$time)


### SUMMARY IN-FRAME
# SUMMARIZE read_count by time, treatment, sample, gene, region_type
(riboseq_summary_30nt_allCond_inFrame <- filter(mapped_reads_by_readlength_wCDS, frame == "Codon pos 1") %>%
    group_by(condition_combn, sample, gene, region_type) %>% 
    summarize(
      mean_depth = mean(read_count)
    )) 

# SUMMARIZE FURTHER by taking the means of the samples just calculated
(riboseq_summary_30nt_allCond_inFrame <- riboseq_summary_30nt_allCond_inFrame %>%
    group_by(condition_combn, gene, region_type) %>%
    summarize(
      mean_depth = mean(mean_depth)
    )) 

riboseq_summary_30nt_allCond_inFrame$frame <- "In-frame"
(riboseq_summary_30nt_allCond_inFrame <- dplyr::select(riboseq_summary_30nt_allCond_inFrame, condition_combn, gene, frame, region_type, mean_depth))


### SUMMARY OUT-OF-FRAME
# SUMMARIZE read_count by time, treatment, sample, gene, region_type
(riboseq_summary_30nt_allCond_outOfFrame <- filter(mapped_reads_by_readlength_wCDS, frame != "Codon pos 1") %>%
    group_by(condition_combn, sample, gene, region_type) %>%
    summarize(
      mean_depth = mean(read_count)
    )) 

# SUMMARIZE FURTHER by taking the means of the samples just calculated
(riboseq_summary_30nt_allCond_outOfFrame <- riboseq_summary_30nt_allCond_outOfFrame %>%
    group_by(condition_combn, gene, region_type) %>%
    summarize(
      mean_depth = mean(mean_depth)
    )) 

riboseq_summary_30nt_allCond_outOfFrame$frame <- "Out-of-frame"
(riboseq_summary_30nt_allCond_outOfFrame <- dplyr::select(riboseq_summary_30nt_allCond_outOfFrame, condition_combn, gene, frame, region_type, mean_depth))

# combine
riboseq_summary_30nt_allCond <- rbind(riboseq_summary_30nt_allCond_inFrame, riboseq_summary_30nt_allCond_outOfFrame)


### Prepare IN-FRAME MEAN of Harr of and LTM for the Ribo-seq
(riboseq_summary_30nt_HarrLTM <- filter(riboseq_summary_30nt_allCond, frame == "In-frame", condition_combn %in% c("Harr / 5 hpi", "LTM / 5 hpi")) %>%
    group_by(gene) %>%
    summarise(
      riboseq = mean(mean_depth)
    ))
riboseq_summary_30nt_HarrLTM$condition_combn <- "In-frame mean"

### Prepare OUT-OF-FRAME MEAN of Harr of and LTM for the Ribo-seq
(riboseq_summary_30nt_HarrLTM_outOfFrame <- filter(riboseq_summary_30nt_allCond, frame == "Out-of-frame", condition_combn %in% c("Harr / 5 hpi", "LTM / 5 hpi")) %>%
    group_by(gene) %>%
    summarise(
      riboseq = mean(mean_depth)
    ))
riboseq_summary_30nt_HarrLTM_outOfFrame$condition_combn <- "Out-of-frame mean"


# combine
riboseq_summary_30nt_HarrLTM <- rbind(riboseq_summary_30nt_HarrLTM, riboseq_summary_30nt_HarrLTM_outOfFrame)


# set aside proteomics for joining
expression_data_PROTEOMICS <- dplyr::select(expression_data, gene, mean_iBAQ_pct)

# prepare for join
riboseq_summary_30nt_HarrLTM[riboseq_summary_30nt_HarrLTM$gene == "ORF1ab_1", ]$gene <- "ORF1ab" # rename first half of ORF1ab

# join
riboseq_summary_30nt_HarrLTM <- left_join(riboseq_summary_30nt_HarrLTM, expression_data_PROTEOMICS, by = "gene")

# remove NA
riboseq_summary_30nt_HarrLTM <- filter(riboseq_summary_30nt_HarrLTM, ! is.na(mean_iBAQ_pct))

# Import UPSTREAM PEAK DATA, which is also based only on 29-31nt reads and Harr/LTM means
(riboseq_upstreamPeak_data <- read_tsv("data/expression_data_by_gene.tsv")) # LTM and Harr mean
riboseq_upstreamPeak_data <- dplyr::select(riboseq_upstreamPeak_data, -method)

# combine
riboseq_summary_30nt_HarrLTM <- rbind(riboseq_summary_30nt_HarrLTM, riboseq_upstreamPeak_data)

### FACTOR GENE BY EXPRESSION (proteomics)
gene_order_proteomics <- arrange(expression_data, desc(mean_iBAQ_pct))$gene
riboseq_summary_30nt_HarrLTM$gene <- factor(riboseq_summary_30nt_HarrLTM$gene, levels = gene_order_proteomics)

### FACTOR CONDITION
#unique(riboseq_summary_30nt_HarrLTM$condition_combn)
riboseq_summary_30nt_HarrLTM$condition_combn <- factor(riboseq_summary_30nt_HarrLTM$condition_combn,
                                                       levels = c("Upstream peak", "In-frame mean", "Out-of-frame mean"))

# Rename Gene
names(riboseq_summary_30nt_HarrLTM)[names(riboseq_summary_30nt_HarrLTM) == "gene"] <- "Gene"

### PLOT main Figure 2B
(riboseq_summary_30nt_HarrLTM_PLOT <- ggplot(data = riboseq_summary_30nt_HarrLTM, mapping = aes(x = log10(mean_iBAQ_pct), y = riboseq, color = Gene, shape = Gene)) +
    geom_point(size = 2, stroke = 1) + 
    geom_smooth(method = 'lm', se = FALSE, size = 0.4, color = "darkgrey") + 
    xlab("Mass spec expression (log10[iBAQ%])") +
    ylab("Ribo-seq coverage (reads per million mapped)") +
    facet_grid(condition_combn ~ ., scales = "free_y") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          legend.title = element_text(size = rel(0.9)),
          legend.key.size = unit(0.75, "line"),
          legend.text = element_text(size = rel(0.9), face = "italic"),
          strip.text = element_text(size = rel(0.9)),
          axis.text = element_text(size = rel(0.9)), 
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(0.9)),
          panel.border = element_rect(),
          strip.background = element_blank()) +
    scale_shape_manual(values = c(1:10, 13)) +
    scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[0:5], brewer.pal(11, 'BrBG')[5], brewer.pal(11, 'RdYlBu')[7:11])))


### CORRELATIONS FOR THESE

cor.test(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak")$mean_iBAQ_pct, 
         filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#
#data:  filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak")$mean_iBAQ_pct and filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak")$riboseq
#S = 24, p-value = 0.0004284
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8909091 


cor.test(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "In-frame mean")$mean_iBAQ_pct, 
         filter(riboseq_summary_30nt_HarrLTM, condition_combn == "In-frame mean")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#
#data:  filter(riboseq_summary_30nt_HarrLTM, condition_combn == "In-frame mean")$mean_iBAQ_pct and filter(riboseq_summary_30nt_HarrLTM, condition_combn == "In-frame mean")$riboseq
#S = 88, p-value = 0.05618
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#rho 
#0.6 

cor.test(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Out-of-frame mean")$mean_iBAQ_pct, 
         filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Out-of-frame mean")$riboseq, 
         method = "spearman")
#Spearman's rank correlation rho
#data:  filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Out-of-frame mean")$mean_iBAQ_pct and filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Out-of-frame mean")$riboseq
#S = 134, p-value = 0.2365
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.3909091  




