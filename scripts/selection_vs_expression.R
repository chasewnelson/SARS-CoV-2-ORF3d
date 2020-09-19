### SARS-CoV-2 selection vs. expression
library(tidyverse)
library(RColorBrewer)
library(scales)

### set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

### SELECTION VS. EXPRESSION

### Import riboseq results shown in Fig. 2B (source data): upstream, in-frame, out-of-frame
(riboseq_summary_30nt_HarrLTM <- read_tsv("data/expression_data_by_gene_frame.tsv"))
# 33 x 4

### IMPORT selection data
(selection_data <- read_tsv("data/selection_three_levels.tsv"))
# 108 x 32

# Process
unique(selection_data$gene_name) # "E"   "M"   "N"   "10"  "1ab" "3a"  "3b"  "3c"  "3d"  "6"   "7a"  "7b"  "8"   "9b"  "9c"  "S"     
selection_gene_IDs <- c("1ab", "S", "3a", "3c", "3d", "3b", "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10")
selection_gene_names <- c("ORF1ab", "S", "ORF3a", "ORF3c", "ORF3d", "ORF3b", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF9c", "ORF10")

selection_data$gene_name <- factor(selection_data$gene_name,
                                   levels = selection_gene_IDs,
                                   labels = selection_gene_names)

names(selection_data)[names(selection_data) == "gene_name"] <- "Gene"


### JOIN expression and selection data

# Proteomics expression
(selection_vs_expression_data <- left_join(selection_data, 
                                           dplyr::select(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak"), Gene, mean_iBAQ_pct), 
                                           by = "Gene"))
selection_vs_expression_data$mean_iBAQ_pct <- log10(selection_vs_expression_data$mean_iBAQ_pct)
names(selection_vs_expression_data)[names(selection_vs_expression_data) == "mean_iBAQ_pct"] <- "log10_iBAQ_pct"

# Upstream peak
(selection_vs_expression_data <- left_join(selection_vs_expression_data, 
                                           dplyr::select(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Upstream peak"), Gene, riboseq), 
                                           by = "Gene"))
names(selection_vs_expression_data)[names(selection_vs_expression_data) == "riboseq"] <- "Upstream peak"

# In-frame mean
(selection_vs_expression_data <- left_join(selection_vs_expression_data, 
                                           dplyr::select(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "In-frame mean"), Gene, riboseq), 
                                           by = "Gene")) 
names(selection_vs_expression_data)[names(selection_vs_expression_data) == "riboseq"] <- "In-frame mean"

# Out-of-frame mean
(selection_vs_expression_data <- left_join(selection_vs_expression_data, 
                                           dplyr::select(filter(riboseq_summary_30nt_HarrLTM, condition_combn == "Out-of-frame mean"), Gene, riboseq), 
                                           by = "Gene"))
names(selection_vs_expression_data)[names(selection_vs_expression_data) == "riboseq"] <- "Out-of-frame mean"


# filter to only those with proteomic data
selection_vs_expression_data <- selection_vs_expression_data[! is.na(selection_vs_expression_data$log10_iBAQ_pct), ] 


### make LONG
(selection_vs_expression_data_LONG <- selection_vs_expression_data %>%
    pivot_longer(cols = c('log10_iBAQ_pct', 'Upstream peak', 'In-frame mean', 'Out-of-frame mean'), names_to = "expression_method", values_to = "expression_value"))


# FACTOR selection level
selection_vs_expression_data_LONG$level = factor(selection_vs_expression_data_LONG$level, levels = c('Interspecies', 'Interhost', 'Intrahost'),
                              labels = c("Between-taxa", "Between-host", "Within-host"))

selection_vs_expression_data_LONG$gene_name <- factor(selection_vs_expression_data_LONG$Gene,
                                   levels = selection_gene_names,
                                   labels = selection_gene_names)

selection_vs_expression_data_LONG$OL_category <- factor(selection_vs_expression_data_LONG$OL_category,
                                     levels = c("Non-OLG regions","OLG regions"))



### Filter to non-OLGs but INCLUDE 9b
selection_data_subset <- filter(selection_vs_expression_data_LONG, OL_category == "Non-OLG regions", d_measure == "dN", ! is.na(dNdS))
selection_data_ORF9b <- filter(selection_vs_expression_data_LONG, Gene == "ORF9b", OL_category == "OLG regions", d_measure == "dN", ! is.na(dNdS))
selection_data_subset <- rbind(selection_data_subset, selection_data_ORF9b)

# Check same number of records for each gene
selection_data_subset %>%
  group_by(Gene) %>%
  summarise(
    count = n()
  ) # 11 genes, 12 records each (4 expression x 3 levels)

# FACTOR expression method
selection_data_subset$expression_method <- factor(selection_data_subset$expression_method,
                                                              levels = c("log10_iBAQ_pct", "Upstream peak", "In-frame mean", "Out-of-frame mean"),
                                                              labels = c("Mass spec", "Upstream peak", "In-frame mean", "Out-of-frame mean"))

### FACTOR GENE BY EXPRESSION (proteomics)
proteomics_data <- unique(dplyr::select(riboseq_summary_30nt_HarrLTM, Gene, mean_iBAQ_pct))
gene_order_proteomics <- arrange(proteomics_data, desc(mean_iBAQ_pct))$Gene
selection_data_subset$Gene <- factor(selection_data_subset$Gene, levels = gene_order_proteomics)


### PLOT Figure 5—figure supplement 1
(expression_vs_selection_PLOT <- ggplot(data = filter(selection_data_subset, expression_method %in% c("Mass spec", "Upstream peak", "In-frame mean")), mapping = aes(x = expression_value, y = dNdS, color = Gene, shape = Gene)) +
    geom_point(alpha = 0.8) +
    facet_grid(level ~ expression_method, scales = "free") +
    xlab("Expression level") +
    ylab("Selection strength") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
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
    scale_y_continuous(expand = expand_scale(mult = c(0.1, 0.1))) +
    scale_shape_manual(values = c(1:11)) +
    scale_color_manual(values = c(brewer.pal(11, 'RdYlBu')[0:5], 'yellow', brewer.pal(11, 'RdYlBu')[7:11])))



###############################################################################
### Correlation

### BETWEEN-TAXA

## Mass spec
cor.test(filter(selection_data_subset, expression_method == "Mass spec", level == "Between-taxa")$expression_value, 
         filter(selection_data_subset, expression_method == "Mass spec", level == "Between-taxa")$dNdS, 
         method = "spearman")
#S = 294, p-value = 0.313
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.3363636 

## Upstream peak
cor.test(filter(selection_data_subset, expression_method == "Upstream peak", level == "Between-taxa")$expression_value, 
         filter(selection_data_subset, expression_method == "Upstream peak", level == "Between-taxa")$dNdS, 
         method = "spearman")
#S = 276, p-value = 0.4512
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.2545455 

## In-frame mean
cor.test(filter(selection_data_subset, expression_method == "In-frame mean", level == "Between-taxa")$expression_value, 
         filter(selection_data_subset, expression_method == "In-frame mean", level == "Between-taxa")$dNdS, 
         method = "spearman")
#S = 218, p-value = 0.9892
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.009090909 


### BETWEEN-HOST

## Mass spec
cor.test(filter(selection_data_subset, expression_method == "Mass spec", level == "Between-host")$expression_value, 
         filter(selection_data_subset, expression_method == "Mass spec", level == "Between-host")$dNdS, 
         method = "spearman")
#S = 232, p-value = 0.8815
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.05454545 

## Upstream peak
cor.test(filter(selection_data_subset, expression_method == "Upstream peak", level == "Between-host")$expression_value, 
         filter(selection_data_subset, expression_method == "Upstream peak", level == "Between-host")$dNdS, 
         method = "spearman")
#S = 256, p-value = 0.6339
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.1636364 

## In-frame mean
cor.test(filter(selection_data_subset, expression_method == "In-frame mean", level == "Between-host")$expression_value, 
         filter(selection_data_subset, expression_method == "In-frame mean", level == "Between-host")$dNdS, 
         method = "spearman")
#S = 278, p-value = 0.4345
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.2636364 


### WITHIN-HOST

## Mass spec
cor.test(filter(selection_data_subset, expression_method == "Mass spec", level == "Within-host")$expression_value, 
         filter(selection_data_subset, expression_method == "Mass spec", level == "Within-host")$dNdS, 
         method = "spearman")
#S = 266, p-value = 0.5391
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.2090909 


## Upstream peak
cor.test(filter(selection_data_subset, expression_method == "Upstream peak", level == "Within-host")$expression_value, 
         filter(selection_data_subset, expression_method == "Upstream peak", level == "Within-host")$dNdS, 
         method = "spearman")
#S = 250, p-value = 0.6935
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.1363636 

## In-frame mean
cor.test(filter(selection_data_subset, expression_method == "In-frame mean", level == "Within-host")$expression_value, 
         filter(selection_data_subset, expression_method == "In-frame mean", level == "Within-host")$dNdS, 
         method = "spearman")
#S = 280, p-value = 0.4182
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.2727273 


