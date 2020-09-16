### SARS-CoV-2 ORF length visualization

library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(patchwork)

# Set the working directory <-- CHANGE THIS
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d")

### INPUT
(frameshift_data <- read_tsv("frameshift_results.txt"))
#140 x 10

### Add columns
frameshift_data$length <- frameshift_data$`genome-co-ord2` - frameshift_data$`genome-co-ord1` + 1
frameshift_data$center <- (frameshift_data$`genome-co-ord2` + frameshift_data$`genome-co-ord1`) / 2
frameshift_data$`reference-frame-gene`

### PLOT exploration
ggplot(data = filter(frameshift_data), mapping = aes(x = center, y = length, color = `pvalue-permutation`)) + 
  geom_point() +
  geom_text_repel(mapping = aes(label = ifelse(`pvalue-permutation` < 0.05, `reference-frame-gene`, "")), color = 'black') +
  xlab("Site") +
  ylab("ORF length (nt)") +
  theme_classic() +
  theme(
    panel.background = element_blank()
  ) +
  scale_color_gradient(low = "red", high = brewer.pal(9, "Greys")[3], na.value = "white", name = "P-value")

### PLOT exploration 2

# FACTOR region type (known vs. not known)
frameshift_data$Type <- factor(frameshift_data$Type, levels = c("Non-OLG", "Overlapping"),
                               labels = c("Known gene", "Overlapping ORF"))

# PLOT
(ORF_length_pvalue_PLOT <- ggplot(data = filter(frameshift_data), mapping = aes(x = center, y = log(1/`pvalue-permutation`), color = Type)) + 
    geom_point() +
    geom_text_repel(mapping = aes(label = ifelse(`pvalue-permutation` < 0.3, `reference-frame-gene`, ""), color = Type), ylim = c(1.5, Inf)) + 
    xlab("Site") +
    ylab("Significance (log[1/P])") +
    theme_classic() +
    theme(legend.position = c(0.25, 0.75),
      legend.justification = c(0.5, 0.5),
      legend.title = element_blank(),
      legend.background = element_rect(color = "black"),
      panel.background = element_blank()
    ) +
    scale_color_manual(values = c(brewer.pal(9, "Greens")[5], 'black')))


### PLOT HEATMAP ("the bee")

# Prepare table
genome_length <- 29903
ORF_length_heatmap_data <- data.frame(site = rep(1:(genome_length), 3))
ORF_length_heatmap_data$frame <-c(rep(0, genome_length), rep(1, genome_length), rep(2, genome_length))
ORF_length_heatmap_data$permutation_pvalue <- NA
ORF_length_heatmap_data$type <- NA
ORF_length_heatmap_data$gene <- NA
ORF_length_heatmap_data$length <- NA
ORF_length_heatmap_data$strand <- NA

# Construct -- takes a LONG TIME (~40 minutes) because it's a loop and this is R; option to simply reload below
for (this_rec_num in 1:nrow(frameshift_data)) {
  this_rec <- frameshift_data[this_rec_num, ]
  this_type <- this_rec$Type
  this_frame <- this_rec$frame
  this_name <- this_rec$`reference-frame-gene`
  this_start <- this_rec$`genome-co-ord1`
  this_end <- this_rec$`genome-co-ord2`
  this_length <- this_rec$length
  this_strand <- this_rec$strand
  this_permut_p <- this_rec$`pvalue-permutation`
  
  for (this_site in (this_start:this_end)) {
    #this_site <- this_start
    ORF_length_heatmap_data[ORF_length_heatmap_data$site == this_site & ORF_length_heatmap_data$frame == this_frame, ]$permutation_pvalue <- this_permut_p
    ORF_length_heatmap_data[ORF_length_heatmap_data$site == this_site & ORF_length_heatmap_data$frame == this_frame, ]$type <- this_type
    ORF_length_heatmap_data[ORF_length_heatmap_data$site == this_site & ORF_length_heatmap_data$frame == this_frame, ]$gene <- this_name
    ORF_length_heatmap_data[ORF_length_heatmap_data$site == this_site & ORF_length_heatmap_data$frame == this_frame, ]$length <- this_length
    ORF_length_heatmap_data[ORF_length_heatmap_data$site == this_site & ORF_length_heatmap_data$frame == this_frame, ]$strand <- this_strand
  }
}

### SAVE <-- CHANGE THIS
write_tsv(ORF_length_heatmap_data, "ORF_length_heatmap_data.tsv")

### RELOAD
(ORF_length_heatmap_data <- read_tsv("ORF_length_heatmap_data.tsv"))

# Prepare gene name data; coordinates from Wuhan-Hu-1
ORF_length_gene_data <- data.frame(
  gene = c("1a", "1ab", "S", "3a", "3c", "3d", "3b", "E", "M", "6", "7a", "7b", "8", "N", "9b", "9c", "10"),
  site = c((266+13483)/2,  (13468+21555)/2,  (21563+25384)/2,  (25393+26220)/2, (25457+25582)/2, (25524+25697)/2, (25814+25882)/2, (26245+26472)/2,  (26523+27191)/2,  (27202+27387)/2,  (27394+27759)/2,  (27756+27887)/2,  (27894+28259)/2,  (28274+29533)/2,  (28284+28577)/2,  (28734+28955)/2,  (29558+29674)/2),
  frame = c((266%%3),  (13468%%3),  (21563%%3),  (25393%%3),  (25457%%3), (25524%%3), (25814%%3), (26245%%3),  (26523%%3),  (27202%%3),  (27394%%3),  (27756%%3),  (27894%%3),  (28274%%3),  (28284%%3),  (28734%%3),  (29558%%3)),
  length = c((13483-266+1), (21555-13468+1), (25384-21563+1), (26220-25393+1), (25582-25457+1), (25697-25524+1), (25882-25814+1), (26472-26245+1), (27191-26523+1), (27387-27202+1), (27759-27394+1), (27887-27756+1), (28259-27894+1), (29533-28274+1), (28577-28284+1), (28955-28734+1), (29674-29558+1))
)

# Order frame
ORF_length_heatmap_data$frame <- factor(ORF_length_heatmap_data$frame, levels = c(1, 2, 0), labels = c('frame 1', 'frame 2', 'frame 3'))
ORF_length_gene_data$frame <- factor(ORF_length_gene_data$frame, levels = c(1, 2, 0), labels = c('frame 1', 'frame 2', 'frame 3'))

# Plot
(ORF_length_heatmap_PLOT <- ggplot(data = filter(ORF_length_heatmap_data), mapping = aes(x = site, y = frame)) + 
    geom_raster(mapping = aes(fill = permutation_pvalue)) +
    geom_text(data = filter(ORF_length_gene_data, gene %in% c("1a", "1ab", "S", "3a", "M", "7a", "8", "N")), 
              mapping = aes(label = gene), size = 4, color = brewer.pal(9, "Greens")[6]) + 
    geom_text(data = filter(ORF_length_gene_data, gene %in% c("E", "6")), # "7a"
              mapping = aes(label = gene), y = 3.75, size = 4, color = brewer.pal(9, "Greens")[6]) + 
    geom_text_repel(data = filter(ORF_length_gene_data, gene %in% c("3d", "7b", "9b", "9c")), # "8"
                    mapping = aes(label = gene), size = 4, color = brewer.pal(9, "Greens")[6], ylim = c(-Inf, 0.5)) + 
    geom_text_repel(data = filter(ORF_length_gene_data, gene %in% c("10")), 
                    mapping = aes(label = gene), size = 4, color = brewer.pal(9, "Greens")[6], xlim = c(30100, Inf)) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.background = element_blank()
    ) + 
    xlab("Site") +
    ylab("") +
    scale_y_discrete(expand = expand_scale(mult = c(0.5, 0.5))) +
    scale_fill_gradient(low = "white", high = "black", na.value = "black", name = "P-value"))


####################################################################################################
### NO LABELS, ONE FRAME AT A TIME

# Plot

#frame 1
(ORF_length_heatmap_PLOT_frame1 <- ggplot(data = filter(ORF_length_heatmap_data, frame == "frame 1"), mapping = aes(x = site, y = frame)) + 
   geom_raster(mapping = aes(fill = (permutation_pvalue))) +
   theme(axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.x= element_blank(),
         axis.ticks.y = element_blank(),
         legend.position = 'none',
         axis.title.x = element_blank(), 
         axis.title.y = element_blank(), 
         panel.background = element_blank()
   ) + 
   xlab("") + 
   ylab("") +
   scale_y_discrete(expand = expand_scale(mult = c(0.9, 0.9))) +
   scale_fill_gradient(low = "yellow", high = "black", name = "P-value"))

# frame 2
(ORF_length_heatmap_PLOT_frame2 <- ggplot(data = filter(ORF_length_heatmap_data, frame == "frame 2"), mapping = aes(x = site, y = frame)) + 
    geom_raster(mapping = aes(fill = (permutation_pvalue))) +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.x= element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          panel.background = element_blank()
    ) + 
    xlab("") + 
    ylab("") +
    scale_y_discrete(expand = expand_scale(mult = c(0.9, 0.9))) +
    scale_fill_gradient(low = "yellow", high = "black", name = "P-value")) 


# frame 3 
(ORF_length_heatmap_PLOT_frame3 <- ggplot(data = filter(ORF_length_heatmap_data, frame == "frame 3"), mapping = aes(x = site, y = frame)) + 
    geom_raster(mapping = aes(fill = (permutation_pvalue))) +
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          panel.background = element_blank()
    ) + 
    xlab("") + 
    ylab("") +
    scale_y_discrete(expand = expand_scale(mult = c(0.9, 0.9))) +
    scale_fill_gradient(low = "yellow", high = "black", name = "P-value")) # na.value = "dark", 


# PLOT TOGETHER USING PATCHWORK
ORF_length_heatmap_PLOT_frame1 / ORF_length_heatmap_PLOT_frame2 / ORF_length_heatmap_PLOT_frame3


