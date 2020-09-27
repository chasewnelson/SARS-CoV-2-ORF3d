### SARS-CoV-2 ORF3a hydrophobicity

library(patchwork)
library(RColorBrewer)
library(scales)
library(tidyverse)


# set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

# Import, for frames 3a(ss11), 3d(ss12), and 3c(ss13) respectively; trailing columns are meaningless
hydroph_data <- read_csv("data_zotero/hydrophobicity_profiles_ORF3a.csv")
names(hydroph_data) <- c('position', 'ss11', 'ss12', 'ss13')
hydroph_data

# Convert to LONG
(hydroph_data_LONG <- hydroph_data %>%
    pivot_longer(cols = c('ss11', 'ss12', 'ss13'), names_to = "frame", values_to = "hydrophobicity"))


# Distribution of x data?
summary(hydroph_data_LONG$position) 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.00   73.75  137.50  137.50  201.25  265.00 

# Dsibtribution of y data?
summary(hydroph_data_LONG$hydrophobicity)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
#-0.77905 -0.34929 -0.13714 -0.16982  0.01571  0.37429        1 

### SLIDING WINDOW FOR INTER-PAIR CORRELATIONS
WIN_SIZE <- 25

# initialize data frame
hydro_corr_results <- data.frame(frame_1 = character(), frame_2 = character(), position = integer(), r = numeric(), p = numeric())

# loop
for(this_frame_1 in unique(hydroph_data_LONG$frame)) {
  #this_frame_1 <- "ss12"
  
  for(this_frame_2 in unique(hydroph_data_LONG$frame)) {
    #this_frame_2 <- "ss13"
    
    if(this_frame_1 != this_frame_2) {
      
      frame_1_data <- filter(hydroph_data_LONG, frame == this_frame_1)
      frame_1_data <- arrange(frame_1_data, position)
      frame_2_data <- filter(hydroph_data_LONG, frame == this_frame_2)
      frame_2_data <- arrange(frame_2_data, position)
      
      if(all(frame_1_data$position == frame_2_data$position)) {
        
        for(i in sort(unique(frame_1_data$position))) {
          #i <- 10
          #cat(i, " ")
          
          max_i <- i + WIN_SIZE - 1
          
          if(max_i <= max(unique(frame_1_data$position))) {
            
            # Extract window; analyze
            frame_1_values <- frame_1_data[frame_1_data$position >= i & frame_1_data$position <= (i + WIN_SIZE - 1), ]$hydrophobicity
            frame_2_values <- frame_2_data[frame_2_data$position >= i & frame_2_data$position <= (i + WIN_SIZE - 1), ]$hydrophobicity
            
            # Determine correlation
            corr_result <- cor.test(frame_1_values, frame_2_values, method = "spearman")
            corr_result_r <- as.vector(corr_result$estimate)
            corr_result_P <- as.vector(corr_result$p.value)
            
            # temp data frame from addition
            hydro_corr_results_NEW <- data.frame(frame_1 = c(this_frame_1), frame_2 = c(this_frame_2), position = c(i), 
                                                 r = c(corr_result_r), p = c(corr_result_P))
            
            # Add to data frame
            hydro_corr_results <- rbind(hydro_corr_results, hydro_corr_results_NEW)
          }
        }
      }
    }
  }
}


# SAVE
#write_tsv(hydro_corr_results, "data/hydrophobicity_profiles_ORF3a_corr.tsv")

############################################################################################################
# Reload
hydro_corr_results <- read_tsv("data/hydrophobicity_profiles_ORF3a_corr.tsv")

# add center
hydro_corr_results$window_center <- ((hydro_corr_results$position) + (hydro_corr_results$position + WIN_SIZE - 1)) / 2

# add comparison
hydro_corr_results$comparison <- paste0(hydro_corr_results$frame_1, '/', hydro_corr_results$frame_2)

# add significance level
hydro_corr_results$significance <- "n.s."
hydro_corr_results[hydro_corr_results$p < 0.05, ]$significance <- "*"
hydro_corr_results[hydro_corr_results$p < 0.01, ]$significance <- "**"
hydro_corr_results[hydro_corr_results$p < 0.001, ]$significance <- "***"
hydro_corr_results$significance <- factor(hydro_corr_results$significance, levels = c('***', '**', '*', "n.s."),
                                          labels = c("p<0.001", "p<0.01", "p<0.05", "n.s."))

############################################################################################################
### PLOT for HYDROPHOBICITY of the ORF3a/ORF3c/ORF3d OLG region frames ###
(hydroph_data_LONG_PLOT <- ggplot(data = filter(hydroph_data_LONG, ! is.na(hydrophobicity)),
                                  mapping = aes(x = position, y = hydrophobicity, color = frame)) +
   
   # ORF3b
   geom_vline(xintercept = 141, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = 141, xmax = 164, ymin = 0.4, ymax = 0.495), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 141, y = (0.4 + 0.495) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3b") + # fontface = 'italic', 
   
   # ORF3a-2
   geom_vline(xintercept = 124, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = 124, xmax = Inf, ymin = 0.5, ymax = 0.595), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 124, y = (0.5 + 0.595) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3a-2") + 
   
   # ORF3d-2
   geom_vline(xintercept = 68, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = 68, xmax = 102, ymin = 0.5, ymax = 0.595), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 68, y = (0.5 + 0.595) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3d-2") + 
   
   # ORF3d
   geom_vline(xintercept = 44, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = 44, xmax = 102, ymin = 0.4, ymax = 0.495), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 44, y = (0.4 + 0.495) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3d") + 
   
   # ORF3c
   geom_vline(xintercept = 22, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = 22, xmax = 64, ymin = 0.5, ymax = 0.595), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 22, y = (0.5 + 0.595) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3c") + 
   
   # ORF3a
   geom_vline(xintercept = -Inf, linetype = "dotted", color = "grey", size = 0.25) +
   geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymin = 0.6, ymax = 0.695), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = -Inf, y = (0.6 + 0.695) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), label = "ORF3a") + 
   
   # DATA
   geom_line(size = rel(0.25)) +
   
   xlab("Residue (ORF3a)") +
   ylab("Hydrophobicity") +
   theme_bw() +
   theme(panel.grid = element_blank(),
         axis.title = element_text(size = rel(0.75)),
         legend.title = element_text(size = rel(0.9)),
         legend.text = element_text(size = rel(0.7)),
         strip.background = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.text.y.right = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x = element_blank()
   ) +
   scale_color_manual(values = c("#00D39B", "#BDB10F", "#D48EB4"), name = "Frame") + 
   scale_fill_manual(values = c("#00D39B", "#BDB10F", "#D48EB4"), name = "Frame") +
   scale_x_continuous(limits = c(min(hydroph_data_LONG$position), max(hydroph_data_LONG$position)), expand = expand_scale(mult = c(0, 0))) +
   scale_y_continuous(limits = c(-0.8, 0.695), breaks = c(-0.8, -0.4, 0, 0.4), expand = expand_scale(mult = c(0.05, 0)))) 


############################################################################################################
### PLOT for CORRELATION

(hydroph_corr_results_PLOT <- ggplot(data = filter(hydro_corr_results, comparison %in% c('ss11/ss12', 'ss11/ss13', 'ss12/ss13')),
                                     mapping = aes(x = window_center, y = r, color = significance, group = 1)) +
   
   # ORF3b
   geom_vline(xintercept = 141, linetype = "dotted", color = "grey", size = 0.25) +
   # ORF3a-2
   geom_vline(xintercept = 124, linetype = "dotted", color = "grey", size = 0.25) +
   # ORF3d-2
   geom_vline(xintercept = 68, linetype = "dotted", color = "grey", size = 0.25) +
   # ORF3d
   geom_vline(xintercept = 44, linetype = "dotted", color = "grey", size = 0.25) +
   # ORF3c
   geom_vline(xintercept = 22, linetype = "dotted", color = "grey", size = 0.25) +
   # ORF3a
   geom_vline(xintercept = -Inf, linetype = "dotted", color = "grey", size = 0.25) +
   # r=0 mark
   geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 0.25) +
   
   # DATA
   geom_line(size = rel(0.25)) +
   
   facet_grid(comparison ~ .) +
   xlab("Residue (ORF3a)") +
   ylab("Correlation") +
   theme_bw() +
   theme(panel.grid = element_blank(),
         axis.title = element_text(size = rel(0.75)),
         legend.title = element_text(size = rel(0.9)),
         legend.text = element_text(size = rel(0.7)),
         strip.background = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.text.y.right = element_blank()
   ) +
   scale_color_manual(values = brewer.pal(4, "Spectral"), name = "Significance") + 
   scale_x_continuous(limits = c(min(hydroph_data_LONG$position), max(hydroph_data_LONG$position)), expand = expand_scale(mult = c(0, 0))) +
   scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.05))))


# Plot both with patchwork
hydroph_data_LONG_PLOT / hydroph_corr_results_PLOT




############################################################################################################
### PLOT for Fig 3C: hydrophobicity zoomed in on the ORF3c/3d region
max_position <- 102+12 # end of ORF3d is 102

(hydroph_data_LONG_ORF3dRegion_PLOT <- ggplot(data = filter(hydroph_data_LONG, ! is.na(hydrophobicity), position < max_position),
                                              mapping = aes(x = position, y = hydrophobicity, color = frame)) +
    
    # ORF3d-2
    geom_vline(xintercept = 68, linetype = "dotted", color = "grey", size = 0.25) +
    geom_rect(mapping = aes(xmin = 68, xmax = 102, ymin = 0.525, ymax = 0.645), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
    geom_text(mapping = aes(x = 68, y = (0.525 + 0.645) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(3), label = "ORF3d-2") + 
    
    # ORF3d
    geom_vline(xintercept = 44, linetype = "dotted", color = "grey", size = 0.25) +
    geom_rect(mapping = aes(xmin = 44, xmax = 102, ymin = 0.4, ymax = 0.52), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
    geom_text(mapping = aes(x = 44, y = (0.4 + 0.52) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(3), label = "ORF3d") + 
    
    # ORF3c
    geom_vline(xintercept = 22, linetype = "dotted", color = "grey", size = 0.25) +
    geom_rect(mapping = aes(xmin = 22, xmax = 64, ymin = 0.525, ymax = 0.645), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
    geom_text(mapping = aes(x = 22, y = (0.525 + 0.645) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(3), label = "ORF3c") + 
    
    # ORF3a
    geom_vline(xintercept = -Inf, linetype = "dotted", color = "grey", size = 0.25) +
    geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymin = 0.65, ymax = 0.77), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
    geom_text(mapping = aes(x = -Inf, y = (0.65 + 0.77) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(3), label = "ORF3a") + 
    
    # DATA
    geom_line(size = rel(0.25)) +
    
    xlab("Residue (ORF3a)") +
    ylab("Hydrophobicity") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = 'none',
          axis.text.x = element_text(size = rel(1.05)), 
          axis.text.y = element_text(size = rel(1.05)),
          axis.title.x = element_text(size = rel(0.9)),
          axis.title.y = element_text(size = rel(0.9)),
          strip.background = element_blank() 
    ) +
    scale_color_manual(values = c("#00D39B", "#BDB10F", "#D48EB4"), name = "Frame") +
    scale_fill_manual(values = c("#00D39B", "#BDB10F", "#D48EB4"), name = "Frame") +
    scale_x_continuous(limits = c(min(hydroph_data_LONG$position), max_position), expand = expand_scale(mult = c(0, 0))) + # <-- CHANGE
    scale_y_continuous(limits = c(-0.8, 0.77), breaks = c(-0.8, -0.4, 0, 0.4), expand = expand_scale(mult = c(0.05, 0)))) 



################################################################################
### FINALLY, calculate correlation in the ORF3c, ORF3c/ORF3d, ORF3d, ORF3b, and remainder

### ORF3c: ORF3a codons 22-63, so 22-43 (43-22+1=22 codons)
# ss11/ss12
ORF3c_ss11_ss12_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 22:43)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss12", position %in% 22:43)$hydrophobicity, 
                                method = "spearman")
#S = 425.1, p-value = 4.069e-05
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.7599663 


# ss11/ss13
ORF3c_ss11_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 22:43)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 22:43)$hydrophobicity, 
                                method = "spearman")
#S = 2488.1, p-value = 0.06157
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.4049236 


# ss12/ss13
ORF3c_ss12_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss12", position %in% 22:43)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 22:43)$hydrophobicity, 
                                method = "spearman")
#S = 1936.7, p-value = 0.6787
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.09358217 


### ORF3c/ORF3d: ORF3a codons 44-102, so 44-63 (63-44+1=20 codons)
# ss11/ss12
ORF3cORF3d_ss11_ss12_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 44:63)$hydrophobicity, 
                                     filter(hydroph_data_LONG, frame == "ss12", position %in% 44:63)$hydrophobicity, 
                                     method = "spearman")
#S = 221.92, p-value = 5.108e-06
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8331455 


# ss11/ss13
ORF3cORF3d_ss11_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 44:63)$hydrophobicity, 
                                     filter(hydroph_data_LONG, frame == "ss13", position %in% 44:63)$hydrophobicity, 
                                     method = "spearman")
#S = 342.03, p-value = 0.0001755
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.7428356 


# ss12/ss13
ORF3cORF3d_ss12_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss12", position %in% 44:63)$hydrophobicity, 
                                     filter(hydroph_data_LONG, frame == "ss13", position %in% 44:63)$hydrophobicity, 
                                     method = "spearman")
#S = 393.74, p-value = 0.0005322
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.7039553 


### ORF3d: ORF3a codons 44-102, so 64-102 (102-64+1=39)
# ss11/ss12
ORF3d_ss11_ss12_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 64:102)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss12", position %in% 64:102)$hydrophobicity, 
                                method = "spearman")
#S = 1276.5, p-value = 5.787e-13
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8707979 


# ss11/ss13
ORF3d_ss11_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 64:102)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 64:102)$hydrophobicity, 
                                method = "spearman")
#S = 1349.3, p-value = 1.512e-12
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8634339 

# ss12/ss13
ORF3d_ss12_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss12", position %in% 64:102)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 64:102)$hydrophobicity, 
                                method = "spearman")
#S = 1657, p-value = 5.125e-11
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.8322871 


### ORF3b: ORF3a codons 141-164 (164-141+1=24)
# ss11/ss12
ORF3b_ss11_ss12_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 141-164)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss12", position %in% 141-164)$hydrophobicity, 
                                method = "spearman")
#S = 1400800, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.4990255 


# ss11/ss13
ORF3b_ss11_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", position %in% 141-164)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 141-164)$hydrophobicity, 
                                method = "spearman")
#S = 1527700, p-value = 6.099e-14
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.4471898 


# ss12/ss13
ORF3b_ss12_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss12", position %in% 141-164)$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", position %in% 141-164)$hydrophobicity, 
                                method = "spearman")
#S = 738080, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.732922 


### ORF3a remainder: 276-((102-22+1)+(164-141+1))=171 codons
# ss11/ss12
ORF3a_ss11_ss12_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss12", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                method = "spearman")
#S = 419910, p-value = 0.0008703
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.2681899 


# ss11/ss13
ORF3a_ss11_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss11", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                method = "spearman")
#S = 370980, p-value = 2.019e-05
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.3404448 


# ss12/ss13
ORF3a_ss12_ss13_cor <- cor.test(filter(hydroph_data_LONG, frame == "ss12", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                filter(hydroph_data_LONG, frame == "ss13", ! position %in% c(22:102, 141:164))$hydrophobicity, 
                                method = "spearman")
#S = 259330, p-value = 1.121e-12
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#0.5389483 


### ORF3a subregion correlation results data frame


# temp data frame from addition
(ORF3a_subregion_hydro_cor <- rbind(
  data.frame(gene = "ORF3c", frame_1 = "ss11", frame_2 = "ss12", r = as.vector(ORF3c_ss11_ss12_cor$estimate), p = as.vector(ORF3c_ss11_ss12_cor$p.value)),
  data.frame(gene = "ORF3c", frame_1 = "ss11", frame_2 = "ss13", r = as.vector(ORF3c_ss11_ss13_cor$estimate), p = as.vector(ORF3c_ss11_ss13_cor$p.value)),
  data.frame(gene = "ORF3c", frame_1 = "ss12", frame_2 = "ss13", r = as.vector(ORF3c_ss12_ss13_cor$estimate), p = as.vector(ORF3c_ss12_ss13_cor$p.value)),
  
  data.frame(gene = "ORF3cORF3d", frame_1 = "ss11", frame_2 = "ss12", r = as.vector(ORF3cORF3d_ss11_ss12_cor$estimate), p = as.vector(ORF3cORF3d_ss11_ss12_cor$p.value)),
  data.frame(gene = "ORF3cORF3d", frame_1 = "ss11", frame_2 = "ss13", r = as.vector(ORF3cORF3d_ss11_ss13_cor$estimate), p = as.vector(ORF3cORF3d_ss11_ss13_cor$p.value)),
  data.frame(gene = "ORF3cORF3d", frame_1 = "ss12", frame_2 = "ss13", r = as.vector(ORF3cORF3d_ss12_ss13_cor$estimate), p = as.vector(ORF3cORF3d_ss12_ss13_cor$p.value)),
  
  data.frame(gene = "ORF3d", frame_1 = "ss11", frame_2 = "ss12", r = as.vector(ORF3d_ss11_ss12_cor$estimate), p = as.vector(ORF3d_ss11_ss12_cor$p.value)),
  data.frame(gene = "ORF3d", frame_1 = "ss11", frame_2 = "ss13", r = as.vector(ORF3d_ss11_ss13_cor$estimate), p = as.vector(ORF3d_ss11_ss13_cor$p.value)),
  data.frame(gene = "ORF3d", frame_1 = "ss12", frame_2 = "ss13", r = as.vector(ORF3d_ss12_ss13_cor$estimate), p = as.vector(ORF3d_ss12_ss13_cor$p.value)),
  
  data.frame(gene = "ORF3b", frame_1 = "ss11", frame_2 = "ss12", r = as.vector(ORF3b_ss11_ss12_cor$estimate), p = as.vector(ORF3b_ss11_ss12_cor$p.value)),
  data.frame(gene = "ORF3b", frame_1 = "ss11", frame_2 = "ss13", r = as.vector(ORF3b_ss11_ss13_cor$estimate), p = as.vector(ORF3b_ss11_ss13_cor$p.value)),
  data.frame(gene = "ORF3b", frame_1 = "ss12", frame_2 = "ss13", r = as.vector(ORF3b_ss12_ss13_cor$estimate), p = as.vector(ORF3b_ss12_ss13_cor$p.value)),
  
  data.frame(gene = "ORF3a", frame_1 = "ss11", frame_2 = "ss12", r = as.vector(ORF3a_ss11_ss12_cor$estimate), p = as.vector(ORF3a_ss11_ss12_cor$p.value)),
  data.frame(gene = "ORF3a", frame_1 = "ss11", frame_2 = "ss13", r = as.vector(ORF3a_ss11_ss13_cor$estimate), p = as.vector(ORF3a_ss11_ss13_cor$p.value)),
  data.frame(gene = "ORF3a", frame_1 = "ss12", frame_2 = "ss13", r = as.vector(ORF3a_ss12_ss13_cor$estimate), p = as.vector(ORF3a_ss12_ss13_cor$p.value))
))

ORF3a_subregion_hydro_cor$frame_comparison <- paste0(ORF3a_subregion_hydro_cor$frame_1, "/", ORF3a_subregion_hydro_cor$frame_2)

# FACTOR gene regions
ORF3a_subregion_hydro_cor$gene <- factor(ORF3a_subregion_hydro_cor$gene, levels = c("ORF3a", "ORF3c", "ORF3cORF3d", "ORF3d", "ORF3b"))


### PLOT CORRELATIONS
(ORF3a_subregion_hydro_cor_PLOT <- ggplot(data = ORF3a_subregion_hydro_cor, mapping = aes(x = frame_comparison, y = r, fill = frame_comparison)) + 
    xlab("Frame comparison") + # xlab("") + # 
    ylab("Correlation (Spearman's rank)") + # ylab("") + 
    geom_bar(stat = "identity") +
    facet_wrap(. ~ gene, nrow = 1) + 
    
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = rel(0.7)),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_text(size = rel(0.7)),
          legend.text = element_text(size = rel(0.7)),
          strip.background = element_blank(),
          strip.text = element_text(size = rel(0.7))
    ) +
    scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.1))) +
    scale_fill_manual(values = c("#BDB10F", "#D48EB4", "gray"), name = "Frame comparison") # "#00D39B"
)


