### SARSCOV2 between-taxa for heatmap/strip
library(tidyverse)
library(boot)
library(RColorBrewer)

# set the working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

# Prepare taxa IDs and order

# LAM ET AL. 2020 ORDERING FOLLOWS
pair_names_SARSCOV2 <- c("SARSCOV2_vs_KY352407.1",
                         "SARSCOV2_vs_GU190215.1", 
                         "SARSCOV2_vs_MG772933.1", 
                         "SARSCOV2_vs_EPI_ISL_410721",
                         "SARSCOV2_vs_MN996532.1",
                         "NC_045512.2",
                         "SARSCOV2_vs_EPI_ISL_410540",
                         "SARSCOV2_vs_GQ153547.1", 
                         "SARSCOV2_vs_GQ153541.1",
                         "SARSCOV2_vs_KJ473814.1",
                         "SARSCOV2_vs_KU182964.1", 
                         "SARSCOV2_vs_MK211377.1", 
                         "SARSCOV2_vs_KY417147.1", 
                         "SARSCOV2_vs_KY417143.1", 
                         "SARSCOV2_vs_FJ588686.1", 
                         "SARSCOV2_vs_MK211376.1", 
                         "SARSCOV2_vs_KC881006.1",
                         "SARSCOV2_vs_KY417151.1",
                         "SARSCOV2_vs_KY417146.1",
                         "SARSCOV2_vs_NC_004718.3",
                         "SARSCOV2_vs_AY502924.1")

pair_names_SARS <- c("SARSCOV2_KY352407.1",
                     "SARSCOV2_GU190215.1", 
                     "SARSCOV2_MG772933.1", 
                     "SARSCOV2_EPI_ISL_410721",
                     "SARSCOV2_MN996532.1",
                     "SARSCOV2_NC_045512.2",
                     "SARSCOV2_EPI_ISL_410540",
                     "SARSCOV2_GQ153547.1", 
                     "SARSCOV2_GQ153541.1", 
                     "SARSCOV2_KJ473814.1",
                     "SARSCOV2_KU182964.1", 
                     "SARSCOV2_MK211377.1", 
                     "SARSCOV2_KY417147.1", 
                     "SARSCOV2_KY417143.1", 
                     "SARSCOV2_FJ588686.1", 
                     "SARSCOV2_MK211376.1", 
                     "SARSCOV2_KC881006.1", 
                     "SARSCOV2_KY417151.1", 
                     "SARSCOV2_KY417146.1",
                     "NC_004718.3",
                     "SARSCOV2_AY502924.1")

pair_name_converter <- c("Bat SARS-like CoV (BtKY72/2007)",
                         "Bat CoV (BM48-31/BGR/2008)", 
                         "Bat SARS-like CoV (ZC45/2017)", 
                         "Pangolin CoV (PD/1/2019)",
                         "Bat CoV (RaTG13/2013)",
                         "SARS-CoV-2 (Wuhan-Hu-1/2019)",
                         "Pangolin CoV (GX/P5L/2017)",
                         "Bat SARS-CoV (HKU3-12/2004-2008)", 
                         "Bat SARS-CoV (HKU3-6/2004-2008)", 
                         "Bat CoV (HuB2013/2013)",
                         "Bat CoV (JTMC15/2013)", 
                         "Bat CoV (YN2018C/2016)", 
                         "Bat SARS-like CoV (Rs4237/2013)", 
                         "Bat SARS-like CoV (Rs4081/2012)", 
                         "Bat SARS-CoV (Rs672/2006)", 
                         "Bat CoV (YN2018B/2016)", 
                         "Bat SARS-like CoV (Rs3367/2012)", 
                         "Bat SARS-like CoV (Rs7327/2014)", 
                         "Bat SARS-like CoV (Rs4231/2013)", 
                         "SARS-CoV (Tor2/2003)",
                         "SARS-CoV (TW11/2003)")


###############################################################################
### LOAD SLIDING WINDOW RESULTS FROM OLGENIE

##############################
### SARSCOV2 vs. each, JC corrected <-- OPTION 1
# Only use (a) 3a or (b) N, NOT BOTH

# (a) 3a/ss12 
data_subset <- read_tsv("data_zotero/SARS-CoV-2-ref_ORF3a_ss12_windows_prepangolin.txt") # 5,480 x 43
data_subset <- filter(data_subset, pair != "SARSCOV2_vs_EPI_ISL_410540") # 5,206 x 43
data_subset_P5L <- read_tsv("data_zotero/SARS-CoV-2-ref_ORF3a_ss12_windows_pangolin.txt") # 274 x 43
data_subset <- rbind(data_subset, data_subset_P5L) # 5,480 x 43
#write_tsv(data_subset, "/Users/cwnelson88/Desktop/SARS-CoV-2/eLife_revision/Figures/Fig 6 - Sliding Window Heatmap/SARS-CoV-2-ref_ORF3a_ss12_windows.txt")


# (b) N/ss13
#data_subset <- read_tsv("data/SARS-CoV-2-ref_N_ss13_windows.tsv")


##############################
### SARS vs. each  <-- OPTION 2
# Only use (a) 3a or (b) N, NOT BOTH

### SARS vs. each, JC corrected
# (a) 3a/ss13
#data_subset <- read_tsv("data/SARS-ref_ORF3a_ss13_windows.tsv")

# (b) N/ss13
#data_subset <- read_tsv("data/SARS-ref_N_ss13_windows.tsv")
##############################


##############################
# tell R the window size we used
win_size <- 50 
step_size <- 1

# add central codon coordinate
data_subset$codon_center <- (2 * data_subset$codon + win_size - 1) / 2


##############################
# NORMALIZE DN/DS: log inverse P of Z
data_subset$dNdS_norm <- NA
data_subset[! is.na(data_subset$sw_dNdS) & data_subset$sw_dN > data_subset$sw_dS, ]$dNdS_norm <- 
  log(1 / data_subset[! is.na(data_subset$sw_dNdS) & data_subset$sw_dN > data_subset$sw_dS, ]$sw_boot_dN_m_dS_P)
data_subset[! is.na(data_subset$sw_dNdS) & data_subset$sw_dN < data_subset$sw_dS, ]$dNdS_norm <- 
  -log(1 / data_subset[! is.na(data_subset$sw_dNdS) & data_subset$sw_dN < data_subset$sw_dS, ]$sw_boot_dN_m_dS_P)


#############################
max_codon <- max(data_subset$codon, na.rm = TRUE)


#############################
# add referent for SARSCOV2 # <-- OPTION 1
names(pair_name_converter) <- pair_names_SARSCOV2
data_subset_referent <- data_subset[data_subset$pair == "SARSCOV2_vs_MN996532.1", ] # some random pair
data_subset_referent$pair <- "NC_045512.2" # Wuhan-Hu-1
data_subset_referent$dNdS_norm <- NA
data_subset <- rbind(data_subset, data_subset_referent)

# reorder names for SARSCOV2 #
data_subset$pair <- factor(data_subset$pair,
                           levels = rev(pair_names_SARSCOV2),
                           labels = rev(pair_name_converter))
#############################

#############################
## add referent for SARS # <-- OPTION 2
#names(pair_name_converter) <- pair_names_SARS
#data_subset_referent <- data_subset[data_subset$pair == "SARSCOV2_AY502924.1", ] # some random pair
#data_subset_referent$pair <- "NC_004718.3" # SARS SARS_CoV_Tor2/2003
#data_subset_referent$dNdS_norm <- NA
#data_subset <- rbind(data_subset, data_subset_referent)
#
## reorder names for SARS#
#data_subset$pair <- factor(data_subset$pair,
#                           levels = rev(pair_names_SARS),
#                           labels = rev(pair_name_converter))
#############################

# MAX SIG ALL
max_sig_value <- 6.907756


### PLOT
(phylogenetic_sw_heatmap_PLOT <- ggplot(filter(data_subset, codon <= (max_codon - win_size + 1)), aes(x = codon_center, y = pair, fill = dNdS_norm)) +
    xlab("Codon (ORF3a)") + # <-- CHANGE
    #xlab("Codon (N)") + # <-- CHANGE
    geom_raster() + # data_subset
    
    ### SARS-CoV-2, ORF3a, ss12 # <-- CHANGE
    # ORF3a / ORF3c gene coordinates
    geom_vline(xintercept = 22, linetype = "dashed", color = "grey") + # it's off the edge given the sliding window size
    geom_vline(xintercept = 64, linetype = "dashed", color = "grey") +
    
    # ORF3a / ORF3d gene coordinates
    geom_vline(xintercept = 44, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 102, linetype = "dashed", color = "black") +
    
    # ORF3a / ORF3b gene coordinates (short isoform)
    geom_vline(xintercept = 141, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 164, linetype = "dashed", color = "grey") +
    
    
    ### SARS-CoV, ORF3a, ss13 # <-- CHANGE
    ## ORF3a / ORF3c gene coordinates
    #geom_vline(xintercept = 64, linetype = "dashed", color = "black") + # start is off the edge given sliding window size
    #
    ## ORF3a / ORF3d gene coordinates
    #geom_vline(xintercept = 44, linetype = "dashed", color = "grey") +
    #geom_vline(xintercept = 102, linetype = "dashed", color = "grey") +
    #
    ## ORF3a / ORF3b gene coordinates
    #geom_vline(xintercept = 141, linetype = "dashed", color = "black") + # then goes further than 164
    
    ## N / 9b and 9c gene coordinates # <-- CHANGE
    #geom_vline(xintercept = 103, linetype = "dashed") +
    #geom_vline(xintercept = 155, linetype = "dashed") +
    #geom_vline(xintercept = 229, linetype = "dashed") +
    
    theme(axis.text.x = element_text(size = 15.5), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          legend.text = element_blank(),
          axis.title.x = element_text(size = 17),
          axis.title.y = element_blank(), 
          legend.title = element_blank(),
          panel.background = element_blank()
    ) + 
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.02))) +
    scale_y_discrete(position = "right") +
    scale_fill_gradient2(limits = c(-max_sig_value, max_sig_value), breaks = c(0),
                         low = brewer.pal(9, "Set1")[2], mid = 'white', high = brewer.pal(9, "Set1")[1], midpoint = 0,  na.value = "white")) 





