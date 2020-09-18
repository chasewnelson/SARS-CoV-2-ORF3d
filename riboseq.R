##############################################################################################################
### Ribo-seq analysis and visualization

library(tidyverse)
library(RColorBrewer)
library(scales)
library(patchwork)

# set working directory
setwd("/Users/cwnelson88/scripts_NGS/SARS-CoV-2-ORF3d/")

##############################################################################################################
### MAPPED READS BY POSITION OF FIRST NUCLEOTIDE IN GENOME

##############################################################################################################
### INPUT
mapped_reads_by_readlength <- read_tsv("data/mapped_reads_by_readlength.tsv")
mapped_reads_by_readlength$position <- mapped_reads_by_readlength$position + 1
mapped_reads_by_readlength
# 2,610,414 x 6
# sample      read_length chromosome  position frame read_count

# input and construct SAMPLE METADATA translator
sample_metadata_joiner <- read_tsv("data/riboseq_sample_metadata.tsv")

# Join
mapped_reads_by_readlength <- left_join(x = mapped_reads_by_readlength, y = sample_metadata_joiner, by = "sample")
mapped_reads_by_readlength
# still 2,610,414 x 6? YES, it's 2,610,414 x 9

# Another condition option
mapped_reads_by_readlength$condition_combn <- paste0(mapped_reads_by_readlength$treatment, " / ", mapped_reads_by_readlength$time)
unique(mapped_reads_by_readlength$condition_combn)
mapped_reads_by_readlength$condition_combn <- factor(mapped_reads_by_readlength$condition_combn,
                                                     levels = c("CHX / 5 hpi", "CHX / 24 hpi",
                                                                "Harr / 5 hpi", "Harr / 24 hpi",
                                                                "LTM / 5 hpi", "LTM / 24 hpi",
                                                                "mRNA / 5 hpi", "mRNA / 24 hpi")) # each of these has two samples

# Factor
mapped_reads_by_readlength$frame <- factor(mapped_reads_by_readlength$frame,
                                           levels = c(0, 1, 2), # was 0-based when determined
                                           labels = c('Frame 1', 'Frame 2', 'Frame 3')) 




####################################################################################################
### codon position x read length while excluding OLGs

# Strip away superfluous columns
(mapped_reads_by_readlength_wCDS <- dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition_combn))

### Add CDS position, knowing position is 1-based ###
mapped_reads_by_readlength_wCDS$region_type <- "UTR"
mapped_reads_by_readlength_wCDS$gene <- NA
mapped_reads_by_readlength_wCDS$CDS <- NA
mapped_reads_by_readlength_wCDS$dist_from_end <- NA

#ORF1ab_1:266-13468
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$gene <- "ORF1ab_1"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$position - 266 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$dist_from_end <- 
  13468 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 266 & mapped_reads_by_readlength_wCDS$position <= 13468, ]$position

#ORF1ab_2:13468-21555
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$gene <- "ORF1ab_2"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$position - 13468 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$dist_from_end <- 
  21555 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 13468 & mapped_reads_by_readlength_wCDS$position <= 21555, ]$position

#S:21563-25384
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$gene <- "S"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$position - 21563 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$dist_from_end <- 
  25384 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 21563 & mapped_reads_by_readlength_wCDS$position <= 25384, ]$position

#ORF3a:25393-26220
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$gene <- "ORF3a"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$position - 25393 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$dist_from_end <- 
  26220 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25393 & mapped_reads_by_readlength_wCDS$position <= 26220, ]$position

#E:26245-26472
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$gene <- "E"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$position - 26245 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$dist_from_end <- 
  26472 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26245 & mapped_reads_by_readlength_wCDS$position <= 26472, ]$position

#M:26523-27191
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$gene <- "M"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$position - 26523 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$dist_from_end <- 
  27191 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 26523 & mapped_reads_by_readlength_wCDS$position <= 27191, ]$position

#ORF6:27202-27387
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$gene <- "ORF6"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$position - 27202 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$dist_from_end <- 
  27387 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27202 & mapped_reads_by_readlength_wCDS$position <= 27387, ]$position

#ORF7a:27394-27759
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$gene <- "ORF7a"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$position - 27394 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$dist_from_end <- 
  27759 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27394 & mapped_reads_by_readlength_wCDS$position <= 27759, ]$position

#ORF7b:27756-27887
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$gene <- "ORF7b"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$position - 27756 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$dist_from_end <- 
  27887 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27756 & mapped_reads_by_readlength_wCDS$position <= 27887, ]$position

#ORF8:27894-28259
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$gene <- "ORF8"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$position - 27894 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$dist_from_end <- 
  28259 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 27894 & mapped_reads_by_readlength_wCDS$position <= 28259, ]$position

#N:28274-29533
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$gene <- "N"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$position - 28274 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$dist_from_end <- 
  29533 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28274 & mapped_reads_by_readlength_wCDS$position <= 29533, ]$position

#ORF10:29558-29674
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$region_type <- "gene"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$gene <- "ORF10"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$position - 29558 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$dist_from_end <- 
  29674 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 29558 & mapped_reads_by_readlength_wCDS$position <= 29674, ]$position


### Exclude/categorize OLG sites ###
#pp1ab/pp1a (266..13468,13468..21555)/(266..13483)
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (13468-2) & mapped_reads_by_readlength_wCDS$position <= (13483+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (13468-2) & mapped_reads_by_readlength_wCDS$position <= (13483+2), ]$CDS <- NA
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (13468-2) & mapped_reads_by_readlength_wCDS$position <= (13483+2), ]$gene <- NA

#ORF3c:25457-25582
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (25457-2) & mapped_reads_by_readlength_wCDS$position <= (25582+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25457 & mapped_reads_by_readlength_wCDS$position <= 25582, ]$gene <- "ORF3c"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25457 & mapped_reads_by_readlength_wCDS$position <= 25582, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25457 & mapped_reads_by_readlength_wCDS$position <= 25582, ]$position - 25457 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25457 & mapped_reads_by_readlength_wCDS$position <= 25582, ]$dist_from_end <- 
  25582 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25457 & mapped_reads_by_readlength_wCDS$position <= 25582, ]$position

#ORF3d:25524-25697
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (25524-2) & mapped_reads_by_readlength_wCDS$position <= (25697+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25524 & mapped_reads_by_readlength_wCDS$position <= 25697, ]$gene <- "ORF3d"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25524 & mapped_reads_by_readlength_wCDS$position <= 25697, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25524 & mapped_reads_by_readlength_wCDS$position <= 25697, ]$position - 25524 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25524 & mapped_reads_by_readlength_wCDS$position <= 25697, ]$dist_from_end <- 
  25697 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25524 & mapped_reads_by_readlength_wCDS$position <= 25697, ]$position

#ORF3b:25814-25882
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (25814-2) & mapped_reads_by_readlength_wCDS$position <= (25882+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25814 & mapped_reads_by_readlength_wCDS$position <= 25882, ]$gene <- "ORF3b"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25814 & mapped_reads_by_readlength_wCDS$position <= 25882, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25814 & mapped_reads_by_readlength_wCDS$position <= 25882, ]$position - 25814 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25814 & mapped_reads_by_readlength_wCDS$position <= 25882, ]$dist_from_end <- 
  25882 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 25814 & mapped_reads_by_readlength_wCDS$position <= 25882, ]$position

#ORF7a/ORF7b (27394-27759)/(27756-27887)
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (27756-2) & mapped_reads_by_readlength_wCDS$position <= (27759+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (27756-2) & mapped_reads_by_readlength_wCDS$position <= (27759+2), ]$CDS <- NA
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (27756-2) & mapped_reads_by_readlength_wCDS$position <= (27759+2), ]$gene <- NA

#ORF9b:28284-28577
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (28284-2) & mapped_reads_by_readlength_wCDS$position <= (28577+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28284 & mapped_reads_by_readlength_wCDS$position <= 28577, ]$gene <- "ORF9b"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28284 & mapped_reads_by_readlength_wCDS$position <= 28577, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28284 & mapped_reads_by_readlength_wCDS$position <= 28577, ]$position - 28284 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28284 & mapped_reads_by_readlength_wCDS$position <= 28577, ]$dist_from_end <- 
  28577 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28284 & mapped_reads_by_readlength_wCDS$position <= 28577, ]$position

#ORF9c:28734-28955
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= (28734-2) & mapped_reads_by_readlength_wCDS$position <= (28955+2), ]$region_type <- "OLG"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28734 & mapped_reads_by_readlength_wCDS$position <= 28955, ]$gene <- "ORF9c"
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28734 & mapped_reads_by_readlength_wCDS$position <= 28955, ]$CDS <- 
  mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28734 & mapped_reads_by_readlength_wCDS$position <= 28955, ]$position - 28734 + 1
mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28734 & mapped_reads_by_readlength_wCDS$position <= 28955, ]$dist_from_end <- 
  28955 - mapped_reads_by_readlength_wCDS[mapped_reads_by_readlength_wCDS$position >= 28734 & mapped_reads_by_readlength_wCDS$position <= 28955, ]$position


## Add frame, 0-based (0, 1, 2)
mapped_reads_by_readlength_wCDS$frame <- NA
mapped_reads_by_readlength_wCDS[! is.na(mapped_reads_by_readlength_wCDS$CDS), ]$frame <- 
  (mapped_reads_by_readlength_wCDS[! is.na(mapped_reads_by_readlength_wCDS$CDS), ]$CDS - 1) %% 3
mapped_reads_by_readlength_wCDS$frame <- factor(mapped_reads_by_readlength_wCDS$frame,
                                                levels = c(0, 1, 2), # 0-based
                                                labels = c('Codon pos 1', 'Codon pos 2', 'Codon pos 3'))

##################################################
### PREPARE FOR PLOTTING, ALL OF SUBSET
### vvv CHANGE THIS vvv ###
genes_included <- c("ORF1ab_1", "ORF1ab_2", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10") # all
#genes_included <- c("S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10") # all except ORF1ab
#genes_included <- "ORF3a"
#genes_included <- "N"
regions_included <- c("gene")
#regions_included <- c("gene", "OLG")
min_read_length <- 26
max_read_length <- 32
min_CDS <- 1
min_dist_from_end <- 0 #50
### ^^^ CHANGE THIS ^^^ ###

# Summarize
(mapped_reads_by_readlength_wCDS_summary <- filter(mapped_reads_by_readlength_wCDS, region_type %in% regions_included, gene %in% genes_included,
                                                   ! is.na(CDS), CDS >= min_CDS,
                                                   ! is.na(dist_from_end), dist_from_end >= min_dist_from_end) %>%
    group_by(read_length, time, treatment, frame) %>%
    summarise(
      read_count = sum(read_count)
    )) # 360 x 5

# Read count sum for each combination (pool frames)
(mapped_reads_by_readlength_wCDS_sums <- filter(mapped_reads_by_readlength_wCDS, region_type %in% regions_included, gene %in% genes_included,
                                                ! is.na(CDS), CDS >= min_CDS,
                                                ! is.na(dist_from_end), dist_from_end >= min_dist_from_end) %>%
    group_by(read_length, time, treatment) %>%
    summarise(
      read_count_sum = sum(read_count) # 
    )) # 120 x 4

# Join
mapped_reads_by_readlength_wCDS_summary <- left_join(mapped_reads_by_readlength_wCDS_summary, mapped_reads_by_readlength_wCDS_sums, 
                                                     by = c("read_length", "time", "treatment"))

# Proportion in each frame
mapped_reads_by_readlength_wCDS_summary$prop <- mapped_reads_by_readlength_wCDS_summary$read_count / mapped_reads_by_readlength_wCDS_summary$read_count_sum

## Calculate maxes -- specific for the range of lengths we want to plot
read_count_sum_max_5hr <- max(mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "5 hpi" & mapped_reads_by_readlength_wCDS_summary$read_length >= 26 & mapped_reads_by_readlength_wCDS_summary$read_length <= 32, ]$read_count_sum)
prop_max_5hr <- max(mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "5 hpi" & mapped_reads_by_readlength_wCDS_summary$read_length >= 26 & mapped_reads_by_readlength_wCDS_summary$read_length <= 32, ]$prop)
read_count_sum_max_24hr <- max(mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "24 hpi" & mapped_reads_by_readlength_wCDS_summary$read_length >= 26 & mapped_reads_by_readlength_wCDS_summary$read_length <= 32, ]$read_count_sum)
prop_max_24hr <- max(mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "24 hpi" & mapped_reads_by_readlength_wCDS_summary$read_length >= 26 & mapped_reads_by_readlength_wCDS_summary$read_length <= 32, ]$prop)

mapped_reads_by_readlength_wCDS_summary$read_count_sum_max <- read_count_sum_max_5hr
mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "24 hpi", ]$read_count_sum_max <- read_count_sum_max_24hr
mapped_reads_by_readlength_wCDS_summary$prop_max <- prop_max_5hr
mapped_reads_by_readlength_wCDS_summary[mapped_reads_by_readlength_wCDS_summary$time == "24 hpi", ]$prop_max <- prop_max_24hr

### PLOT Figure 2—figure supplement 3
(mapped_reads_by_readlength_wCDS_summary_PLOT <- ggplot(data = filter(mapped_reads_by_readlength_wCDS_summary, read_length >= min_read_length, read_length <= max_read_length),
                                                        mapping = aes(x = read_length, y = 100*prop, color = frame)) +
    geom_bar(data = filter(mapped_reads_by_readlength_wCDS_summary, read_length >= min_read_length, read_length <= max_read_length, frame == "Codon pos 1"), # otherwise triplicate plotting
             mapping = aes(x = read_length, y = read_count_sum / read_count_sum_max * 100 * prop_max), fill = "#E0E0E0", stat = "identity", inherit.aes = FALSE) +
    geom_hline(yintercept = 100*1/3, linetype = "dashed", color = "lightgrey") + 
    geom_point(size = rel(0.8)) +
    geom_line(size = rel(0.25)) +
    facet_grid(time ~ treatment) +
    xlab("Read length") +
    ylab("Fraction reads mapping (%)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = rel(0.75)),
          legend.title = element_blank(),
          legend.text = element_text(size = rel(0.7)), 
          axis.text = element_text(size = rel(0.7)), 
          strip.background = element_blank(),
          strip.text.y = element_text(angle = 0),
          axis.ticks.y.right = element_blank(),
          axis.text.y.right = element_blank()
    ) +
    scale_color_manual(values = c("#00D39B", "#D48EB4", "#D0C311")) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))))



####################################################################################################
### SLIDING WINDOWS prop frame analysis, pooling the samples for a given condition

### USE PYTHON with desired parameter combination, e.g.,
#DO: Rscript riboseq_sliding_window.R riboseq_sample_metadata.tsv mapped_reads_by_readlength_ALL.tsv 30 30

############################################################################
######### vv CHANGE THIS: SELECTION WINDOW SIZE AND READ LENGTH vv #########
############################################################################
# SPECIFY THE WINDOW SIZE YOU USED (CRITICAL)
WIN_SIZE <- 30 # 19 29 30 # also defined with you run the bacth R script for sliding windows; MAKE SURE THEY MATCH

### WIN_SIZE=19 ###
#READ_LEN=29
#mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win19read29.tsv")
#READ_LEN=30
#mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win19read30.tsv")

### WIN_SIZE=29 ###
#READ_LEN=29
#mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win29read29.tsv")
#READ_LEN=30
#mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win29read30.tsv")

### WIN_SIZE=30 ###
#READ_LEN=29
#mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win30read29.tsv")
#READ_LEN=30
mapped_reads_by_readlength_sumSamples <- read_tsv("data/mapped_reads_by_readlength_sumSamples_win30read30.tsv")

############################################################################
######### ^^ CHANGE THIS: SELECT WINDOW SIZE AND READ LENGTH ^^ #########
############################################################################


############################################################################
### Wrangle data
# Sum of all read counts mapping to window
mapped_reads_by_readlength_sumSamples$window_sum <- mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_1 +
  mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_2 +
  mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_3

# Center of window, if important
mapped_reads_by_readlength_sumSamples$window_center <- ((2 * mapped_reads_by_readlength_sumSamples$position) + WIN_SIZE - 1) / 2

### MAKE LONG
(mapped_reads_by_readlength_sumSamples_LONG <- dplyr::select(mapped_reads_by_readlength_sumSamples, 
                                                             condition_combn, position, window_center, window_sum, 
                                                             window_sum_codon_pos_1, window_sum_codon_pos_2, window_sum_codon_pos_3) %>%
    pivot_longer(cols = c('window_sum_codon_pos_1', 'window_sum_codon_pos_2', 'window_sum_codon_pos_3'), names_to = "frame", values_to = "read_count"))

# Factor frame
mapped_reads_by_readlength_sumSamples_LONG$frame <- factor(mapped_reads_by_readlength_sumSamples_LONG$frame,
                                                           levels = c('window_sum_codon_pos_1', 'window_sum_codon_pos_2', 'window_sum_codon_pos_3'),
                                                           labels = c("Frame 1", "Frame 2", "Frame 3"))

# Factor condition
mapped_reads_by_readlength_sumSamples_LONG$condition_combn <- factor(mapped_reads_by_readlength_sumSamples_LONG$condition_combn,
                                                                     levels = c("CHX / 5 hpi", "CHX / 24 hpi", 
                                                                                "Harr / 5 hpi", "LTM / 5 hpi", 
                                                                                "Harr / 24 hpi", "LTM / 24 hpi", 
                                                                                "mRNA / 5 hpi", "mRNA / 24 hpi"))


############################################################################################################################
### PLOT

### GENE

#ORF1ab_1:266-13468
#ORF1ab_2:13468-21555
#S:21563-25384

#ORF3a:25393-26220
#ORF3c:25457-25582
#ORF3d:25524-25697
#ORF3d_short:25596-25697
#ORF3a_short:25765-26220
#ORF3b:25814-25882

#ORF3b_1:25814-25882
#ORF3b_2:25910-25984
#ORF3b_3:26072-26170
#ORF3b_4:26183-26281

#E:26245-26472
#M:26523-27191
#ORF6:27202-27387
#ORF7a:27394-27759
#ORF7b:27756-27887
#ORF8:27894-28259
#N:28274-29533
#ORF9b:28284-28577
#ORF9c:28734-28955
#ORF10:29558-29674

### vvv CHANGE THIS vvv ###
#curr_conditions <- c("CHX / 5 hpi", "Harr / 5 hpi", "LTM / 5 hpi", "mRNA / 5 hpi")
curr_conditions <- c("Harr / 5 hpi", "LTM / 5 hpi")
#curr_conditions <- c("CHX / 5 hpi") 
# "Harr / 5 hpi", "LTM / 5 hpi", "Harr / 24 hpi", "LTM / 24 hpi", "CHX / 5 hpi", "CHX / 24 hpi", "mRNA / 5 hpi", "mRNA / 24 hpi"
### ^^^ CHANGE THIS ^^^ ###


(mapped_reads_by_readlength_sumSamples_LONG_SUBSET <- dplyr::filter(mapped_reads_by_readlength_sumSamples_LONG, 
                                                                    condition_combn %in% curr_conditions, 
                                                                    ### CHANGE THIS ###
                                                                    #window_center >= 28274 - 200, # for 8-N-9b-9c-10: window_center >= 28274 - 200
                                                                    #window_center <= 29533 + 200)) # for 8-N-9b-9c-10: window_center <= 29533 + 200
                                                                    #window_center >= 25393 - 200, # for S - 3a - E: window_center >= 25393 - 200
                                                                    #window_center <= 26220 + 200)) # for S - 3a - E: window_center <= 26220 + 200))
                                                                    #window_center >= 13468 - 250, # for ORF1ab ribosomal slip: window_center >= 13468 - 250
                                                                    #window_center <= 13468 + 250)) # for ORF1ab ribosomal slip: window_center <= 13468 + 250
                                                                    window_center >= 25457-50, # for just 3a OLGs: 25457-50 # for just 3a OLGs, wider: 25400 
                                                                    window_center <= 25697+50)) # for just 3a OLGs, alt: 25697+50 # for just 3a OLGs, wider: 25800 

###########################
### vvv CHANGE THIS vvv ###
# IFF you want to POOL all selected treatments together, as for Fig 2
(mapped_reads_by_readlength_sumSamples_LONG_SUBSET <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET %>%
   group_by(position, window_center, frame) %>%
   summarise(
     window_sum = sum(window_sum),
     read_count = sum(read_count)
   )) 
## ^^^ CHANGE THIS ^^^ ###
###########################

### Add binomial confidence intervals
mapped_reads_by_readlength_sumSamples_LONG_SUBSET$prop_reads <- NA
mapped_reads_by_readlength_sumSamples_LONG_SUBSET$binom_CI_point <- NA
mapped_reads_by_readlength_sumSamples_LONG_SUBSET$binom_CI_low <- NA
mapped_reads_by_readlength_sumSamples_LONG_SUBSET$binom_CI_high <- NA

# loop for CI
for (i in 1:nrow(mapped_reads_by_readlength_sumSamples_LONG_SUBSET)) {
  #i <- 1
  mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$prop_reads <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$read_count / mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$window_sum
  
  binom_CI <- prop.test(x = mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$read_count, n = mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$window_sum, conf.level = 0.95)
  mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$binom_CI_point <- as.vector(binom_CI$estimate)
  mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$binom_CI_low <- as.vector(binom_CI$conf.int[1])
  mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$binom_CI_high <- as.vector(binom_CI$conf.int[2])
  mapped_reads_by_readlength_sumSamples_LONG_SUBSET[i, ]$binom_CI_point <- as.vector(binom_CI$estimate)
}

#View(mapped_reads_by_readlength_sumSamples_LONG_SUBSET)
(read_count_sum_max <- max(mapped_reads_by_readlength_sumSamples_LONG_SUBSET$window_sum))
(prop_max <- max(mapped_reads_by_readlength_sumSamples_LONG_SUBSET$prop_reads))


# REGIONAL CONTEXT <-- cannot easily be generalized ### CHANGE for region ###
(mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT <- dplyr::filter(mapped_reads_by_readlength_sumSamples_LONG, 
                                                                            condition_combn %in% curr_conditions,
                                                                            ### CHANGE THIS ###
                                                                            #(window_center >= 28274 & window_center < 28284) | (window_center > 28577 & window_center < 28734) | (window_center > 28955 & window_center <= 29533))) # N excluding OLGs
                                                                            (window_center >= 25393 & window_center < 25457) | (window_center > 25697 & window_center <= 26220))) # ORF3a excluding OLGs
#(window_center >= 266 & window_center < 13468) | (window_center > 13483 & window_center <= 21555))) # ORF1ab excluding OL region
##window_center >= 25393, # ORF3a start; not used because includes OLGs
##window_center <= 26220)) # ORF3a end; not used because includes OLGs


###########################
### OPTION 1: POOL TREATMENTS (main Fig. 2)
### vvv CHANGE THIS vvv ###
# IFF you want to POOL all selected treatments together, as for Fig 2B
(mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT %>%
   group_by(position, window_center, frame) %>%
   summarise(
     window_sum = sum(window_sum),
     read_count = sum(read_count)
   )) 

mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$prop <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$read_count / mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$window_sum

(mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT_summary <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT %>%
    group_by(frame) %>%
    summarise(
      region_mean_prop = mean(prop)
    ))
### ^^^ CHANGE THIS ^^^ ###
###########################

###########################
### OPTION 2: DON'T POOL TREATMENTS (supplement)
### vvv CHANGE THIS vvv ###
#mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$prop <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$read_count / mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT$window_sum
#(mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT_summary <- mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT %>%
#  group_by(condition_combn, frame) %>%
#  summarise(
#    region_mean_prop = mean(prop)
#  ))
### ^^^ CHANGE THIS ^^^ ###
###########################


############################################################################################################
### PLOT Figure 2C: ZOOM-IN of ORF3a/ORF3c/ORF3d OLG region ###
(mapped_reads_by_readlength_sumSamples_PLOT <- ggplot(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET,
                                                      mapping = aes(x = window_center, y = 100 * prop_reads, color = frame)) +
   
   # LINES
   geom_vline(xintercept = 25697, linetype = "dotted", color = "grey", size = 0.25) +# linetype = "dashed", # ORF3d-2
   geom_vline(xintercept = 25596, linetype = "dotted", color = "grey", size = 0.25) +# linetype = "dashed", # ORF3d-2
   geom_vline(xintercept = 25697, linetype = "dotted", color = "grey", size = 0.25) +# linetype = "dashed", # ORF3d
   geom_vline(xintercept = 25524, linetype = "dotted", color = "grey", size = 0.25) +# linetype = "dashed", # ORF3d
   geom_vline(xintercept = 25582, linetype = "dotted", color = "grey", size = 0.25) +# linetype = "dashed", # ORF3c
   geom_vline(xintercept = 25457, linetype = "dotted", color = "grey", size = 0.25) + # linetype = "dashed", # ORF3c
   
   # ORF3d-2
   geom_rect(mapping = aes(xmin = 25596, xmax = 25697, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 25596, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF3d-2") + 
   
   # ORF3d
   geom_rect(mapping = aes(xmin = 25524, xmax = 25697, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 25524, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF3d") + 
   
   # ORF3c
   geom_rect(mapping = aes(xmin = 25457, xmax = 25582, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 25457, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF3c") + 
   
   # ORF3a
   geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
   geom_text(mapping = aes(x = 25457-50, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF3a") +
   
   geom_bar(data = filter(mapped_reads_by_readlength_sumSamples_LONG_SUBSET, frame == "Frame 1"), # redundant; otherwise triplicate plotting
            mapping = aes(x = window_center, y = window_sum / read_count_sum_max * 100 * prop_max), fill = "#E0E0E0", stat = "identity", inherit.aes = FALSE) + 
   
   # Region mean prop
   geom_hline(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT_summary, 
              mapping = aes(yintercept = 100 * region_mean_prop, color = frame), linetype = "dashed", size = rel(0.25)) +
   
   # ERROR RIBBON
   geom_ribbon(mapping = aes(ymin = 100 * binom_CI_low, ymax = 100 * binom_CI_high, fill = frame), alpha = 0.25, linetype = 0) +
   geom_line(size = rel(0.25)) +
   
   # -->CHANGE THIS<-- FACET if more than one treatment; COMMENT OUT otherwise
   #facet_grid(condition_combn ~ .) +
   xlab("Genome position") +
   ylab("Fraction reads (%)") +
   theme_bw() +
   theme(panel.grid = element_blank(),
         axis.title = element_text(size = rel(0.75)),
         legend.title = element_blank(),
         legend.text = element_text(size = rel(0.7)),
         strip.background = element_blank(),
         axis.ticks.y.right = element_blank(),
         axis.text.y.right = element_blank()
   ) +
   scale_color_manual(values = c("#00D39B", "#D48EB4", "#BDB10F")) +
   scale_fill_manual(values = c("#00D39B", "#D48EB4", "#F0E442")) +
   scale_x_continuous(expand = expand_scale(mult = c(0, 0))) +
   scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75), expand = expand_scale(mult = c(0, 0))))


###########################
### PLOT Figure 2—figure supplement 8, for S - 3a - E ###
#(mapped_reads_by_readlength_sumSamples_PLOT <- ggplot(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET,
#                                                      mapping = aes(x = window_center, y = 100 * prop_reads, color = frame)) +
#   
#   # E
#   geom_vline(xintercept = 26245, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 26245, xmax = Inf, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "white", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 26245, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "E") + 
#   
#   # ORF3b-4
#   geom_vline(xintercept = 26183, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 26183, xmax = 26281, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 26183, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3b-4") + 
#   
#   # ORF3b-3
#   geom_vline(xintercept = 26072, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 26072, xmax = 26170, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 26072, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3b-3") + 
#   
#   # ORF3b-2
#   geom_vline(xintercept = 25910, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25910, xmax = 25984, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25910, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3b-2") + 
#   
#   # ORF3b-1
#   geom_vline(xintercept = 25814, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25814, xmax = 25882, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25814, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3b") + 
#   
#   # ORF3a-2
#   geom_vline(xintercept = 25765, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25765, xmax = 26220, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25765, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3a-2") + 
#   
#   # ORF3d-2
#   geom_vline(xintercept = 25596, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25596, xmax = 25697, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25596, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3d-2") + 
#   
#   # ORF3d
#   geom_vline(xintercept = 25524, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25524, xmax = 25697, ymin = 85, ymax = 92), color = 'black', size = rel(0.25), fill = "#F0E442", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25524, y = (85 + 92) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3d") + 
#   
#   # ORF3c
#   geom_vline(xintercept = 25457, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25457, xmax = 25582, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25457, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3c") + 
#   
#   # ORF3a
#   geom_vline(xintercept = 25393, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 25393, xmax = 26220, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "#00D39B", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 25393, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "3a") + 
#   
#   # S
#   geom_rect(mapping = aes(xmin = -Inf, xmax = 25384, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "white", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = -Inf, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "S") + 
#   
#   geom_bar(data = filter(mapped_reads_by_readlength_sumSamples_LONG_SUBSET, frame == "Frame 1"), # otherwise triplicate plotting
#            mapping = aes(x = window_center, y = window_sum / read_count_sum_max * 100 * prop_max), fill = "#E0E0E0", stat = "identity", inherit.aes = FALSE) + 
#   
#   # Region mean prop
#   geom_hline(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT_summary, 
#              mapping = aes(yintercept = 100 * region_mean_prop, color = frame), linetype = "dashed", size = rel(0.25)) +
#   
#   # ERROR RIBBON
#   geom_ribbon(mapping = aes(ymin = 100 * binom_CI_low, ymax = 100 * binom_CI_high, fill = frame), alpha = 0.25, linetype = 0) +
#   
#   geom_line(size = rel(0.25)) +
#   facet_grid(condition_combn ~ .) +
#   xlab("Genome position") +
#   ylab("Fraction reads (%)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = rel(0.9)),
#         axis.title = element_text(size = rel(0.75)),
#         legend.title = element_blank(),
#         legend.text = element_text(size = rel(0.7)), 
#         axis.text = element_text(size = rel(0.7)), 
#         strip.background = element_blank(),
#         strip.text.y = element_text(angle = 0), 
#         axis.ticks.y.right = element_blank(),
#         axis.text.y.right = element_blank()
#   ) +
#   scale_color_manual(values = c("#00D39B", "#D48EB4", "#BDB10F")) +
#   scale_fill_manual(values = c("#00D39B", "#D48EB4", "#F0E442")) +
#   scale_x_continuous(expand = expand_scale(mult = c(0, 0))) +
#   scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75), expand = expand_scale(mult = c(0, 0)))) 



###########################
### PLOT sliding window for ORF8 - N - 9b - 9c - 10 ###

#ORF8:27894-28259, so frame 27894 %% 3 = 0 (Codon pos 1) # white
#N:28274-29533, so frame 28274 %% 3 = 2 (Codon pos 3), # structural gene blue: #71BFED
#ORF9b:28284-28577, so frame 28284 %% 3 = 0 (Codon pos 1) # overlapping gene burgundy: #D48EB4
#ORF9c:28734-28955, so frame 28734 %% 3 = 0 (Codon pos 1) # overlapping gene burgundy: #D48EB4
#ORF10:29558-29674, so frame 29558 %% 3 = 2 (Codon pos 3) # white

#(mapped_reads_by_readlength_sumSamples_PLOT <- ggplot(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET,
#                                                      mapping = aes(x = window_center, y = 100 * prop_reads, color = frame)) +
#   # ORF10
#   geom_vline(xintercept = 29558, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 29558, xmax = 29674, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "white", alpha = 0.1, inherit.aes = FALSE) +
#   geom_text(mapping = aes(x = 29558, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF10") + 
#   
#   # ORF9c
#   geom_vline(xintercept = 28734, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 28734, xmax = 28955, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 28734, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF9c") + 
#   
#   # ORF9b
#   geom_vline(xintercept = 28284, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 28284, xmax = 28577, ymin = 93, ymax = 100), color = 'black', size = rel(0.25), fill = "#D48EB4", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 28284, y = (93 + 100) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF9b") + 
#   
#   # N
#   geom_vline(xintercept = 28274, linetype = "dotted", color = "grey", size = 0.25) +
#   geom_rect(mapping = aes(xmin = 28274, xmax = 29533, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "#71BFED", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = 28274, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "N") + 
#   
#   # ORF8
#   geom_rect(mapping = aes(xmin = -Inf, xmax = 28259, ymin = 101, ymax = 108), color = 'black', size = rel(0.25), fill = "white", alpha = 0.1, inherit.aes = FALSE) + 
#   geom_text(mapping = aes(x = -Inf, y = (101 + 108) / 2), color = 'black', hjust = -0.02, vjust = 0.5, size = rel(2.5), fontface = 'italic', label = "ORF8") + 
#   
#   geom_bar(data = filter(mapped_reads_by_readlength_sumSamples_LONG_SUBSET, frame == "Frame 1"), # otherwise triplicate plotting
#            mapping = aes(x = window_center, y = window_sum / read_count_sum_max * 100 * prop_max), fill = "#E0E0E0", stat = "identity", inherit.aes = FALSE) + 
#   
#   # Region mean prop
#   geom_hline(data = mapped_reads_by_readlength_sumSamples_LONG_SUBSET_CONTEXT_summary, 
#              mapping = aes(yintercept = 100 * region_mean_prop, color = frame), linetype = "dashed", size = rel(0.25)) +
#   
#   # ERROR RIBBON
#   geom_ribbon(mapping = aes(ymin = 100 * binom_CI_low, ymax = 100 * binom_CI_high, fill = frame), alpha = 0.25, linetype = 0) +
#   
#   geom_line(size = rel(0.25)) +
#   facet_grid(condition_combn ~ .) +
#   xlab("Genome position") +
#   ylab("Fraction reads (%)") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         plot.title = element_text(hjust = 0.5, size = rel(0.9)),
#         axis.title = element_text(size = rel(0.75)),
#         legend.title = element_blank(),
#         legend.text = element_text(size = rel(0.7)), 
#         axis.text = element_text(size = rel(0.7)), 
#         strip.background = element_blank(),
#         strip.text.y = element_text(angle = 0),
#         axis.ticks.y.right = element_blank(),
#         axis.text.y.right = element_blank()
#   ) +
#   scale_color_manual(values = c("#00D39B", "#71BFED", "#D48EB4")) +
#   scale_fill_manual(values = c("#00D39B", "#71BFED", "#D48EB4")) +
#   scale_x_continuous(expand = expand_scale(mult = c(0, 0))) +
#   scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75), expand = expand_scale(mult = c(0, 0)))) 



########################################################################################################################
########################################################################################################################
### GET GENE START SITES AND ANNOTATE
(Wuhan_Hu_1_genes <- read_tsv("data/Wuhan_Hu_1_genes.tsv"))
unique(Wuhan_Hu_1_genes$gene)

genes_to_include <- c("ORF1ab_1", "ORF1ab_2", "S", 
                      "ORF3a", "ORF3c", "ORF3d", "ORF3d_short", "ORF3a_short",
                      "ORF3b_1", "ORF3b_2", "ORF3b_3", "ORF3b_4", 
                      "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                      "N", "ORF9b", "ORF9c", "ORF10")

gene_position_joiner <- tibble(position = unique(mapped_reads_by_readlength$position)) # now 1-based
gene_position_joiner # 29,884 x 1

## DID THIS ONCE OVERNIGHT AND SAVED
## Actually LOOPING the 2-million record file would take about 5 hours per gene, because it's R, so just loop the joiner
## STILL takes a LONG TIME; do it once and SAVE for reloading
#for (this_gene in unique(Wuhan_Hu_1_genes$gene)) {
#  if(this_gene %in% genes_to_include) {
#    this_gene_start <- Wuhan_Hu_1_genes[Wuhan_Hu_1_genes$gene == this_gene, ]$start
#    this_gene_end <- Wuhan_Hu_1_genes[Wuhan_Hu_1_genes$gene == this_gene, ]$end
#    
#    # Initialize a new column by the name of this gene
#    gene_position_joiner[[this_gene]] <- NA
#    
#    # Populate this column with genome's position relative to the start site of this gene, or NA if after its end
#    for (this_site in sort(unique(gene_position_joiner$position))) {
#      #this_gene_start is 29558
#      #this_site <- 29558
#      #this_site <- 29557
#      #this_site <- 29559
#      
#      if(this_site <= this_gene_end) {
#        gene_position_joiner[gene_position_joiner$position == this_site, ][[this_gene]] <- this_site - this_gene_start
#      } # else it's NA
#    }
#  }
#}

### SAVE
#write_tsv(gene_position_joiner, "data/gene_position_joiner.tsv")

########################################
# RELOAD
gene_position_joiner <- read_tsv("data/gene_position_joiner.tsv")

### JOIN --> SUPER FAST
mapped_reads_by_readlength <- left_join(x = mapped_reads_by_readlength, y = gene_position_joiner, by = "position")
#2,610,414 x 32 = SUCCESS!


################################################################################
### POOL THE TWO SAMPLES FOR EACH TREATMENT; CUT DOWN SIZE; SPECIFY READ LENGTH RANGE; LIMIT TO ONE UPSTREAM

### GENE

#ORF1ab_1:266-13468
#ORF1ab_2:13468-21555
#S:21563-25384

#ORF3a:25393-26220
#ORF3c:25457-25582
#ORF3d:25524-25697
#ORF3d_short:25596-25697
#ORF3a_short:25765-26220
#ORF3b:25814-25882

#E:26245-26472
#M:26523-27191
#ORF6:27202-27387
#ORF7a:27394-27759
#ORF7b:27756-27887
#ORF8:27894-28259
#N:28274-29533
#ORF9b:28284-28577
#ORF9c:28734-28955
#ORF10:29558-29674


### CHANGE THIS ### replace with name and start site of desired gene
this_gene_name <- "ORF3b" # ORF3a ORF3c ORF3d ORF3d-2 ORF3a-2 ORF3b
this_gene_start <- 25814 # 25393 25457 25524 25596
curr_conditions <- c("Harr / 5 hpi", "LTM / 5 hpi")
### CHANGE THIS ###

# Adds the two samples together
(riboseq_sumSamplesUPSTREAM <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), 
                                      read_length >= 27, read_length <= 32,
                                      
                                      # ALL PURPOSE
                                      #position >= this_gene_start - 25, position <= this_gene_start + 25 # centered
                                      position >= this_gene_start - 25, position <= this_gene_start # upstream only
) %>% 
    group_by(condition_combn, read_length, position) %>%
    summarise(
      sample_count = n(),
      read_count_sum = sum(read_count)
    )) 

# sum ALL READS IN THE START UPSTREAM for each condition
(riboseq_readSumPerCombnUPSTREAM <- riboseq_sumSamplesUPSTREAM %>%
    group_by(condition_combn) %>% 
    summarise(
      read_sum_UPSTREAMxCONDITION = sum(read_count_sum)
    ))

# LEFT JOIN TOTAL READS FOR EACH COMBN; calculate reads per million mapped reads
riboseq_sumSamplesUPSTREAM <- left_join(x = riboseq_sumSamplesUPSTREAM, y = riboseq_readSumPerCombnUPSTREAM, by = c("condition_combn")) 
riboseq_sumSamplesUPSTREAM

# prop of reads in UPSTREAMxCONDITION
riboseq_sumSamplesUPSTREAM$prop_reads_UPSTREAMxCONDITION <- riboseq_sumSamplesUPSTREAM$read_count_sum / riboseq_sumSamplesUPSTREAM$read_sum_UPSTREAMxCONDITION
riboseq_sumSamplesUPSTREAM

### PREP BARS DUMMY DATAFRAME
(riboseq_sumSamplesUPSTREAM_summary <- riboseq_sumSamplesUPSTREAM %>%
    group_by(condition_combn, position) %>%
    summarise(
      position_total_depth = sum(read_count_sum)
    ))

(position_total_depth_max <- max(filter(riboseq_sumSamplesUPSTREAM_summary, condition_combn %in% curr_conditions)$position_total_depth))
(prop_reads_UPSTREAMxCONDITION_max <- max(filter(riboseq_sumSamplesUPSTREAM, condition_combn %in% curr_conditions)$prop_reads_UPSTREAMxCONDITION))

### PLOT Figure 2—figure supplement 4
### CHANGE THIS ###
(riboseq_sumSamplesUPSTREAM_PLOT <- 
    ggplot(data = filter(riboseq_sumSamplesUPSTREAM, condition_combn %in% curr_conditions), 
           mapping = aes(x = position, y = 100 * prop_reads_UPSTREAMxCONDITION, color = factor(read_length))) + 
    
    # P SITE OFFSET at -13, -12, -11 upstream start sites
    geom_rect(mapping = aes(xmin = this_gene_start - 13.5, xmax = this_gene_start - 10.5, ymin = -Inf, ymax = Inf), fill = "#F7FFD6", inherit.aes = FALSE) + # FBFFEB F7FFD6 F3FFC2
    
    # ORF9b ONLY: P SITE OFFSET of N at -23, -22, -21 upstream start sites <-- CHANGE THIS
    #geom_rect(mapping = aes(xmin = this_gene_start - 23.5, xmax = this_gene_start - 20.5, ymin = -Inf, ymax = Inf), fill = "#F7FFD6", inherit.aes = FALSE) + # FBFFEB F7FFD6 F3FFC2
    
    # GENE START
    geom_vline(xintercept = this_gene_start, linetype = "dashed", color = "grey", size = 0.3) +
    
    # BAR version -- worry about bars later  
    geom_bar(data = filter(riboseq_sumSamplesUPSTREAM_summary, condition_combn %in% curr_conditions),
             mapping = aes(x = position, y = position_total_depth / position_total_depth_max * 100 * prop_reads_UPSTREAMxCONDITION_max), fill = "#E0E0E0", stat = "identity", inherit.aes = FALSE) + 
    
    # POINT/LINE version
    geom_point(size = rel(0.8)) +
    geom_line(size = rel(0.25)) +
    
    ggtitle(this_gene_name) + # "ORF3c" "ORF3d" "ORF3d-2" ORF3a etc.
    
    facet_grid(condition_combn ~ ., scales = "free_y") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "italic", size = rel(1)),
          legend.text = element_text(size = rel(1.15)),
          axis.text.x = element_text(size = rel(1.15)), 
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = rel(1.15)),
          axis.title.y = element_blank(),
          panel.border = element_rect(),
          strip.background = element_blank(),
          strip.text.x = element_text(size = rel(1.15)),
          strip.text.y = element_text(size = rel(1.15), angle = 0)) +
    scale_color_manual(values = rev(hue_pal()(6)), name = "Read length (bp)") +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0))) +
    scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.1)))
)



####################
### POOL THE TWO SAMPLES FOR EACH TREATMENT; CUT DOWN SIZE
(riboseq_sumSamples <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition)) %>% # , read_length >=26, read_length <= 32
   group_by(condition_combn, read_length, position) %>%
   summarise(
     sample_count = n(),
     read_count_sum = sum(read_count)
   )) # should cut in half (two samples per condition)

(riboseq_readSumPerCombn <- riboseq_sumSamples %>%
    group_by(condition_combn, read_length) %>%
    summarise(
      read_sum_genome = sum(read_count_sum)
    ))

# LEFT JOIN TOTAL READS FOR EACH COMBN AND READ LENGTH; calculate reads per million mapped reads
riboseq_sumSamples <- left_join(x = riboseq_sumSamples, y = riboseq_readSumPerCombn, by = c("condition_combn", "read_length"))
riboseq_sumSamples$reads_per_million <- riboseq_sumSamples$read_count_sum / riboseq_sumSamples$read_sum_genome * 1e6
riboseq_sumSamples

# LEFT JOIN RELATIVE POSITIONS OF GENES
riboseq_sumSamples <- left_join(x = riboseq_sumSamples, y = gene_position_joiner, by = "position")
riboseq_sumSamples # 1,725,599 x 29

# make long
riboseq_sumSamples_long <- riboseq_sumSamples %>%
  pivot_longer(cols = genes_to_include, names_to = "gene", values_to = "relative_position")
riboseq_sumSamples_long

# filter to within 40bp of SOME start site, so it's not too big
(riboseq_sumSamples_long <- filter(riboseq_sumSamples_long, relative_position >= -40, relative_position <= 40))
# 159,391 x 9

# Factor genes 5' to 3'
genes_to_include_labels <- c("1a", "1b", "S", 
                             "3a", "3c", "3d", "3d-2", "3a-2",
                             "3b-1", "3b-2", "3b-3", "3b-4", 
                             "E", "M", "6", "7a", "7b", "8",
                             "N", "9b", "9c", "10")
riboseq_sumSamples_long$gene <- factor(riboseq_sumSamples_long$gene, levels = genes_to_include, labels = genes_to_include_labels)

genes_to_include_labels_subset <- c("1a", "1b", "S", 
                                    "3a", "3c", "3d", "3d-2", "3a-2",
                                    "3b-1", #"3b-2", "3b-3", "3b-4", 
                                    "E", "M", "6", "7a", "7b", "8",
                                    "N", "9b", "9c", "10")

####################
# PLOT Harr at 5 and 24 hpi, or whatever combination is of interest

# CHX / 5 hpi -- CHX / 24 hpi -- Harr / 5 hpi -- Harr / 24 hpi -- LTM / 5 hpi -- LTM / 24 hpi -- mRNA / 5 hpi -- mRNA / 24 hpi
(riboseq_sumSamples_long_PLOT <- ggplot(data = filter(riboseq_sumSamples_long, 
                                                      gene %in% genes_to_include_labels_subset,
                                                      relative_position >= -20, relative_position <= 20,
                                                      read_length == 30,
                                                      condition_combn %in% c("Harr / 5 hpi", "Harr / 24 hpi")), 
                                        mapping = aes(x = relative_position, y = reads_per_million)) + # read_count_sum
   
   # 12 bp upstream
   geom_rect(mapping = aes(xmin = -12.5, xmax = -11.5, ymin = -Inf, ymax = Inf), fill = "#F3FFC2", alpha = 0.1, inherit.aes = FALSE) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.3) +
   
   # BAR version
   geom_bar(stat = "identity") +
   
   ggtitle("Read length 30 bp - free y axis") + # "Harr / 5 hpi"
   
   facet_grid(gene ~ condition_combn, scales = "free_y") + 
   xlab("Distance from start (nt)") +
   ylab("Reads per million mapped") +
   theme_bw() +
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5),
         strip.text = element_text(size = 10),
         axis.text.x = element_text(size = 9.5), 
         axis.title.x = element_text(size = 10),
         axis.text.y = element_text(size = 9.5),
         axis.title.y = element_text(size = 10),
         panel.border = element_rect(),
         strip.background = element_blank(),
         strip.text.y = element_text(size = rel(0.9), angle = 0)) +
   scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.05))) 
)



####################
# plot ALL conditions for 29-31bp pooled lengths
(riboseq_sumSamples2 <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), read_length >= 29, read_length <= 31) %>%
   group_by(condition_combn, position) %>%
   summarise(
     record_count = n(),
     read_count_sum = sum(read_count)
   )) 

(riboseq_readSumPerCombn2 <- riboseq_sumSamples2 %>%
    group_by(condition_combn) %>%
    summarise(
      read_sum_genome = sum(read_count_sum)
    ))

# LEFT JOIN TOTAL READS FOR EACH COMBN AND READ LENGTH; calculate reads per million mapped reads
riboseq_sumSamples2 <- left_join(x = riboseq_sumSamples2, y = riboseq_readSumPerCombn2, by = c("condition_combn"))
riboseq_sumSamples2$reads_per_million <- riboseq_sumSamples2$read_count_sum / riboseq_sumSamples2$read_sum_genome * 1e6
riboseq_sumSamples2

# LEFT JOIN RELATIVE POSITIONS OF GENES
riboseq_sumSamples2 <- left_join(x = riboseq_sumSamples2, y = gene_position_joiner, by = "position")
riboseq_sumSamples2 

# make long
riboseq_sumSamples2_long <- riboseq_sumSamples2 %>%
  pivot_longer(cols = genes_to_include, names_to = "gene", values_to = "relative_position")
riboseq_sumSamples2_long

# filter to within 40bp of SOME start site, so it's not too big
(riboseq_sumSamples2_long <- filter(riboseq_sumSamples2_long, relative_position >= -40, relative_position <= 40))

# Factor genes 5' to 3'
genes_to_include_labels <- c("1ab", "1b", "S", 
                             "3a", "3c", "3d", "3d-2", "3a-2",
                             "3b-1", "3b-2", "3b-3", "3b-4", 
                             "E", "M", "6", "7a", "7b", "8",
                             "N", "9b", "9c", "10")
riboseq_sumSamples2_long$gene <- factor(riboseq_sumSamples2_long$gene, levels = genes_to_include, labels = genes_to_include_labels)

genes_to_include_labels_subset <- c("1ab", "1b", "S", 
                                    "3a", "3c", "3d", "3d-2", "3a-2",
                                    "3b-1", #"3b-2", "3b-3", "3b-4", 
                                    "E", "M", "6", "7a", "7b", "8",
                                    "N", "9b", "9c", "10")

####################
# plot ALL condition combinations for SPECIFIED READ LENGTH RANGES

# CHX / 5 hpi -- CHX / 24 hpi -- Harr / 5 hpi -- Harr / 24 hpi -- LTM / 5 hpi -- LTM / 24 hpi -- mRNA / 5 hpi -- mRNA / 24 hpi
(riboseq_sumSamples2_long_PLOT <- ggplot(data = filter(riboseq_sumSamples2_long, 
                                                       gene %in% genes_to_include_labels_subset,
                                                       gene != "1b",
                                                       relative_position >= -20, relative_position <= 20,
                                                       condition_combn %in% c("Harr / 5 hpi", "LTM / 5 hpi")),
                                         mapping = aes(x = relative_position, y = reads_per_million)) +
   
   # 12 bp upstream
   geom_rect(mapping = aes(xmin = -13.5, xmax = -10.5, ymin = -Inf, ymax = Inf), fill = "#F3FFC2", alpha = 0.1, inherit.aes = FALSE) +
   geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.3) +
   geom_bar(stat = "identity") + # , position = position_dodge()
   
   facet_grid(gene ~ condition_combn, scales = "free_y") +
   xlab("Distance from start (nt)") +
   ylab("Reads per million mapped") + 
   theme_bw() +
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5),
         strip.text = element_text(size = 10),
         axis.text.x = element_text(size = 9.5), 
         axis.title.x = element_text(size = 10),
         axis.text.y = element_text(size = 9.5),
         axis.title.y = element_text(size = 10),
         panel.border = element_rect(),
         strip.background = element_blank(),
         strip.text.y = element_text(size = rel(0.9), angle = 0)) +
   scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.05))) 
)




########################################
### BOOK-KEEPING: WHAT ARE THE SAMPLES FOR EACH TREATMENT?
########################################
mapped_reads_by_readlength %>%
  group_by(condition_combn, sample) %>%
  summarise(
    count = n()
  )

#condition_combn sample       count
#1 CHX / 5 hpi     SRR11713356  76388
#2 CHX / 5 hpi     SRR11713357  78242
#3 CHX / 24 hpi    SRR11713364  70230
#4 CHX / 24 hpi    SRR11713365 327831
#5 Harr / 5 hpi    SRR11713360  57931
#6 Harr / 5 hpi    SRR11713361  50150
#7 Harr / 24 hpi   SRR11713368 313066
#8 Harr / 24 hpi   SRR11713369 258834
#9 LTM / 5 hpi     SRR11713358  64767
#10 LTM / 5 hpi     SRR11713359  75772
#11 LTM / 24 hpi    SRR11713366 160926
#12 LTM / 24 hpi    SRR11713367 229330
#13 mRNA / 5 hpi    SRR11713354 122305
#14 mRNA / 5 hpi    SRR11713355 124381
#15 mRNA / 24 hpi   SRR11713362 292810
#16 mRNA / 24 hpi   SRR11713363 307451


########################################
### FOUR-SAMPLE VERSION
########################################

####################
# plot ALL for 29-31bp pooled lengths, "Harr / 5 hpi" and "LTM / 5 hpi" only
(riboseq_sumSamples4 <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), read_length >= 29, read_length <= 31,
                               condition_combn %in% c("Harr / 5 hpi", "LTM / 5 hpi")) %>%
   group_by(sample, position) %>%
   summarise(
     record_count = n(),
     read_count_sum = sum(read_count)
   )) 

(riboseq_readSumPerCombn4 <- riboseq_sumSamples4 %>%
    group_by(sample) %>%
    summarise(
      read_sum_genome = sum(read_count_sum)
    ))
#sample      read_sum_genome
#SRR11713358          344461
#SRR11713359          569576
#SRR11713360          223651
#SRR11713361          122341

# LEFT JOIN TOTAL READS FOR EACH COMBN AND READ LENGTH; calculate reads per million mapped reads
riboseq_sumSamples4 <- left_join(x = riboseq_sumSamples4, y = riboseq_readSumPerCombn4, by = c("sample"))
riboseq_sumSamples4$reads_per_million <- riboseq_sumSamples4$read_count_sum / riboseq_sumSamples4$read_sum_genome * 1e6
riboseq_sumSamples4

# LEFT JOIN RELATIVE POSITIONS OF GENES
riboseq_sumSamples4 <- left_join(x = riboseq_sumSamples4, y = gene_position_joiner, by = "position")
riboseq_sumSamples4 # 47,129 x 28

# make long
riboseq_sumSamples4_long <- riboseq_sumSamples4 %>%
  pivot_longer(cols = genes_to_include, names_to = "gene", values_to = "relative_position")
riboseq_sumSamples4_long

# filter to within 40bp of SOME start site, so it's not too big
(riboseq_sumSamples4_long <- filter(riboseq_sumSamples4_long, relative_position >= -20, relative_position <= 20))

# FACTOR genes 5' to 3'
genes_to_include_labels2 <- c("ORF1ab", "ORF1b", "S", 
                              "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                              "ORF3b", "ORF3b-2", "ORF3b-3", "ORF3b-4", 
                              "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                              "N", "ORF9b", "ORF9c", "ORF10")
riboseq_sumSamples4_long$gene <- factor(riboseq_sumSamples4_long$gene, levels = genes_to_include, labels = genes_to_include_labels2)

genes_to_include_labels2_subset <- c("ORF1ab", "ORF1b", "S", 
                                     "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                                     "ORF3b", #"ORF3b-2", "ORF3b-3", "ORF3b-4", 
                                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                     "N", "ORF9b", "ORF9c", "ORF10")

# FIND MAX BAR (red in Figure 2A, Figure 2—figure supplement 1)
riboseq_sumSamples4_long$peak <- "Non-Peak"
for (this_gene in unique(riboseq_sumSamples4_long$gene)) {
  #this_gene <- '1b'
  this_gene_data <- filter(riboseq_sumSamples4_long, gene == this_gene)
  
  for (this_sample in unique(this_gene_data$sample)) {
    #this_sample <- 'SRR11713358'
    this_gene_sample_data <- filter(this_gene_data, sample == this_sample)
    this_gene_max_depth <- max(this_gene_sample_data$reads_per_million)
    this_gene_max_depth_position <- this_gene_sample_data[this_gene_sample_data$reads_per_million == this_gene_max_depth, ]$relative_position
    
    # Cat
    cat(this_gene, this_sample, this_gene_max_depth_position, "\n")
    
    # Add peak position height
    riboseq_sumSamples4_long[riboseq_sumSamples4_long$gene == this_gene & 
                               riboseq_sumSamples4_long$sample == this_sample &
                               riboseq_sumSamples4_long$relative_position == this_gene_max_depth_position, ]$peak = "Peak" 
  }
}

riboseq_sumSamples4_long$peak <- factor(riboseq_sumSamples4_long$peak, levels = c("Peak", "Non-Peak"))

# Factor sample order Harr first
riboseq_sumSamples4_long$sample <- factor(riboseq_sumSamples4_long$sample, levels = c("SRR11713360", "SRR11713361", "SRR11713358", "SRR11713359"))


####################
# PLOT Figure 2—figure supplement 1, ALL FOUR SAMPLES with the SPECIFIED READ LENGTH RANGES

(riboseq_sumSamples4_long_PLOT <- ggplot(data = filter(riboseq_sumSamples4_long, 
                                                       gene %in% genes_to_include_labels2_subset,
                                                       gene != "1b", gene != "ORF1b",
                                                       relative_position >= -20, relative_position <= 20),
                                         mapping = aes(x = relative_position, y = reads_per_million, fill = peak)) + 
   
   geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey", size = 0.5) + 
   geom_bar(stat = "identity") + 
   facet_grid(gene ~ sample, scales = "free_y") + 
   xlab("Distance from start (nt)") +
   ylab("Reads per million mapped") + 
   theme_bw() +
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5),
         legend.position = 'none',
         strip.text = element_text(size = 10),
         panel.border = element_rect(),
         strip.background = element_blank(),
         strip.text.y = element_text(size = rel(0.9), angle = 0, face = "italic")) +
   scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.05))) + 
   scale_fill_manual(values = c("red", "grey35"))
)



########################################
### FOUR-SAMPLE VERSION, **24 HPI** to see how it compares
########################################

####################
# plot ALL for 29-31bp pooled lengths, "Harr / 5 hpi" and "LTM / 5 hpi" only
(riboseq_sumSamples24hpi <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), read_length >= 29, read_length <= 31,
                                   condition_combn %in% c("Harr / 24 hpi", "LTM / 24 hpi")) %>%
   group_by(sample, position) %>%
   summarise(
     record_count = n(),
     read_count_sum = sum(read_count)
   )) 

(riboseq_readSumPerCombn24hpi <- riboseq_sumSamples24hpi %>%
    group_by(sample) %>%
    summarise(
      read_sum_genome = sum(read_count_sum)
    ))
#sample      read_sum_genome
#SRR11713366          833006
#SRR11713367         2463553
#SRR11713368        10245313
#SRR11713369         3828880

# LEFT JOIN TOTAL READS FOR EACH COMBN AND READ LENGTH; calculate reads per million mapped reads
riboseq_sumSamples24hpi <- left_join(x = riboseq_sumSamples24hpi, y = riboseq_readSumPerCombn24hpi, by = c("sample"))
riboseq_sumSamples24hpi$reads_per_million <- riboseq_sumSamples24hpi$read_count_sum / riboseq_sumSamples24hpi$read_sum_genome * 1e6
riboseq_sumSamples24hpi

# LEFT JOIN RELATIVE POSITIONS OF GENES
riboseq_sumSamples24hpi <- left_join(x = riboseq_sumSamples24hpi, y = gene_position_joiner, by = "position")
riboseq_sumSamples24hpi # 94,780 x 28

# make long
riboseq_sumSamples24hpi_long <- riboseq_sumSamples24hpi %>%
  pivot_longer(cols = genes_to_include, names_to = "gene", values_to = "relative_position")
riboseq_sumSamples24hpi_long
# 2,085,160 x 8

# filter to within 40bp of SOME start site, so it's not too big
(riboseq_sumSamples24hpi_long <- filter(riboseq_sumSamples24hpi_long, relative_position >= -20, relative_position <= 20))
# 3,415 x 8

# FACTOR genes 5' to 3'
genes_to_include_labels2 <- c("ORF1ab", "ORF1b", "S", 
                              "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                              "ORF3b", "ORF3b-2", "ORF3b-3", "ORF3b-4", 
                              "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                              "N", "ORF9b", "ORF9c", "ORF10")
riboseq_sumSamples24hpi_long$gene <- factor(riboseq_sumSamples24hpi_long$gene, levels = genes_to_include, labels = genes_to_include_labels2)

genes_to_include_labels2_subset <- c("ORF1ab", "ORF1b", "S", 
                                     "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                                     "ORF3b", #"ORF3b-2", "ORF3b-3", "ORF3b-4", 
                                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                     "N", "ORF9b", "ORF9c", "ORF10")

# FIND MAX BAR
riboseq_sumSamples24hpi_long$peak <- "Non-Peak"
for (this_gene in unique(riboseq_sumSamples24hpi_long$gene)) {
  #this_gene <- '1b'
  this_gene_data <- filter(riboseq_sumSamples24hpi_long, gene == this_gene)
  
  for (this_sample in unique(this_gene_data$sample)) {
    #this_sample <- 'SRR11713358'
    this_gene_sample_data <- filter(this_gene_data, sample == this_sample)
    this_gene_max_depth <- max(this_gene_sample_data$reads_per_million)
    this_gene_max_depth_position <- this_gene_sample_data[this_gene_sample_data$reads_per_million == this_gene_max_depth, ]$relative_position
    
    # Cat
    cat(this_gene, this_sample, this_gene_max_depth_position, "\n")
    
    # Add peak position height
    riboseq_sumSamples24hpi_long[riboseq_sumSamples24hpi_long$gene == this_gene & 
                                   riboseq_sumSamples24hpi_long$sample == this_sample &
                                   riboseq_sumSamples24hpi_long$relative_position == this_gene_max_depth_position, ]$peak = "Peak" 
  }
}

riboseq_sumSamples24hpi_long$peak <- factor(riboseq_sumSamples24hpi_long$peak, levels = c("Peak", "Non-Peak"))

# Factor sample order Harr first
riboseq_sumSamples24hpi_long$sample <- factor(riboseq_sumSamples24hpi_long$sample, levels = c("SRR11713368", "SRR11713369", "SRR11713366", "SRR11713367"))


####################
# PLOT a version of Figure 2—figure supplement 1 at 24 hpi

(riboseq_sumSamples24hpi_long_PLOT <- ggplot(data = filter(riboseq_sumSamples24hpi_long, 
                                                           gene %in% genes_to_include_labels2_subset,
                                                           gene != "1b", gene != "ORF1b",
                                                           relative_position >= -20, relative_position <= 20),
                                             mapping = aes(x = relative_position, y = reads_per_million, fill = peak)) + 
   
   geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey", size = 0.5) + 
   geom_bar(stat = "identity") + 
   facet_grid(gene ~ sample, scales = "free_y") +
   xlab("Distance from start (nt)") +
   ylab("Reads per million mapped") + 
   theme_bw() +
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5),
         legend.position = 'none',
         strip.text = element_text(size = 10),
         panel.border = element_rect(),
         strip.background = element_blank(),
         strip.text.y = element_text(size = rel(0.9), angle = 0, face = "italic")) +
   scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.05))) +
   scale_fill_manual(values = c("red", "grey35"))
)



########################################
### Plot ALL for 29-31bp pooled lengths, BUT NOW POOL Harr/5hpi and LTM/5hpi for a summary Figure 2

### POOL HARR SAMPLES
(riboseq_sumSamples3_Harr <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), read_length >= 29, read_length <= 31,
                                    condition_combn == "Harr / 5 hpi") %>%
   group_by(position) %>%
   summarise(
     record_count = n(),
     read_count_sum = sum(read_count)
   )) 

# calculate total number of mapped reads across genome
read_sum_genome_POOLED_Harr <- sum(riboseq_sumSamples3_Harr$read_count_sum)

# calculate reads per million mapped reads
riboseq_sumSamples3_Harr$reads_per_million <- riboseq_sumSamples3_Harr$read_count_sum / read_sum_genome_POOLED_Harr * 1e6
riboseq_sumSamples3_Harr$condition <- "Harr_5hpi"

### POOL LTM SAMPLES
(riboseq_sumSamples3_LTM <- filter(dplyr::select(mapped_reads_by_readlength, -chromosome, -frame, -condition), read_length >= 29, read_length <= 31,
                                   condition_combn == "LTM / 5 hpi") %>%
    group_by(position) %>%
    summarise(
      record_count = n(),
      read_count_sum = sum(read_count)
    )) 

# calculate total number of mapped reads across genome
read_sum_genome_POOLED_LTM <- sum(riboseq_sumSamples3_LTM$read_count_sum)

# calculate reads per million mapped reads
riboseq_sumSamples3_LTM$reads_per_million <- riboseq_sumSamples3_LTM$read_count_sum / read_sum_genome_POOLED_LTM * 1e6
riboseq_sumSamples3_LTM$condition <- "LTM_5hpi"

### MEAN of Harr & LTM
riboseq_sumSamples3 <- rbind(riboseq_sumSamples3_Harr, riboseq_sumSamples3_LTM) 
(riboseq_sumSamples3 <- riboseq_sumSamples3 %>%
    group_by(position) %>%
    summarise(
      record_count = sum(record_count),
      read_count_sum = sum(read_count_sum),
      reads_per_million = mean(reads_per_million)
    )) # 18,960 x 4

# LEFT JOIN RELATIVE POSITIONS OF GENES
riboseq_sumSamples3 <- left_join(x = riboseq_sumSamples3, y = gene_position_joiner, by = "position")
riboseq_sumSamples3

# make long
riboseq_sumSamples3_long <- riboseq_sumSamples3 %>%
  pivot_longer(cols = genes_to_include, names_to = "gene", values_to = "relative_position")
riboseq_sumSamples3_long

# filter to EXACT start site region
(riboseq_sumSamples3_long <- filter(riboseq_sumSamples3_long, relative_position >= -20, relative_position <= 20))

# Factor genes 5' to 3'
genes_to_include_groupByOLG <- c("ORF1ab_1", "ORF1ab_2", "S", 
                                 "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                 "ORF10",
                                 "ORF3a", "ORF3c", "ORF3d", "ORF3d_short", "ORF3a_short",
                                 "ORF3b_1", "ORF3b_2", "ORF3b_3", "ORF3b_4", 
                                 "N", "ORF9b", "ORF9c")
genes_to_include_groupByOLG_labels <- c("1ab", "1b", "S", 
                                        "E", "M", "6", "7a", "7b", "8",
                                        "10",
                                        "3a", "3c", "3d", "3d-2", "3a-2",
                                        "3b", "3b-2", "3b-3", "3b-4", 
                                        "N", "9b", "9c")
genes_to_include_groupByOLG_labels2 <- c("ORF1ab", "ORF1b", "S", 
                                         "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                         "ORF10",
                                         "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                                         "ORF3b", "ORF3b-2", "ORF3b-3", "ORF3b-4", 
                                         "N", "ORF9b", "ORF9c")
riboseq_sumSamples3_long$gene <- factor(riboseq_sumSamples3_long$gene, levels = genes_to_include_groupByOLG, labels = genes_to_include_groupByOLG_labels2)

genes_to_include_labels_subset2 <- c("ORF1ab", "ORF1b", "S", 
                                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                     "ORF10",
                                     "ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                                     "ORF3b", #"ORF3b-2", "ORF3b-3", "ORF3b-4", 
                                     "N", "ORF9b", "ORF9c")

# Factor for gridding
genes_to_include_labels_OLG2 <- c("ORF3a", "ORF3c", "ORF3d", "ORF3d-2", "ORF3a-2",
                                  "ORF3b", "ORF3b-2", "ORF3b-3", "ORF3b-4", 
                                  "N", "ORF9b", "ORF9c")
genes_to_include_labels_NonOLG2 <- c("ORF1ab", "ORF1b", "S", 
                                     "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8",
                                     "ORF10")

riboseq_sumSamples3_long$region_type <- NA
riboseq_sumSamples3_long[riboseq_sumSamples3_long$gene %in% genes_to_include_labels_OLG2, ]$region_type <- "OLG loci"
riboseq_sumSamples3_long[riboseq_sumSamples3_long$gene %in% genes_to_include_labels_NonOLG2, ]$region_type <- "Non-OLG loci"
riboseq_sumSamples3_long$region_type <- factor(riboseq_sumSamples3_long$region_type, levels = c("Non-OLG loci", "OLG loci"))

# DUMMY DATA FRAME to highlight the site in this region with the max depth
riboseq_sumSamples3_long$peak <- "Non-Peak"
for (this_gene in unique(riboseq_sumSamples3_long$gene)) {
  #this_gene <- '1b'
  this_gene_data <- filter(riboseq_sumSamples3_long, gene == this_gene)
  this_gene_max_depth <- max(this_gene_data$reads_per_million)
  this_gene_max_depth_position <- this_gene_data[this_gene_data$reads_per_million == this_gene_max_depth, ]$relative_position
  
  # Cat
  cat(this_gene, this_gene_max_depth_position, "\n")
  
  # Add peak position height
  riboseq_sumSamples3_long[riboseq_sumSamples3_long$gene == this_gene & riboseq_sumSamples3_long$relative_position == this_gene_max_depth_position, ]$peak = "Peak"
}

riboseq_sumSamples3_long$peak <- factor(riboseq_sumSamples3_long$peak, levels = c("Peak", "Non-Peak"))


####################
# PLOT Figure 2C, POOLED: Harr / 5 hpi AND LTM / 5 hpi
(riboseq_sumSamples3_long_PLOT <- ggplot(data = filter(riboseq_sumSamples3_long, 
                                                       gene %in% genes_to_include_labels_subset2,
                                                       gene != "1b", gene != "ORF1b"),
                                         mapping = aes(x = relative_position, y = reads_per_million, fill = peak)) +
   
   # gene start site
   geom_vline(xintercept = 0, linetype = "dotted", color = "darkgrey", size = 0.5) +
   geom_bar(stat = "identity") + 
   
   facet_wrap(. ~ gene, ncol = 2, scales = "free_y", dir = "v", drop = TRUE) + 
   xlab("Distance from start (nt)") + 
   ylab("Reads per million mapped") + 
   theme_bw() +
   theme(panel.grid = element_blank(),
         legend.position = 'none',
         strip.background = element_blank(),
         strip.text = element_text(size = rel(0.9), margin = margin(0, 0, 1, 0, "pt"),  face = 'italic'), # size = 10, 
         axis.title = element_text(size = rel(0.9)),
         axis.text = element_text(size = rel(0.9)),
         axis.ticks = element_line(size = rel(0.9)),
         panel.border = element_rect(size = rel(0.9))
   ) +
   scale_y_continuous(breaks = pretty_breaks(2), expand = expand_scale(mult = c(0, 0.1))) +
   scale_fill_manual(values = c("red", "grey35"))
)


