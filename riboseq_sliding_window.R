##############################################################################################################
### Ribo-seq sliding windows

# Mac
suppressMessages(library(package = readr))
suppressMessages(library(package = stringr))
suppressMessages(library(package = dplyr))
suppressMessages(library(package = tidyr))


############################################################################################################
############################################################################################################
### GATHER PARAMETERS FROM COMMAND LINE
ARGV <- commandArgs(trailingOnly = T)

kill_script <- FALSE

if(! (length(ARGV) == 4)) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[1], pattern = "\\w")) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[2], pattern = "\\w")) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[3], pattern = "\\d")) {
  kill_script <- TRUE
} else if(! str_detect(string = ARGV[4], pattern = "\\d")) {
  kill_script <- TRUE
} # could add more conditions later

if(kill_script) {
  cat("\n\n### WARNING: there must be 4 command line arguments, in this order:\n")
  cat("    (1) INFILE_FRAMES: file with condition metadata for samples\n")
  cat("    (2) INFILE_READS: file with read data\n")
  cat("    (3) WIN_SIZE: size of the sliding window\n")
  cat("    (4) READ_LEN: length of reads to consider\n\n")
  quit(save = 'no', status = 1, runLast = TRUE)
}

### PARAMETERS ###
INFILE_FRAMES <- as.character(ARGV[1])
INFILE_READS <- as.character(ARGV[2])
WIN_SIZE <- as.integer(ARGV[3]) 
READ_LEN <- as.integer(ARGV[4]) 
##################

# Produce some helpful warning messages
if(! (WIN_SIZE >= 1 && WIN_SIZE <= 1000)) {
  cat("### WARNING: WIN_SIZE must be in the range [1,1000]. Using: 29.\n")
  WIN_SIZE <- 29
} else {
  
}

if(! (READ_LEN >= 1 && READ_LEN <= 100)) {
  cat("### WARNING: READ_LEN must be in the range [1,100]. Using: 29.\n")
  READ_LEN <- 29
}

# Product output file name
OUTFILE_NAME <- str_replace(string = INFILE_FRAMES, pattern = "\\w+.\\w+$", replacement = "")
OUTFILE_NAME <- paste0(OUTFILE_NAME, "mapped_reads_by_readlength_sumSamples_win", WIN_SIZE, "read", READ_LEN, ".tsv")

cat("\n\nOUTPUT will be written to: ", OUTFILE_NAME, "\n")


##############################################################################################################
### INPUT FILES
## File 1
suppressMessages(frames_table <- read_tsv(INFILE_FRAMES,
                                          col_names = c("prop_frame0", "prop_frame1", "prop_frame2", "read_length", 
                                                        "read_count_frame0", "frame_frame0", 
                                                        "read_count_frame1", "frame_frame1", 
                                                        "read_count_frame2", "frame_frame2", 
                                                        "sample", "condition")))

## File 2
suppressMessages(mapped_reads_by_readlength <- read_tsv(INFILE_READS,
                                       col_names = c("sample", "read_length", "chromosome", "position", "frame", "read_count")))

# Convert to 1-based
mapped_reads_by_readlength$position <- mapped_reads_by_readlength$position + 1


##############################################################################################################
### JOIN sample data

# construct sample metadata translator
sample_metadata_joiner <- dplyr::select(frames_table, sample, condition)

# Extract time
sample_metadata_joiner$time <- sample_metadata_joiner$condition
sample_metadata_joiner$time <- str_replace(string = sample_metadata_joiner$time,
                                            pattern = "\\w+(\\d+)hr", replacement = "\\1")
sample_metadata_joiner$time <- factor(sample_metadata_joiner$time,
                                       levels = c("5", "4"),
                                       labels = c("5 hpi", "24 hpi"))

# Extract treatment
sample_metadata_joiner$treatment <- as.character(sample_metadata_joiner$condition)
sample_metadata_joiner[str_detect(string = sample_metadata_joiner$treatment, "chx"), ]$treatment <- "CHX"
sample_metadata_joiner[str_detect(string = sample_metadata_joiner$treatment, "harr"), ]$treatment <- "Harr"
sample_metadata_joiner[str_detect(string = sample_metadata_joiner$treatment, "ltm"), ]$treatment <- "LTM"
sample_metadata_joiner[str_detect(string = sample_metadata_joiner$treatment, "mrna"), ]$treatment <- "mRNA"

# unique only
sample_metadata_joiner 
sample_metadata_joiner <- unique(sample_metadata_joiner)

# Join
mapped_reads_by_readlength <- left_join(x = mapped_reads_by_readlength, y = sample_metadata_joiner, by = "sample")

# Another condition option
mapped_reads_by_readlength$condition_combn <- paste0(mapped_reads_by_readlength$treatment, " / ", mapped_reads_by_readlength$time)
mapped_reads_by_readlength$condition_combn <- factor(mapped_reads_by_readlength$condition_combn,
                                                     levels = c("CHX / 5 hpi", "CHX / 24 hpi",
                                                                "Harr / 5 hpi", "Harr / 24 hpi",
                                                                "LTM / 5 hpi", "LTM / 24 hpi",
                                                                "mRNA / 5 hpi", "mRNA / 24 hpi")) # each of these has two samples

# Factor
mapped_reads_by_readlength$frame <- factor(mapped_reads_by_readlength$frame,
                                           levels = c(0, 1, 2),
                                           labels = c('Codon pos 1', 'Codon pos 2', 'Codon pos 3'))


####################################################################################################
### SLIDING WINDOWS prop frame analysis, pooling the samples for a given condition

# Summarize read counts by condition, position, frame (pool samples)
mapped_reads_by_readlength_sumSamples <- filter(mapped_reads_by_readlength, read_length == READ_LEN) %>%
    group_by(condition_combn, position, frame) %>%
    summarise(
      read_count_sumSamples = sum(read_count)
    )

### PERFORM SLIDING WINDOWS,  summing each window's total reads for each codon position

# Prepare new columns
mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_1 <- NA
mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_2 <- NA
mapped_reads_by_readlength_sumSamples$window_sum_codon_pos_3 <- NA

### TAKES FOREVER ###
# specific to read lengths specified above # TODO: vectorize the approach
for(this_condition in unique(mapped_reads_by_readlength_sumSamples$condition_combn)) {
  #this_condition <- "CHX / 5 hpi"
  condition_data <- mapped_reads_by_readlength_sumSamples[mapped_reads_by_readlength_sumSamples$condition_combn == this_condition, ]
  
  for(i in sort(unique(condition_data$position))) { # no genes at tail end of genome, so we're good
    #i <- 1
    
    # Extract window; analyze
    window_data <- condition_data[condition_data$position >= i & condition_data$position <= (i + WIN_SIZE - 1), ]
    
    ### Add to data frame
    # Codon pos 1
    mapped_reads_by_readlength_sumSamples[mapped_reads_by_readlength_sumSamples$condition_combn == this_condition &
                                            mapped_reads_by_readlength_sumSamples$position == i, ]$window_sum_codon_pos_1 <- 
      sum(window_data[window_data$frame == "Codon pos 1", ]$read_count_sumSamples)
    
    # Codon pos 2
    mapped_reads_by_readlength_sumSamples[mapped_reads_by_readlength_sumSamples$condition_combn == this_condition &
                                            mapped_reads_by_readlength_sumSamples$position == i, ]$window_sum_codon_pos_2 <- 
      sum(window_data[window_data$frame == "Codon pos 2", ]$read_count_sumSamples)
    
    # Codon pos 3
    mapped_reads_by_readlength_sumSamples[mapped_reads_by_readlength_sumSamples$condition_combn == this_condition &
                                            mapped_reads_by_readlength_sumSamples$position == i, ]$window_sum_codon_pos_3 <- 
      sum(window_data[window_data$frame == "Codon pos 3", ]$read_count_sumSamples)
    
  } # end last window
} # end last condition_combn

#View(mapped_reads_by_readlength_sumSamples)

### SAVE ###
write_tsv(mapped_reads_by_readlength_sumSamples, OUTFILE_NAME)


