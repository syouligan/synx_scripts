#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of transcript lengths in digested datasets
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/tim_projects/synx/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/tim/mer/scott/synx/")
  place <- "timmer"
}

# Load libraries
library('stringr')
library('dplyr')
library('ggplot2')
library('rtracklayer')
library('Gviz')
library('ggsci')
library('runner')
library('matrixStats')
library('spgs')
library('ggformula')
library('seqinr')
library('plyr')
library('reshape2')

# Set colour palette
npg_cols <- pal_npg("nrc")(7)

# Load feature data
# --------------------------------------------------------------------------

digested_data_sets <- c('ONT_DNA_BamHI', 'ONT_DNA_EcoRI', 'ONT_DNA_HindIII', 'ONT_DNA_Fragmentase')

length_data <- data.frame()
for (ds in digested_data_sets) {
  tmp <- read.table(paste0("project_results/", ds, "/", ds, ".read_length.tsv"), header = FALSE)
  colnames(tmp) <- c("Length", "Count")
  tmp$Dataset <- ds
  length_data <- rbind(length_data, tmp)
  }

# Remove transcripts with length = 1
length_data <- length_data[length_data$Length != 1, ]

# Plot fragment size distribution
ggplot(length_data, aes(x = log(Length), y = Count)) +
  geom_line(aes(color = Dataset, fill = Dataset)) +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  scale_fill_npg(palette = c("nrc"), alpha = 0.3) +
  xlab("Transcript Length (log)") +
  ylab("Count") +
  theme_classic() +
  facet_wrap(~Dataset, scales = 'free_y', nrow = 4, ncol = 1) +
  theme(legend.position = "none") +
  ggsave(paste0("project_results/Digested_fragment_lengths.pdf"))
