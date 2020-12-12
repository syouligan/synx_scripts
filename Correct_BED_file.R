#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Correct bed file annotation of SynX elements
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

# Correct missing E5 restriction site (Spel) from start of sequence and export bam
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/synx_annotations.bed", genome = 'SynX')

# Add 6 base-pairs (restriction site length) to feature start
synx_ranges@ranges@start <-as.integer(synx_ranges@ranges@start + 6)
synx_ranges$thick <- synx_ranges@ranges

# Create granges for E5 site
E5_ranges <- synx_ranges[synx_ranges$name == 'E5', ]
E5_ranges@ranges@start <- as.integer(1)
E5_ranges$thick <- E5_ranges@ranges

# Prepend corrected bed file and save
corrected_ranges <- c(E5_ranges, synx_ranges)
export.bed(corrected_ranges, "data/reference_files/corrected_synx_annotations.bed")
