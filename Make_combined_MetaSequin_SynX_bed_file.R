#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make combined MetaSequin SynX GTF
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
library('rtracklayer')
library('GenomicAlignments')

# Load MetaGenome Sequins
# --------------------------------------------------------------------------

seqins_ranges <- import.bed("data/reference_files/metasequin_corrected.bed", genome = 'MetaSequin_Synx')

# Load SynX bed file
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations_CN_mask.bed", genome = 'MetaSequin_Synx')
synx_ranges@elementMetadata <- synx_ranges@elementMetadata[,c('name'), drop = FALSE]

# Combine bed files
combined_ranges <- c(seqins_ranges, synx_ranges)

# Export to gtf
export(combined_ranges, "data/reference_files/MetaSequin_synx_CN_mask.bed")
