#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make single-feature bed files for all performance elements in Synx
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
require('scales')

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')

# Subset to performance elements
synx_ranges_perf <- synx_ranges[grepl("^P[0-9]|^PH[0-9]", synx_ranges$name),]

# Make bed file for each performance element
for(feature in synx_ranges_perf$name) {
  export(synx_ranges_perf[synx_ranges_perf$name == feature, ], paste0("data/reference_files/", feature, "_one_line.bed"))
}
