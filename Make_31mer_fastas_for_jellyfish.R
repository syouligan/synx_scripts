#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make fasta files for kmers in copy number elements for use with Jellyfish
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

# Set colour palette
npg_cols <- pal_npg("nrc")(7)

# Load feature data
# --------------------------------------------------------------------------

# Load synx bed file
synx_bed <- read.table("data/reference_files/corrected_synx_annotations.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')

# Load synx feature annotations and annotate ranges
synx_annotations <- read.table("design/custom_features.annotations.tab", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")

idx <- match(synx_ranges$name, synx_annotations$NAME)
synx_ranges$TYPE <- synx_annotations$TYPE [idx]
synx_ranges$TYPE <- gsub(" ", "_", synx_ranges$TYPE)

synx_ranges$ORIGIN <- synx_annotations$ORIGIN [idx]
synx_ranges$COPY <- synx_annotations$COPY [idx]
synx_ranges$SIZE <- synx_annotations$SIZE [idx]
synx_ranges$GC_pct <- synx_annotations$GC... [idx]
synx_ranges$SEQUENCE <- synx_annotations$SEQUENCE [idx]

# Create unique feature idenitfiers
synx_ranges$Identifier <- paste(synx_ranges$name, synx_ranges@ranges@start, sep = '_')

# Identify all 31mers in SynX sequence
# --------------------------------------------------------------------------

ONTDNA_pileup <- read.table(paste0("project_results/ONT_DNA/barcode06.sorted.bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

# Calculate all unique 31mers in Synx
KD_k <- function(x){paste(x, collapse = "")}
synx_seq_F <- ONTDNA_pileup$REF_NT
synx_seq_F_31mers <- runner(x = synx_seq_F, k = 31, f = KD_k)
synx_seq_total_31mers <- synx_seq_F_31mers[str_length(synx_seq_F_31mers) == 31]
synx_seq_total_31mers <- unique(synx_seq_total_31mers)
write.fasta(sequences = as.list(synx_seq_total_31mers), names = synx_seq_total_31mers, as.string = TRUE, file.out = paste0('project_results/ONT_DNA/jellyfish/All_31mers.fasta'))

# Identify 31mers in quantitative ladder for quantification using jellyfish
# --------------------------------------------------------------------------
quantitative_ladder <- synx_ranges[grepl('cn', synx_ranges$name), ]
KD_k <- function(x){paste(x, collapse = "")}

for(cn in c('1cn', '2cn', '3cn', '4cn')){
  cnseq <- quantitative_ladder[quantitative_ladder$name == cn,]$SEQUENCE[1]
  cn_31mers <- runner(x = unlist(strsplit(cnseq, split="")), k = 31, f = KD_k)
  cn_31mers <- cn_31mers[str_length(cn_31mers) == 31]
  write.fasta(sequences = as.list(cn_31mers), names = paste0(cn, '_', cn_31mers), as.string = TRUE, file.out = paste0('project_results/ONT_DNA/jellyfish/', cn, '_31mers.fasta'))
}




