#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Compare PhiX sequence composition
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
library('GenomicAlignments')
library('sarlacc')
library('Biostrings')

# Set colour palette
npg_cols <- pal_npg("nrc")(7)

phix_seq <- unlist(str_split(read.csv("data/reference_files/phix.fa", skip = 1, header = FALSE)[1,1], pattern = ""))

# Calculate kmer representation in PhiX sequence
# --------------------------------------------------------------------------

# PhiX sequences (5' - 3', forward and reverse)
phix_seq_F <- phix_seq
phix_seq_R <- rev(phix_seq)
phix_seq_C <- toupper(spgs::complement(phix_seq))
phix_seq_RC <- toupper(spgs::reverseComplement(phix_seq))

# Plot kmer representation in PhiX sequence
KD_k <- function(x){paste(x, collapse = "")}

phix_proportion <- data.frame()
for(i in 1:10) {
  # kmers in forward sequence
  seq_F <- runner(x = phix_seq_F, k = i, f = KD_k)
  seq_F <- seq_F[str_length(seq_F) == i]
  seq_F_unique <- unique(seq_F)
  
  # kmers in reverse sequence
  seq_R <- runner(x = phix_seq_R, k = i, f = KD_k)
  seq_R <- seq_R[str_length(seq_R) == i]
  seq_R_unique <- unique(seq_R)
  
  # kmers in complement sequence
  seq_C <- runner(x = phix_seq_C, k = i, f = KD_k)
  seq_C <- seq_C[str_length(seq_C) == i]
  seq_C_unique <- unique(seq_C)
  
  # kmers in reverse complement sequence
  seq_RC <- runner(x = phix_seq_RC, k = i, f = KD_k)
  seq_RC <- seq_RC[str_length(seq_RC) == i]
  seq_RC_unique <- unique(seq_RC)
  
  # Total unique kmers in forward and reverse complement
  seq_FRC_unique <- unique(c(seq_F_unique, seq_RC_unique))
  
  # Total unique kmers in all directions
  seq_total_unique <- unique(c(seq_F_unique, seq_R_unique, seq_C_unique, seq_RC_unique))
  
  # Kmer data table
  kmer_counts <- data.frame('Name' = c('Possible', 'Total', 'Forward + RC', 'Forward', 'Reverse', 'Complement', 'Reverse complement'),
                            'Number' = c(4^i, length(seq_total_unique), length(seq_FRC_unique), length(seq_F_unique), length(seq_R_unique), length(seq_C_unique), length(seq_RC_unique)))
  kmer_counts$Name <- factor(kmer_counts$Name, levels = kmer_counts$Name)
  
  # Plot
  ggplot(kmer_counts) +
    geom_bar(aes(x = Name, y = Number, fill = Name, color = Name), stat = 'identity') +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Kmer_representation_", i, "mer_PhiX.pdf"))
  
  phix_proportion <- rbind(phix_proportion, c(i, length(seq_FRC_unique)/4^i))
}

# Line plot of proportion by kmer length
colnames(phix_proportion) <- c("Length", "Proportion")
phix_proportion$Genome <- "PhiX"
ggplot(phix_proportion, aes(x=Length, y=Proportion)) +
  geom_line(group = 1) +
  xlab("Kmer length") +
  ylab("Proportion of total possible kmers") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Kmer_representation_vs_length_PhiX.pdf")

# Calculate kmer representation in synx sequence
# --------------------------------------------------------------------------

ONTDNA_pileup <- read.table("project_results/ONT_DNA/barcode06.sorted.bam.bed.tsv", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
# Synx sequences (5' - 3', forward and reverse)
synx_seq_F <- ONTDNA_pileup$REF_NT
synx_seq_R <- rev(ONTDNA_pileup$REF_NT)
synx_seq_C <- toupper(spgs::complement(ONTDNA_pileup$REF_NT))
synx_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup$REF_NT))

# Plot kmer representation in synx sequence
KD_k <- function(x){paste(x, collapse = "")}

synx_proportion <- data.frame()
for(i in 1:10) {
  # kmers in forward sequence
  seq_F <- runner(x = synx_seq_F, k = i, f = KD_k)
  seq_F <- seq_F[str_length(seq_F) == i]
  seq_F_unique <- unique(seq_F)
  
  # kmers in reverse sequence
  seq_R <- runner(x = synx_seq_R, k = i, f = KD_k)
  seq_R <- seq_R[str_length(seq_R) == i]
  seq_R_unique <- unique(seq_R)
  
  # kmers in complement sequence
  seq_C <- runner(x = synx_seq_C, k = i, f = KD_k)
  seq_C <- seq_C[str_length(seq_C) == i]
  seq_C_unique <- unique(seq_C)
  
  # kmers in reverse complement sequence
  seq_RC <- runner(x = synx_seq_RC, k = i, f = KD_k)
  seq_RC <- seq_RC[str_length(seq_RC) == i]
  seq_RC_unique <- unique(seq_RC)
  
  # Total unique kmers in forward and reverse complement
  seq_FRC_unique <- unique(c(seq_F_unique, seq_RC_unique))
  
  # Total unique kmers in all directions
  seq_total_unique <- unique(c(seq_F_unique, seq_R_unique, seq_C_unique, seq_RC_unique))
  
  # Kmer data table
  kmer_counts <- data.frame('Name' = c('Possible', 'Total', 'Forward + RC', 'Forward', 'Reverse', 'Complement', 'Reverse complement'),
                            'Number' = c(4^i, length(seq_total_unique), length(seq_FRC_unique), length(seq_F_unique), length(seq_R_unique), length(seq_C_unique), length(seq_RC_unique)))
  kmer_counts$Name <- factor(kmer_counts$Name, levels = kmer_counts$Name)
  
  # Plot
  ggplot(kmer_counts) +
    geom_bar(aes(x = Name, y = Number, fill = Name, color = Name), stat = 'identity') +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Kmer_representation_", i, "mer.pdf"))
  
  synx_proportion <- rbind(synx_proportion, c(i, length(seq_FRC_unique)/4^i))
}

# Line plot of proportion by kmer length
colnames(synx_proportion) <- c("Length", "Proportion")
synx_proportion$Genome <- "SynX"
ggplot(synx_proportion, aes(x=Length, y=Proportion)) +
  geom_line(group = 1) +
  xlab("Kmer length") +
  ylab("Proportion of total possible kmers") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Kmer_representation_vs_length_SynX.pdf")

# SynX vs PhiX kmer representation
proportion <- rbind(synx_proportion, phix_proportion)
ggplot(proportion, aes(x=Length, y=Proportion)) +
  geom_point(aes(color = Genome)) +
  geom_line(aes(color = Genome)) +
  scale_x_continuous(breaks = 1:10) +
  xlab("Kmer length") +
  ylab("Proportion of total possible kmers") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Kmer_representation_vs_length_SynXvsPhiX.pdf")

# Plot frequency of homopolymers by length in PhiX
# --------------------------------------------------------------------------

# Import bed file of homopolymer sites
synx_phix_hp_ranges <- import.bed("project_results/ONT_DNA/synx_phix_homopolymers.bed", genome = 'Synx_PhiX')
synx_phix_hp_ranges$genome <- as.character(synx_phix_hp_ranges@seqnames)

# Plot homopolymer frequency over synx vs phix synx
phix_hp_freq <- data.frame(table(synx_phix_hp_ranges[synx_phix_hp_ranges$genome == 'PhiX',]$score))
phix_hp_freq$Genome <- "PhiX"
synx_hp_freq <- data.frame(table(synx_phix_hp_ranges[synx_phix_hp_ranges$genome == 'SynX',]$score))
synx_hp_freq$Genome <- "SynX"
hp_freq <- rbind(phix_hp_freq, synx_hp_freq)

ggplot(hp_freq, aes(x=Var1, y=log2(Freq), color = Genome)) +
  geom_point(aes(color = Genome)) +
  geom_line(aes(group = Genome)) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency (log2)") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_phix_vs_synx_total_frequency.pdf")





