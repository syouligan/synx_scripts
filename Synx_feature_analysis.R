#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX features ONT-DNA
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

# Calculate the feature characteristics
# --------------------------------------------------------------------------

# Calculate GATC content of each feature
synx_ranges$feature_length <- str_length(synx_ranges$SEQUENCE)
for(Base in c('A', 'T', 'G', 'C')){
  synx_ranges@elementMetadata[,paste0(Base, '_base_count')] <- list(str_count(synx_ranges$SEQUENCE, pattern = Base))
  
  synx_ranges@elementMetadata[,paste0(Base, '_base_percent')] <- list(synx_ranges@elementMetadata[,paste0(Base, '_base_count')] / synx_ranges$feature_length * 100)
}  

# Make GRanges object for each feature group
synx_ranges_by_types <- split(synx_ranges, synx_ranges$TYPE)
promoters <- synx_ranges_by_types$Promoter
restriction <- synx_ranges_by_types$Restriction_site
polyA <- synx_ranges_by_types$`Poly-A_tail`
performance <- synx_ranges_by_types$Sequencing_performance_unit
performanceGC <- performance[! grepl("PH", performance$name),]
performanceHP <- performance[grepl("PH", performance$name),]
quantitative_ladder <- synx_ranges_by_types$Quantitative_ladder_unit
size_ladder <- synx_ranges_by_types$Size_ladder_unit

# Calculate error rates for each base (SUBS and INDELS)
# --------------------------------------------------------------------------

# Load sequencing error rates
ONTDNA_pileup <- read.table("project_results/ONT_DNA/barcode06.sorted.bam.bed.tsv", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

# Find reference mapping rates
get_correct <- function(df){
  bases <- c('A', 'T', 'G', 'C')
  ref <- df['REF_NT']
  correct <- bases[bases == ref]
  sum(as.numeric(df[correct]))
}

ONTDNA_pileup$Depth <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del + ONTDNA_pileup$Coverage
ONTDNA_pileup$REF_counts <- apply(ONTDNA_pileup, 1, get_correct)
ONTDNA_pileup$REF_rate <- ONTDNA_pileup$REF_counts / ONTDNA_pileup$Depth

# Find substitution rates
get_subs <- function(df){
  bases <- c('A', 'T', 'G', 'C')
  ref <- df['REF_NT']
  subs <- bases[!(bases == ref)]
  sum(as.numeric(df[subs]))
}

ONTDNA_pileup$SUB_counts <- apply(ONTDNA_pileup, 1, get_subs)
ONTDNA_pileup$SUB_rate <- ONTDNA_pileup$SUB_counts / ONTDNA_pileup$Depth

# Find INDEL rates
ONTDNA_pileup$Ins_rate <- ONTDNA_pileup$Ins / ONTDNA_pileup$Depth
ONTDNA_pileup$Del_rate <- ONTDNA_pileup$Del / ONTDNA_pileup$Depth
ONTDNA_pileup$INDEL_counts <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del
ONTDNA_pileup$INDEL_rate <- ONTDNA_pileup$INDEL_counts / ONTDNA_pileup$Depth

# Find error rates
ONTDNA_pileup$Error_counts <- ONTDNA_pileup$SUB_counts + ONTDNA_pileup$INDEL_counts
ONTDNA_pileup$Error_rate <- ONTDNA_pileup$Error_counts / ONTDNA_pileup$Depth

# Plot error rates across GC performance features
# --------------------------------------------------------------------------

performanceGC$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, performanceGC@ranges))
performanceGC$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, performanceGC@ranges))
performanceGC$Error_upper <- performanceGC$Error_mean + performanceGC$Error_sd
performanceGC$Error_lower <- performanceGC$Error_mean - performanceGC$Error_sd

performanceGC$SUB_mean <- mean(extractList(ONTDNA_pileup$SUB_rate, performanceGC@ranges))
performanceGC$SUB_sd <- sd(extractList(ONTDNA_pileup$SUB_rate, performanceGC@ranges))
performanceGC$SUB_upper <- performanceGC$SUB_mean + performanceGC$SUB_sd
performanceGC$SUB_lower <- performanceGC$SUB_mean - performanceGC$SUB_sd

performanceGC$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_rate, performanceGC@ranges))
performanceGC$INDEL_sd <- sd(extractList(ONTDNA_pileup$INDEL_rate, performanceGC@ranges))
performanceGC$INDEL_upper <- performanceGC$INDEL_mean + performanceGC$INDEL_sd
performanceGC$INDEL_lower <- performanceGC$INDEL_mean - performanceGC$INDEL_sd

performanceGC$INS_mean <- mean(extractList(ONTDNA_pileup$Ins_rate, performanceGC@ranges))
performanceGC$INS_sd <- sd(extractList(ONTDNA_pileup$Ins_rate, performanceGC@ranges))
performanceGC$INS_upper <- performanceGC$INS_mean + performanceGC$INS_sd
performanceGC$INS_lower <- performanceGC$INS_mean - performanceGC$INS_sd

performanceGC$DEL_mean <- mean(extractList(ONTDNA_pileup$Del_rate, performanceGC@ranges))
performanceGC$DEL_sd <- sd(extractList(ONTDNA_pileup$Del_rate, performanceGC@ranges))
performanceGC$DEL_upper <- performanceGC$DEL_mean + performanceGC$DEL_sd
performanceGC$DEL_lower <- performanceGC$DEL_mean - performanceGC$DEL_sd

errors <- c('Error', 'SUB', 'INS', 'DEL')

# Plot sequencing errors across performance units
performanceGC <- performanceGC[order(performanceGC$name), ]

# Plot GC content of GC performance units
ggplot(data.frame(performanceGC@elementMetadata)) +
  geom_bar(aes(y = GC_pct, x = name), stat = 'identity') +
  theme_minimal() +
  ggsave(paste0("project_results/ONT_DNA/SynX_performance_GC_content.pdf"))

# Error rates (individual) across GC-performance units
for(i in errors) {
  Mean <- paste0(i, "_mean")
  Lower <- paste0(i, "_lower")
  Upper <- paste0(i, "_upper")
  
  ggplot(data.frame(performanceGC@elementMetadata)) +
    geom_line(aes_string(x = 'GC_pct', y = Mean)) +
    geom_ribbon(aes_string(x = 'GC_pct', y = Mean, ymin=Lower, ymax=Upper), fill="blue", alpha=0.2) +
    xlab("GC content") +
    ylab(i) +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/SynX_performance_GC_", i, ".pdf"))
}

# Error rates (all) across GC-performance units
GC_mean_long <- melt(data.frame(performanceGC@elementMetadata[,c(paste0(errors, "_mean"), 'GC_pct')]), id = 'GC_pct')
GC_lower_long <- melt(data.frame(performanceGC@elementMetadata[,c(paste0(errors, "_lower"), 'GC_pct')]), id = 'GC_pct')
GC_upper_long <- melt(data.frame(performanceGC@elementMetadata[,c(paste0(errors, "_upper"), 'GC_pct')]), id = 'GC_pct')

GC_error_long <- data.frame("GC_pct" = GC_mean_long$GC_pct, "Type" = GC_mean_long$variable, "Mean" = GC_mean_long$value, "Upper" = GC_upper_long$value, "Lower" = GC_lower_long$value)

ggplot(GC_error_long) +
  geom_line(aes(x=GC_pct, y = Mean, colour = Type)) +
  geom_ribbon(aes(x = GC_pct, y = Mean, ymin=Lower, ymax=Upper, fill = Type), alpha=0.2) +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  scale_fill_npg(palette = c("nrc"), alpha = 1) +
  xlab("GC content") +
  ylab('Error rate') +
  theme_classic() +
  ggsave("project_results/ONT_DNA/SynX_performance_GC_all_errors.pdf")

# Plot error rates across homopolymers performance features
# --------------------------------------------------------------------------

# Import bed file of homopolymer sites
synx_hp_ranges <- import.bed("project_results/ONT_DNA/synx_homopolymers.bed", genome = 'SynX')

synx_hp_ranges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, synx_hp_ranges@ranges))
synx_hp_ranges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, synx_hp_ranges@ranges))
synx_hp_ranges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, synx_hp_ranges@ranges))
synx_hp_ranges$Error_upper <- synx_hp_ranges$Error_mean + synx_hp_ranges$Error_sd
synx_hp_ranges$Error_lower <- synx_hp_ranges$Error_mean - synx_hp_ranges$Error_sd

synx_hp_ranges$SUB_freq <- sum(extractList(ONTDNA_pileup$SUB_counts, synx_hp_ranges@ranges))
synx_hp_ranges$SUB_mean <- mean(extractList(ONTDNA_pileup$SUB_rate, synx_hp_ranges@ranges))
synx_hp_ranges$SUB_sd <- sd(extractList(ONTDNA_pileup$SUB_rate, synx_hp_ranges@ranges))
synx_hp_ranges$SUB_upper <- synx_hp_ranges$SUB_mean + synx_hp_ranges$SUB_sd
synx_hp_ranges$SUB_lower <- synx_hp_ranges$SUB_mean - synx_hp_ranges$SUB_sd

synx_hp_ranges$INDEL_freq <- sum(extractList(ONTDNA_pileup$INDEL_counts, synx_hp_ranges@ranges))
synx_hp_ranges$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_rate, synx_hp_ranges@ranges))
synx_hp_ranges$INDEL_sd <- sd(extractList(ONTDNA_pileup$INDEL_rate, synx_hp_ranges@ranges))
synx_hp_ranges$INDEL_upper <- synx_hp_ranges$INDEL_mean + synx_hp_ranges$INDEL_sd
synx_hp_ranges$INDEL_lower <- synx_hp_ranges$INDEL_mean - synx_hp_ranges$INDEL_sd

synx_hp_ranges$INS_freq <- sum(extractList(ONTDNA_pileup$Ins, synx_hp_ranges@ranges))
synx_hp_ranges$INS_mean <- mean(extractList(ONTDNA_pileup$Ins_rate, synx_hp_ranges@ranges))
synx_hp_ranges$INS_sd <- sd(extractList(ONTDNA_pileup$Ins_rate, synx_hp_ranges@ranges))
synx_hp_ranges$INS_upper <- synx_hp_ranges$INS_mean + synx_hp_ranges$INS_sd
synx_hp_ranges$INS_lower <- synx_hp_ranges$INS_mean - synx_hp_ranges$INS_sd

synx_hp_ranges$DEL_freq <- sum(extractList(ONTDNA_pileup$Del, synx_hp_ranges@ranges))
synx_hp_ranges$DEL_mean <- mean(extractList(ONTDNA_pileup$Del_rate, synx_hp_ranges@ranges))
synx_hp_ranges$DEL_sd <- sd(extractList(ONTDNA_pileup$Del_rate, synx_hp_ranges@ranges))
synx_hp_ranges$DEL_upper <- synx_hp_ranges$DEL_mean + synx_hp_ranges$DEL_sd
synx_hp_ranges$DEL_lower <- synx_hp_ranges$DEL_mean - synx_hp_ranges$DEL_sd

# Plot homopolymer frequency over whole synx
ggplot(data.frame(table(synx_hp_ranges$score)), aes(x=Var1, y=log2(Freq))) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency (log2)") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_frequency.pdf")

# Subset to just HP performance elements
performanceHP1_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP[performanceHP$name == 'PH1',])
performanceHP2_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP[performanceHP$name == 'PH1',])
performanceHPtot_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP)

# Plot homopolymer frequency over performance element
ggplot(data.frame(table(performanceHP1_ranges$score)), aes(x=Var1, y=Freq)) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element1_frequency_total.pdf")

ggplot(data.frame(table(performanceHP2_ranges$score)), aes(x=Var1, y=Freq)) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element2_frequency_total.pdf")

ggplot(data.frame(table(performanceHPtot_ranges$score)), aes(x=Var1, y=Freq)) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_total_frequency_total.pdf")

# Plot error rate by homopolymer length
for( i in paste0(errors, '_freq')) {
  ggplot(data.frame(performanceHPtot_ranges@elementMetadata), ) +
    geom_boxplot(aes_string(x='score', group='score', y = i), color = "grey80") +
    geom_point(aes_string(x='score', y = i, colour = 'name')) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    xlab("Homopolymer (bp)") +
    ylab(i) +
    theme_classic() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_error_rate_",  i,".pdf"))
}

# Error rates (all) across HP-performance units
HP_fraction_long <- melt(data.frame(performanceHPtot_ranges@elementMetadata[,c(paste0(errors, "_freq"), 'score')]), id = 'score')

ggplot(HP_fraction_long) +
  geom_boxplot(aes(x=score, y = value, colour = variable, group = score)) +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  scale_fill_npg(palette = c("nrc"), alpha = 1) +
  xlab("HP content") +
  ylab('Error rate') +
  theme_classic() +
  facet_wrap(~variable) +
  ggsave("project_results/ONT_DNA/SynX_performance_HP_all_errors.pdf")


# Identify 31mers in quantitative ladder for quantification using jellyfish
# --------------------------------------------------------------------------

# Import pileup stats for read data aligned to copy-number festures only
ONTDNA_CN_pileup <- read.table("project_results/ONT_DNA/barcode06.synx_cn_only.bam.bed.tsv", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
synx_cn_ranges <- import.bed("data/reference_files/synx_cn_only.bed", genome = 'SynX_cn_only')

ONTDNA_CN_pileup$Copy_number <- c(rep("cn1", synx_cn_ranges@ranges@width[1]), rep("cn2", synx_cn_ranges@ranges@width[2]), rep("cn3", synx_cn_ranges@ranges@width[3]), rep("cn4", synx_cn_ranges@ranges@width[4]))

ggplot(ONTDNA_CN_pileup) +
  geom_density(aes(x = Coverage, color = Copy_number, fill = Copy_number, group = Copy_number)) +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  scale_fill_npg(palette = c("nrc"), alpha = 0.3) +
  xlab("Copy number") +
  ylab("Coverage (density)") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Quantitative_ladder_coverage_density_with_median.pdf")

ggplot(ONTDNA_CN_pileup) +
  geom_boxplot(aes(y = Coverage, group = Copy_number, x = Copy_number, fill = Copy_number)) +
  scale_fill_npg(palette = c("nrc"), alpha = 1) +
  xlab("Copy number") +
  ylab("Per-base coverage") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Quantitative_ladder_coverage_boxplot.pdf")

