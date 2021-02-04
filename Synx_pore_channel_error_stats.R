#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX each read in ONT-DNA
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


# Set colour palette
npg_cols <- pal_npg("nrc")(7)

# Calculate error rates for for each read and each pore
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')

# Subset to performance elements
synx_ranges_perf <- synx_ranges[grepl("^P[0-9]|^PH[0-9]", synx_ranges$name),]

# Calculate error rates with time and pore over each performance element
# --------------------------------------------------------------------------

for(feature in c(synx_ranges_perf$name, "SynX")) {
  # Load sequencing error rates
  reads_1linebed <- read.table(paste0("project_results/ONT_DNA/", feature, "_perRead/", feature, "_one_line.barcode06.sorted.INTERSECT.perRead.tsv"), header = TRUE, sep=";", stringsAsFactors=FALSE, quote="")
  reads_1linebed$start_time <- as.POSIXct(reads_1linebed$start_time, format="%Y-%m-%dT%H:%M:%SZ")
  reads_1linebed$ch <- factor(as.character(reads_1linebed$ch))
  
  # Plot length distribution
  ggplot(reads_1linebed) +
    geom_density(aes(x=Length))
  
  # Calculate per base error frequency
  reads_1linebed$GC_pct <- (reads_1linebed$G + reads_1linebed$C)/reads_1linebed$Length*100
  reads_1linebed$Error_freq <- reads_1linebed$Error/reads_1linebed$Length
  reads_1linebed$SUB_freq <- reads_1linebed$Mismatch/reads_1linebed$Length
  reads_1linebed$INS_freq <- reads_1linebed$Ins/reads_1linebed$Length
  reads_1linebed$DEL_freq <- reads_1linebed$Del/reads_1linebed$Length
  
  # Subset to reads those around the mean (+/- 1 SD)
  if(feature != 'SynX') {
    feat_mean <- mean(reads_1linebed$Length)
    feat_sd <- sd(reads_1linebed$Length)
    reads_1linebed <- reads_1linebed[reads_1linebed$Length < (feat_mean + feat_sd) & reads_1linebed$Length > (feat_mean - feat_sd),]
  } else {
    print('Feature = SynX')
  }
  
  # Order channels by means per base error frequency across all transcripts
  sum_stats <- reads_1linebed %>%
    group_by(ch) %>%
    dplyr::summarize(Mean_error_freq = mean(Error_freq), count(ch)) %>%
    mutate(ch = as.character(ch)) %>%
    arrange(Mean_error_freq)
  reads_1linebed$ch <- factor(reads_1linebed$ch, levels = sum_stats$ch)
  
  # PLot error rate over time
  error_time <- melt(reads_1linebed %>%
                       select(start_time, Error_freq, SUB_freq, INS_freq, DEL_freq), id.vars = 'start_time')
  
  ggplot(error_time, aes(x = start_time, y = value)) +
    geom_line(aes(color = variable), alpha = 1) +
    geom_smooth(alpha = 0.5, color = 'grey50', fill = 'grey90', size = 0.5) +
    scale_colour_npg(palette = c("nrc")) +
    scale_fill_npg(palette = c("nrc")) +
    scale_x_datetime(labels = date_format("%Y-%m-%d %H"), date_breaks = "4 hours") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_classic() +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/ONT_DNA/", feature, "_perRead/Synx_error_rates_over_time_", feature, ".pdf"))
  
  # Plot error rate by channel over time
  error_ch_time <- melt(reads_1linebed %>%
                          select(ch, start_time, Error_freq, SUB_freq, INS_freq, DEL_freq),
                        id.vars = c("ch", "start_time"),
                        measure.vars = c('Error_freq', 'SUB_freq', 'INS_freq', 'DEL_freq'))
  
  ggplot(error_ch_time, aes(x = start_time, y = ch)) +
    geom_point(aes(color = value, size = value), shape = 0) +
    scale_color_gradientn(colours = hcl.colors(11, palette = 'Blue-Yellow', rev = TRUE)) +
    scale_size_area(max_size = 1) +
    scale_x_datetime(labels = date_format("%Y-%m-%d %H"), date_breaks = "4 hours") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_classic() +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/ONT_DNA/", feature, "_perRead/Synx_error_rates_over_time_by_channel_", feature, ".pdf"))
}

# Calculate error rates without hp elements
# --------------------------------------------------------------------------

for(feature in c(synx_ranges_perf$name, "SynX")) {
  # Load sequencing error rates
  reads_1linebed <- read.table(paste0("project_results/ONT_DNA/", feature, "_perRead/", feature, "_one_line.barcode06.sorted.INTERSECT.perRead.no_HP.tsv"), header = TRUE, sep=";", stringsAsFactors=FALSE, quote="")
  reads_1linebed$start_time <- as.POSIXct(reads_1linebed$start_time, format="%Y-%m-%dT%H:%M:%SZ")
  reads_1linebed$ch <- factor(as.character(reads_1linebed$ch))
  
  # Plot length distribution
  ggplot(reads_1linebed) +
    geom_density(aes(x=Length))
  
  # Calculate per base error frequency
  reads_1linebed$GC_pct <- (reads_1linebed$G + reads_1linebed$C)/reads_1linebed$Length*100
  reads_1linebed$Error_freq <- reads_1linebed$Error/reads_1linebed$Length
  reads_1linebed$SUB_freq <- reads_1linebed$Mismatch/reads_1linebed$Length
  reads_1linebed$INS_freq <- reads_1linebed$Ins/reads_1linebed$Length
  reads_1linebed$DEL_freq <- reads_1linebed$Del/reads_1linebed$Length
  
  # Subset to reads those around the mean (+/- 1 SD)
  if(feature != 'SynX') {
    feat_mean <- mean(reads_1linebed$Length)
    feat_sd <- sd(reads_1linebed$Length)
    reads_1linebed <- reads_1linebed[reads_1linebed$Length < (feat_mean + feat_sd) & reads_1linebed$Length > (feat_mean - feat_sd),]
  } else {
    print('Feature = SynX')
  }
  
  # Order channels by means per base error frequency across all transcripts
  sum_stats <- reads_1linebed %>%
    group_by(ch) %>%
    dplyr::summarize(Mean_error_freq = mean(Error_freq), count(ch)) %>%
    mutate(ch = as.character(ch)) %>%
    arrange(Mean_error_freq)
  reads_1linebed$ch <- factor(reads_1linebed$ch, levels = sum_stats$ch)
  
  # PLot error rate over time
  error_time <- melt(reads_1linebed %>%
                       select(start_time, Error_freq, SUB_freq, INS_freq, DEL_freq), id.vars = 'start_time')
  
  ggplot(error_time, aes(x = start_time, y = value)) +
    geom_line(aes(color = variable), alpha = 1) +
    geom_smooth(alpha = 0.5, color = 'grey50', fill = 'grey90', size = 0.5) +
    scale_colour_npg(palette = c("nrc")) +
    scale_fill_npg(palette = c("nrc")) +
    scale_x_datetime(labels = date_format("%Y-%m-%d %H"), date_breaks = "4 hours") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_classic() +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/ONT_DNA/", feature, "_perRead/Synx_error_rates_over_time_HP_masked_", feature, ".pdf"))
  
  # Plot error rate by channel over time
  error_ch_time <- melt(reads_1linebed %>%
                          select(ch, start_time, Error_freq, SUB_freq, INS_freq, DEL_freq),
                        id.vars = c("ch", "start_time"),
                        measure.vars = c('Error_freq', 'SUB_freq', 'INS_freq', 'DEL_freq'))
  
  ggplot(error_ch_time, aes(x = start_time, y = ch)) +
    geom_point(aes(color = value, size = value), shape = 0) +
    scale_color_gradientn(colours = hcl.colors(11, palette = 'Blue-Yellow', rev = TRUE)) +
    scale_size_area(max_size = 1) +
    scale_x_datetime(labels = date_format("%Y-%m-%d %H"), date_breaks = "4 hours") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme_classic() +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/ONT_DNA/", feature, "_perRead/Synx_error_rates_over_time_by_channel_HP_masked_", feature, ".pdf"))
}

# Plot error rates with time for individual pores over each performance element
# --------------------------------------------------------------------------

for(feature in c(synx_ranges_perf$name)) {
  # Load sequencing error rates
  reads_1linebed <- read.table(paste0("project_results/ONT_DNA/", feature, "_perRead/", feature, "_one_line.barcode06.sorted.INTERSECT.perRead.tsv"), header = TRUE, sep=";", stringsAsFactors=FALSE, quote="")
  reads_1linebed$start_time <- as.POSIXct(reads_1linebed$start_time, format="%Y-%m-%dT%H:%M:%SZ")
  reads_1linebed$start_order <- order(reads_1linebed$start_time)
  reads_1linebed$ch <- factor(as.character(reads_1linebed$ch))
  
  # Plot length distribution
  ggplot(reads_1linebed) +
    geom_density(aes(x=Length))
  
  # Calculate per base error frequency
  reads_1linebed$GC_pct <- (reads_1linebed$G + reads_1linebed$C)/reads_1linebed$Length*100
  reads_1linebed$Error_freq <- reads_1linebed$Error/reads_1linebed$Length

  # Subset to reads those around the mean (+/- 1 SD)
  if(feature != 'SynX') {
    feat_mean <- mean(reads_1linebed$Length)
    feat_sd <- sd(reads_1linebed$Length)
    reads_1linebed <- reads_1linebed[reads_1linebed$Length < (feat_mean + feat_sd) & reads_1linebed$Length > (feat_mean - feat_sd),]
  } else {
    print('Feature = SynX')
  }

  # Calculate pore summary statistics
  sum_stats <- reads_1linebed %>%
    group_by(ch) %>%
    dplyr::summarize(Mean_error_freq = mean(Error_freq), count(ch)) %>%
    mutate(ch = as.character(ch)) %>%
    arrange(Mean_error_freq)
  reads_1linebed$ch <- factor(reads_1linebed$ch, levels = sum_stats$ch)

  # Plot error rate by channel over time
  error_ch_time <- melt(reads_1linebed %>%
                          select(ch, start_time, Error_freq),
                        id.vars = c("ch", "start_time"),
                        measure.vars = c('Error_freq'))
    
  # Subset to pore with >= 20 transcripts
  gt20 <- sum_stats[sum_stats$freq >= 20, 'ch']
  error_ch_time_gt20 <- error_ch_time[is.element(as.character(error_ch_time$ch), pull(gt20, ch)),]
  
  # Plot pore performance over time for pore with >20 transcripts
  dir.create(paste0("project_results/ONT_DNA/", feature, "_perRead/Pore_over_time"))
  for(pore in unique(error_ch_time_gt20$ch)) {
    error_ch_time_pore <- error_ch_time[error_ch_time$ch == pore, ]
    ggplot(error_ch_time_pore, aes(x = start_time, y = value)) +
      geom_point(alpha = 1) +
      geom_smooth(alpha = 0.5, color = 'grey50', fill = 'grey90', size = 0.5) +
      xlab("Time") +
      ylab("Error frequency (per base)") +
      ylim(0, NA) +
      scale_x_datetime(labels = date_format("%Y-%m-%d %H"), date_breaks = "4 hours") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme_classic() +
      ggsave(paste0("project_results/ONT_DNA/", feature, "_perRead/Pore_over_time/Synx_error_rates_over_time_", feature, "_", pore, ".pdf"))
  }
  
}

# Plot error rates for individual pores over each performance element
# --------------------------------------------------------------------------

for(feature in c(synx_ranges_perf$name, "SynX")) {
  # Load sequencing error rates
  reads_1linebed <- read.table(paste0("project_results/ONT_DNA/", feature, "_perRead/", feature, "_one_line.barcode06.sorted.INTERSECT.perRead.tsv"), header = TRUE, sep=";", stringsAsFactors=FALSE, quote="")
  reads_1linebed$start_time <- as.POSIXct(reads_1linebed$start_time, format="%Y-%m-%dT%H:%M:%SZ")
  reads_1linebed$ch <- factor(as.character(reads_1linebed$ch))
  
  # Plot length distribution
  ggplot(reads_1linebed) +
    geom_density(aes(x=Length))
  
  # Calculate per base error frequency
  reads_1linebed$GC_pct <- (reads_1linebed$G + reads_1linebed$C)/reads_1linebed$Length*100
  reads_1linebed$Error_freq <- reads_1linebed$Error/reads_1linebed$Length

  # Subset to reads those around the mean (+/- 1 SD)
  if(feature != 'SynX') {
    feat_mean <- mean(reads_1linebed$Length)
    feat_sd <- sd(reads_1linebed$Length)
    reads_1linebed <- reads_1linebed[reads_1linebed$Length < (feat_mean + feat_sd) & reads_1linebed$Length > (feat_mean - feat_sd),]
  } else {
    print('Feature = SynX')
  }
  
  # Order channels by means per base error frequency across all transcripts
  sum_stats <- reads_1linebed %>%
    group_by(ch) %>%
    dplyr::summarize(Mean_error_freq = mean(Error_freq), count(ch)) %>%
    mutate(ch = as.numeric(as.character(ch))) %>%
    arrange(ch)
  reads_1linebed$ch <- factor(reads_1linebed$ch, levels = sum_stats$ch)
  
  # Annotate with position in flow cell (row = A:G, column = 1:18)
  
  ONT_array <- data.frame("Pore" = 1:126)
  idx <- match(ONT_array$Pore, as.numeric(sum_stats$ch))
  ONT_array$ch <- sum_stats$ch [idx]
  ONT_array$Mean_error_freq <- sum_stats$Mean_error_freq  [idx]
  ONT_array$freq <- sum_stats$freq  [idx]
  ONT_array$Row <- rep(LETTERS[1:7], 18)
  ONT_array$Column <- rep(1:18, 7)
  
  dir.create(paste0("project_results/ONT_DNA/Error_perPore"))
  ggplot(ONT_array, aes(x = Column, y = Row)) +
    geom_point(aes(fill = Mean_error_freq), shape = 21, size = 6) +
    scale_fill_gradientn(colours = hcl.colors(3, palette = 'Zissou 1', rev = FALSE), na.value = "white") +
    theme_classic() +
    scale_x_continuous(breaks=c(1:18)) +
    ggsave(paste0("project_results/ONT_DNA/Error_perPore/Mean_error_rates_per_pore_", feature, ".pdf"))
}

