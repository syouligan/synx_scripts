#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX HP performance features
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/tim_projects/synx/")
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

# List datasets to loop script over
data_sets <- c('Illumina_DNASynX', 'Illumina_SP6', 'Illumina_T7', 'ONT_DNA', 'ONT_DNA_BamHI', 'ONT_DNA_EcoRI', 'ONT_DNA_HindIII', 'ONT_DNA_Fragmentase', 'ONT_cDNA_SP6', 'ONT_cDNA_SP6_2', 'ONT_cDNA_T7', 'ONT_cDNA_T7_2', 'ONT_cDNA_SP6_3', 'ONT_cDNA_T7_3', 'ONT_RNA_SP6')

for (ds in data_sets) {
  # Calculate error rates for each homopolymer in the homopolymer performance elements
  # --------------------------------------------------------------------------
  
  # Load sequencing error rates
  ONTDNA_pileup <- read.table(paste0("project_results/", ds, "/", ds, ".bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  
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
  
  # Calculate error rates across all features in SynX
  # --------------------------------------------------------------------------
  
  synx_ranges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, synx_ranges@ranges))
  synx_ranges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, synx_ranges@ranges))
  synx_ranges$Error_upper <- synx_ranges$Error_mean + synx_ranges$Error_sd
  synx_ranges$Error_lower <- synx_ranges$Error_mean - synx_ranges$Error_sd
  
  synx_ranges$SUB_mean <- mean(extractList(ONTDNA_pileup$SUB_rate, synx_ranges@ranges))
  synx_ranges$SUB_sd <- sd(extractList(ONTDNA_pileup$SUB_rate, synx_ranges@ranges))
  synx_ranges$SUB_upper <- synx_ranges$SUB_mean + synx_ranges$SUB_sd
  synx_ranges$SUB_lower <- synx_ranges$SUB_mean - synx_ranges$SUB_sd
  
  synx_ranges$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_rate, synx_ranges@ranges))
  synx_ranges$INDEL_sd <- sd(extractList(ONTDNA_pileup$INDEL_rate, synx_ranges@ranges))
  synx_ranges$INDEL_upper <- synx_ranges$INDEL_mean + synx_ranges$INDEL_sd
  synx_ranges$INDEL_lower <- synx_ranges$INDEL_mean - synx_ranges$INDEL_sd
  
  synx_ranges$INS_mean <- mean(extractList(ONTDNA_pileup$Ins_rate, synx_ranges@ranges))
  synx_ranges$INS_sd <- sd(extractList(ONTDNA_pileup$Ins_rate, synx_ranges@ranges))
  synx_ranges$INS_upper <- synx_ranges$INS_mean + synx_ranges$INS_sd
  synx_ranges$INS_lower <- synx_ranges$INS_mean - synx_ranges$INS_sd
  
  synx_ranges$DEL_mean <- mean(extractList(ONTDNA_pileup$Del_rate, synx_ranges@ranges))
  synx_ranges$DEL_sd <- sd(extractList(ONTDNA_pileup$Del_rate, synx_ranges@ranges))
  synx_ranges$DEL_upper <- synx_ranges$DEL_mean + synx_ranges$DEL_sd
  synx_ranges$DEL_lower <- synx_ranges$DEL_mean - synx_ranges$DEL_sd
  
  # Make GRanges object for each feature group
  synx_ranges_by_types <- split(synx_ranges, synx_ranges$TYPE)
  promoters <- synx_ranges_by_types$Promoter
  restriction <- synx_ranges_by_types$Restriction_site
  polyA <- synx_ranges_by_types$`Poly-A_tail`
  performance <- synx_ranges_by_types$Sequencing_performance_unit
  performanceGC <- performance[grepl("^P[1-6]", performance$name),]
  performanceHP <- performance[grepl("PH", performance$name),]
  quantitative_ladder <- synx_ranges_by_types$Quantitative_ladder_unit
  size_ladder <- synx_ranges_by_types$Size_ladder_unit
  
  # Plot summary error rates across HP performance features
  # --------------------------------------------------------------------------
  errors <- c('Error', 'SUB', 'INS', 'DEL')
  
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
    ggsave(paste0("project_results/HP_performance/Homopolymer_synx_total_frequency.pdf"))
  
  # Subset to just HP performance elements
  performanceHP1_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP[performanceHP$name == 'PH1',])
  performanceHP2_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP[performanceHP$name == 'PH2',])
  performanceHPtot_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP)
  
  # Plot homopolymer frequency over performance element
  ggplot(data.frame(table(performanceHP1_ranges$score)), aes(x=Var1, y=Freq)) +
    geom_line(group = 1) +
    xlab("Homopolymer (bp)") +
    ylab("Frequency") +
    theme_classic() +
    ggsave(paste0("project_results/HP_performance/HPperformance_element1_homopolymer_frequency_", ds, ".pdf"))
  
  ggplot(data.frame(table(performanceHP2_ranges$score)), aes(x=Var1, y=Freq)) +
    geom_line(group = 1) +
    xlab("Homopolymer (bp)") +
    ylab("Frequency") +
    theme_classic() +
    ggsave(paste0("project_results/HP_performance/HPperformance_element2_homopolymer_frequency_", ds, ".pdf"))
  
  ggplot(data.frame(table(performanceHPtot_ranges$score)), aes(x=Var1, y=Freq)) +
    geom_line(group = 1) +
    xlab("Homopolymer (bp)") +
    ylab("Frequency") +
    theme_classic() +
    ggsave(paste0("project_results/HP_performance/HPperformance_total_homopolymer_frequency_", ds, ".pdf"))
  
  # Plot error rate by homopolymer length
  for( i in paste0(errors, '_mean')) {
    ggplot(data.frame(performanceHPtot_ranges@elementMetadata), ) +
      geom_boxplot(aes_string(x='score', group='score', y = i), color = "grey80") +
      geom_point(aes_string(x='score', y = i, colour = 'name')) +
      scale_color_npg(palette = c("nrc"), alpha = 1) +
      xlab("Homopolymer (bp)") +
      ylab(i) +
      theme_classic() +
      ggsave(paste0("project_results/HP_performance/HPperformance_homopolymer_error_rate_",  i,"_", ds, ".pdf"))
  }
  
  # Error rates (all) across HP-performance units
  HP_fraction_long <- melt(data.frame(performanceHPtot_ranges@elementMetadata[,c(paste0(errors, "_mean"), 'score')]), id = 'score')
  
  ggplot(HP_fraction_long) +
    geom_boxplot(aes(x=score, y = value, colour = variable, group = score)) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    xlab("HP content") +
    ylab('Error rate') +
    theme_classic() +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/HP_performance/HPperformance_all_errors_", ds, ".pdf"))
}
