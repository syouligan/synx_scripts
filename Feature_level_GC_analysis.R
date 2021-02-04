#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX GC performance features
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
  # Calculate error rates for each base (SUBS and INDELS)
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
  
  # Plot summary error rates across GC performance features
  # --------------------------------------------------------------------------
  errors <- c('Error', 'SUB', 'INS', 'DEL')
  
  # Plot sequencing errors across performance units
  performanceGC <- performanceGC[order(performanceGC$name), ]
  
  # Plot GC content of GC performance units
  ggplot(data.frame(performanceGC@elementMetadata)) +
    geom_bar(aes(y = GC_pct, x = name), stat = 'identity') +
    theme_minimal() +
    ggsave(paste0("project_results/GC_performance/SynX_performance_GC_content_", ds, ".pdf"))
  
  # # Error rates (individual) across GC-performance units
  # for(i in errors) {
  #   Mean <- paste0(i, "_mean")
  #   Lower <- paste0(i, "_lower")
  #   Upper <- paste0(i, "_upper")
  #   
  #   ggplot(data.frame(performanceGC@elementMetadata)) +
  #     geom_line(aes_string(x = 'GC_pct', y = Mean)) +
  #     geom_ribbon(aes_string(x = 'GC_pct', y = Mean, ymin=Lower, ymax=Upper), fill="blue", alpha=0.2) +
  #     xlab("GC content") +
  #     ylab(i) +
  #     theme_minimal() +
  #     ggsave(paste0("project_results/", ds, "/SynX_performance_GC_", i, "_", ds, ".pdf"))
  # }

  # Error rates across GC-performance units (mean only)
  GC_fraction_long_mean <- melt(data.frame(performanceGC@elementMetadata[,c(paste0(errors, "_mean"), 'name')]),  id.vars = 'name')
  GC_fraction_long_mean$Measure <- gsub("_.*","", GC_fraction_long_mean$variable)
  GC_fraction_long_mean$Measure <- factor(GC_fraction_long_mean$Measure, levels = errors)
  
  GC_fraction_long_mean$Mean <- GC_fraction_long_mean$value
  GC_fraction_long_sd <- melt(data.frame(performanceGC@elementMetadata[,c(paste0(errors, "_sd"), 'name')]),  id.vars = 'name')
  GC_fraction_long_sd$Measure <- gsub("_.*","", GC_fraction_long_sd$variable)
  GC_fraction_long_mean$SD <- GC_fraction_long_sd$value 
  
  ggplot(GC_fraction_long_mean, aes(x=name, y=Mean, fill=Measure)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    # geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.2, position=position_dodge(.9)) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    xlab("GC content") +
    ylab('Error rate') +
    theme_classic() +
    facet_wrap(~Measure) +
    ggsave(paste0("project_results/GC_performance/SynX_performance_GC_mean_error_per_feature_", ds, ".pdf"))

  # Plot basewise error rates across GC performance features
  # --------------------------------------------------------------------------
  
  # Make Granges to plot single base resolution data tracks
  pileup_ranges <- GRanges(seqnames = "SynX",
                           ranges = IRanges(start = ONTDNA_pileup$Number, end = ONTDNA_pileup$Number, width = 1),
                           strand = "*")
  
  for(col in colnames(ONTDNA_pileup)){
    pileup_ranges@elementMetadata[,col] <- list(ONTDNA_pileup[,col])
  }
  
  # Subset to bases in GC performance features
  for(pe in performanceGC$name){
    tmp_range <- performanceGC[performanceGC$name == pe,]
    st <- tmp_range@ranges@start
    ed <- tmp_range@ranges@start + tmp_range@ranges@width
    pe_range <- pileup_ranges[pileup_ranges$Number >= st & pileup_ranges$Number < ed,]
    pe_range$name <- pe
    pe_range$GC_pct <- tmp_range$GC_pct
    
    if(exists('GCpe_ranges')){
      GCpe_ranges <- c(GCpe_ranges, pe_range)
    } else {
      GCpe_ranges <- pe_range
    }
  }
  
  # Make box plot of basewise error rate in each GC performance feature
  melt(id.vars = 'name', data = data.frame(GCpe_ranges@elementMetadata) %>% 
         select(name, Error_rate, SUB_rate, Ins_rate, Del_rate)) %>%
    ggplot() +
    geom_boxplot(aes(y = value, x = name, group = name, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_discrete(breaks = 1:6) +
    xlab("Kmer Homopolymer Length") +
    ylab("Error rate") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/GC_performance/SynX_performance_GC_basewise_error_per_feature_", ds, ".pdf"))
}
