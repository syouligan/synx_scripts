#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX kmers
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
library('tidyr')
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

# List datasets to loop script over
data_sets <- c('Illumina_DNASynX', 'Illumina_SP6', 'Illumina_T7', 'ONT_DNA', 'ONT_DNA_Fragmentase', 'ONT_cDNA_SP6', 'ONT_cDNA_T7')

# Store kmer analysis for each dataset in lists
dslist.names <- data_sets
dslist1 <- vector("list", length(dslist.names))
names(dslist1) <- dslist.names

dslist2 <- vector("list", length(dslist.names))
names(dslist2) <- dslist.names

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
  
  # Calculate GC content indepenedent and Homopolymer content of kmers
  # --------------------------------------------------------------------------
  
  # Calculate all unique 6mers in Synx
  synx_seq_F <- ONTDNA_pileup$REF_NT
  
  KD_k <- function(x){paste(x, collapse = "")}
  synx_seq_F_6mers <- runner(x = synx_seq_F, k = 6, f = KD_k)
  synx_seq_total_6mers <- synx_seq_F_6mers[str_length(synx_seq_F_6mers) == 6]
  
  # Make ranges of 6mers in synx sequence
  k_total_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_total_6mers), width = 6), strand = "+", SEQUENCE = synx_seq_total_6mers)
  
  # Calculate GC content
  GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G')) / 6 *100}
  k_total_granges$GC_pct <- GC_k(k_total_granges$SEQUENCE)
  
  # Calculate homopolymers in kmers
  hp_max <- data.frame(k_total_granges$SEQUENCE)
  hp_max$HP2 <- grepl("AA|TT|GG|CC", k_total_granges$SEQUENCE)
  hp_max$HP3 <- grepl("AAA|TTT|GGG|CCC", k_total_granges$SEQUENCE)
  hp_max$HP4 <- grepl("AAAA|TTTT|GGGG|CCCC", k_total_granges$SEQUENCE)
  hp_max$HP5 <- grepl("AAAAA|TTTTT|GGGGG|CCCCC", k_total_granges$SEQUENCE)
  hp_max$HP6 <- grepl("AAAAAA|TTTTTT|GGGGGG|CCCCCC", k_total_granges$SEQUENCE)
  k_total_granges$HP_length <- rowSums(hp_max[,2:6]) + 1
  
  k_total_granges$Depth_mean <- mean(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges))
  
  # Calculate error rates
  k_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, k_total_granges@ranges))
  k_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
  k_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
  k_total_granges$Error_upper <- k_total_granges$Error_mean + k_total_granges$Error_sd
  k_total_granges$Error_lower <- k_total_granges$Error_mean - k_total_granges$Error_sd
  
  k_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup$SUB_counts, k_total_granges@ranges))
  k_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup$SUB_rate, k_total_granges@ranges))
  k_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup$SUB_rate, k_total_granges@ranges))
  k_total_granges$SUB_upper <- k_total_granges$SUB_mean + k_total_granges$SUB_sd
  k_total_granges$SUB_lower <- k_total_granges$SUB_mean - k_total_granges$SUB_sd
  
  k_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup$INDEL_counts, k_total_granges@ranges))
  k_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_rate, k_total_granges@ranges))
  k_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup$INDEL_rate, k_total_granges@ranges))
  k_total_granges$INDEL_upper <- k_total_granges$INDEL_mean + k_total_granges$INDEL_sd
  k_total_granges$INDEL_lower <- k_total_granges$INDEL_mean - k_total_granges$INDEL_sd
  
  k_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup$Ins, k_total_granges@ranges))
  k_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup$Ins_rate, k_total_granges@ranges))
  k_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup$Ins_rate, k_total_granges@ranges))
  k_total_granges$INS_upper <- k_total_granges$INS_mean + k_total_granges$INS_sd
  k_total_granges$INS_lower <- k_total_granges$INS_mean - k_total_granges$INS_sd
  
  k_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup$Del, k_total_granges@ranges))
  k_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup$Del_rate, k_total_granges@ranges))
  k_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup$Del_rate, k_total_granges@ranges))
  k_total_granges$DEL_upper <- k_total_granges$DEL_mean + k_total_granges$DEL_sd
  k_total_granges$DEL_lower <- k_total_granges$DEL_mean - k_total_granges$DEL_sd
  
  # Add to ds1 lists
  dslist1[[ds]] <- k_total_granges
  
  
  # Calculate kmer frequency
  kmer_frequency <- data.frame(table(k_total_granges@elementMetadata$SEQUENCE))
  kmer_frequency <- kmer_frequency[order(kmer_frequency$Freq, decreasing = TRUE),]
  kmer_frequency$Var1 <- factor(kmer_frequency$Var1, levels = kmer_frequency$Var1)
  ggplot(kmer_frequency) +
    geom_point(aes(x = Var1, y = log2(Freq)))+
    theme_classic() +
    xlab("Kmer") +
    ylab("Kmer Frequency (log2)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_frequency.pdf"))
  
  # Calculate error rates across kmers
  # --------------------------------------------------------------------------
  
  # Calculate error rate per 6mer
  kmer_error_freq <- data.frame(k_total_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Kmer_count = n(),
                     Error_rate_mean = mean(Error_mean, na.rm = TRUE),
                     SUB_rate_mean = mean(SUB_mean, na.rm = TRUE),
                     INS_rate_mean = mean(INS_mean, na.rm = TRUE),
                     DEL_rate_mean = mean(DEL_mean, na.rm = TRUE),
                     GC_pct = mean(GC_pct, na.rm = TRUE),
                     HP_length = mean(HP_length, na.rm = TRUE),
                     Depth = max(Depth_mean, na.rm = TRUE))
  
  # Add to ds2 lists
  dslist2[[ds]] <- kmer_error_freq
  
  # Plot error frequency by GC content
  melt(id.vars = 'GC_pct', data = kmer_error_freq %>% 
         select(GC_pct, Error_rate_mean, SUB_rate_mean, INS_rate_mean, DEL_rate_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = value, x = GC_pct, group = GC_pct, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    xlab("Kmer GC content") +
    ylab("Error rate") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_error_frequency_w_GC_content.pdf"))
  
  # Plot error frequency by homopolymer length
  melt(id.vars = 'HP_length', data = kmer_error_freq %>% 
         select(HP_length, Error_rate_mean, SUB_rate_mean, INS_rate_mean, DEL_rate_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = value, x = HP_length, group = HP_length, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = 1:6) +
    xlab("Kmer Homopolymer Length") +
    ylab("Error rate") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_error_frequency_w_HP_length.pdf"))
  
  # Plot error frequency by kmer
  kmer_error_freq %>%
    arrange(desc(Error_rate_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = Error_rate_mean)) +
    theme_classic() +
    xlab("Kmer") +
    ylab("Error rate") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_Error_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(INS_rate_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = INS_rate_mean)) +
    theme_classic() +
    xlab("Kmer") +
    ylab("INS rate") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_INS_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(DEL_rate_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = DEL_rate_mean)) +
    theme_classic() +
    xlab("Kmer") +
    ylab("DEL rate") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_DEL_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(SUB_rate_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = SUB_rate_mean)) +
    theme_classic() +
    xlab("Kmer") +
    ylab("SUB rate") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_SUB_frequency.pdf"))
}

# Plot error rate of individual kmers in each dataset
# --------------------------------------------------------------------------

for (ds in data_sets) {
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
  
  # Make Granges to plot single base resolution data tracks
  pileup_ranges <- GRanges(seqnames = "SynX",
                           ranges = IRanges(start = ONTDNA_pileup$Number, end = ONTDNA_pileup$Number, width = 1),
                           strand = "*")
  
  for(col in colnames(ONTDNA_pileup)){
    pileup_ranges@elementMetadata[,col] <- list(ONTDNA_pileup[,col])
  }
  
  kmer_error_freq <- dslist2[[ds]]
  k_total_granges <- dslist1[[ds]]
  
  # List kmers represented at least three times in the SynX genome
  kmerGT2 <- kmer_error_freq[kmer_error_freq$Kmer_count > 2, ]$SEQUENCE
  
  # Plot kmers with mean and SD
  for (kmer in c(kmerGT2)) {
    dir.create(paste0('project_results/SynX_kmers/kmer_plots/', kmer))
    KOI <- k_total_granges[k_total_granges$SEQUENCE == kmer, ]
    
    ggplot(data.frame(KOI@elementMetadata), aes(x=SEQUENCE, y=Error_mean, color=SEQUENCE)) +
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
      stat_summary(fun=mean, geom="crossbar", color="black", width=0.4) +
      geom_point(size = 3, alpha = 0.5) +
      scale_color_npg(palette = c("nrc"), alpha = 1) +
      scale_fill_npg(palette = c("nrc"), alpha = 1) +
      ylim(0, max(k_total_granges$Error_mean)) +
      xlab("Kmer sequence") +
      ylab('Error rate') +
      theme_classic() +
      ggsave(paste0('project_results/SynX_kmers/kmer_plots/', kmer, '/Error_rates', '_', kmer, '_', ds ,'.pdf'))
    
    # Plot error rate of individual bases in each kmer
    # --------------------------------------------------------------------------
    
    tmp_range <- k_total_granges[k_total_granges$SEQUENCE == kmer,]
    kmer_tot_df <- data.frame()
    for(i in 1:length(tmp_range@ranges)){
      tmp <- tmp_range@ranges[i]
      kmer_df <- data.frame('Error_rate' = as.numeric(ONTDNA_pileup$Error_rate[tmp@start:(tmp@start + tmp@width - 1)]))
      kmer_df$REF_NT <- c(ONTDNA_pileup$REF_NT[tmp@start:(tmp@start + tmp@width - 1)])
      kmer_df$Kmer_NT <- factor(paste(1:6, ONTDNA_pileup$REF_NT[tmp@start:(tmp@start + tmp@width - 1)], sep = '_'))
      kmer_df$Occurence <- i
      kmer_tot_df <- rbind(kmer_tot_df, kmer_df)
      }

    melt(id.vars = c('Kmer_NT', 'Occurence', 'REF_NT'), data = kmer_tot_df %>% 
           select(Kmer_NT, Error_rate, Occurence, REF_NT)) %>%
      ggplot(aes(x=Kmer_NT, y=value, color=Kmer_NT)) +
      geom_point(size = 3, alpha = 0.5, position = position_dodge2(w = 0.4)) +
      stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="black", width=0.2) +
      stat_summary(fun=mean, geom="crossbar", color="black", width=0.4) +
      scale_color_npg(palette = c("nrc"), alpha = 1) +
      scale_fill_npg(palette = c("nrc"), alpha = 1) +
      scale_x_discrete(expand=c(0.1, 0.5)) +
      ylim(0, max(pileup_ranges$Error_rate)) +
      xlab("Kmer sequence") +
      ylab("Error rate") +
      theme_classic() +
      ggsave(paste0('project_results/SynX_kmers/kmer_plots/', kmer, '/Base_wise_error_rates_', '_', kmer, '_', ds ,'.pdf'))
    
    # Make kmer motif plot
    data.frame(kmer_tot_df) %>%
      group_by(Kmer_NT) %>%
      select(Kmer_NT, Error_rate, Occurence , REF_NT) %>%
      dplyr::summarize(Error_rate_mean = mean(Error_rate, na.rm = TRUE),
                       REF_NT = unique(REF_NT)) %>%
      ggplot(aes(x=Kmer_NT, color = REF_NT, size = Error_rate_mean, y = 1)) +
      geom_text(aes(label = REF_NT), fontface = "bold") +
      scale_color_npg(palette = c("nrc"), alpha = 1, guide = NULL) +
      scale_size(range = c(3, 15), limits = c(0, 1)) +
      xlab("Kmer sequence") +
      theme_void() +
      theme(panel.border = element_rect(colour = "black", fill=NA)) +
      ggsave(paste0('project_results/SynX_kmers/kmer_plots/', kmer, '/Motif_error_rates_', '_', kmer, '_', ds ,'.pdf'))
  }
}
