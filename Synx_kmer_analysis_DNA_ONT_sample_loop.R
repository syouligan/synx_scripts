#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX kmers ONT-DNA
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
digested_data_sets <- c('ONT_DNA_BamHI', 'ONT_DNA_EcoRI', 'ONT_DNA_HindIII', 'ONT_DNA_Fragmentase')
digested_data_sets <- c('ONT_cDNA_SP6', 'ONT_cDNA_SP6_2', 'ONT_cDNA_T7', 'ONT_cDNA_T7_2', 'ONT_RNA_SP6')
digested_data_sets <- c('ONT_cDNA_SP6_3', 'ONT_cDNA_T7_3')

for (ds in digested_data_sets) {
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
  synx_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup$REF_NT))
  
  KD_k <- function(x){paste(x, collapse = "")}
  synx_seq_F_6mers <- runner(x = synx_seq_F, k = 6, f = KD_k)
  synx_seq_F_6mers <- synx_seq_F_6mers[str_length(synx_seq_F_6mers) == 6]
  synx_seq_RC_6mers <- runner(x = synx_seq_RC, k = 6, f = KD_k)
  synx_seq_RC_6mers <- synx_seq_RC_6mers[str_length(synx_seq_RC_6mers) == 6]
  synx_seq_total_6mers <- unique(synx_seq_F_6mers, synx_seq_RC_6mers)
  
  # Make ranges of 10mers in synx sequence
  k_F_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_F_6mers), width = 6), strand = "+", SEQUENCE = synx_seq_F_6mers)
  k_RC_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_RC_6mers), width = 6), strand = "+", SEQUENCE = rev(synx_seq_RC_6mers))
  k_total_granges <- c(k_F_granges, k_RC_granges)
  
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
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_frequency.pdf"))
  
  # Plot Error vs GC content vs homopolymer length
  ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(Error_mean)) +
    geom_jitter(aes(x = GC_pct, y = HP_length, fill = Error_mean), alpha = 0.3, shape = 21) +
    theme_classic() +
    scale_fill_viridis_c(option = 'B') +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    scale_y_continuous(breaks = c(1:6)) +
    xlab("GC percent") +
    ylab("HP length") +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_GC_vs_HP_Error.pdf"))
  
  ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(INS_mean)) +
    geom_jitter(aes(x = GC_pct, y = HP_length, fill = INS_mean), alpha = 0.3, shape = 21) +
    theme_classic() +
    scale_fill_viridis_c(option = 'B') +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    scale_y_continuous(breaks = c(1:6)) +
    xlab("GC percent") +
    ylab("HP length") +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_GC_vs_HP_INS.pdf"))
  
  ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(SUB_mean)) +
    geom_jitter(aes(x = GC_pct, y = HP_length, fill = SUB_mean), alpha = 0.3, shape = 21) +
    theme_classic() +
    scale_fill_viridis_c(option = 'B') +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    scale_y_continuous(breaks = c(1:6)) +
    xlab("GC percent") +
    ylab("HP length") +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_GC_vs_HP_SUB.pdf"))
  
  ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(DEL_mean)) +
    geom_jitter(aes(x = GC_pct, y = HP_length, fill = DEL_mean), alpha = 0.3, shape = 21) +
    theme_classic() +
    scale_fill_viridis_c(option = 'B') +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    scale_y_continuous(breaks = c(1:6)) +
    xlab("GC percent") +
    ylab("HP length") +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_GC_vs_HP_DEL.pdf"))
  
  # Calculate error rates across kmers
  # --------------------------------------------------------------------------
  
  # Calculate error rate per 6mer
  kmer_error_freq <- data.frame(k_total_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Error_freq_mean = mean(Error_freq, na.rm = TRUE),
                     SUB_freq_mean = mean(SUB_freq, na.rm = TRUE),
                     INS_freq_mean = mean(INS_freq, na.rm = TRUE),
                     DEL_freq_mean = mean(DEL_freq, na.rm = TRUE),
                     GC_pct = mean(GC_pct, na.rm = TRUE),
                     HP_length = mean(HP_length, na.rm = TRUE))
  
  # Plot error frequency by GC content
  melt(id.vars = 'GC_pct', data = kmer_error_freq %>% 
         select(GC_pct, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = log2(value), x = GC_pct, group = GC_pct, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    xlab("Kmer GC content") +
    ylab("Error Frequency (log2)") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_error_frequency_w_GC_content.pdf"))
  
  # Plot error frequency by homopolymer length
  melt(id.vars = 'HP_length', data = kmer_error_freq %>% 
         select(HP_length, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = log2(value), x = HP_length, group = HP_length, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = 1:6) +
    xlab("Kmer Homopolymer Length") +
    ylab("Error Frequency (log2)") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_error_frequency_w_HP_length.pdf"))
  
  # Plot error frequency by kmer
  kmer_error_freq %>%
    arrange(desc(Error_freq_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = log2(Error_freq_mean)))+
    theme_classic() +
    xlab("Kmer") +
    ylab("Error Frequency (log2)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_Error_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(INS_freq_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = log2(INS_freq_mean))) +
    theme_classic() +
    xlab("Kmer") +
    ylab("INS Frequency (log2)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_INS_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(DEL_freq_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = log2(DEL_freq_mean))) +
    theme_classic() +
    xlab("Kmer") +
    ylab("DEL Frequency (log2)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_DEL_frequency.pdf"))
  
  kmer_error_freq %>%
    arrange(desc(SUB_freq_mean)) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(x = SEQUENCE, y = log2(SUB_freq_mean))) +
    theme_classic() +
    xlab("Kmer") +
    ylab("SUB Frequency (log2)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    ggsave(paste0("project_results/", ds, "/SynX_kmer_", ds,"_SUB_frequency.pdf"))
}
