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
data_sets <- c('ONT_DNA_phix_synx', 'ONT_DNA_phix_synx_2')

for (ds in data_sets) {
  # Calculate error rates for each base (SUBS and INDELS) for each vector individually
  # --------------------------------------------------------------------------
  
  # Load sequencing error rates
  ONTDNA_pileup_total <- read.table(paste0("project_results/", ds, "/", ds, ".bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  
  # Find reference mapping rates
  get_correct <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    correct <- bases[bases == ref]
    sum(as.numeric(df[correct]))
  }
  
  ONTDNA_pileup_total$Depth <- ONTDNA_pileup_total$Ins + ONTDNA_pileup_total$Del + ONTDNA_pileup_total$Coverage
  ONTDNA_pileup_total$REF_counts <- apply(ONTDNA_pileup_total, 1, get_correct)
  ONTDNA_pileup_total$REF_rate <- ONTDNA_pileup_total$REF_counts / ONTDNA_pileup_total$Depth
  
  # Find substitution rates
  get_subs <- function(df){
    bases <- c('A', 'T', 'G', 'C')
    ref <- df['REF_NT']
    subs <- bases[!(bases == ref)]
    sum(as.numeric(df[subs]))
  }
  
  ONTDNA_pileup_total$SUB_counts <- apply(ONTDNA_pileup_total, 1, get_subs)
  ONTDNA_pileup_total$SUB_rate <- ONTDNA_pileup_total$SUB_counts / ONTDNA_pileup_total$Depth
  
  # Find INDEL rates
  ONTDNA_pileup_total$Ins_rate <- ONTDNA_pileup_total$Ins / ONTDNA_pileup_total$Depth
  ONTDNA_pileup_total$Del_rate <- ONTDNA_pileup_total$Del / ONTDNA_pileup_total$Depth
  ONTDNA_pileup_total$INDEL_counts <- ONTDNA_pileup_total$Ins + ONTDNA_pileup_total$Del
  ONTDNA_pileup_total$INDEL_rate <- ONTDNA_pileup_total$INDEL_counts / ONTDNA_pileup_total$Depth
  
  # Find error rates
  ONTDNA_pileup_total$Error_counts <- ONTDNA_pileup_total$SUB_counts + ONTDNA_pileup_total$INDEL_counts
  ONTDNA_pileup_total$Error_rate <- ONTDNA_pileup_total$Error_counts / ONTDNA_pileup_total$Depth
  
  # Calculate GC content indepenedent and Homopolymer content of 6mers
  # --------------------------------------------------------------------------
  
  for(vec in unique(ONTDNA_pileup_total$Chrom)) {
    
    # Subset to each genome
    ONTDNA_pileup <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == vec, ]
    
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
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_frequency.pdf"))
    
    # Plot Error vs GC content vs homopolymer length
    ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(Error_mean)) +
      geom_jitter(aes(x = GC_pct, y = HP_length, fill = Error_mean), alpha = 0.3, shape = 21) +
      theme_classic() +
      scale_fill_viridis_c(option = 'B') +
      scale_x_continuous(breaks = c(0, 50, 100)) +
      scale_y_continuous(breaks = c(1:6)) +
      xlab("GC percent") +
      ylab("HP length") +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_GC_vs_HP_Error.pdf"))
    
    ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(INS_mean)) +
      geom_jitter(aes(x = GC_pct, y = HP_length, fill = INS_mean), alpha = 0.3, shape = 21) +
      theme_classic() +
      scale_fill_viridis_c(option = 'B') +
      scale_x_continuous(breaks = c(0, 50, 100)) +
      scale_y_continuous(breaks = c(1:6)) +
      xlab("GC percent") +
      ylab("HP length") +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_GC_vs_HP_INS.pdf"))
    
    ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(SUB_mean)) +
      geom_jitter(aes(x = GC_pct, y = HP_length, fill = SUB_mean), alpha = 0.3, shape = 21) +
      theme_classic() +
      scale_fill_viridis_c(option = 'B') +
      scale_x_continuous(breaks = c(0, 50, 100)) +
      scale_y_continuous(breaks = c(1:6)) +
      xlab("GC percent") +
      ylab("HP length") +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_GC_vs_HP_SUB.pdf"))
    
    ggplot(data.frame(k_total_granges@elementMetadata) %>% arrange(DEL_mean)) +
      geom_jitter(aes(x = GC_pct, y = HP_length, fill = DEL_mean), alpha = 0.3, shape = 21) +
      theme_classic() +
      scale_fill_viridis_c(option = 'B') +
      scale_x_continuous(breaks = c(0, 50, 100)) +
      scale_y_continuous(breaks = c(1:6)) +
      xlab("GC percent") +
      ylab("HP length") +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_GC_vs_HP_DEL.pdf"))
    
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
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_error_frequency_w_GC_content.pdf"))
    
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
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_error_frequency_w_HP_length.pdf"))
    
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
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_Error_frequency.pdf"))
    
    kmer_error_freq %>%
      arrange(desc(INS_freq_mean)) %>%
      mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
      ggplot( ) +
      geom_point(aes(x = SEQUENCE, y = log2(INS_freq_mean))) +
      theme_classic() +
      xlab("Kmer") +
      ylab("INS Frequency (log2)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_INS_frequency.pdf"))
    
    kmer_error_freq %>%
      arrange(desc(DEL_freq_mean)) %>%
      mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
      ggplot( ) +
      geom_point(aes(x = SEQUENCE, y = log2(DEL_freq_mean))) +
      theme_classic() +
      xlab("Kmer") +
      ylab("DEL Frequency (log2)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_DEL_frequency.pdf"))
    
    kmer_error_freq %>%
      arrange(desc(SUB_freq_mean)) %>%
      mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
      ggplot( ) +
      geom_point(aes(x = SEQUENCE, y = log2(SUB_freq_mean))) +
      theme_classic() +
      xlab("Kmer") +
      ylab("SUB Frequency (log2)") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      ggsave(paste0("project_results/", ds, "/", vec,"_kmer_", ds,"_SUB_frequency.pdf")) 
  }
  # Compare GC content of 100mers between SynX and PhiX
  # --------------------------------------------------------------------------
  
  ### Calculate all unique 100mers in SynX
  
  ONTDNA_pileup_SynX <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == "SynX", ]
  synx_seq_F <- ONTDNA_pileup_SynX$REF_NT
  synx_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup_SynX$REF_NT))
  
  KD_k <- function(x){paste(x, collapse = "")}
  synx_seq_F_100mers <- runner(x = synx_seq_F, k = 100, f = KD_k)
  synx_seq_F_100mers <- synx_seq_F_100mers[str_length(synx_seq_F_100mers) == 100]
  synx_seq_RC_100mers <- runner(x = synx_seq_RC, k = 100, f = KD_k)
  synx_seq_RC_100mers <- synx_seq_RC_100mers[str_length(synx_seq_RC_100mers) == 100]
  synx_seq_total_100mers <- unique(synx_seq_F_100mers, synx_seq_RC_100mers)
  
  # Make ranges of 100mers in synx sequence
  synx100_F_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_F_100mers), width = 100), strand = "+", SEQUENCE = synx_seq_F_100mers)
  synx100_RC_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_RC_100mers), width = 100), strand = "+", SEQUENCE = rev(synx_seq_RC_100mers))
  synx100_total_granges <- c(synx100_F_granges, synx100_RC_granges)
  
  # Calculate GC content in SynX 100mers
  GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G'))}
  synx100_total_granges$GC_pct <- GC_k(synx100_total_granges$SEQUENCE)
  
  # Calculate error rates in SynX 100mers
  synx100_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup_SynX$Error_counts, synx100_total_granges@ranges))
  synx100_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup_SynX$Error_rate, synx100_total_granges@ranges))
  synx100_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup_SynX$Error_rate, synx100_total_granges@ranges))
  synx100_total_granges$Error_upper <- synx100_total_granges$Error_mean + synx100_total_granges$Error_sd
  synx100_total_granges$Error_lower <- synx100_total_granges$Error_mean - synx100_total_granges$Error_sd
  
  synx100_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup_SynX$SUB_counts, synx100_total_granges@ranges))
  synx100_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup_SynX$SUB_rate, synx100_total_granges@ranges))
  synx100_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup_SynX$SUB_rate, synx100_total_granges@ranges))
  synx100_total_granges$SUB_upper <- synx100_total_granges$SUB_mean + synx100_total_granges$SUB_sd
  synx100_total_granges$SUB_lower <- synx100_total_granges$SUB_mean - synx100_total_granges$SUB_sd
  
  synx100_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup_SynX$INDEL_counts, synx100_total_granges@ranges))
  synx100_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup_SynX$INDEL_rate, synx100_total_granges@ranges))
  synx100_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup_SynX$INDEL_rate, synx100_total_granges@ranges))
  synx100_total_granges$INDEL_upper <- synx100_total_granges$INDEL_mean + synx100_total_granges$INDEL_sd
  synx100_total_granges$INDEL_lower <- synx100_total_granges$INDEL_mean - synx100_total_granges$INDEL_sd
  
  synx100_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup_SynX$Ins, synx100_total_granges@ranges))
  synx100_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup_SynX$Ins_rate, synx100_total_granges@ranges))
  synx100_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup_SynX$Ins_rate, synx100_total_granges@ranges))
  synx100_total_granges$INS_upper <- synx100_total_granges$INS_mean + synx100_total_granges$INS_sd
  synx100_total_granges$INS_lower <- synx100_total_granges$INS_mean - synx100_total_granges$INS_sd
  
  synx100_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup_SynX$Del, synx100_total_granges@ranges))
  synx100_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup_SynX$Del_rate, synx100_total_granges@ranges))
  synx100_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup_SynX$Del_rate, synx100_total_granges@ranges))
  synx100_total_granges$DEL_upper <- synx100_total_granges$DEL_mean + synx100_total_granges$DEL_sd
  synx100_total_granges$DEL_lower <- synx100_total_granges$DEL_mean - synx100_total_granges$DEL_sd
  
  ### Calculate all unique 100mers in PhiX
  
  ONTDNA_pileup_PhiX <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == "PhiX", ]
  phix_seq_F <- ONTDNA_pileup_PhiX$REF_NT
  phix_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup_PhiX$REF_NT))
  
  KD_k <- function(x){paste(x, collapse = "")}
  phix_seq_F_100mers <- runner(x = phix_seq_F, k = 100, f = KD_k)
  phix_seq_F_100mers <- phix_seq_F_100mers[str_length(phix_seq_F_100mers) == 100]
  phix_seq_RC_100mers <- runner(x = phix_seq_RC, k = 100, f = KD_k)
  phix_seq_RC_100mers <- phix_seq_RC_100mers[str_length(phix_seq_RC_100mers) == 100]
  phix_seq_total_100mers <- unique(phix_seq_F_100mers, phix_seq_RC_100mers)
  
  # Make ranges of 100mers in phix sequence
  phix100_F_granges <- GRanges(seqnames = 'PhiX', ranges = IRanges(start = 1:length(phix_seq_F_100mers), width = 100), strand = "+", SEQUENCE = phix_seq_F_100mers)
  phix100_RC_granges <- GRanges(seqnames = 'PhiX', ranges = IRanges(start = 1:length(phix_seq_RC_100mers), width = 100), strand = "+", SEQUENCE = rev(phix_seq_RC_100mers))
  phix100_total_granges <- c(phix100_F_granges, phix100_RC_granges)
  
  # Calculate GC content
  GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G'))}
  phix100_total_granges$GC_pct <- GC_k(phix100_total_granges$SEQUENCE)
  
  # Calculate error rates in PhiX 100mers
  phix100_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup_PhiX$Error_counts, phix100_total_granges@ranges))
  phix100_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup_PhiX$Error_rate, phix100_total_granges@ranges))
  phix100_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup_PhiX$Error_rate, phix100_total_granges@ranges))
  phix100_total_granges$Error_upper <- phix100_total_granges$Error_mean + phix100_total_granges$Error_sd
  phix100_total_granges$Error_lower <- phix100_total_granges$Error_mean - phix100_total_granges$Error_sd
  
  phix100_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup_PhiX$SUB_counts, phix100_total_granges@ranges))
  phix100_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup_PhiX$SUB_rate, phix100_total_granges@ranges))
  phix100_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup_PhiX$SUB_rate, phix100_total_granges@ranges))
  phix100_total_granges$SUB_upper <- phix100_total_granges$SUB_mean + phix100_total_granges$SUB_sd
  phix100_total_granges$SUB_lower <- phix100_total_granges$SUB_mean - phix100_total_granges$SUB_sd
  
  phix100_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup_PhiX$INDEL_counts, phix100_total_granges@ranges))
  phix100_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup_PhiX$INDEL_rate, phix100_total_granges@ranges))
  phix100_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup_PhiX$INDEL_rate, phix100_total_granges@ranges))
  phix100_total_granges$INDEL_upper <- phix100_total_granges$INDEL_mean + phix100_total_granges$INDEL_sd
  phix100_total_granges$INDEL_lower <- phix100_total_granges$INDEL_mean - phix100_total_granges$INDEL_sd
  
  phix100_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup_PhiX$Ins, phix100_total_granges@ranges))
  phix100_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup_PhiX$Ins_rate, phix100_total_granges@ranges))
  phix100_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup_PhiX$Ins_rate, phix100_total_granges@ranges))
  phix100_total_granges$INS_upper <- phix100_total_granges$INS_mean + phix100_total_granges$INS_sd
  phix100_total_granges$INS_lower <- phix100_total_granges$INS_mean - phix100_total_granges$INS_sd
  
  phix100_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup_PhiX$Del, phix100_total_granges@ranges))
  phix100_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup_PhiX$Del_rate, phix100_total_granges@ranges))
  phix100_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup_PhiX$Del_rate, phix100_total_granges@ranges))
  phix100_total_granges$DEL_upper <- phix100_total_granges$DEL_mean + phix100_total_granges$DEL_sd
  phix100_total_granges$DEL_lower <- phix100_total_granges$DEL_mean - phix100_total_granges$DEL_sd
  
  # Calculate error rate per 100mer in each genome
  kmer_error_freq_synx <- data.frame(synx100_total_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Error_freq_mean = mean(Error_freq, na.rm = TRUE),
                     SUB_freq_mean = mean(SUB_freq, na.rm = TRUE),
                     INS_freq_mean = mean(INS_freq, na.rm = TRUE),
                     DEL_freq_mean = mean(DEL_freq, na.rm = TRUE),
                     GC_pct = mean(GC_pct, na.rm = TRUE),
                     Genome = 'SynX')
  
  kmer_error_freq_phix <- data.frame(phix100_total_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Error_freq_mean = mean(Error_freq, na.rm = TRUE),
                     SUB_freq_mean = mean(SUB_freq, na.rm = TRUE),
                     INS_freq_mean = mean(INS_freq, na.rm = TRUE),
                     DEL_freq_mean = mean(DEL_freq, na.rm = TRUE),
                     GC_pct = mean(GC_pct, na.rm = TRUE),
                     Genome = 'PhiX')
  
  # Plot 100mer GC content distribution
  kmer_error_freq_synx %>%
    arrange(GC_pct) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(y = SEQUENCE, x = GC_pct)) +
    theme_classic() +
    xlim(0, 100) +
    xlab("GC (%)") +
    ylab("100mer") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ggsave(paste0("project_results/", ds, "/SynX_100mer_", ds,"_GC_percent.pdf")) 
  
  kmer_error_freq_phix %>%
    arrange(GC_pct) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(y = SEQUENCE, x = GC_pct)) +
    theme_classic() +
    xlim(0, 100) +
    xlab("GC (%)") +
    ylab("100mer") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ggsave(paste0("project_results/", ds, "/PhiX_100mer_", ds,"_GC_percent.pdf")) 
  
  rbind(kmer_error_freq_synx, kmer_error_freq_phix) %>%
    arrange(GC_pct) %>%
    mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
    ggplot( ) +
    geom_point(aes(y = SEQUENCE, x = GC_pct, color = Genome)) +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    scale_x_continuous(breaks = c(25, 50, 75)) +
    theme_classic() +
    xlab("GC (%)") +
    ylab("100mer") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 0.5) +
    facet_wrap(~Genome, nrow = 2, ncol = 1) +
    ggsave(paste0("project_results/", ds, "/Synx_PhiX_100mer_", ds,"_GC_percent.pdf"))
  
  # Plot error frequency by GC content
  melt(id.vars = 'GC_pct', data = kmer_error_freq_synx %>% 
         select(GC_pct, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = log2(value), x = GC_pct, group = GC_pct, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    xlab("100mer GC content") +
    ylab("Error Frequency (log2)") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/", ds, "/SynX_100mer_", ds,"_error_frequency_w_GC_content.pdf"))
  
  melt(id.vars = 'GC_pct', data = kmer_error_freq_phix %>% 
         select(GC_pct, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
    ggplot() +
    geom_boxplot(aes(y = log2(value), x = GC_pct, group = GC_pct, fill = variable)) +
    scale_fill_npg(palette = c("nrc"), alpha = 1) +
    theme_classic() +
    scale_x_continuous(breaks = c(0, 50, 100)) +
    xlab("100mer GC content") +
    ylab("Error Frequency (log2)") +
    facet_wrap(~variable) +
    ggsave(paste0("project_results/", ds, "/PhiX_100mer_", ds,"_error_frequency_w_GC_content.pdf"))
  
  for(i in c('Error_freq_mean', 'SUB_freq_mean', 'INS_freq_mean', 'DEL_freq_mean')){
    rbind(kmer_error_freq_synx, kmer_error_freq_phix) %>%
      arrange(GC_pct) %>%
      ggplot() +
      geom_boxplot(aes_string(y = i, x = 'GC_pct', group = 'GC_pct', fill = 'Genome'), outlier.shape = NA) +
      scale_fill_npg(palette = c("nrc"), alpha = 1) +
      scale_x_continuous(breaks = c(25, 50, 75)) +
      scale_y_continuous(trans='log2') +
      theme_classic() +
      xlab("GC (%)") +
      ylab(i) +
      theme(aspect.ratio = 0.5) +
      facet_wrap(~Genome, nrow = 2, ncol = 1) +
      ggsave(paste0("project_results/", ds, "/Synx_PhiX_100mer_", ds,"_", i, "_frequency.pdf"))
  }
  
  # # Plot error rate correlation by 6mer in both genomes
  # # --------------------------------------------------------------------------
  # 
  # ### Calculate all unique 6mers in SynX
  # 
  # ONTDNA_pileup_SynX <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == "SynX", ]
  # synx_seq_F <- ONTDNA_pileup_SynX$REF_NT
  # synx_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup_SynX$REF_NT))
  # 
  # KD_k <- function(x){paste(x, collapse = "")}
  # synx_seq_F_6mers <- runner(x = synx_seq_F, k = 6, f = KD_k)
  # synx_seq_F_6mers <- synx_seq_F_6mers[str_length(synx_seq_F_6mers) == 6]
  # synx_seq_RC_6mers <- runner(x = synx_seq_RC, k = 6, f = KD_k)
  # synx_seq_RC_6mers <- synx_seq_RC_6mers[str_length(synx_seq_RC_6mers) == 6]
  # synx_seq_total_6mers <- unique(synx_seq_F_6mers, synx_seq_RC_6mers)
  # 
  # # Make ranges of 6mers in synx sequence
  # synx6_F_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_F_6mers), width = 6), strand = "+", SEQUENCE = synx_seq_F_6mers)
  # synx6_RC_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_RC_6mers), width = 6), strand = "+", SEQUENCE = rev(synx_seq_RC_6mers))
  # synx6_total_granges <- c(synx6_F_granges, synx6_RC_granges)
  # 
  # # Calculate GC content in SynX 6mers
  # GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G'))}
  # synx6_total_granges$GC_pct <- GC_k(synx6_total_granges$SEQUENCE)
  # 
  # # Calculate error rates in SynX 6mers
  # synx6_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup_SynX$Error_counts, synx6_total_granges@ranges))
  # synx6_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup_SynX$Error_rate, synx6_total_granges@ranges))
  # synx6_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup_SynX$Error_rate, synx6_total_granges@ranges))
  # synx6_total_granges$Error_upper <- synx6_total_granges$Error_mean + synx6_total_granges$Error_sd
  # synx6_total_granges$Error_lower <- synx6_total_granges$Error_mean - synx6_total_granges$Error_sd
  # 
  # synx6_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup_SynX$SUB_counts, synx6_total_granges@ranges))
  # synx6_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup_SynX$SUB_rate, synx6_total_granges@ranges))
  # synx6_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup_SynX$SUB_rate, synx6_total_granges@ranges))
  # synx6_total_granges$SUB_upper <- synx6_total_granges$SUB_mean + synx6_total_granges$SUB_sd
  # synx6_total_granges$SUB_lower <- synx6_total_granges$SUB_mean - synx6_total_granges$SUB_sd
  # 
  # synx6_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup_SynX$INDEL_counts, synx6_total_granges@ranges))
  # synx6_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup_SynX$INDEL_rate, synx6_total_granges@ranges))
  # synx6_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup_SynX$INDEL_rate, synx6_total_granges@ranges))
  # synx6_total_granges$INDEL_upper <- synx6_total_granges$INDEL_mean + synx6_total_granges$INDEL_sd
  # synx6_total_granges$INDEL_lower <- synx6_total_granges$INDEL_mean - synx6_total_granges$INDEL_sd
  # 
  # synx6_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup_SynX$Ins, synx6_total_granges@ranges))
  # synx6_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup_SynX$Ins_rate, synx6_total_granges@ranges))
  # synx6_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup_SynX$Ins_rate, synx6_total_granges@ranges))
  # synx6_total_granges$INS_upper <- synx6_total_granges$INS_mean + synx6_total_granges$INS_sd
  # synx6_total_granges$INS_lower <- synx6_total_granges$INS_mean - synx6_total_granges$INS_sd
  # 
  # synx6_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup_SynX$Del, synx6_total_granges@ranges))
  # synx6_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup_SynX$Del_rate, synx6_total_granges@ranges))
  # synx6_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup_SynX$Del_rate, synx6_total_granges@ranges))
  # synx6_total_granges$DEL_upper <- synx6_total_granges$DEL_mean + synx6_total_granges$DEL_sd
  # synx6_total_granges$DEL_lower <- synx6_total_granges$DEL_mean - synx6_total_granges$DEL_sd
  # 
  # ### Calculate all unique 6mers in PhiX
  # 
  # ONTDNA_pileup_PhiX <- ONTDNA_pileup_total[ONTDNA_pileup_total$Chrom == "PhiX", ]
  # phix_seq_F <- ONTDNA_pileup_PhiX$REF_NT
  # phix_seq_RC <- toupper(spgs::reverseComplement(ONTDNA_pileup_PhiX$REF_NT))
  # 
  # KD_k <- function(x){paste(x, collapse = "")}
  # phix_seq_F_6mers <- runner(x = phix_seq_F, k = 6, f = KD_k)
  # phix_seq_F_6mers <- phix_seq_F_6mers[str_length(phix_seq_F_6mers) == 6]
  # phix_seq_RC_6mers <- runner(x = phix_seq_RC, k = 6, f = KD_k)
  # phix_seq_RC_6mers <- phix_seq_RC_6mers[str_length(phix_seq_RC_6mers) == 6]
  # phix_seq_total_6mers <- unique(phix_seq_F_6mers, phix_seq_RC_6mers)
  # 
  # # Make ranges of 6mers in phix sequence
  # phix6_F_granges <- GRanges(seqnames = 'PhiX', ranges = IRanges(start = 1:length(phix_seq_F_6mers), width = 6), strand = "+", SEQUENCE = phix_seq_F_6mers)
  # phix6_RC_granges <- GRanges(seqnames = 'PhiX', ranges = IRanges(start = 1:length(phix_seq_RC_6mers), width = 6), strand = "+", SEQUENCE = rev(phix_seq_RC_6mers))
  # phix6_total_granges <- c(phix6_F_granges, phix6_RC_granges)
  # 
  # # Calculate GC content
  # GC_k <- function(x){(str_count(x, pattern = 'C') + str_count(x, pattern = 'G'))}
  # phix6_total_granges$GC_pct <- GC_k(phix6_total_granges$SEQUENCE)
  # 
  # # Calculate error rates in PhiX 6mers
  # phix6_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup_PhiX$Error_counts, phix6_total_granges@ranges))
  # phix6_total_granges$Error_mean <- mean(extractList(ONTDNA_pileup_PhiX$Error_rate, phix6_total_granges@ranges))
  # phix6_total_granges$Error_sd <- sd(extractList(ONTDNA_pileup_PhiX$Error_rate, phix6_total_granges@ranges))
  # phix6_total_granges$Error_upper <- phix6_total_granges$Error_mean + phix6_total_granges$Error_sd
  # phix6_total_granges$Error_lower <- phix6_total_granges$Error_mean - phix6_total_granges$Error_sd
  # 
  # phix6_total_granges$SUB_freq <- sum(extractList(ONTDNA_pileup_PhiX$SUB_counts, phix6_total_granges@ranges))
  # phix6_total_granges$SUB_mean <- mean(extractList(ONTDNA_pileup_PhiX$SUB_rate, phix6_total_granges@ranges))
  # phix6_total_granges$SUB_sd <- sd(extractList(ONTDNA_pileup_PhiX$SUB_rate, phix6_total_granges@ranges))
  # phix6_total_granges$SUB_upper <- phix6_total_granges$SUB_mean + phix6_total_granges$SUB_sd
  # phix6_total_granges$SUB_lower <- phix6_total_granges$SUB_mean - phix6_total_granges$SUB_sd
  # 
  # phix6_total_granges$INDEL_freq <- sum(extractList(ONTDNA_pileup_PhiX$INDEL_counts, phix6_total_granges@ranges))
  # phix6_total_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup_PhiX$INDEL_rate, phix6_total_granges@ranges))
  # phix6_total_granges$INDEL_sd <- sd(extractList(ONTDNA_pileup_PhiX$INDEL_rate, phix6_total_granges@ranges))
  # phix6_total_granges$INDEL_upper <- phix6_total_granges$INDEL_mean + phix6_total_granges$INDEL_sd
  # phix6_total_granges$INDEL_lower <- phix6_total_granges$INDEL_mean - phix6_total_granges$INDEL_sd
  # 
  # phix6_total_granges$INS_freq <- sum(extractList(ONTDNA_pileup_PhiX$Ins, phix6_total_granges@ranges))
  # phix6_total_granges$INS_mean <- mean(extractList(ONTDNA_pileup_PhiX$Ins_rate, phix6_total_granges@ranges))
  # phix6_total_granges$INS_sd <- sd(extractList(ONTDNA_pileup_PhiX$Ins_rate, phix6_total_granges@ranges))
  # phix6_total_granges$INS_upper <- phix6_total_granges$INS_mean + phix6_total_granges$INS_sd
  # phix6_total_granges$INS_lower <- phix6_total_granges$INS_mean - phix6_total_granges$INS_sd
  # 
  # phix6_total_granges$DEL_freq <- sum(extractList(ONTDNA_pileup_PhiX$Del, phix6_total_granges@ranges))
  # phix6_total_granges$DEL_mean <- mean(extractList(ONTDNA_pileup_PhiX$Del_rate, phix6_total_granges@ranges))
  # phix6_total_granges$DEL_sd <- sd(extractList(ONTDNA_pileup_PhiX$Del_rate, phix6_total_granges@ranges))
  # phix6_total_granges$DEL_upper <- phix6_total_granges$DEL_mean + phix6_total_granges$DEL_sd
  # phix6_total_granges$DEL_lower <- phix6_total_granges$DEL_mean - phix6_total_granges$DEL_sd
  # 
  # # Calculate error rate per 6mer in each genome
  # kmer_error_freq_synx <- data.frame(synx6_total_granges@elementMetadata) %>%
  #   group_by(SEQUENCE) %>%
  #   dplyr::summarize(Error_freq_mean = mean(Error_freq, na.rm = TRUE),
  #                    SUB_freq_mean = mean(SUB_freq, na.rm = TRUE),
  #                    INS_freq_mean = mean(INS_freq, na.rm = TRUE),
  #                    DEL_freq_mean = mean(DEL_freq, na.rm = TRUE),
  #                    GC_pct = mean(GC_pct, na.rm = TRUE),
  #                    Genome = 'SynX')
  # 
  # kmer_error_freq_phix <- data.frame(phix6_total_granges@elementMetadata) %>%
  #   group_by(SEQUENCE) %>%
  #   dplyr::summarize(Error_freq_mean = mean(Error_freq, na.rm = TRUE),
  #                    SUB_freq_mean = mean(SUB_freq, na.rm = TRUE),
  #                    INS_freq_mean = mean(INS_freq, na.rm = TRUE),
  #                    DEL_freq_mean = mean(DEL_freq, na.rm = TRUE),
  #                    GC_pct = mean(GC_pct, na.rm = TRUE),
  #                    Genome = 'PhiX')
  # 
  # # Plot 6mer GC content distribution
  # kmer_error_freq_synx %>%
  #   arrange(GC_pct) %>%
  #   mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
  #   ggplot( ) +
  #   geom_point(aes(y = SEQUENCE, x = GC_pct)) +
  #   theme_classic() +
  #   xlim(0, 6) +
  #   xlab("GC (%)") +
  #   ylab("6mer") +
  #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #   ggsave(paste0("project_results/", ds, "/SynX_6mer_", ds,"_GC_percent.pdf")) 
  # 
  # kmer_error_freq_phix %>%
  #   arrange(GC_pct) %>%
  #   mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
  #   ggplot( ) +
  #   geom_point(aes(y = SEQUENCE, x = GC_pct)) +
  #   theme_classic() +
  #   xlim(0, 6) +
  #   xlab("GC (%)") +
  #   ylab("6mer") +
  #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  #   ggsave(paste0("project_results/", ds, "/PhiX_6mer_", ds,"_GC_percent.pdf")) 
  # 
  # rbind(kmer_error_freq_synx, kmer_error_freq_phix) %>%
  #   arrange(GC_pct) %>%
  #   mutate(SEQUENCE = factor(SEQUENCE, levels = SEQUENCE)) %>%
  #   ggplot( ) +
  #   geom_point(aes(y = SEQUENCE, x = GC_pct, color = Genome)) +
  #   scale_color_npg(palette = c("nrc"), alpha = 1) +
  #   scale_x_continuous(breaks = c(25, 50, 75)) +
  #   theme_classic() +
  #   xlab("GC (%)") +
  #   ylab("6mer") +
  #   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), aspect.ratio = 0.5) +
  #   facet_wrap(~Genome, nrow = 2, ncol = 1) +
  #   ggsave(paste0("project_results/", ds, "/Synx_PhiX_6mer_", ds,"_GC_percent.pdf"))
  # 
  # # Plot error frequency by GC content
  # melt(id.vars = 'GC_pct', data = kmer_error_freq_synx %>% 
  #        select(GC_pct, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
  #   ggplot() +
  #   geom_boxplot(aes(y = log2(value), x = GC_pct, group = GC_pct, fill = variable)) +
  #   scale_fill_npg(palette = c("nrc"), alpha = 1) +
  #   theme_classic() +
  #   scale_x_continuous(breaks = c(0, 50, 6)) +
  #   xlab("6mer GC content") +
  #   ylab("Error Frequency (log2)") +
  #   facet_wrap(~variable) +
  #   ggsave(paste0("project_results/", ds, "/SynX_6mer_", ds,"_error_frequency_w_GC_content.pdf"))
  # 
  # melt(id.vars = 'GC_pct', data = kmer_error_freq_phix %>% 
  #        select(GC_pct, Error_freq_mean, SUB_freq_mean, INS_freq_mean, DEL_freq_mean)) %>%
  #   ggplot() +
  #   geom_boxplot(aes(y = log2(value), x = GC_pct, group = GC_pct, fill = variable)) +
  #   scale_fill_npg(palette = c("nrc"), alpha = 1) +
  #   theme_classic() +
  #   scale_x_continuous(breaks = c(0, 50, 6)) +
  #   xlab("6mer GC content") +
  #   ylab("Error Frequency (log2)") +
  #   facet_wrap(~variable) +
  #   ggsave(paste0("project_results/", ds, "/PhiX_6mer_", ds,"_error_frequency_w_GC_content.pdf"))
  # 
  # for(i in c('Error_freq_mean', 'SUB_freq_mean', 'INS_freq_mean', 'DEL_freq_mean')){
  #   rbind(kmer_error_freq_synx, kmer_error_freq_phix) %>%
  #     arrange(GC_pct) %>%
  #     ggplot() +
  #     geom_boxplot(aes_string(y = i, x = 'GC_pct', group = 'GC_pct', fill = 'Genome'), outlier.shape = NA) +
  #     scale_fill_npg(palette = c("nrc"), alpha = 1) +
  #     scale_x_continuous(breaks = c(25, 50, 75)) +
  #     scale_y_continuous(trans='log2') +
  #     theme_classic() +
  #     xlab("GC (%)") +
  #     ylab(i) +
  #     theme(aspect.ratio = 0.5) +
  #     facet_wrap(~Genome, nrow = 2, ncol = 1) +
  #     ggsave(paste0("project_results/", ds, "/Synx_PhiX_6mer_", ds,"_", i, "_frequency.pdf"))
  # }
  # synx_errors <- kmer_error_freq_synx %>%
  #   select(SEQUENCE, DEL_freq_mean)
  # 
  # phix_errors <- kmer_error_freq_phix %>%
  #   select(SEQUENCE, DEL_freq_mean)
  # 
  # idx <- match(synx_errors$SEQUENCE, phix_errors$SEQUENCE)
  # synx_errors$DEL_freq_mean_synx <- synx_errors$DEL_freq_mean
  # synx_errors$DEL_freq_mean_phix <- phix_errors$DEL_freq_mean [idx]
  # 
  # synx_errors <- synx_errors[! is.na(synx_errors$DEL_freq_mean_phix),]
  # 
  # ggplot(synx_errors) +
  #   geom_point(aes(x=log2(DEL_freq_mean_synx), y= log2(DEL_freq_mean_phix))) +
  #   geom_smooth(aes(x=log2(DEL_freq_mean_synx), y= log2(DEL_freq_mean_phix)), method = 'lm') +
  #   theme_classic() +
  #   theme(aspect.ratio = 1)
  # 
  # cor(synx_errors$DEL_freq_mean_synx, synx_errors$DEL_freq_mean_phix, method = 'pearson')
}
