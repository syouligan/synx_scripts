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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
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
    theme(aspect.ratio = 1) +
    ggsave(paste0("project_results/SynX_kmers/SynX_kmer_", ds,"_SUB_frequency.pdf"))
}

# Plot kmer error distribution for ONT datasets across expressed regions
# --------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')

# Load synx feature annotations and annotate ranges
synx_annotations <- read.table("design/custom_features.annotations.tab", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
idx <- match(synx_ranges$name, synx_annotations$NAME)
synx_ranges$TYPE <- synx_annotations$TYPE [idx]
synx_ranges$TYPE <- gsub(" ", "_", synx_ranges$TYPE)
synx_ranges$SEQUENCE <- synx_annotations$SEQUENCE [idx]

# Subset to expressed features in SynX sequence
synx_ranges_expressed <- synx_ranges[synx_ranges@ranges@width > 200]

# List datasets to plot script over
data_sets <- c('ONT_DNA_Fragmentase', 'ONT_cDNA_SP6', 'ONT_cDNA_T7')

# Subset kmers to those expressed in SynX sequence and caluclate mean error rate per unique kmer
dslist.names <- data_sets
dslist3 <- vector("list", length(dslist.names))
names(dslist3) <- dslist.names

for(ds in data_sets) {
  expressed_kmers <- subsetByOverlaps(dslist1[[ds]], synx_ranges_expressed)
  dslist3[[ds]] <- data.frame(expressed_kmers@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Kmer_count = n(),
                     Error_rate_mean = mean(Error_mean, na.rm = TRUE),
                     Depth = max(Depth_mean, na.rm = TRUE),
                     Dataset = unique(ds))
}

# Create dataframe with error rates for each dataset
kmer_exp_df <- data.frame()
for(ds in data_sets) {
  tmp <- data.frame(dslist3[[ds]])
  kmer_exp_df <- rbind(kmer_exp_df, tmp)
}

# Order points based on ONT-DNA-Fragmentase
kmer_order <- kmer_exp_df[kmer_exp_df$Dataset == 'ONT_DNA_Fragmentase',] %>%
  arrange(desc(Error_rate_mean))
kmer_order$SEQUENCE <- factor(kmer_order$SEQUENCE, kmer_order$SEQUENCE)

# Plot distribution of kmers in each sample in datasets
kmer_exp_df$SEQUENCE <- factor(kmer_exp_df$SEQUENCE, kmer_order$SEQUENCE)
kmer_exp_df$Dataset <- factor(kmer_exp_df$Dataset, data_sets)

ggplot(kmer_exp_df) +
  geom_point(aes(x = SEQUENCE, y = Error_rate_mean, color = Dataset)) +
  theme_classic() +
  xlab("Kmer") +
  ylab("Error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(aspect.ratio = 0.4) +
  facet_wrap(~ Dataset, nrow = length(data_sets), ncol = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_ONT_sample_comparison_Error_frequency.pdf"))

# Plot kmer correlation in DNA and RNA in ONT data
SP6 <- kmer_exp_df[kmer_exp_df$Dataset == 'ONT_cDNA_SP6',]
idx <- match(kmer_order$SEQUENCE, SP6$SEQUENCE)
kmer_order$SP6_error <- SP6$Error_rate_mean [idx]
mod <- lm(Error_rate_mean ~ SP6_error, kmer_order)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

ggplot(kmer_order, aes(y = Error_rate_mean, x = SP6_error)) +
  annotate(geom = "text", y = 0.5, x = 0.1, label = eq_r2) +
  geom_point(alpha = 0.3) +
  geom_smooth(method='lm', se=TRUE, color = 'black') +
  theme_classic() +
  xlab("ONT_cDNA_SP6_error") +
  ylab("ONT_DNA_error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(aspect.ratio = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_ONT_DNA_cDNA_SP6_Error_frequency.pdf"))

T7 <- kmer_exp_df[kmer_exp_df$Dataset == 'ONT_cDNA_T7',]
idx <- match(kmer_order$SEQUENCE, T7$SEQUENCE)
kmer_order$T7_error <- T7$Error_rate_mean [idx]
mod <- lm(Error_rate_mean ~ T7_error, kmer_order)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

ggplot(kmer_order, aes(y = Error_rate_mean, x = T7_error)) +
  annotate(geom = "text", y = 0.5, x = 0.1, label = eq_r2) +
  geom_point(alpha = 0.3) +
  geom_smooth(method='lm', se=TRUE, color = 'black') +
  theme_classic() +
  xlab("ONT_cDNA_T7_error") +
  ylab("ONT_DNA_error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(aspect.ratio = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_ONT_DNA_cDNA_T7_Error_frequency.pdf"))

# Plot kmer error distribution for Illumina datasets across expressed regions
# --------------------------------------------------

# List datasets to plot script over
data_sets <- c('Illumina_DNASynX', 'Illumina_SP6', 'Illumina_T7')

# Subset kmers to those expressed in SynX sequence and caluclate mean error rate per unique kmer
dslist.names <- data_sets
dslist3 <- vector("list", length(dslist.names))
names(dslist3) <- dslist.names

for(ds in data_sets) {
  expressed_kmers <- subsetByOverlaps(dslist1[[ds]], synx_ranges_expressed)
  dslist3[[ds]] <- data.frame(expressed_kmers@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    dplyr::summarize(Kmer_count = n(),
                     Error_rate_mean = mean(Error_mean, na.rm = TRUE),
                     Depth = max(Depth_mean, na.rm = TRUE),
                     Dataset = unique(ds))
}

# Create dataframe with error rates for each dataset
kmer_exp_df <- data.frame()
for(ds in data_sets) {
  tmp <- data.frame(dslist3[[ds]])
  kmer_exp_df <- rbind(kmer_exp_df, tmp)
}

# Order points based on ONT-DNA-Fragmentase
kmer_order <- kmer_exp_df[kmer_exp_df$Dataset == 'Illumina_DNASynX',] %>%
  arrange(desc(Error_rate_mean))
kmer_order$SEQUENCE <- factor(kmer_order$SEQUENCE, kmer_order$SEQUENCE)

# Plot distribution of kmers in each sample in datasets
kmer_exp_df$SEQUENCE <- factor(kmer_exp_df$SEQUENCE, kmer_order$SEQUENCE)
kmer_exp_df$Dataset <- factor(kmer_exp_df$Dataset, data_sets)

ggplot(kmer_exp_df) +
  geom_point(aes(x = SEQUENCE, y = Error_rate_mean, color = Dataset)) +
  theme_classic() +
  xlab("Kmer") +
  ylab("Error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(aspect.ratio = 0.4) +
  facet_wrap(~ Dataset, nrow = length(data_sets), ncol = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_Illumina_sample_comparison_Error_frequency.pdf"))

# Plot kmer correlation in DNA and RNA in Illumina data
SP6 <- kmer_exp_df[kmer_exp_df$Dataset == 'Illumina_SP6',]
idx <- match(kmer_order$SEQUENCE, SP6$SEQUENCE)
kmer_order$SP6_error <- SP6$Error_rate_mean [idx]
mod <- lm(Error_rate_mean ~ SP6_error, kmer_order)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

ggplot(kmer_order, aes(y = Error_rate_mean, x = SP6_error)) +
  annotate(geom = "text", y = 0.15, x = 0.1, label = eq_r2) +
  geom_point(alpha = 0.3) +
  geom_smooth(method='lm', se=TRUE, color = 'black') +
  theme_classic() +
  xlab("Illumina_cDNA_SP6_error") +
  ylab("Illumina_DNA_error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(aspect.ratio = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_Illumina_DNA_cDNA_SP6_Error_frequency.pdf"))

T7 <- kmer_exp_df[kmer_exp_df$Dataset == 'Illumina_T7',]
idx <- match(kmer_order$SEQUENCE, T7$SEQUENCE)
kmer_order$T7_error <- T7$Error_rate_mean [idx]
mod <- lm(Error_rate_mean ~ T7_error, kmer_order)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

ggplot(kmer_order, aes(y = Error_rate_mean, x = T7_error)) +
  annotate(geom = "text", y = 0.15, x = 0.1, label = eq_r2) +
  geom_point(alpha = 0.3) +
  geom_smooth(method='lm', se=TRUE, color = 'black') +
  theme_classic() +
  xlab("Illumina_cDNA_T7_error") +
  ylab("Illumina_DNA_error rate") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  theme(aspect.ratio = 1) +
  ggsave(paste0("project_results/SynX_kmers/SynX_kmer_Illumina_DNA_cDNA_T7_Error_frequency.pdf"))


