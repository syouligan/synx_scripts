#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX copy number ladder
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

# Identify unique 31mers in copy number features
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')

# Load synx feature annotations and annotate ranges
synx_annotations <- read.table("design/custom_features.annotations.tab", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
idx <- match(synx_ranges$name, synx_annotations$NAME)
synx_ranges$TYPE <- synx_annotations$TYPE [idx]
synx_ranges$TYPE <- gsub(" ", "_", synx_ranges$TYPE)
synx_ranges$SEQUENCE <- synx_annotations$SEQUENCE [idx]

# Make list of 31mers for each cn element
KD_k <- function(x){paste(x, collapse = "")}

cnlist.names <- c('1cn', '2cn', '3cn', '4cn')
cnlist <- vector("list", length(cnlist.names))
names(cnlist) <- cnlist.names

for (cn in c('1cn', '2cn', '3cn', '4cn')) {
  cn_seq <- unlist(strsplit(unique(synx_ranges[synx_ranges$name == cn]$SEQUENCE), split = ""))
  cn_31mers <- runner(x = cn_seq, k = 31, f = KD_k)
  cn_31mers <- cn_31mers[str_length(cn_31mers) == 31]
  cnlist[[cn]] <- cn_31mers
}

# Calculate coverage across copy number features
# --------------------------------------------------------------------------

# List datasets to loop script over
data_sets <- c('Illumina_DNASynX', 'Illumina_SP6', 'Illumina_T7', 'ONT_DNA', 'ONT_DNA_BamHI', 'ONT_DNA_EcoRI', 'ONT_DNA_HindIII', 'ONT_DNA_Fragmentase', 'ONT_cDNA_SP6', 'ONT_cDNA_T7', 'badread_simulator', 'Illumina_DNASynX_simulation')

for (ds in data_sets) {
  # # Load sequencing error rates
  # ONTDNA_pileup <- read.table(paste0("project_results/", ds, "/", ds, ".bam.bed.tsv"), header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
  # 
  # # Find reference mapping rates
  # get_correct <- function(df){
  #   bases <- c('A', 'T', 'G', 'C')
  #   ref <- df['REF_NT']
  #   correct <- bases[bases == ref]
  #   sum(as.numeric(df[correct]))
  # }
  # 
  # ONTDNA_pileup$Depth <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del + ONTDNA_pileup$Coverage
  # ONTDNA_pileup$REF_counts <- apply(ONTDNA_pileup, 1, get_correct)
  # ONTDNA_pileup$REF_rate <- ONTDNA_pileup$REF_counts / ONTDNA_pileup$Depth
  # 
  # # Find substitution rates
  # get_subs <- function(df){
  #   bases <- c('A', 'T', 'G', 'C')
  #   ref <- df['REF_NT']
  #   subs <- bases[!(bases == ref)]
  #   sum(as.numeric(df[subs]))
  # }
  # 
  # ONTDNA_pileup$SUB_counts <- apply(ONTDNA_pileup, 1, get_subs)
  # ONTDNA_pileup$SUB_rate <- ONTDNA_pileup$SUB_counts / ONTDNA_pileup$Depth
  # 
  # # Find INDEL rates
  # ONTDNA_pileup$Ins_rate <- ONTDNA_pileup$Ins / ONTDNA_pileup$Depth
  # ONTDNA_pileup$Del_rate <- ONTDNA_pileup$Del / ONTDNA_pileup$Depth
  # ONTDNA_pileup$INDEL_counts <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del
  # ONTDNA_pileup$INDEL_rate <- ONTDNA_pileup$INDEL_counts / ONTDNA_pileup$Depth
  # 
  # # Find error rates
  # ONTDNA_pileup$Error_counts <- ONTDNA_pileup$SUB_counts + ONTDNA_pileup$INDEL_counts
  # ONTDNA_pileup$Error_rate <- ONTDNA_pileup$Error_counts / ONTDNA_pileup$Depth
  # 
  # # Find mean coverage
  # mean_depth <- mean(ONTDNA_pileup$Depth)
  # 
  # # Calculate error rates across 31mers
  # # --------------------------------------------------------------------------
  # 
  # # Calculate all unique 31mers in Synx
  # synx_seq_F <- ONTDNA_pileup$REF_NT
  # synx_seq_F_31mers <- runner(x = synx_seq_F, k = 31, f = KD_k)
  # synx_seq_total_31mers <- synx_seq_F_31mers[str_length(synx_seq_F_31mers) == 31]
  # 
  # # Make ranges of 10mers in synx sequence
  # k_total_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(synx_seq_F_31mers), width = 31), strand = "+", SEQUENCE = synx_seq_F_31mers)
  # 
  # # Calculate error rates
  # k_total_granges$Error_freq <- sum(extractList(ONTDNA_pileup$Error_counts, k_total_granges@ranges))
  # k_total_granges$Error_rate <- mean(extractList(ONTDNA_pileup$Error_rate, k_total_granges@ranges))
  # k_total_granges$Coverage_mean <- mean(extractList(ONTDNA_pileup$Depth, k_total_granges@ranges)) / mean_depth

  # # Calculate mean mean error rates and coverage across kmers
  # # --------------------------------------------------------------------------
  # 
  # # Calculate error rate per 31mer
  # kmer_error_freq <- data.frame(k_total_granges@elementMetadata) %>%
  #   group_by(SEQUENCE) %>%
  #   dplyr::summarize(Error_rate_mean = mean(Error_rate, na.rm = TRUE),
  #                    Coverage_total = sum(Coverage_mean))
  # 
  # # Identify 31mers found in the copy number elements
  # for (cn in c('1cn', '2cn', '3cn', '4cn')) {
  #   cn_kmers <- cnlist[[cn]]
  #   kmer_error_freq[,paste0('CN', gsub("cn", "", cn))] <- is.element(kmer_error_freq$SEQUENCE, cn_kmers)
  # }
  # 
  # kmer_error_freq <- kmer_error_freq %>%
  #   mutate(
  #     Copy_number = case_when(
  #       CN1 == TRUE ~ "CN1",
  #       CN2 == TRUE ~ "CN2",
  #       CN3 == TRUE ~ "CN3",
  #       CN4 == TRUE ~ "CN4"
  #     )
  #   )
  # 
  # kmer_error_freq$Coverage <- as.numeric(kmer_error_freq$Coverage_total)
  # kmer_error_freq$Copies <- as.numeric(gsub('CN', '', kmer_error_freq$Copy_number))
  # 
  # # Subset to copy number 31mers and plot coverage
  # cn_elements <- kmer_error_freq[!is.na(kmer_error_freq$Copy_number),]
  # 
  # # GET EQUATION AND R-SQUARED AS STRING
  # # SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
  # lm_eqn <- function(df, xx, yy){
  #   m <- lm(Coverage ~ Copies, df);
  #   eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
  #                    list(a = format(unname(coef(m)[1]), digits = 2),
  #                         b = format(unname(coef(m)[2]), digits = 2),
  #                         r2 = format(summary(m)$r.squared, digits = 3)))
  #   as.character(as.expression(eq));
  # }
  # 
  # mod <- lm(Coverage ~ Copies, cn_elements)
  # 
  # ggplot(mod, aes(x = .resid)) +
  #   geom_histogram() +
  #   xlim(-max(abs(mod$resid)), max(abs(mod$resid))) +
  #   theme_classic() +
  #   ggsave(paste0("project_results/copy_number/Quantitative_ladder_residuals_", ds, ".pdf"))
  # 
  # ggplot(cn_elements, aes(y = Coverage, x = Copies)) +
  #   geom_point(aes(color=Copy_number)) +
  #   geom_smooth(method='lm', se=TRUE, color = 'black') +
  #   # geom_text(x = 1.5, y = 8000, label = lm_eqn(cn_elements), parse = TRUE) +
  #   scale_color_npg(palette = c("nrc"), alpha = 1) +
  #   xlab("Copy number") +
  #   ylab("Per-base coverage") +
  #   # ylim(1000, 8500) +
  #   theme_classic() +
  #   ggsave(paste0("project_results/copy_number/Quantitative_ladder_lm_", ds, ".pdf"))

  # Plot counts of 31mers based on jellyfish
  # --------------------------------------------------------------------------
  
  all_kmer_counts <- read.table(paste0("project_results/", ds, "/jellyfish/All_31mers_counts.tsv"), header = TRUE)
  colnames(all_kmer_counts) <- c("Sequence", "Count")
  
  # Identify 31mers found in the copy number elements
  for (cn in c('1cn', '2cn', '3cn', '4cn')) {
    cn_kmers <- cnlist[[cn]]
    all_kmer_counts[,paste0('CN', gsub("cn", "", cn))] <- is.element(all_kmer_counts$Sequence, cn_kmers)
  }
  
  all_kmer_counts <- all_kmer_counts %>%
    mutate(
      Copy_number = case_when(
        CN1 == TRUE ~ "CN1",
        CN2 == TRUE ~ "CN2",
        CN3 == TRUE ~ "CN3",
        CN4 == TRUE ~ "CN4"
      )
    )
  
  # Normalise counts by mean across whole SynX
  all_kmer_counts$Count <- as.numeric(all_kmer_counts$Count)
  all_kmer_counts$Normalised_count <- as.numeric(all_kmer_counts$Count/mean(all_kmer_counts$Count))
  all_kmer_counts$Copies <- as.numeric(gsub('CN', '', all_kmer_counts$Copy_number))
  
  # Subset to kmers in copy number elements
  cn_elements <- all_kmer_counts[!is.na(all_kmer_counts$Copy_number),]
  
  # Calculate linear model based on CN counts, find equation and R2
  mod <- lm(Normalised_count ~ Copies, cn_elements)
  eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
  r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
  eq_r2 <- paste0(eq, " ", r2)
  
  # Plot jely fish kmer counts across copy number elements with linear model
  ggplot(cn_elements, aes(y = Normalised_count, x = Copies)) +
    annotate(geom = "text", x = 2, y = 3, label = eq_r2) +
    geom_point(aes(color=Copy_number)) +
    geom_smooth(method='lm', se=TRUE, color = 'black') +
    scale_color_npg(palette = c("nrc"), alpha = 1) +
    xlab("Copy number") +
    ylab("Per-base coverage") +
    # ylim(1000, 8500) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    ggsave(paste0("project_results/copy_number/Quantitative_ladder_jellyfish_lm_", ds, ".pdf"))
  
  # Plot jely fish kmer counts across copy number elements linear model residuals
  ggplot(mod, aes(x = .resid)) +
    geom_density(adjust = 2) +
    xlim(-max(abs(mod$resid)), max(abs(mod$resid))) +
    xlim(-2.5, 2.5) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    ggsave(paste0("project_results/copy_number/Quantitative_ladder_jellyfish_residuals_", ds, ".pdf"))
  
  }
 