#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Error rate (subs, indels) across synX elements ONT-DNA
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
ONTDNA_pileup$REF_pct <- ONTDNA_pileup$REF_counts / ONTDNA_pileup$Depth

# Find substitution rates
get_subs <- function(df){
  bases <- c('A', 'T', 'G', 'C')
  ref <- df['REF_NT']
  subs <- bases[!(bases == ref)]
  sum(as.numeric(df[subs]))
}

ONTDNA_pileup$SUBS_counts <- apply(ONTDNA_pileup, 1, get_subs)
ONTDNA_pileup$SUBS_pct <- ONTDNA_pileup$SUBS_counts / ONTDNA_pileup$Depth

# Find INDEL rates
ONTDNA_pileup$Ins_pct <- ONTDNA_pileup$Ins / ONTDNA_pileup$Depth
ONTDNA_pileup$Del_pct <- ONTDNA_pileup$Del / ONTDNA_pileup$Depth
ONTDNA_pileup$INDEL_counts <- ONTDNA_pileup$Ins + ONTDNA_pileup$Del
ONTDNA_pileup$INDEL_pct <- ONTDNA_pileup$INDEL_counts / ONTDNA_pileup$Depth

# Find error rates
ONTDNA_pileup$Error_counts <- ONTDNA_pileup$SUBS_counts + ONTDNA_pileup$INDEL_counts
ONTDNA_pileup$Error_pct <- ONTDNA_pileup$Error_counts / ONTDNA_pileup$Depth

# Caluculate GCscore for each base
GC_k <- function(x){str_count(paste(x, collapse = ""), pattern = 'C') + str_count(paste(x, collapse = ""), pattern = 'G')}
GC_score <- data.frame(row.names = 1:length(ONTDNA_pileup$REF_NT))
for(kn in 0:5){
  GC_score[,paste0('GC_k', kn)]<- runner(x = ONTDNA_pileup$REF_NT, k = 6, f = GC_k, lag = -kn)
}
ONTDNA_pileup$GC_score <- rowMeans(as.matrix(GC_score))

# Caluculate homopolymer score for each base
homo_score <- function(x){length(unique(x)) == 1}
HP_score <- data.frame(row.names = 1:length(ONTDNA_pileup$REF_NT))
for(kn in 2:6){
  tmp1 <- data.frame(row.names = 1:length(ONTDNA_pileup$REF_NT))
  for(k in 0:(kn-1)){
    tmp1[,paste0('HP_k', kn, '_lag_', k)] <- runner(x = ONTDNA_pileup$REF_NT, k = kn, f = homo_score, lag = -k)
    }
  HP_score[,paste0('HP_k', kn)] <- rowSums(tmp1) > 0
  }
ONTDNA_pileup$HP_score <- rowSums(HP_score) + 1

# Caluculate k-mer diversity score for each base
KD_k <- function(x){paste(x, collapse = "")}
KD_score <- data.frame(row.names = 1:length(ONTDNA_pileup$REF_NT))
for(kn in 0:5){
  KD_score[,paste0('KD_k', kn)]<- runner(x = ONTDNA_pileup$REF_NT, k = 6, f = KD_k, lag = -kn)
}

div_k <- function(x){length(unique(x))}
ONTDNA_pileup$KD_score <- (6 - apply(X = KD_score, 1, FUN = div_k))

# Plot error rates across features
# --------------------------------------------------------------------------

synx_ranges$Error_fraction <- sum(relist(ONTDNA_pileup$Error_counts, synx_ranges@ranges)) / sum(relist(ONTDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$SUB_fraction <- sum(relist(ONTDNA_pileup$SUBS_counts, synx_ranges@ranges)) / sum(relist(ONTDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$INDEL_fraction <- sum(relist(ONTDNA_pileup$INDEL_counts, synx_ranges@ranges)) / sum(relist(ONTDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$INS_fraction <- sum(relist(ONTDNA_pileup$Ins, synx_ranges@ranges)) / sum(relist(ONTDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$DEL_fraction <- sum(relist(ONTDNA_pileup$Del, synx_ranges@ranges)) / sum(relist(ONTDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$HP_score <- mean(relist(ONTDNA_pileup$HP_score, synx_ranges@ranges))

# Plot sequencing errors across performance units
synx_pus <- synx_ranges[synx_ranges$TYPE == "Sequencing_performance_unit", ]
synx_pus <- synx_pus[order(synx_pus$name), ]

# Plot GC content of GC performance units
synx_pus_GCs <- synx_pus[!grepl('PH', synx_pus$name), ]

ggplot(data.frame(synx_pus_GCs@elementMetadata)) +
  geom_bar(aes(y = GC_pct, x = name), stat = 'identity') +
  theme_minimal() +
  ggsave(paste0("project_results/ONT_DNA/SynX_performance_GC_content.pdf"))

# Error rates across GC-performance units
for(quant in c("Error_fraction", "SUB_fraction", "INDEL_fraction", "INS_fraction", "DEL_fraction")){
  ggplot(data.frame(synx_pus_GCs@elementMetadata)) +
    geom_point(aes_string(x = 'GC_pct', y = quant)) +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/SynX_performance_GC_", quant, ".pdf"))
}

# Error rates across homopolymer-performance units
synx_pus_HPs <- synx_pus[grepl('PH', synx_pus$name), ]

for(quant in c("Error_fraction", "SUB_fraction", "INDEL_fraction")){
  ggplot(data.frame(synx_pus_HPs@elementMetadata)) +
    geom_point(aes_string(x = 'HP_score', y = quant)) +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/SynX_performance_HomoPolymer_", quant, ".pdf"))
}

# Plot error rates across kmers
# --------------------------------------------------------------------------

# List all kmers in SynX genome and make Granges object
all_k <- runner(x = ONTDNA_pileup$REF_NT, k = 6, f = KD_k)
all_k <- all_k[str_length(all_k) == 6]
k_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(all_k), width = 6), strand = "*", SEQUENCE = all_k)

# Reverse k-mers
synx_rev <- toupper(reverseComplement(ONTDNA_pileup$REF_NT))
rev_k <- runner(x = synx_rev, k = 6, f = KD_k)
rev_k <- rev_k[str_length(rev_k) == 6]

# Calculate kmer GC content
k_granges$GC_pct <- unlist(lapply(k_granges$SEQUENCE, GC_k))/6*100

# Calculate longest Homopolymer string
HP_score <- data.frame(row.names = 1:length(k_granges$SEQUENCE))
for(kn in 2:6){
  matches <- paste(c(paste(rep('A', kn), collapse = ''), paste(rep('T', kn), collapse = ''), paste(rep('G', kn), collapse = ''), paste(rep('C', kn), collapse = '')), collapse = "|")
  HP_score[,paste0('HP_k', kn)] <- unlist(lapply(k_granges$SEQUENCE, function(x){grepl(matches, x)}))
  }
k_granges$HP_score <- rowSums(HP_score) + 1

# Calculate mean error rate across kmers
k_granges$Error_fraction <- sum(extractList(ONTDNA_pileup$Error_counts, k_granges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, k_granges@ranges))
k_granges$SUB_fraction <- sum(extractList(ONTDNA_pileup$SUBS_counts, k_granges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, k_granges@ranges))
k_granges$INDEL_fraction <- sum(extractList(ONTDNA_pileup$INDEL_counts, k_granges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, k_granges@ranges))
k_granges$INS_fraction <- sum(extractList(ONTDNA_pileup$Ins, k_granges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, k_granges@ranges))
k_granges$DEL_fraction <- sum(extractList(ONTDNA_pileup$Del, k_granges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, k_granges@ranges))

k_granges$Error_mean <- mean(extractList(ONTDNA_pileup$Error_pct, k_granges@ranges))
k_granges$SUB_mean <- mean(extractList(ONTDNA_pileup$SUBS_pct, k_granges@ranges))
k_granges$INDEL_mean <- mean(extractList(ONTDNA_pileup$INDEL_pct, k_granges@ranges))

for(quant in c("Error_fraction", "SUB_fraction", "INDEL_fraction", "INS_fraction", "DEL_fraction", "GC_pct", "HP_score")){
  tmp <- data.frame(k_granges@elementMetadata) %>%
    group_by(SEQUENCE) %>%
    summarise('Error' = mean(!!as.name(quant))) %>%
    arrange(desc(Error))
  tmp$SEQUENCE <- factor(tmp$SEQUENCE, levels = tmp$SEQUENCE)
  ggplot(tmp) +
    geom_point(aes_string(x = 'SEQUENCE', y = 'Error')) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggsave(paste0("project_results/ONT_DNA/SynX_kmer_", quant, ".pdf"))
}

for(quant in c("Error_fraction", "SUB_fraction", "INDEL_fraction", "INS_fraction", "DEL_fraction")){
  ggplot(data.frame(k_granges@elementMetadata)) +
    geom_point(aes_string(x = 'GC_pct', y = quant)) +
    geom_spline(aes_string(x = 'GC_pct', y = quant), colour = 'red') +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/SynX_kmer_GCpct_", quant, ".pdf"))
  
  ggplot(data.frame(k_granges@elementMetadata)) +
    geom_point(aes_string(x = 'HP_score', y = quant)) +
    geom_spline(aes_string(x = 'HP_score', y = quant), colour = 'red') +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/SynX_kmer_HP_score_", quant, ".pdf"))
  }

# Visualisation of performance units
# --------------------------------------------------------------------------
options(ucscChromosomeNames=FALSE)

# Make genome position track
axisTrack <- GenomeAxisTrack(add53 = TRUE)

# Make sequence data track
sTrack <- SequenceTrack("data/reference_files/synx.fa", genome = 'SynX', chromosome = 'SynX')

# Make coverage track
depthTrack <- DataTrack(range = "data/ONT_DNA/barcode06.sorted.bam", type = "l", name = "Coverage", window = -1, genome = 'SynX', chromosome = 'SynX', transformation = function(x){ log10(x) }, size=1, ylim = c(0, 4), yTicksAt = seq(0, 4))

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

# Make annotation track
aTrack <- AnnotationTrack(range = synx_ranges, shape = "box", Feature = 'name', collapse = FALSE, name = 'Features')
aProm_track <- AnnotationTrack(range = promoters, shape = "box", stacking = 'squish', fill = npg_cols[1], name = NULL)
aRest_track <- AnnotationTrack(range = restriction, shape = "box", stacking = 'squish', fill = npg_cols[2], name = NULL)
aPAT_track <- AnnotationTrack(range = polyA, shape = "box", collapse = FALSE, fill = npg_cols[3], name = NULL)
aPerfGC_track <- AnnotationTrack(range = performanceGC, feature = , showFeatureId = TRUE, shape = "box", collapse = FALSE, fill = npg_cols[4], name = NULL)
feature(aPerfGC_track) <- performanceGC$name
aPerfHP_track <- AnnotationTrack(range = performanceHP, shape = "box", collapse = FALSE, fill = npg_cols[5], name = NULL)
feature(aPerfHP_track) <- performanceHP$name
aQuant_track <- AnnotationTrack(range = quantitative_ladder, shape = "box", collapse = FALSE, fill = npg_cols[6], name = NULL)
aSize_track <- AnnotationTrack(range = size_ladder, shape = "box", collapse = FALSE, fill = npg_cols[7], name = NULL)

# Make Granges to plot single base resolution data tracks
pileup_ranges <- GRanges(seqnames = "SynX",
                         ranges = IRanges(start = ONTDNA_pileup$Number, end = ONTDNA_pileup$Number, width = 1),
                         strand = "*")

for(col in colnames(ONTDNA_pileup)){
  pileup_ranges@elementMetadata[,col] <- list(ONTDNA_pileup[,col])
}

# Make coverage track
covTrack <- DataTrack(pileup_ranges, data = 'Coverage', name = "Coverage", type = c("l"), size = 1, col = 'grey80')

# Make GC% tracks
gcPctTrack <- DataTrack(performanceGC, data = 'GC_pct', name = "GC_pct", type = c("gradient"), size = 1, gradient = hcl.colors(15, 'Blues', rev = TRUE))
gcscoreTrack <- DataTrack(pileup_ranges, data = 'GC_score', name = "GC_score", type = c("gradient"), size = 2, gradient = hcl.colors(15, 'Blues', rev = TRUE))
gcFeaturesTrack <- DataTrack(subsetByOverlaps(pileup_ranges, performanceGC), data = 'GC_score', name = "GC_score", type = c("gradient"), size = 2, gradient = hcl.colors(15, 'Blues', rev = TRUE), ylim = c(0, 6), yTicksAt = seq(0, 6))

# Make HP% track
hpscoreTrack <- DataTrack(pileup_ranges, data = 'HP_score', name = "HP_score", type = c("gradient"), size = 2, gradient = hcl.colors(15, 'Reds', rev = TRUE))
hpFeaturesTrack <- DataTrack(subsetByOverlaps(pileup_ranges, performanceHP), data = 'HP_score', name = "HP_score", type = c("gradient"), size = 2, gradient = hcl.colors(15, 'Reds', rev = TRUE), ylim = c(1, 6), yTicksAt = seq(1, 6))

# Make K-mer diversity track
kdscoreTrack <- DataTrack(pileup_ranges, data = 'KD_score', name = "KD_score", type = c("gradient"), size = 2, gradient = hcl.colors(15, 'Greens', rev = TRUE), ylim = c(1, 6), yTicksAt = seq(1, 6))

# Make data track showing substitution rate
subTrack <- DataTrack(pileup_ranges, data = 'SUBS_pct', name = "Subs", type = c("l"), size = 2, col = npg_cols[1])
insTrack <- DataTrack(pileup_ranges, data = 'Ins_pct', name = "Ins", type = c("l"), size = 2, col = npg_cols[2])
delTrack <- DataTrack(pileup_ranges, data = 'Del_pct', name = "Del", type = c("l"), size = 2, col = npg_cols[4])
indelTrack <- DataTrack(pileup_ranges, data = 'INDEL_pct', name = "InDel", type = c("l"), size = 2, col = npg_cols[2])
errorTrack <- DataTrack(pileup_ranges, data = 'Error_pct', name = "Error", type = c("l"), size = 2, col = npg_cols[3])

# Plot tracks across whole SynX
pdf("project_results/ONT_DNA/Synx_error_profile.pdf")
plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack))
dev.off()

# Plots for HP performance elements
for(feat in performanceHP$name) {
  tmp <- performanceHP[performanceHP$name == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/ONT_DNA/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

# Plots for GC performance elements
for(feat in performanceGC$name) {
  tmp <- performanceGC[performanceGC$name == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/ONT_DNA/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfGC_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

# Plots for quantification elements
for(feat in quantitative_ladder$Identifier) {
  tmp <- quantitative_ladder[quantitative_ladder$Identifier == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/ONT_DNA/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, aQuant_track, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

pdf(paste0("project_results/ONT_DNA/Synx_error_profile_polyA.pdf"))
plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = 5350, to = 5400)
dev.off()

plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = 5390, to = 5450)

# Identify 31mers in quantitative ladder for quantification using jellyfish
# --------------------------------------------------------------------------

cn4_sequence <- quantitative_ladder[quantitative_ladder$name == '4cn',]$SEQUENCE[1]
cn4_31mers <- runner(x = unlist(strsplit(cn4_sequence, split="")), k = 31, f = KD_k)
cn4_31mers <- cn4_31mers[str_length(cn4_31mers) == 31]
write.fasta(sequences = as.list(cn4_31mers), names = paste0('cn4_', 1:length(cn4_31mers)), as.string = TRUE, file.out = 'project_results/ONT_DNA/jellyfish/cn4_31mers.fasta')

cn3_sequence <- quantitative_ladder[quantitative_ladder$name == '3cn',]$SEQUENCE[1]
cn3_31mers <- runner(x = unlist(strsplit(cn3_sequence, split="")), k = 31, f = KD_k)
cn3_31mers <- cn3_31mers[str_length(cn3_31mers) == 31]
write.fasta(sequences = as.list(cn3_31mers), names = paste0('cn4_', 1:length(cn3_31mers)), as.string = TRUE, file.out = 'project_results/ONT_DNA/jellyfish/cn3_31mers.fasta')

cn2_sequence <- quantitative_ladder[quantitative_ladder$name == '2cn',]$SEQUENCE[1]
cn2_31mers <- runner(x = unlist(strsplit(cn2_sequence, split="")), k = 31, f = KD_k)
cn2_31mers <- cn2_31mers[str_length(cn2_31mers) == 31]
write.fasta(sequences = as.list(cn2_31mers), names = paste0('cn4_', 1:length(cn2_31mers)), as.string = TRUE, file.out = 'project_results/ONT_DNA/jellyfish/cn2_31mers.fasta')

cn1_sequence <- quantitative_ladder[quantitative_ladder$name == '1cn',]$SEQUENCE[1]
cn1_31mers <- runner(x = unlist(strsplit(cn1_sequence, split="")), k = 31, f = KD_k)
cn1_31mers <- cn1_31mers[str_length(cn1_31mers) == 31]
write.fasta(sequences = as.list(cn1_31mers), names = paste0('cn4_', 1:length(cn1_31mers)), as.string = TRUE, file.out = 'project_results/ONT_DNA/jellyfish/cn1_31mers.fasta')

cn1_counts <- read.table("project_results/ONT_DNA/jellyfish/cn1_31mers_counts.tsv", col.names = c("Sequence", "Counts"))
cn_stats <- data.frame("cn1" = fivenum(cn1_counts$Counts))

cn2_counts <- read.table("project_results/ONT_DNA/jellyfish/cn2_31mers_counts.tsv", col.names = c("Sequence", "Counts"))
cn_stats$cn2 <- fivenum(cn2_counts$Counts)

cn3_counts <- read.table("project_results/ONT_DNA/jellyfish/cn3_31mers_counts.tsv", col.names = c("Sequence", "Counts"))
cn_stats$cn3 <- fivenum(cn3_counts$Counts)

cn4_counts <- read.table("project_results/ONT_DNA/jellyfish/cn4_31mers_counts.tsv", col.names = c("Sequence", "Counts"))
cn_stats$cn4 <- fivenum(cn4_counts$Counts)

rownames(cn_stats) <- c("min", "low", "mid", "top", "max")
cn_stats <- t(cn_stats)
cn_stats <- data.frame(cn_stats)
cn_stats$Copy_number <- rownames(cn_stats)

cn_stats$cn <- 1:4

ggplot(cn_stats, aes(x = cn, y = mid)) +
         geom_line()


ggplot(cn_stats, aes(x=Copy_number, ymin = min, lower = low, middle = mid, upper = top, ymax = max)) +
  geom_boxplot(stat = "identity", aes(color = Copy_number)) +
  theme_classic() +
  xlab("Copy number") +
  ylab("Counts") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  ggsave("project_results/ONT_DNA/Quantitative_ladder_counts_boxplot.pdf")


ggplot() +
  geom_density(data = cn1_counts, aes( x = Counts), color = npg_cols[1], fill = npg_cols[1], alpha = 0.1) +
  geom_histogram(data = cn2_counts, aes( x = Counts), color = npg_cols[2], fill = npg_cols[2], alpha = 0.1) +
  geom_histogram(data = cn3_counts, aes( x = Counts), color = npg_cols[3], fill = npg_cols[3], alpha = 0.1) +
  geom_histogram(data = cn4_counts, aes( x = Counts), color = npg_cols[4], fill = npg_cols[4], alpha = 0.1) +
  theme_classic()
  theme_classic() +
  xlab("Copy number") +
  ylab("Counts") +
  scale_color_npg(palette = c("nrc"), alpha = 1) +
  ggsave("Quantitative_ladder_counts_boxplot.pdf")

 
# Homopolymer analysis
# --------------------------------------------------------------------------

# Import bed file of homopolymer sites
synx_hp_ranges <- import.bed("project_results/ONT_DNA/synx_homopolymers.bed", genome = 'SynX')

synx_hp_ranges$Error_fraction <- sum(extractList(ONTDNA_pileup$Error_counts, synx_hp_ranges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, synx_hp_ranges@ranges))
synx_hp_ranges$SUB_fraction <- sum(extractList(ONTDNA_pileup$SUBS_counts, synx_hp_ranges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, synx_hp_ranges@ranges))
synx_hp_ranges$INDEL_fraction <- sum(extractList(ONTDNA_pileup$INDEL_counts, synx_hp_ranges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, synx_hp_ranges@ranges))
synx_hp_ranges$INS_fraction <- sum(extractList(ONTDNA_pileup$Ins, synx_hp_ranges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, synx_hp_ranges@ranges))
synx_hp_ranges$DEL_fraction <- sum(extractList(ONTDNA_pileup$Del, synx_hp_ranges@ranges)) / sum(extractList(ONTDNA_pileup$Depth, synx_hp_ranges@ranges))

# Plot homopolymer frequency over performance element
ggplot(data.frame(table(synx_hp_ranges$score)), aes(x=Var1, y=log2(Freq))) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency (log2)") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_frequency.pdf")

# Plot error statistics over whole synx
st_emd <- data.frame(synx_hp_ranges@elementMetadata)
st_emd_mean <- ddply(st_emd, .(score), colwise(mean))
st_emd_sd <- ddply(st_emd, .(score), colwise(sd))

st_emd_Error <- data.frame("Length" = st_emd_mean$score, "mean" = st_emd_mean$Error_fraction, "sd" = st_emd_sd$Error_fraction, "upper" = st_emd_mean$Error_fraction + st_emd_sd$Error_fraction, "lower" = st_emd_mean$Error_fraction - st_emd_sd$Error_fraction)
ggplot(st_emd_Error, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Error (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_Error_fraction_total.pdf")

st_emd_DEL <- data.frame("Length" = st_emd_mean$score, "mean" = st_emd_mean$DEL_fraction, "sd" = st_emd_sd$DEL_fraction, "upper" = st_emd_mean$DEL_fraction + st_emd_sd$DEL_fraction, "lower" = st_emd_mean$DEL_fraction - st_emd_sd$DEL_fraction)
ggplot(st_emd_DEL, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Deletion (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_deletion_fraction_total.pdf")

st_emd_INS <- data.frame("Length" = st_emd_mean$score, "mean" = st_emd_mean$INS_fraction, "sd" = st_emd_sd$INS_fraction, "upper" = st_emd_mean$INS_fraction + st_emd_sd$INS_fraction, "lower" = st_emd_mean$INS_fraction - st_emd_sd$INS_fraction)
ggplot(st_emd_INS, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Insertion (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_Insertion_fraction_total.pdf")

st_emd_SUB <- data.frame("Length" = st_emd_mean$score, "mean" = st_emd_mean$SUB_fraction, "sd" = st_emd_sd$SUB_fraction, "upper" = st_emd_mean$SUB_fraction + st_emd_sd$SUB_fraction, "lower" = st_emd_mean$SUB_fraction - st_emd_sd$SUB_fraction)
ggplot(st_emd_SUB, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Substitution (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_synx_total_Substitution_fraction_total.pdf")


# Subset to just HP performance elements
performanceHP_ranges <- subsetByOverlaps(synx_hp_ranges, performanceHP)

# Plot homopolymer frequency over performance element
ggplot(data.frame(table(performanceHP_ranges$score)), aes(x=Var1, y=Freq)) +
  geom_line(group = 1) +
  xlab("Homopolymer (bp)") +
  ylab("Frequency") +
  theme_classic() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_frequency_total.pdf")

# Plot error statistics over performance elements
php_emd <- data.frame(performanceHP_ranges@elementMetadata)
php_emd_mean <- ddply(php_emd, .(score), colwise(mean))
php_emd_sd <- ddply(php_emd, .(score), colwise(sd))

php_emd_Error <- data.frame("Length" = php_emd_mean$score, "mean" = php_emd_mean$Error_fraction, "sd" = php_emd_sd$Error_fraction, "upper" = php_emd_mean$Error_fraction + php_emd_sd$Error_fraction, "lower" = php_emd_mean$Error_fraction - php_emd_sd$Error_fraction)
ggplot(php_emd_Error, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Error (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_Error_fraction_total.pdf")

php_emd_DEL <- data.frame("Length" = php_emd_mean$score, "mean" = php_emd_mean$DEL_fraction, "sd" = php_emd_sd$DEL_fraction, "upper" = php_emd_mean$DEL_fraction + php_emd_sd$DEL_fraction, "lower" = php_emd_mean$DEL_fraction - php_emd_sd$DEL_fraction)
ggplot(php_emd_DEL, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Deletion (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_deletion_fraction_total.pdf")

php_emd_INS <- data.frame("Length" = php_emd_mean$score, "mean" = php_emd_mean$INS_fraction, "sd" = php_emd_sd$INS_fraction, "upper" = php_emd_mean$INS_fraction + php_emd_sd$INS_fraction, "lower" = php_emd_mean$INS_fraction - php_emd_sd$INS_fraction)
ggplot(php_emd_INS, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Insertion (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_Insertion_fraction_total.pdf")

php_emd_SUB <- data.frame("Length" = php_emd_mean$score, "mean" = php_emd_mean$SUB_fraction, "sd" = php_emd_sd$SUB_fraction, "upper" = php_emd_mean$SUB_fraction + php_emd_sd$SUB_fraction, "lower" = php_emd_mean$SUB_fraction - php_emd_sd$SUB_fraction)
ggplot(php_emd_SUB, aes(x=Length, y=mean)) +
  geom_line(group = 1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
  xlab("Homopolymer (bp)") +
  ylab("Substitution (fraction)") +
  theme_minimal() +
  ggsave("project_results/ONT_DNA/Homopolymer_performance_element_Substitution_fraction_total.pdf")

# Make plots for Homopolymers of each base individually
for (base in c('A', 'T', 'G', 'C')) {
  ggplot(data.frame(table(performanceHP_ranges[performanceHP_ranges$name == base,]$score)), aes(x=Var1, y=Freq)) +
    geom_line(group = 1) +
    xlab("Homopolymer (bp)") +
    ylab("Frequency") +
    theme_classic() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_performance_element_frequency_", base, ".pdf"))
  
  base_emd <- data.frame(performanceHP_ranges@elementMetadata)
  base_emd <- base_emd[base_emd$name == base, ]
  base_emd_mean <- ddply(base_emd, .(score), colwise(mean))
  base_emd_sd <- ddply(base_emd, .(score), colwise(sd))
  base_emd_Error <- data.frame("Length" = base_emd_mean$score, "mean" = base_emd_mean$Error_fraction, "sd" = base_emd_sd$Error_fraction, "upper" = base_emd_mean$Error_fraction + base_emd_sd$Error_fraction, "lower" = base_emd_mean$Error_fraction - base_emd_sd$Error_fraction)
  ggplot(base_emd_Error, aes(x=Length, y=mean)) +
    geom_line(group = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
    xlab("Homopolymer (bp)") +
    ylab("Error (fraction)") +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_performance_element_Error_fraction_", base, "_nt.pdf"))
  
  base_emd_DEL <- data.frame("Length" = base_emd_mean$score, "mean" = base_emd_mean$DEL_fraction, "sd" = base_emd_sd$DEL_fraction, "upper" = base_emd_mean$DEL_fraction + base_emd_sd$DEL_fraction, "lower" = base_emd_mean$DEL_fraction - base_emd_sd$DEL_fraction)
  ggplot(base_emd_DEL, aes(x=Length, y=mean)) +
    geom_line(group = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
    xlab("Homopolymer (bp)") +
    ylab("Deletion (fraction)") +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_performance_element_deletion_fraction_", base, "_nt.pdf"))
  
  base_emd_INS <- data.frame("Length" = base_emd_mean$score, "mean" = base_emd_mean$INS_fraction, "sd" = base_emd_sd$INS_fraction, "upper" = base_emd_mean$INS_fraction + base_emd_sd$INS_fraction, "lower" = base_emd_mean$INS_fraction - base_emd_sd$INS_fraction)
  ggplot(base_emd_INS, aes(x=Length, y=mean)) +
    geom_line(group = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
    xlab("Homopolymer (bp)") +
    ylab("Insertion (fraction)") +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_performance_element_Insertion_fraction_", base, "_nt.pdf"))
  
  base_emd_SUB <- data.frame("Length" = base_emd_mean$score, "mean" = base_emd_mean$SUB_fraction, "sd" = base_emd_sd$SUB_fraction, "upper" = base_emd_mean$SUB_fraction + base_emd_sd$SUB_fraction, "lower" = base_emd_mean$SUB_fraction - base_emd_sd$SUB_fraction)
  ggplot(base_emd_SUB, aes(x=Length, y=mean)) +
    geom_line(group = 1) +
    geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue", alpha=0.2) +
    xlab("Homopolymer (bp)") +
    ylab("Substitution (fraction)") +
    theme_minimal() +
    ggsave(paste0("project_results/ONT_DNA/Homopolymer_performance_element_Substitution_fraction_", base, "_nt.pdf"))
}





