#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make Gviz alignment across synX elements Illumina
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
IlluminaDNA_pileup <- read.table("project_results/Illumina_DNASynX/Illumina_DNASynX.bam.bed.tsv", header = TRUE, sep="\t", stringsAsFactors=FALSE, quote="")

# Find reference mapping rates
get_correct <- function(df){
  bases <- c('A', 'T', 'G', 'C')
  ref <- df['REF_NT']
  correct <- bases[bases == ref]
  sum(as.numeric(df[correct]))
}

IlluminaDNA_pileup$Depth <- IlluminaDNA_pileup$Ins + IlluminaDNA_pileup$Del + IlluminaDNA_pileup$Coverage
IlluminaDNA_pileup$REF_counts <- apply(IlluminaDNA_pileup, 1, get_correct)
IlluminaDNA_pileup$REF_pct <- IlluminaDNA_pileup$REF_counts / IlluminaDNA_pileup$Depth

# Find substitution rates
get_subs <- function(df){
  bases <- c('A', 'T', 'G', 'C')
  ref <- df['REF_NT']
  subs <- bases[!(bases == ref)]
  sum(as.numeric(df[subs]))
}

IlluminaDNA_pileup$SUBS_counts <- apply(IlluminaDNA_pileup, 1, get_subs)
IlluminaDNA_pileup$SUBS_pct <- IlluminaDNA_pileup$SUBS_counts / IlluminaDNA_pileup$Depth

# Find INDEL rates
IlluminaDNA_pileup$Ins_pct <- IlluminaDNA_pileup$Ins / IlluminaDNA_pileup$Depth
IlluminaDNA_pileup$Del_pct <- IlluminaDNA_pileup$Del / IlluminaDNA_pileup$Depth
IlluminaDNA_pileup$INDEL_counts <- IlluminaDNA_pileup$Ins + IlluminaDNA_pileup$Del
IlluminaDNA_pileup$INDEL_pct <- IlluminaDNA_pileup$INDEL_counts / IlluminaDNA_pileup$Depth

# Find error rates
IlluminaDNA_pileup$Error_counts <- IlluminaDNA_pileup$SUBS_counts + IlluminaDNA_pileup$INDEL_counts
IlluminaDNA_pileup$Error_pct <- IlluminaDNA_pileup$Error_counts / IlluminaDNA_pileup$Depth

# Caluculate GCscore for each base
GC_k <- function(x){str_count(paste(x, collapse = ""), pattern = 'C') + str_count(paste(x, collapse = ""), pattern = 'G')}
GC_score <- data.frame(row.names = 1:length(IlluminaDNA_pileup$REF_NT))
for(kn in 0:5){
  GC_score[,paste0('GC_k', kn)]<- runner(x = IlluminaDNA_pileup$REF_NT, k = 6, f = GC_k, lag = -kn)
}
IlluminaDNA_pileup$GC_score <- rowMeans(as.matrix(GC_score))

# Caluculate homopolymer score for each base
homo_score <- function(x){length(unique(x)) == 1}
HP_score <- data.frame(row.names = 1:length(IlluminaDNA_pileup$REF_NT))
for(kn in 2:6){
  tmp1 <- data.frame(row.names = 1:length(IlluminaDNA_pileup$REF_NT))
  for(k in 0:(kn-1)){
    tmp1[,paste0('HP_k', kn, '_lag_', k)] <- runner(x = IlluminaDNA_pileup$REF_NT, k = kn, f = homo_score, lag = -k)
  }
  HP_score[,paste0('HP_k', kn)] <- rowSums(tmp1) > 0
}
IlluminaDNA_pileup$HP_score <- rowSums(HP_score) + 1

# Caluculate k-mer diversity score for each base
KD_k <- function(x){paste(x, collapse = "")}
KD_score <- data.frame(row.names = 1:length(IlluminaDNA_pileup$REF_NT))
for(kn in 0:5){
  KD_score[,paste0('KD_k', kn)]<- runner(x = IlluminaDNA_pileup$REF_NT, k = 6, f = KD_k, lag = -kn)
}

div_k <- function(x){length(unique(x))}
IlluminaDNA_pileup$KD_score <- (6 - apply(X = KD_score, 1, FUN = div_k))

# Plot error rates across features
# --------------------------------------------------------------------------

synx_ranges$Error_fraction <- sum(relist(IlluminaDNA_pileup$Error_counts, synx_ranges@ranges)) / sum(relist(IlluminaDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$SUB_fraction <- sum(relist(IlluminaDNA_pileup$SUBS_counts, synx_ranges@ranges)) / sum(relist(IlluminaDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$INDEL_fraction <- sum(relist(IlluminaDNA_pileup$INDEL_counts, synx_ranges@ranges)) / sum(relist(IlluminaDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$INS_fraction <- sum(relist(IlluminaDNA_pileup$Ins, synx_ranges@ranges)) / sum(relist(IlluminaDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$DEL_fraction <- sum(relist(IlluminaDNA_pileup$Del, synx_ranges@ranges)) / sum(relist(IlluminaDNA_pileup$Depth, synx_ranges@ranges))
synx_ranges$HP_score <- mean(relist(IlluminaDNA_pileup$HP_score, synx_ranges@ranges))

# Plot error rates across kmers
# --------------------------------------------------------------------------

# List all kmers in SynX genome and make Granges object
all_k <- runner(x = IlluminaDNA_pileup$REF_NT, k = 6, f = KD_k)
all_k <- all_k[str_length(all_k) == 6]
k_granges <- GRanges(seqnames = 'SynX', ranges = IRanges(start = 1:length(all_k), width = 6), strand = "*", SEQUENCE = all_k)

# Reverse k-mers
synx_rev <- toupper(reverseComplement(IlluminaDNA_pileup$REF_NT))
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
k_granges$Error_fraction <- sum(extractList(IlluminaDNA_pileup$Error_counts, k_granges@ranges)) / sum(extractList(IlluminaDNA_pileup$Depth, k_granges@ranges))
k_granges$SUB_fraction <- sum(extractList(IlluminaDNA_pileup$SUBS_counts, k_granges@ranges)) / sum(extractList(IlluminaDNA_pileup$Depth, k_granges@ranges))
k_granges$INDEL_fraction <- sum(extractList(IlluminaDNA_pileup$INDEL_counts, k_granges@ranges)) / sum(extractList(IlluminaDNA_pileup$Depth, k_granges@ranges))
k_granges$INS_fraction <- sum(extractList(IlluminaDNA_pileup$Ins, k_granges@ranges)) / sum(extractList(IlluminaDNA_pileup$Depth, k_granges@ranges))
k_granges$DEL_fraction <- sum(extractList(IlluminaDNA_pileup$Del, k_granges@ranges)) / sum(extractList(IlluminaDNA_pileup$Depth, k_granges@ranges))

k_granges$Error_mean <- mean(extractList(IlluminaDNA_pileup$Error_pct, k_granges@ranges))
k_granges$SUB_mean <- mean(extractList(IlluminaDNA_pileup$SUBS_pct, k_granges@ranges))
k_granges$INDEL_mean <- mean(extractList(IlluminaDNA_pileup$INDEL_pct, k_granges@ranges))

# Visualisation of performance units
# --------------------------------------------------------------------------
options(ucscChromosomeNames=FALSE)

# Make genome position track
axisTrack <- GenomeAxisTrack(add53 = TRUE)

# Make sequence data track
sTrack <- SequenceTrack("data/reference_files/synx.fa", genome = 'SynX', chromosome = 'SynX')

# Make coverage track
depthTrack <- DataTrack(range = "project_results/Illumina_DNASynX/Illumina_DNASynX_1mil.sorted.bam", type = "l", name = "Coverage", window = -1, genome = 'SynX', chromosome = 'SynX', transformation = function(x){ log10(x) }, size=1, ylim = c(0, 4), yTicksAt = seq(0, 4))

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
                         ranges = IRanges(start = IlluminaDNA_pileup$Number, end = IlluminaDNA_pileup$Number, width = 1),
                         strand = "*")

for(col in colnames(IlluminaDNA_pileup)){
  pileup_ranges@elementMetadata[,col] <- list(IlluminaDNA_pileup[,col])
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
pdf("project_results/Illumina_DNASynX/Synx_error_profile.pdf")
plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack))
dev.off()

# Plots for HP performance elements
for(feat in performanceHP$name) {
  tmp <- performanceHP[performanceHP$name == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/Illumina_DNASynX/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

# Plots for GC performance elements
for(feat in performanceGC$name) {
  tmp <- performanceGC[performanceGC$name == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/Illumina_DNASynX/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfGC_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

# Plots for quantification elements
for(feat in quantitative_ladder$Identifier) {
  tmp <- quantitative_ladder[quantitative_ladder$Identifier == feat, ]
  startbp <- tmp@ranges@start - 50
  endbp <- tmp@ranges@start + tmp@ranges@width + 50
  pdf(paste0("project_results/Illumina_DNASynX/Synx_error_profile_feature_", feat, ".pdf"))
  plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, aQuant_track, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = startbp, to = endbp)
  dev.off()
}

pdf(paste0("project_results/Illumina_DNASynX/Synx_error_profile_polyA.pdf"))
plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = 5350, to = 5400)
dev.off()

plotTracks(list(axisTrack, aTrack, errorTrack, subTrack, insTrack, delTrack, hpFeaturesTrack, aPerfHP_track, gcFeaturesTrack, gcPctTrack, gcscoreTrack, hpscoreTrack, kdscoreTrack, depthTrack, covTrack, sTrack), from = 5390, to = 5450)
