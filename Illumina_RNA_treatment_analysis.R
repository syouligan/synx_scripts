#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX in Illumina cancer samples
# --------------------------------------------------------------------------

# Set adaptive working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/tim_projects/synx/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/tim/mer/scott/synx/")
  place <- "timmer"
}

# Libraries
library(edgeR)
library(limma)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
library(ggsci)
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(subSeq)

# Prepare data for analysis
# --------------------------------------------------------------------------

# Load GTF of feature annotations
GENCODE_gtf <- import.gff("data/reference_files/gencode.v36.primary_assembly.annotation.gtf")
synx_gtf <- import.gff("data/reference_files/SynX_cn_masked_transcripts.gtf")

feature_annot <- data.frame(c(GENCODE_gtf, synx_gtf))
feature_annot <- feature_annot[,c("seqnames", "gene_id", "gene_type", "gene_name")]
feature_annot <- feature_annot %>% distinct()

# Generate transcript to gene database for gencode
gencodetx2gene <- data.frame(GENCODE_gtf@elementMetadata[GENCODE_gtf$type == 'transcript', c("transcript_id", "gene_id")])
colnames(gencodetx2gene) <- c("TXNAME", "GENEID")
gencodetx2gene <- gencodetx2gene %>%
  distinct()
gencodetx2gene <- as_tibble(gencodetx2gene)

# Generate transcript to gene database for synx_cn_masked
synxtx2gene <- data.frame(synx_gtf@elementMetadata[, c("gene_name", "gene_id")])
colnames(synxtx2gene) <- c("TXNAME", "GENEID")
synxtx2gene <- synxtx2gene %>%
  distinct()
synxtx2gene <- as_tibble(synxtx2gene)

# Make combined tx2gene object

tx2gene <- rbind(gencodetx2gene, synxtx2gene)

# List files
samples <- c("Illumina_T7untreated1", "Illumina_T7untreated2", "Illumina_T7torkinib1", "Illumina_T7torkinib2")
files <- file.path("project_results", samples, "abundance.tsv")
names(files) <- samples
all(file.exists(files))

# Import counts from Kallisto
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# Remove unaligned read info and unexpressed features from gene counts and gene annotation
gene_counts <- round(txi.kallisto.tsv$counts)
gene_counts <- gene_counts[! rowSums(gene_counts == 0) == ncol(gene_counts),]

feature_annot <- feature_annot[feature_annot$gene_id %in% rownames(gene_counts),]
rownames(feature_annot) <- feature_annot$gene_id
feature_annot <- feature_annot[rownames(gene_counts), ]
feature_annot_subseq <- feature_annot #store for later

# Split gene counts into GENCODE and SynX genes
synx_counts <- gene_counts[feature_annot$seqnames == "SynX",]
gencode_counts <- gene_counts[feature_annot$seqnames != "SynX",]

# Scale SynX counts to 1% of gencode library size (NOTE: something weird going on with ratio of SynX to GENCODE counts)
synx_size <- colSums(synx_counts)
gencode_size <- colSums(gencode_counts)
lib_size <- colSums(gene_counts)
scaling_pct <- (lib_size*0.01)/synx_size
synx_counts_scaled <- data.frame(round(t(t(synx_counts)*scaling_pct)))
gene_counts_scaled <- rbind(gencode_counts, synx_counts_scaled)
feature_annot <- feature_annot[rownames(gene_counts_scaled), ]

# Add copy number inversions for SynX conu genes
cn_counts <- synx_counts_scaled[grepl("cn", rownames(synx_counts_scaled)),]
cn_counts_untreated <- cn_counts[ ,grepl("untreated", colnames(cn_counts))]
cn_counts_torkinib <- cn_counts[ ,grepl("torkinib", colnames(cn_counts))]

cn_counts_1 <- cbind(cn_counts_untreated[1,], cn_counts_torkinib[-1, ])
rownames(cn_counts_1) <- c("2cnVS1cn", "3cnVS1cn", "4cnVS1cn")
cn_counts_2 <- cbind(cn_counts_untreated[-1,], cn_counts_torkinib[1, ])
rownames(cn_counts_2) <- c("1cnVS2cn", "1cnVS3cn", "1cnVS4cn")

cn_inverted <- rbind(cn_counts_1, cn_counts_2)

# Add to features
cn_inverted_feats <- data.frame("seqnames" = 'conu_inverted', "gene_id" = rownames(cn_inverted), gene_type = "Inverted_cn_ladder", gene_name = rownames(cn_inverted), row.names = rownames(cn_inverted))
feature_annot <- rbind(feature_annot, cn_inverted_feats)

feature_annot <- feature_annot %>%
  mutate('Feature' = case_when(seqnames == 'conu_inverted' ~ 'Copy number',
            seqnames == 'SynX' ~ 'SynX',
            TRUE ~ 'Gencode'))

# Add onverted copy number genes to scaled counts
gene_counts_scaled <- rbind(gene_counts_scaled, cn_inverted)

# Create sample design matrix
sampleName <- gsub("Illumina_T7", "", colnames(gene_counts_scaled))
sampleType <- factor(c(gsub( "[0-9]$", "", sampleName)))
sampleRep <- factor(gsub("[^0-9]", "", sampleName))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Find differentially expressed genes (without TMM) normalisation
# --------------------------------------------------------------------------

# Plot dataset before normalisation
set <- newSeqExpressionSet(as.matrix(gene_counts_scaled), phenoData = data.frame(sampleType, row.names=colnames(gene_counts_scaled)))

# Plot unnormalised date
pdf("project_results/Illumina_RNAseq/RLE_unnormalised.pdf")
plotRLE(set, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_unnormalised.pdf")
plotPCA(set, cex=1)
dev.off()

# Plot Upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")

pdf("project_results/Illumina_RNAseq/RLE_Upper_Quartile_normalised.pdf")
plotRLE(set, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_Upper_Quartile_normalised.pdf")
plotPCA(set, cex=1)
dev.off()

# Find differentially expressed genes with TMM normalisation
design <- model.matrix(~sampleType, data=pData(set))
y <- DGEList(counts=counts(set), group=sampleType, genes = feature_annot)
keep <- rowSums(counts(set)[,1:2] > 10) == 2 | rowSums(counts(set)[,3:4] > 10) == 2
y <- y[keep,]
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
res <- data.frame(topTags(lrt, n = nrow(counts(set))))

# Volcano plot of DEGs after TMM
top20 <- res[1:10, "gene_name"]
res$Top20 <- res$gene_name %in% top20
res$synx <- ! grepl("^ENS", rownames(res))
res$Copy_number <- grepl("cnVS", rownames(res))
res <- res[order(res$synx),]
res$Significant <- abs(res$logFC) > 3 & res$FDR < 0.05

ggplot(res, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Feature), alpha = 0.3) +
  geom_text_repel(aes(label = ifelse(Top20, as.character(gene_name), ""))) +
  geom_text_repel(aes(label = ifelse(Copy_number, as.character(gene_name), ""))) +
  theme_classic() +
  scale_color_manual(values = c('Gencode' = "grey50", 'SynX' = "Blue", 'Copy number' = "Red")) +
  xlab("Fold change (log2)") +
  ylab("FDR (-log10)") +
  xlim(-max(res$logFC), max(res$logFC)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/Illumina_RNAseq/Treated_untreated_volcano_plot_TMM.pdf")

# Find differentially expressed genes (with RUVg) normalisation
# --------------------------------------------------------------------------

# SynX control features for normalisation
spikes <- rownames(gene_counts_scaled)[!grepl("^ENS", rownames(gene_counts_scaled))]
spikes <- spikes[!grepl("cn", spikes)]
spikes <- spikes[grepl("P[0-9]", spikes)]

# Perform RUVg unwanted variation removal
set1 <- RUVg(set, spikes, k=1)

# Plot TMM normalisation
pdf("project_results/Illumina_RNAseq/RLE_RUVg_normalised.pdf")
plotRLE(set1, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_RUVg_normalised.pdf")
plotPCA(set1, cex=1)
dev.off()

# Find differentially expressed genes with RUVg normalisation
design <- model.matrix(~sampleType + W_1, data=pData(set1))
y1 <- DGEList(counts=counts(set), group=sampleType, genes = feature_annot)
keep <- rowSums(counts(set)[,1:2] > 10) == 2 | rowSums(counts(set)[,3:4] > 10) == 2
y1 <- y1[keep,]
y1 <- calcNormFactors(y1, method="TMM")
y1 <- estimateGLMCommonDisp(y1, design)
y1 <- estimateGLMTagwiseDisp(y1, design)
fit1 <- glmFit(y1, design)
lrt1 <- glmLRT(fit1, coef=2)
topTags(lrt1)
res1 <- data.frame(topTags(lrt1, n = nrow(counts(set1))))

# Volcano plot of DEGs after RUVg
top20 <- res1[1:10, "gene_name"]
res1$Top20 <- res1$gene_name %in% top20
res1$synx <- ! grepl("^ENS", rownames(res1))
res1$Copy_number <- grepl("cnVS", rownames(res1))
res1 <- res1[order(res1$synx),]
res1$Significant <- abs(res1$logFC) > 3 & res1$FDR < 0.05

ggplot(res1, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Feature), alpha = 0.3) +
  geom_text_repel(aes(label = ifelse(Top20, as.character(gene_name), ""))) +
  geom_text_repel(aes(label = ifelse(Copy_number, as.character(gene_name), ""))) +
  theme_classic() +
  scale_color_manual(values = c('Gencode' = "grey50", 'SynX' = "Blue", 'Copy number' = "Red")) +
  xlab("Fold change (log2)") +
  ylab("FDR (-log10)") +
  xlim(-max(res1$logFC), max(res1$logFC)) +
  theme_classic() +
  theme(aspect.ratio = 1)  +
  ggsave("project_results/Illumina_RNAseq/Treated_untreated_volcano_plot_RUVg.pdf")

# Create object with gene counts for Torkinib subsampled to 10%
# --------------------------------------------------------------------------

# Subsample Torkinib samples
gene_counts_subseq <- cbind(gene_counts[,1:2], generateSubsampledMatrix(gene_counts[,3:4], 0.1, 42))

# Split gene counts into GENCODE and SynX genes
synx_counts_subseq <- gene_counts_subseq[feature_annot_subseq$seqnames == "SynX",]
gencode_counts_subseq <- gene_counts_subseq[feature_annot_subseq$seqnames != "SynX",]

# Scale SynX counts to 1% of gencode library size (NOTE: something weird going on with ratio of SynX to GENCODE counts)
synx_size <- colSums(synx_counts_subseq)
gencode_size <- colSums(gencode_counts_subseq)
lib_size <- colSums(gene_counts_subseq)
scaling_pct <- (lib_size*0.01)/synx_size
synx_counts_subseq_scaled <- data.frame(round(t(t(synx_counts_subseq)*scaling_pct)))
gene_counts_subseq_scaled <- rbind(gencode_counts_subseq, synx_counts_subseq_scaled)
feature_annot_subseq <- feature_annot_subseq[rownames(gene_counts_subseq_scaled), ]

# Add copy number inversions for SynX conu genes
cn_counts_subseq <- synx_counts_subseq_scaled[grepl("cn", rownames(synx_counts_subseq_scaled)),]
cn_counts_subseq_untreated <- cn_counts_subseq[ ,grepl("untreated", colnames(cn_counts_subseq))]
cn_counts_subseq_torkinib <- cn_counts_subseq[ ,grepl("torkinib", colnames(cn_counts_subseq))]

cn_counts_subseq_1 <- cbind(cn_counts_subseq_untreated[1,], cn_counts_subseq_torkinib[-1, ])
rownames(cn_counts_subseq_1) <- c("2cnVS1cn", "3cnVS1cn", "4cnVS1cn")
cn_counts_subseq_2 <- cbind(cn_counts_subseq_untreated[-1,], cn_counts_subseq_torkinib[1, ])
rownames(cn_counts_subseq_2) <- c("1cnVS2cn", "1cnVS3cn", "1cnVS4cn")

cn_inverted <- rbind(cn_counts_subseq_1, cn_counts_subseq_2)

# Add to features
cn_inverted_feats <- data.frame("seqnames" = 'conu_inverted', "gene_id" = rownames(cn_inverted), gene_type = "Inverted_cn_ladder", gene_name = rownames(cn_inverted), row.names = rownames(cn_inverted))
feature_annot_subseq <- rbind(feature_annot_subseq, cn_inverted_feats)

feature_annot_subseq <- feature_annot_subseq %>%
  mutate('Feature' = case_when(seqnames == 'conu_inverted' ~ 'Copy number',
                               seqnames == 'SynX' ~ 'SynX',
                               TRUE ~ 'Gencode'))

# Add onverted copy number genes to scaled counts
gene_counts_subseq_scaled <- rbind(gene_counts_subseq_scaled, cn_inverted)

# Create sample design matrix
sampleName <- gsub("Illumina_T7", "", colnames(gene_counts_subseq_scaled))
sampleType <- factor(c(gsub( "[0-9]$", "", sampleName)))
sampleRep <- factor(gsub("[^0-9]", "", sampleName))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Find differentially expressed genes (without TMM) normalisation
# --------------------------------------------------------------------------

# Plot dataset before normalisation
set <- newSeqExpressionSet(as.matrix(gene_counts_subseq_scaled), phenoData = data.frame(sampleType, row.names=colnames(gene_counts_subseq_scaled)))

# Plot unnormalised data
pdf("project_results/Illumina_RNAseq/RLE_unnormalised_subseq.pdf")
plotRLE(set, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_unnormalised_subseq.pdf")
plotPCA(set, cex=1)
dev.off()

# Plot Upper quartile normalisation
set <- betweenLaneNormalization(set, which="upper")

pdf("project_results/Illumina_RNAseq/RLE_Upper_Quartile_normalised_subseq.pdf")
plotRLE(set, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_Upper_Quartile_normalised_subseq.pdf")
plotPCA(set)
dev.off()

# Find differentially expressed genes with TMM normalisation
design <- model.matrix(~sampleType, data=pData(set))
y <- DGEList(counts=counts(set), group=sampleType, genes = feature_annot_subseq)
keep <- rowSums(counts(set)[,1:2] > 10) == 2 | rowSums(counts(set)[,3:4] > 10) == 2
y <- y[keep,]
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
res <- data.frame(topTags(lrt, n = nrow(counts(set))))

# Volcano plot of DEGs after TMM
top20 <- res[1:10, "gene_name"]
res$Top20 <- res$gene_name %in% top20
res$synx <- ! grepl("^ENS", rownames(res))
res$Copy_number <- grepl("cnVS", rownames(res))
res <- res[order(res$synx),]
res$Significant <- abs(res$logFC) > 3 & res$FDR < 0.05

ggplot(res, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Feature), alpha = 0.3) +
  geom_text_repel(aes(label = ifelse(Top20, as.character(gene_name), ""))) +
  geom_text_repel(aes(label = ifelse(Copy_number, as.character(gene_name), ""))) +
  theme_classic() +
  scale_color_manual(values = c('Gencode' = "grey50", 'SynX' = "Blue", 'Copy number' = "Red")) +
  xlab("Fold change (log2)") +
  ylab("FDR (-log10)") +
  xlim(-max(res$logFC), max(res$logFC)) +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/Illumina_RNAseq/Treated_untreated_volcano_plot_TMM_subseq.pdf")

# Find differentially expressed genes (with RUVg) normalisation
# --------------------------------------------------------------------------

# SynX control features for normalisation
spikes <- rownames(gene_counts_subseq_scaled)[!grepl("^ENS", rownames(gene_counts_subseq_scaled))]
spikes <- spikes[!grepl("cn", spikes)]
spikes <- spikes[grepl("P[0-9]", spikes)]

# Perform RUVg unwanted variation removal
set1 <- RUVg(set, spikes, k=1)

# Plot RUVg normalisation
pdf("project_results/Illumina_RNAseq/RLE_RUVg_normalised_subseq.pdf")
plotRLE(set1, outline=FALSE)
dev.off()

pdf("project_results/Illumina_RNAseq/PCA_RUVg_normalised_subseq.pdf")
plotPCA(set1, cex=1)
dev.off()

# Find differentially expressed genes with RUVg normalisation
design <- model.matrix(~sampleType + W_1, data=pData(set1))
y1 <- DGEList(counts=counts(set), group=sampleType, genes = feature_annot_subseq)
keep <- rowSums(counts(set)[,1:2] > 10) == 2 | rowSums(counts(set)[,3:4] > 10) == 2
y1 <- y1[keep,]
y1 <- calcNormFactors(y1, method="TMM")
y1 <- estimateGLMCommonDisp(y1, design)
y1 <- estimateGLMTagwiseDisp(y1, design)
fit1 <- glmFit(y1, design)
lrt1 <- glmLRT(fit1, coef=2)
topTags(lrt1)
res1 <- data.frame(topTags(lrt1, n = nrow(counts(set1))))

# Volcano plot of DEGs after RUVg
top20 <- res1[1:10, "gene_name"]
res1$Top20 <- res1$gene_name %in% top20
res1$synx <- ! grepl("^ENS", rownames(res1))
res1$Copy_number <- grepl("cnVS", rownames(res1))
res1 <- res1[order(res1$synx),]
res1$Significant <- abs(res1$logFC) > 3 & res1$FDR < 0.05

ggplot(res1, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = Feature), alpha = 0.3) +
  geom_text_repel(aes(label = ifelse(Top20, as.character(gene_name), ""))) +
  geom_text_repel(aes(label = ifelse(Copy_number, as.character(gene_name), ""))) +
  theme_classic() +
  scale_color_manual(values = c('Gencode' = "grey50", 'SynX' = "Blue", 'Copy number' = "Red")) +
  xlab("Fold change (log2)") +
  ylab("FDR (-log10)") +
  xlim(-max(res1$logFC), max(res1$logFC)) +
  theme_classic() +
  theme(aspect.ratio = 1)  +
  ggsave("project_results/Illumina_RNAseq/Treated_untreated_volcano_plot_RUVg_subseq.pdf")


