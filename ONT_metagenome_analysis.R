#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Analysis of synX in ONT MetaSequin samples
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
library(reshape2)

samples <- c("ONT_SeqMetGenA1", "ONT_SeqMetGenA2", "ONT_SeqMetGenA3", "ONT_SeqMetGenB1", "ONT_SeqMetGenB2", "ONT_SeqMetGenB3")

# Prepare data for analysis
# --------------------------------------------------------------------------

# Load sequin feature concentrations
sequin_conc <- read.table("data/reference_files/metasequin_features_corrected.tsv", header = TRUE)

# Make dataframe of HTSeq counts
rm(gene_counts)
for (samp in samples) {
  tmp <-  read.table(paste0("project_results/", samp,"/", samp, ".MetaSequin_SynX.coverage"), header = FALSE)
  
  tmp <- tmp[, c(4, 7), drop = FALSE]
  colnames(tmp) <- c("Features", samp)
  tmp <- data.frame(tapply(tmp[,2], tmp[,1], FUN=sum))
  colnames(tmp) <- samp
  
  if(exists("gene_counts")) {
    gene_counts <- cbind(gene_counts, tmp)
  } else {
    gene_counts <- tmp
  }
}

# Plot observed vs expected concentration
# --------------------------------------------------------------------------

# Subset to dynamic sequins
dyn_seq <- sequin_conc[!is.na(sequin_conc$MIX_A), "ID"]

# Split gene counts into GENCODE and SynX genes
sequin_counts <- gene_counts[rownames(gene_counts) %in% dyn_seq,]

# Split into mix A and mix B concentrations
sequinMixA_counts <- sequin_counts[,grepl("GenA", colnames(sequin_counts))]
sequinMixA_counts$Features <- rownames(sequinMixA_counts)
sequinMixB_counts <- sequin_counts[,grepl("GenB", colnames(sequin_counts))]
sequinMixB_counts$Features <- rownames(sequinMixB_counts)

# Annotate with known concentration
sequinMixA <- melt(sequinMixA_counts, id.vars = 'Features')
colnames(sequinMixA) <- c("Features", "Sample", "Measured")

idx <- match(sequinMixA$Features, sequin_conc$ID)
sequinMixA$Expected <- sequin_conc$MIX_A [idx]

# Calculate mean and SD per feature
MixA <- sequinMixA %>%
  group_by(Features) %>%
  mutate("Mean" = mean(log10(Measured + 1)),
         "SD" = sd(log10(Measured + 1))) %>%
  select(Features, Expected, Mean, SD) %>%
  distinct()

# Calculate linear model of Observed vs Expected
mod <- lm(Mean ~ log10(Expected), MixA)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

# Plot observed vs Expected for MixA
ggplot(MixA, aes(x=log10(Expected), y=Mean)) + 
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD), shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  xlab("Expected conc (10^x)") +
  ylab("Coverage (10^x)") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/ONT_metagenome/MixA_quantitative_accuracy.pdf")

# Annotate with known concentration
sequinMixB <- melt(sequinMixB_counts, id.vars = 'Features')
colnames(sequinMixB) <- c("Features", "Sample", "Measured")

idx <- match(sequinMixB$Features, sequin_conc$ID)
sequinMixB$Expected <- sequin_conc$MIX_B [idx]

# Calculate mean and SD per feature
MixB <- sequinMixB %>%
  group_by(Features) %>%
  mutate("Mean" = mean(log10(Measured + 1)),
         "SD" = sd(log10(Measured + 1))) %>%
  select(Features, Expected, Mean, SD) %>%
  distinct()

# Calculate linear model of Observed vs Expected
mod <- lm(Mean ~ log10(Expected), MixB)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

# Plot observed vs Expected for MixB
ggplot(MixB, aes(x=log10(Expected), y=Mean)) + 
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_pointrange(aes(ymin=Mean-SD, ymax=Mean+SD), shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  xlab("Expected conc (10^x)") +
  ylab("Coverage (10^x)") +
  theme_classic() +
  theme(aspect.ratio = 1) +
  ggsave("project_results/ONT_metagenome/MixB_quantitative_accuracy.pdf")

# Create sample design matrix
sampleName <- gsub("ONT_SeqMetGen", "Mix", colnames(sequin_counts))
sampleType <- factor(c(gsub( "[0-9]$", "", sampleName)))
sampleRep <- factor(gsub("[^0-9]", "", sampleName))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

colnames(sequin_counts) <- sampleName

genemap <- sequin_conc[sequin_conc$ID %in% rownames(sequin_counts),]
rownames(genemap) <- genemap$ID
genemap <- genemap[rownames(sequin_counts),]
genemap$Expected_LFC <- log2(genemap$MIX_B/genemap$MIX_A)

design <- model.matrix(~sampleType)
y <- DGEList(counts=sequin_counts, group=sampleType, genes = genemap)
y <- calcNormFactors(y, method="TMM")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
res <- data.frame(topTags(lrt, n = nrow(sequin_counts)))

mod <- lm(logFC ~ Expected_LFC, res)
eq <- paste0("y = ", format(unname(coef(mod)[1]), digits = 2), " + ", format(unname(coef(mod)[2]), digits = 2), "x")
r2 <- paste0("r2 ", format(summary(mod)$r.squared, digits = 3))
eq_r2 <- paste0(eq, " ", r2)

ggplot(res, aes(x=Expected_LFC, y=logFC)) +
  annotate(geom = "text", y = 3, x = -3, label = eq_r2) +
  geom_point(shape = 1) +
  geom_smooth(method='lm', se=TRUE, color = 'blue', linetype = 'dashed') +
  theme_classic() +
  theme(aspect.ratio = 1) +
  xlab("Expected LFC") +
  ylab("Observed LFC") +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ggsave("project_results/ONT_metagenome/MixA_vs_MixB_LFC.pdf")

plotRLE(log2(fit$fitted.values + 1), outline=FALSE)
plotPCA(log2(fit$fitted.values + 1), cex=1)

# Scale SynX counts to 1% of gencode library size (NOTE: something weird going on with ratio of SynX to GENCODE counts)
synx_counts <- gene_counts[!grepl("_", rownames(gene_counts)),]

synx_size <- colSums(synx_counts)
lib_size <- colSums(gene_counts)
scaling_pct <- (lib_size*0.01)/synx_size
synx_counts_scaled <- data.frame(round(t(t(synx_counts)*scaling_pct)))

# Subset to dynamic sequins
dyn_seq <- sequin_conc[!is.na(sequin_conc$MIX_A), "ID"]

# Split gene counts into GENCODE and SynX genes
sequin_counts <- gene_counts[rownames(gene_counts) %in% dyn_seq,]

gene_counts_scaled <- rbind(synx_counts_scaled, sequin_counts)
colnames(gene_counts_scaled) <- sampleName


genemap <- genemap[rownames(gene_counts_scaled), ]

# SynX control features for normalisation
spikes <- rownames(gene_counts_scaled)[!grepl("_", rownames(gene_counts_scaled))]
# spikes <- spikes[grepl("P[0-9]", spikes)]

# Perform RUVg unwanted variation removal
set <- newSeqExpressionSet(as.matrix(gene_counts_scaled), phenoData = data.frame(sampleType, row.names=colnames(gene_counts_scaled)))
set1 <- RUVg(set, spikes, k=1)

# Find differentially expressed genes with RUVg normalisation
design <- model.matrix(~sampleType + W_1, data=pData(set1))
y1 <- DGEList(counts=counts(set), group=sampleType, genes = genemap)
y1 <- calcNormFactors(y1, method="TMM")
y1 <- estimateGLMCommonDisp(y1, design)
y1 <- estimateGLMTagwiseDisp(y1, design)
fit1 <- glmFit(y1, design)
lrt1 <- glmLRT(fit1, coef=2)
res <- data.frame(topTags(lrt1, n = nrow(counts(set))))

# Plot RUVg normalisation
plotRLE(set1, outline=FALSE)
plotPCA(set1, cex=1)




