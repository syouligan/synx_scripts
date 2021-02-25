#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Make SynX GTF
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
library('rtracklayer')
library('GenomicAlignments')
library('seqQTL')

# Load GECODE and SEQUIN GTF files
# --------------------------------------------------------------------------

GENCODE_ranges <- import.gff("data/reference_files/gencode.v36.primary_assembly.annotation.gtf", genome = 'GENCODE_RNASequin_SynX')

# Load SynX bed file
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations_CN_mask.bed", genome = 'GENCODE_RNASequin_SynX')

# Load synx feature annotations and annotate ranges
synx_annotations <- read.table("design/custom_features.annotations.tab", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
idx <- match(synx_ranges$name, synx_annotations$NAME)
synx_ranges$source <- synx_annotations$ORIGIN [idx]
synx_ranges$type <- 'exon'
synx_ranges$score <- NA
synx_ranges$phase <- NA
synx_ranges$gene_id <- paste(synx_ranges$name, synx_ranges@ranges@start, sep = '_')
synx_ranges$gene_type <- synx_annotations$TYPE [idx]
synx_ranges$gene_type <- gsub(" ", "_", synx_ranges$gene_type)
synx_ranges$gene_name <- synx_ranges$name
synx_ranges$transcript_id <- paste0(synx_ranges$gene_id, "_1")
synx_ranges$transcript_type <- synx_ranges$gene_type
synx_ranges$transcript_name <- paste0(synx_ranges$gene_name, "_1")
synx_ranges$exon_number <- "1"
synx_ranges$exon_id <- paste0(synx_ranges$gene_id, "_1", ".1")

synx_ranges@elementMetadata <- synx_ranges@elementMetadata[,c("source", "type", "score", "phase", "gene_id", "gene_type", "gene_name", "transcript_id", "transcript_type", "transcript_name", "exon_number", "exon_id")]

# Make synx of exons
synx_ranges_exon <- synx_ranges
export(synx_ranges_exon, "data/reference_files/SynX_cn_masked_exon.gtf")

# Make synx of transcripts
synx_ranges_transcripts <- synx_ranges_exon
synx_ranges_transcripts$type <- 'transcript'
synx_ranges_transcripts <- synx_ranges_transcripts[,!grepl("exon", colnames(synx_ranges_transcripts@elementMetadata))]
export(synx_ranges_transcripts, "data/reference_files/SynX_cn_masked_transcript.gtf")

# Make synx of genes
synx_ranges_gene <- synx_ranges_transcripts
synx_ranges_gene$type <- 'gene'
synx_ranges_gene <- synx_ranges_gene[,!grepl("transcript", colnames(synx_ranges_gene@elementMetadata))]
export(synx_ranges_gene, "data/reference_files/SynX_cn_masked_gene.gtf")

# Combine Granges objects and save GTF
# --------------------------------------------------------------------------

# Concatenate Granges objects into a single object
combined_ranges <- c(synx_ranges_gene, synx_ranges_transcripts, synx_ranges_exon)
combined_ranges <- combined_ranges[order(combined_ranges@ranges@start, decreasing = FALSE),]

# Subset to common metadata columns
GR2gtf(combined_ranges, "data/reference_files/SynX_cn_masked.gtf")

# Combine with GENCODE genome
GENCODE_combined_ranges <- c(GENCODE_ranges, combined_ranges)
GR2gtf(GENCODE_combined_ranges, "data/reference_files/GENCODE_SynX_cn_masked.gtf")

