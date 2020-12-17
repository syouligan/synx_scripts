#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# List position of bases in a homopolymer sequence in the Synx genome
# --------------------------------------------------------------------------

if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/tim_projects/synx/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/tim/mer/scott/synx/")
  place <- "timmer"
}

synx_hp_ranges <- import.bed("project_results/ONT_DNA/synx_homopolymers.bed", genome = 'SynX')
synx_hp_ranges <- synx_hp_ranges[synx_hp_ranges$score > 2,]

hp_bases <- c()
for(i in 1:length(synx_hp_ranges)){
  st <- synx_hp_ranges[i, ]@ranges@start
  fi <- st + synx_hp_ranges[i, ]@ranges@width
  hp_bases <- c(hp_bases, seq(st, fi))
  }

write.table(hp_bases, "Homopolymer_positions.txt", col.names = FALSE, row.names = FALSE)

# List sites that are not homopolymer regions
# --------------------------------------------------------------------------

# Load synx bed file
synx_ranges <- import.bed("data/reference_files/corrected_synx_annotations.bed", genome = 'SynX')


synx_total <- seq(1, 7824)
non_hp_bases <- synx_total[!synx_total %in% hp_bases]

write.table(non_hp_bases, "Non_homopolymer_positions.txt", col.names = FALSE, row.names = FALSE)



