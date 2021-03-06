# --------------------------------------------------------------------------
# Command line for ONT-DNA samples
# --------------------------------------------------------------------------

# Pile up per base across all reads
# --------------------------------------------------------------------------

pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/data/ONT_DNA/barcode06.sorted.bam > /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed
python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed > /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed.tsv

# Extracts kmer sequences on forward strand
# --------------------------------------------------------------------------

python /tim/mer/scott/synx/scripts/buildKM.py /tim/mer/scott/synx/data/reference_files/synx.fa 31 > /tim/mer/scott/synx/project_results/ONT_DNA/synx_31mer.bed

# Make bed file of homopolymer sequences (from https://www.biostars.org/p/379454/)
# --------------------------------------------------------------------------

# Get phix sequence from here https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta
HOMOPOLYMER_LENGTH=2
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/synx.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/synx_homopolymers.bed
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/phix.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/phix_homopolymers.bed
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/synx_phix.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/synx_phix_homopolymers.bed

# Basewise error for each read single entry bed file (change back hardcoded analyzePile.py location, added to PATH)
# --------------------------------------------------------------------------

# For flongle
feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'
  samtools view -b -L /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/data/ONT_DNA/barcode06.sorted.bam > /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.bam
  samtools index /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.bam
  python3.6 /tim/mer/scott/synx/scripts/synx/script.py /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.bam  > /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perBase.tsv
  python3.6 /tim/mer/scott/synx/scripts/synx/merge.py /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perBase.tsv > /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perRead.tsv
  done

feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/ONT_DNA/$feature'_perRead'
  cp ~/Desktop/ClusterHome/synx/project_results/ONT_DNA/$feature'_perRead'/*perRead.tsv ~/cloudstor/tim_projects/synx/project_results/ONT_DNA/$feature'_perRead'
  done

# For minion, subsample fragmentase library first
seqtk sample -s100 /tim/mer/scott/synx/data/ONT_DNA_Fragmentase/DNA-Fragmentase.fastq 20000 > /tim/mer/scott/synx/data/ONT_DNA_Fragmentase/ONT_DNA_Fragmentase_20000.fastq 
minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/ONT_DNA_Fragmentase/ONT_DNA_Fragmentase_20000.fastq | samtools sort - > /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/ONT_DNA_Fragmentase_20000.bam
samtools index /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/ONT_DNA_Fragmentase_20000.bam

feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'
  samtools view -b -L /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/ONT_DNA_Fragmentase_20000.bam > /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.bam
  samtools index /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.bam
  python3.6 /tim/mer/scott/synx/scripts/synx/script.py /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/ONT_DNA_Fragmentase/'DNA-Fragmentase.fastq' /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.bam  > /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.perBase.tsv
  python3.6 /tim/mer/scott/synx/scripts/synx/merge.py /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.perBase.tsv > /tim/mer/scott/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/$feature'_one_line'.ONT_DNA_Fragmentase_20000.INTERSECT.perRead.tsv
  done

feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'
  cp ~/Desktop/ClusterHome/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'/*20000*perRead.tsv ~/cloudstor/tim_projects/synx/project_results/ONT_DNA_Fragmentase/$feature'_perRead'
  done

# Remove homopolymer regions from error counts for each read
# --------------------------------------------------------------------------

feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  awk -F '\t' 'NR==FNR {id[$1]; next} $11 in id' /tim/mer/scott/synx/data/reference_files/Non_homopolymer_positions.txt /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perBase.tsv > /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perBase.no_HP.tsv
  python3.6 /tim/mer/scott/synx/scripts/synx/merge.py /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perBase.no_HP.tsv > /tim/mer/scott/synx/project_results/ONT_DNA/$feature'_perRead'/$feature'_one_line'.barcode06.sorted.INTERSECT.perRead.no_HP.tsv
  done

feature_arr=(P1 P2 P3 P4 P5 P6 PH1 PH2 SynX)
for feature in ${feature_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/ONT_DNA/$feature'_perRead'
  cp ~/Desktop/ClusterHome/synx/project_results/ONT_DNA/$feature'_perRead'/*perRead.no_HP.tsv ~/cloudstor/tim_projects/synx/project_results/ONT_DNA/$feature'_perRead'
  done

# Align read data and find base wise error rates
# --------------------------------------------------------------------------

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done


# Use jelly fish to count 31mers aligning to cn ladder. Now redundant to R-script.
# --------------------------------------------------------------------------

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII ONT_DNA ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6 ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample/jellyfish/
  jellyfish count -m 31 -t 30 -s 100M -C /tim/mer/scott/synx/data/$sample/*.fastq -o /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf
  jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/All_31mers.fasta /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf > /tim/mer/scott/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv
  done

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII ONT_DNA ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6 ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish/
  cp -f ~/Desktop/ClusterHome/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish/
  done

# Align to copy number information. Alternative to using jellyfish above.  Now redundant to R-script.
# --------------------------------------------------------------------------

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx_cn_only.fa /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx_cn_only.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed.tsv
  done

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample/
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.synx_cn_only.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample/
  done


# Calculate per read statistics for samples in sample_arr across features in feature_arr
# --------------------------------------------------------------------------

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  feature_arr=(SynX)
  for feature in ${feature_arr[@]}; do
    mkdir /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'
    samtools view -b -L /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.bam
    samtools index /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.bam
    python3.6 /tim/mer/scott/synx/scripts/synx/script.py /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/$sample/*.fastq /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.bam  > /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.perBase.tsv
    python3.6 /tim/mer/scott/synx/scripts/synx/merge.py /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.perBase.tsv > /tim/mer/scott/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.perRead.tsv
    done
  done

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample/$feature'_perRead'/
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$feature'_perRead'/$feature'_one_line'.$sample.INTERSECT.perRead.tsv ~/cloudstor/tim_projects/synx/project_results/$sample/$feature'_perRead'/
  done

# Find distribution of read lengths in digested samples
# --------------------------------------------------------------------------

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  feature_arr=(SynX)
  for feature in ${feature_arr[@]}; do
    mkdir /tim/mer/scott/synx/project_results/$sample
    samtools view -L /tim/mer/scott/synx/data/reference_files/$feature'_one_line'.bed /tim/mer/scott/synx/project_results/$sample/$sample.bam \
    | mawk '{hist[length($10)]++} END {for (l in hist) print l"\t"hist[l]}' \
    | sort -n -k1 > /tim/mer/scott/synx/project_results/$sample/$sample.read_length.tsv
    done
  done

sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII)
sample_arr=(ONT_cDNA_SP6 ONT_cDNA_SP6_2 ONT_cDNA_T7 ONT_cDNA_T7_2 ONT_RNA_SP6)
sample_arr=(ONT_cDNA_SP6_3 ONT_cDNA_T7_3)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.read_length.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done

# Assemble DCS sequence (didnt find what this sequence was in the end...)
# --------------------------------------------------------------------------

# First establish estimated genome size and QC characteristics (flye tells N50/N90, mean depth, etc)
/home/tedwon/Sources/Flye/bin/flye --nano-raw /tim/mer/scott/synx/data/ONT_DNA_DCS/ONT_DNA_Control_Strand.fastq --out-dir /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly --threads 20 --meta
mv /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/assembly.fasta /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/meta_assembly.fasta

# # QC stats
# Total length:   16488
# Fragments:      7
# Fragments N50:  3247
# Largest frg:    5225
# Scaffolds:      0
# Mean coverage:  339
# 
# # Top blast hit for each contig
# Contig 1 - Return many bacterial genomes, Most likely Enterobacteria phage phiX174, complete genome.
# ! Perhaps an interesting observation - could be really important for metagenomes etc or even assessing contamination
# Contig 2 - Probably E.coli, hits many strains. Escherichia coli BL21(DE3) chromosome, complete genome
# Contig 3 - Another E.coli fragment
# Contig 4 - Short lmabda phage fragment
# Contig 5 - (longest) Matches many bacterial genome and phages, including phiX174
# Contig 6 - Another phix174 fragment
# Contig 7 - Bacteriophage sp. isolate 181, complete genome


# Final assembly (with estmiated genome length, min overlap of 1000bp and min depth of coverage of 100x)
/home/tedwon/Sources/Flye/bin/flye --nano-raw /tim/mer/scott/synx/data/ONT_DNA_DCS/ONT_DNA_Control_Strand.fastq --out-dir /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly --threads 20  -m 1000 --genome-size 5k --asm-coverage 100

# # QC stats
# Total length:   3819
# Fragments:      2
# Fragments N50:  2688
# Largest frg:    2688
# Scaffolds:      0
# Mean coverage:  867
# 
# # # Top blast hit for each contig
# Contig 2 - Top hit Bacteriophage sp. isolate 181, complete genome
# Contig 3 - Top hit Escherichia coli BL21(DE3) chromosome, complete genome

# Align DCS sample and view bam with IGV (both individual and combined contigs to examine contig spanning reads)
sample_arr=(ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/assembly.fasta /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.bam
  pysamstats --fasta /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/assembly.fasta --type variation /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

sample_arr=(ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.bam* ~/cloudstor/tim_projects/synx/project_results/$sample
  done
  
sample_arr=(ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/assembly_combined.fasta /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam
  pysamstats --fasta /tim/mer/scott/synx/project_results/ONT_DNA_DCS/assembly/assembly_combined.fasta --type variation /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam > /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample_combined.bam.bed.tsv
  done

sample_arr=(ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample_combined.bam* ~/cloudstor/tim_projects/synx/project_results/$sample
  done  

# Analylse phiX/synX error rates
# --------------------------------------------------------------------------

# Pile up error rates per base phix and synx
sample_arr=(ONT_DNA_phix_synx ONT_DNA_phix_synx_2)

for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx_phix.fa /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx_phix.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

sample_arr=(ONT_DNA_phix_synx ONT_DNA_phix_synx_2)

for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done

# Simulate Synx Reads using NanoSim
# --------------------------------------------------------------------------

pip3.6 --user install six numpy HTSeq==0.9.1 Pysam scipy scikit-learn==0.20.0
python3.6 /tim/mer/scott/synx/scripts/NanoSim3.0.0/src/read_analysis.py genome --read /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq --ref_g /tim/mer/scott/synx/data/reference_files/synx.fa --aligner minimap2 --num_threads 16 --output /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/NanoSim

# With errors simulated - giving different read/quality lengths...
python3.6 /tim/mer/scott/synx/scripts/NanoSim3.0.0/src/simulator.py genome --ref_g /tim/mer/scott/synx/data/reference_files/synx.fa --model_prefix /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/NanoSim --fastq --basecaller guppy --num_threads 16 --output /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/Simulated

# Without error simulated
python3.6 /tim/mer/scott/synx/scripts/NanoSim3.0.0/src/simulator.py genome --ref_g /tim/mer/scott/synx/data/reference_files/synx.fa --model_prefix /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/NanoSim --perfect --fastq --basecaller guppy --num_threads 16 --output /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/Perfect

# Align libraries
sample_arr=(Perfect Simulated)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample'_aligned_reads.fastq' | samtools sort - > /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam > /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam.bed > /tim/mer/scott/synx/project_results/ONT_DNA/NanoSim/$sample/$sample.bam.bed.tsv
done
  
# Simulate ONT reads using BadRead
minimap2 -t 32 -c -x map-ont /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq > /tim/mer/scott/synx/project_results/badread_simulator/alignments.paf
python3.6 /tim/mer/scott/synx/scripts/Badread/badread-runner.py error_model --reference /tim/mer/scott/synx/data/reference_files/synx.fa --reads /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq --alignment /tim/mer/scott/synx/project_results/badread_simulator/alignments.paf > /tim/mer/scott/synx/project_results/badread_simulator/badread_errors
python3.6 /tim/mer/scott/synx/scripts/Badread/badread-runner.py qscore_model --reference /tim/mer/scott/synx/data/reference_files/synx.fa --reads /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq --alignment /tim/mer/scott/synx/project_results/badread_simulator/alignments.paf > /tim/mer/scott/synx/project_results/badread_simulator/badread_qscores
python3.6 /tim/mer/scott/synx/scripts/Badread/badread-runner.py simulate --quantity 14500x --reference /tim/mer/scott/synx/data/reference_files/synx.fa --error_model random \
    --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 \
    --chimeras 0 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" > /tim/mer/scott/synx/project_results/badread_simulator/badread_reads.fastq

jellyfish count -m 31 -t 30 -s 100M -C /tim/mer/scott/synx/project_results/badread_simulator/*.fastq -o /tim/mer/scott/synx/project_results/badread_simulator/jellyfish/mer_counts.jf
jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/All_31mers.fasta /tim/mer/scott/synx/project_results/badread_simulator/jellyfish/mer_counts.jf > /tim/mer/scott/synx/project_results/badread_simulator/jellyfish/All_31mers_counts.tsv

mkdir ~/cloudstor/tim_projects/synx/project_results/badread_simulator/jellyfish/
cp -f ~/Desktop/ClusterHome/synx/project_results/badread_simulator/jellyfish/All_31mers_counts.tsv ~/cloudstor/tim_projects/synx/project_results/badread_simulator/jellyfish/

# Align ONT metagenome samples
# --------------------------------------------------------------------------

# Make genome containing SynX, Metagenome sequins, PhiX sequence
cat /tim/mer/scott/synx/data/reference_files/Metagenome_sequin_resources/metasequin_sequences_3.0.fa /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/reference_files/phix.fa > /tim/mer/scott/synx/data/reference_files/MetaSequin_SynX_PhiX.fa

# Concatenate pass fastq files
sample_arr=(barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08 barcode09 barcode10 barcode11 barcode12)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/data/ONT_metagenome/$sample
  cat /tim/mer/scott/synx/data/ONT_metagenome/GLFL278165/custom_ont_adaptor/20201130_0632_X4_FAO14644_036d5344/fastq_pass/$sample/*.fastq > /tim/mer/scott/synx/data/ONT_metagenome/$sample/$sample'_combined.fastq'
  done

cp -r barcode05 ../ONT_SeqMetGenA1
cp -r barcode06 ../ONT_SeqMetGenB1
cp -r barcode01 ../ONT_SeqMetGenA2
cp -r barcode02 ../ONT_SeqMetGenA3
cp -r barcode03 ../ONT_SeqMetGenB2
cp -r barcode04 ../ONT_SeqMetGenB3

# Align samples using minimap2
sample_arr=(ONT_SeqMetGenA1 ONT_SeqMetGenA2 ONT_SeqMetGenA3 ONT_SeqMetGenB1 ONT_SeqMetGenB2 ONT_SeqMetGenB3)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/MetaSequin_SynX_PhiX.fa /tim/mer/scott/synx/data/$sample/*'_combined.fastq' | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.bam
done

# Quantify covergae across features using bedtools (minimum 50% overlap with feature)
sample_arr=(ONT_SeqMetGenA1 ONT_SeqMetGenA2 ONT_SeqMetGenA3 ONT_SeqMetGenB1 ONT_SeqMetGenB2 ONT_SeqMetGenB3)
for sample in ${sample_arr[@]}; do
bedtools coverage -counts -f 0.5 -a /tim/mer/scott/synx/data/reference_files/MetaSequin_SynX.bed -b /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$sample.MetaSequin_SynX.coverage
done

# Transfer coverage estimates to local drive
sample_arr=(ONT_SeqMetGenA1 ONT_SeqMetGenA2 ONT_SeqMetGenA3 ONT_SeqMetGenB1 ONT_SeqMetGenB2 ONT_SeqMetGenB3)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.MetaSequin_SynX.coverage ~/cloudstor/tim_projects/synx/project_results/$sample
  done



