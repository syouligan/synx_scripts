# --------------------------------------------------------------------------
# Run analyzePile.py on ONT-DNA samples
# --------------------------------------------------------------------------

# Pile up per base across all reads
pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/data/ONT_DNA/barcode06.sorted.bam > /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed
python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed > /tim/mer/scott/synx/project_results/ONT_DNA/barcode06.sorted.bam.bed.tsv

# Extracts kmer sequences on forward strand
python /tim/mer/scott/synx/scripts/buildKM.py /tim/mer/scott/synx/data/reference_files/synx.fa 6 > /tim/mer/scott/synx/project_results/ONT_DNA/synx_6mer.bed

# Make bed file of homopolymer sequences (from https://www.biostars.org/p/379454/)
# Get phix sequence from here https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta
HOMOPOLYMER_LENGTH=2
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/synx.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/synx_homopolymers.bed
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/phix.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/phix_homopolymers.bed
awk '/^>/{if (NR==1) {print $0} else if (NR>1) {print "\n"$0}}; !/^>/ {printf toupper($0)}' /tim/mer/scott/synx/data/reference_files/synx_phix.fa | awk -v len="${HOMOPOLYMER_LENGTH}" -F '' '/^>/{gsub("^>", "", $0); chr=$0}; !/^>/ {base=$1; start=1; end=1; for (i=2; i<=NF+1; i++) {if (base==$(i)) {end=i} else if (base!=$(i)) {if (end-start>=len-1) {print chr"\t"start"\t"end"\t"base"\t"(end-start)+1}; start=i; end=i; base=$(i)}}}' > /tim/mer/scott/synx/project_results/ONT_DNA/synx_phix_homopolymers.bed

# Basewise error for each read single entry bed file (change back hardcoded analyzePile.py location, added to PATH)
samtools view -b -L /tim/mer/scott/synx/scripts/synx/test_bed_1_interval.bed /tim/mer/scott/synx/data/ONT_DNA/barcode06.sorted.bam > /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.bam
samtools index /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.bam
python3.6 /tim/mer/scott/synx/scripts/synx/script.py /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq /tim/mer/scott/synx/data/reference_files/test_bed_1_interval.bed /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.bam  > /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.perBase.tsv
python3.6 /tim/mer/scott/synx/scripts/synx/merge.py /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.perBase.tsv > /tim/mer/scott/synx/project_results/ONT_DNA/test_bed_1_interval.barcode06.sorted.INTERSECT.perRead.tsv

# Align read data and find base wise error rates
sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.bam > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

# Copy output to local directory
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/*.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done


# Use jelly fish to count 31mers aligning to cn ladder. Didnt use in the end, aligned specifically to cn sequences instead.
jellyfish count -m 31 -s 100M -t 10 -C /tim/mer/scott/synx/data/ONT_DNA/barcode06.fastq
jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn1_31mers.fasta mer_counts.jf > /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn1_31mers_counts.tsv
jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn2_31mers.fasta mer_counts.jf > /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn2_31mers_counts.tsv
jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn3_31mers.fasta mer_counts.jf > /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn3_31mers_counts.tsv
jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn4_31mers.fasta mer_counts.jf > /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/cn4_31mers_counts.tsv

# Align to copy number information
sample_arr=(ONT_DNA_BamHI ONT_DNA_EcoRI ONT_DNA_Fragmentase ONT_DNA_HindIII ONT_DNA_DCS)
for sample in ${sample_arr[@]}; do
  minimap2 -ax map-ont -t 8 /tim/mer/scott/synx/data/reference_files/synx_cn_only.fa /tim/mer/scott/synx/data/$sample/*.fastq | samtools sort - > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx_cn_only.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.synx_cn_only.bam.bed.tsv
done
