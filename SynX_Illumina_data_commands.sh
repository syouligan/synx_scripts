# --------------------------------------------------------------------------
# Bash commands used to analyse Illumna SynX samples
# --------------------------------------------------------------------------

# Align short read RNA sequencing with Bowtie (splice aware)
# --------------------------------------------------------------------------

# Copy RNA files to working directories
cd /tim/mer/scott/synx/data/
mkdir Illumina_SP6 Illumina_T7 Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2

cp /tim/mer/scott/data/illumina/*FD02702647*.fastq.gz /tim/mer/scott/synx/data/Illumina_SP6
cp /tim/mer/scott/data/illumina/*FD02702706*.fastq.gz /tim/mer/scott/synx/data/Illumina_T7
cp /tim/mer/scott/data/illumina/*FD02702797*.fastq.gz /tim/mer/scott/synx/data/Illumina_T7untreated1
cp /tim/mer/scott/data/illumina/*FD02702998*.fastq.gz /tim/mer/scott/synx/data/Illumina_T7untreated2
cp /tim/mer/scott/data/illumina/*FD02704376*.fastq.gz /tim/mer/scott/synx/data/Illumina_T7torkinib1
cp /tim/mer/scott/data/illumina/*FD02704436*.fastq.gz /tim/mer/scott/synx/data/Illumina_T7torkinib2

# Make custom genome including SynX + GRCh38
cd  /tim/mer/scott/synx/data/reference_files/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz
cat  /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.fa /tim/mer/scott/synx/data/reference_files/synx.fa >  /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.SynX.fa
bowtie2-build /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.SynX.fa /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.SynX

# Align and quantify errors for each sample (SynX + GRCh38)
sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  bowtie2 -x /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.SynX -1 /tim/mer/scott/synx/data/$sample/*_R1.fastq.gz -2 /tim/mer/scott/synx/data/$sample/*_R2.fastq.gz --threads 32 -S /tim/mer/scott/synx/project_results/$sample/$sample.sam
  samtools view -bS /tim/mer/scott/synx/project_results/$sample/$sample.sam > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools sort /tim/mer/scott/synx/project_results/$sample/$sample.bam -o /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  samtools index /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/GRCh38.primary_assembly.genome.SynX.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam > /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed.tsv
done

# Subsample to 10 million reads
sample_arr=(Illumina_SP6 Illumina_T7)
for samp in ${sample_arr[@]}; do
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R1.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq' 
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq'
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R2.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  done

# Use jelly fish to count 31mers aligning to cn ladder.
sample_arr=(Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample/jellyfish/
  jellyfish count -m 31 -t 30 -s 1G -C -o /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz') <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz')
  jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/All_31mers.fasta /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf > /tim/mer/scott/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv
  done

# Make bowtie index for SynX genome
bowtie2-build /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/reference_files/synx

# Align and quantify errors for each sample (SynX only)
sample_arr=(Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  bowtie2 -x /tim/mer/scott/synx/data/reference_files/synx -1 /tim/mer/scott/synx/data/$sample/*_10mil_R1.fastq.gz -2 /tim/mer/scott/synx/data/$sample/*_10mil_R2.fastq.gz --threads 32 -S /tim/mer/scott/synx/project_results/$sample/$sample.sam
  samtools view -@12 -bS /tim/mer/scott/synx/project_results/$sample/$sample.sam > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools sort -@12 /tim/mer/scott/synx/project_results/$sample/$sample.bam -o /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  samtools index -@12 /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam > /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

sample_arr=(Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done

# Align short read DNA sequencing with BWA-MEM2 - continue once we have DNA data
# --------------------------------------------------------------------------

# Get latest version of BWA-MEM2
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.0pre2/bwa-mem2-2.0pre2_x64-linux.tar.bz2 \
>   | tar jxf -

# Make index
/tim/mer/scott/synx/scripts/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index /tim/mer/scott/synx/data/reference_files/synx.fa

# Copy reads to working directory
cd /tim/mer/scott/synx/data/
mkdir Illumina_DNASynX Illumina_NadiaAburntgrassland Illumina_Naturalpineforest Illumina_Naturalgrassland
cp /tim/mer/scott/data/illumina_DNA/*FD03061374*.fastq.gz /tim/mer/scott/synx/data/Illumina_DNASynX
cp /tim/mer/scott/data/illumina_DNA/*FD09251653*.fastq.gz /tim/mer/scott/synx/data/Illumina_NadiaAburntgrassland
cp /tim/mer/scott/data/illumina_DNA/*FD09251654*.fastq.gz /tim/mer/scott/synx/data/Illumina_Naturalpineforest
cp /tim/mer/scott/data/illumina_DNA/*FD09251655*.fastq.gz /tim/mer/scott/synx/data/Illumina_Naturalgrassland

# Subsample to 10 million reads
sample_arr=(Illumina_DNASynX)
for samp in ${sample_arr[@]}; do
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R1.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq' 
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq'
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R2.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  done

# Use jelly fish to count 31mers aligning to cn ladder.
sample_arr=(Illumina_DNASynX)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample/jellyfish
  jellyfish count -m 31 -t 30 -s 1G -C -o /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz') <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz')
  jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/All_31mers.fasta /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf > /tim/mer/scott/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv
  done

sample_arr=(Illumina_DNASynX)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish/
  done

# Align reads and quanitfy errors basewise
sample_arr=(Illumina_DNASynX)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  /tim/mer/scott/synx/scripts/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 30 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz' /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz' -o /tim/mer/scott/synx/project_results/$sample/$sample.sam
  samtools view -@30 -bS /tim/mer/scott/synx/project_results/$sample/$sample.sam > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools sort -@30 /tim/mer/scott/synx/project_results/$sample/$sample.bam -o /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  samtools index -@30 /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam > /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done

sample_arr=(Illumina_DNASynX)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample.bam.bed.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done
