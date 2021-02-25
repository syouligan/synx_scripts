# --------------------------------------------------------------------------
# Bash commands used to analyse Illumna SynX samples
# --------------------------------------------------------------------------

# Align short read RNA sequencing with Bowtie2
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

# Align short read DNA sequencing with BWA-MEM2
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

# Subsample to 1 million reads for visualisation
samtools view -s 0.1 -b /tim/mer/scott/synx/project_results/Illumina_DNASynX/Illumina_DNASynX.sorted.bam -O BAM --threads 32 | samtools sort > /tim/mer/scott/synx/project_results/Illumina_DNASynX/Illumina_DNASynX_1mil.sorted.bam

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

# Simulate reads using Wgsim
mkdir /tim/mer/scott/synx/data/Illumina_DNASynX_simulation/
wgsim -N 10000000 -1 150 -2 150 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/Illumina_DNASynX_simulation/'Illumina_DNASynX_simulation_10mil_R1.fastq' /tim/mer/scott/synx/data/Illumina_DNASynX_simulation/'Illumina_DNASynX_simulation_10mil_R2.fastq'

# Use jelly fish to count 31mers aligning to cn ladder.
sample_arr=(Illumina_DNASynX_simulation)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample/jellyfish
  jellyfish count -m 31 -t 30 -s 1G -C -o /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz') <(zcat /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz')
  jellyfish query -s /tim/mer/scott/synx/project_results/ONT_DNA/jellyfish/All_31mers.fasta /tim/mer/scott/synx/project_results/$sample/jellyfish/mer_counts.jf > /tim/mer/scott/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv
  done

sample_arr=(Illumina_DNASynX_simulation)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/jellyfish/All_31mers_counts.tsv ~/cloudstor/tim_projects/synx/project_results/$sample/jellyfish/
  done

# Align short read RNA sequencing with bowtie2 and count with HTSeq
# --------------------------------------------------------------------------

# Make custom genome including SynX  and GRCh38
cd  /tim/mer/scott/synx/data/reference_files/
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz

# Make SynX GTF
sed -i '/^#/d' /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_gene.gtf
sed -i '/^#/d' /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_transcript.gtf
sed -i '/^#/d' /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_exon.gtf
cat /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_gene.gtf /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_transcript.gtf /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_exon.gtf > /tim/mer/scott/synx/data/reference_files/SynX_cn_masked.gtf

# Need ; after transcript id for RSEM and sort by position
sed -i 's/$/;/' /tim/mer/scott/synx/data/reference_files/SynX_cn_masked.gtf
# sort -t$'\t' -g -k4,4 /tim/mer/scott/synx/data/reference_files/SynX_cn_masked.gtf > /tim/mer/scott/synx/data/reference_files/SynX_cn_masked_sorted.gtf

# Generate genome reference containing GENCODE v36 and SynX sequences
cat /tim/mer/scott/synx/data/reference_files/GENCODE_v36/GRCh38.primary_assembly.genome.fa /tim/mer/scott/synx/data/reference_files/synx_CN_mask.fa > /tim/mer/scott/synx/data/reference_files/GENCODE_SynX_cn_masked.fa
cat /tim/mer/scott/synx/data/reference_files/GENCODE_v36/gencode.v36.primary_assembly.annotation.gtf /tim/mer/scott/synx/data/reference_files/SynX_cn_masked.gtf > /tim/mer/scott/synx/data/reference_files/GENCODE_SynX_cn_masked.gtf

# Trim low uality reads and subsample to 10 million reads
sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2)
for samp in ${sample_arr[@]}; do
  /tim/mer/scott/synx/scripts/fastp -i /tim/mer/scott/synx/data/$samp/*_R1.fastq.gz -I /tim/mer/scott/synx/data/$samp/*_R2.fastq.gz -o /tim/mer/scott/synx/data/$samp/$samp'_trimmed_R1.fastq.gz' -O /tim/mer/scott/synx/data/$samp/$samp'_trimmed_R2.fastq.gz'
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/$samp'_trimmed_R1.fastq.gz' 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq' 
  gzip -f /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq'
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/$samp'_trimmed_R2.fastq.gz' 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  gzip -f /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  done

# # Generate STAR and RSEM genome references
/tim/mer/scott/synx/scripts/STAR-2.7.7a/bin/Linux_x86_64_static/STAR --runMode genomeGenerate --genomeDir /tim/mer/scott/synx/data/reference_files/STAR_150 --genomeFastaFiles /tim/mer/scott/synx/data/reference_files/STAR_150/GENCODE_SynX_cn_masked.fa --outFileNamePrefix /tim/mer/scott/synx/data/reference_files/STAR_150/ --runThreadN 32 --sjdbGTFfile /tim/mer/scott/synx/data/reference_files/STAR_150/GENCODE_SynX_cn_masked.gtf --sjdbOverhang 149

# Align using STAR and quantify using HTseq
sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2 Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[0]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample/STAR
  /tim/mer/scott/synx/scripts/STAR-2.7.7a/bin/Linux_x86_64_static/STAR \
  --genomeDir /tim/mer/scott/synx/data/reference_files/STAR_150 \
  --readFilesIn /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz' /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz' \
  --readFilesCommand zcat \
  --outSAMtype BAM Unsorted \
  --quantMode GeneCounts \
  --runThreadN 20 \
  --outFileNamePrefix /tim/mer/scott/synx/project_results/$sample/STAR/$sample.
  done
  
htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --order=name /tim/mer/scott/synx/project_results/Illumina_T7untreated1/STAR/Illumina_T7untreated1.Aligned.out.bam /tim/mer/scott/synx/data/reference_files/GENCODE_SynX_cn_masked.gtf > /tim/mer/scott/synx/project_results/Illumina_T7untreated1/STAR/Illumina_T7untreated1'_HTSeq.txt'

#
#   # samtools index -@30 /tim/mer/scott/synx/project_results/$sample/$sample.sorted.Aligned.out.bam
#   # htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --order=pos /tim/mer/scott/synx/project_results/$sample/$sample.sorted.Aligned.out.bam /tim/mer/scott/synx/data/reference_files/GENCODE_RNASequin_SynX.gtf > /tim/mer/scott/synx/project_results/$sample/$sample.txt
#
#     mkdir /tim/mer/scott/synx/project_results/$sample/RSEM
#   rsem-calculate-expression --alignments /tim/mer/scott/synx/project_results/$sample/*Aligned.toTranscriptome.out.bam \
#   --no-bam-output \
#   --num-threads 20 \
#   --paired-end \
#   --strandedness reverse \
#   /tim/mer/scott/synx/data/reference_files/GENCODE_v36/RSEM/RSEM \
#   /tim/mer/scott/synx/project_results/$sample/RSEM/$sample.
#   done

# Make Kallisto index
samtools faidx /tim/mer/scott/synx/data/reference_files/synx_CN_mask.fa
bedtools getfasta -name -fi /tim/mer/scott/synx/data/reference_files/synx_CN_mask.fa -bed /tim/mer/scott/synx/data/reference_files/corrected_synx_annotations_CN_mask.bed -fo /tim/mer/scott/synx/data/reference_files/synx_cn_mask_transcripts.fa
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.pc_transcripts.fa.gz
gunzip /tim/mer/scott/synx/data/reference_files/GENCODE_v36/gencode.v36.pc_transcripts.fa.gz
cat /tim/mer/scott/synx/data/reference_files/GENCODE_v36/gencode.v36.pc_transcripts.fa /tim/mer/scott/synx/data/reference_files/synx_cn_mask_transcripts.fa > /tim/mer/scott/synx/data/reference_files/GENCODEpc_SynX_cn_masked_transcripts.fa
/tim/mer/scott/synx/scripts/kallisto/kallisto index --index=/tim/mer/scott/synx/data/reference_files/GENCODEpc_SynX_cn_masked_transcripts.idx --make-unique /tim/mer/scott/synx/data/reference_files/GENCODEpc_SynX_cn_masked_transcripts.fa

# Calculate abundance with Kallisto
sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2 Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  /tim/mer/scott/synx/scripts/kallisto/kallisto quant --index=/tim/mer/scott/synx/data/reference_files/GENCODEpc_SynX_cn_masked_transcripts.idx --genomebam --gtf /tim/mer/scott/synx/data/reference_files/GENCODE_SynX_cn_masked.gtf --threads 30 --rf-stranded --plaintext -o /tim/mer/scott/synx/project_results/$sample/ /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz' /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz'
  done

sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2 Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/abundance.tsv ~/cloudstor/tim_projects/synx/project_results/$sample
  done


# Align reads using bowtie and count using HTseq - note HTSeq need counts sorted by name otherwise needs alot of memory.
bowtie2-build /tim/mer/scott/synx/data/reference_files/GENCODE_RNASequin_SynX.fa /tim/mer/scott/synx/data/reference_files/GENCODE_RNASequin_SynX

sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2 Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  bowtie2 -x /tim/mer/scott/synx/data/reference_files/GENCODE_RNASequin_SynX -1 /tim/mer/scott/synx/data/$sample/HTCHGDSXY*_R1.fastq.gz -2 /tim/mer/scott/synx/data/$sample/HTCHGDSXY*_R2.fastq.gz --threads 32 | samtools view -@32 -bS - > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools sort -@32 -n /tim/mer/scott/synx/project_results/$sample/$sample.bam -o /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --order=name /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam /tim/mer/scott/synx/data/reference_files/GENCODE_RNASequin_SynX_cn_masked.gtf > /tim/mer/scott/synx/project_results/$sample/$sample'_HTSeq.txt'
  done

sample_arr=(Illumina_T7untreated1 Illumina_T7untreated2 Illumina_T7torkinib1 Illumina_T7torkinib2 Illumina_SP6 Illumina_T7)
for sample in ${sample_arr[@]}; do
  mkdir ~/cloudstor/tim_projects/synx/project_results/$sample
  cp ~/Desktop/ClusterHome/synx/project_results/$sample/$sample'_HTSeq.txt' ~/cloudstor/tim_projects/synx/project_results/$sample
  done

# Count kmers in Illumina Metagenome data
# --------------------------------------------------------------------------

# Subsample to 10 million reads
sample_arr=(Illumina_Burntgrassland Illumina_Naturalgrassland Illumina_Naturalpineforest)
for samp in ${sample_arr[@]}; do
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R1.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq' 
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R1.fastq'
  seqtk sample -s100 /tim/mer/scott/synx/data/$samp/*_R2.fastq.gz 10000000 > /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  gzip /tim/mer/scott/synx/data/$samp/$samp'_10mil_R2.fastq'
  done

sample_arr=(Illumina_Burntgrassland Illumina_Naturalgrassland Illumina_Naturalpineforest)
for sample in ${sample_arr[@]}; do
  mkdir /tim/mer/scott/synx/project_results/$sample
  /tim/mer/scott/synx/scripts/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 30 /tim/mer/scott/synx/data/reference_files/synx.fa /tim/mer/scott/synx/data/$sample/$sample'_10mil_R1.fastq.gz' /tim/mer/scott/synx/data/$sample/$sample'_10mil_R2.fastq.gz' -o /tim/mer/scott/synx/project_results/$sample/$sample.sam
  samtools view -@30 -bS /tim/mer/scott/synx/project_results/$sample/$sample.sam > /tim/mer/scott/synx/project_results/$sample/$sample.bam
  samtools sort -@30 /tim/mer/scott/synx/project_results/$sample/$sample.bam -o /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  samtools index -@30 /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam
  pysamstats --fasta /tim/mer/scott/synx/data/reference_files/synx.fa --type variation /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam > /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed
  python3 /tim/mer/scott/synx/scripts/analyzePile.py /tim/mer/scott/synx/project_results/$sample/$sample.sorted.bam.bed > /tim/mer/scott/synx/project_results/$sample/$sample.bam.bed.tsv
  done


