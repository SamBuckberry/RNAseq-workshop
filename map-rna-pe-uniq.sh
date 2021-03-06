#!/bin/bash

# Fastq input
readR1="$1"
readR2="$2"
cores="$3"

# Set base name for output
baseName=$(basename "$readR1" .fastq.gz)

# Argument 5: Path to the index folder
indexPath="/home/buckberrydata_gmail_com/working/grch37_tran/genome_tran"

# Trim the reads
fastp -i "$readR1" -I "$readR2" -o "$baseName"_trimmed_R1.fastq.gz -O "$baseName"_trimmed_R2.fastq.gz

# Align the reads using hisat2
hisat2 --time --dta \
--rna-strandness RF \
--threads "$cores" \
-k 1 \
--no-mixed \
--met-file "$baseName"_hisat_metrics.txt \
-x "$indexPath" \
-1 "$baseName"_trimmed_R1.fastq.gz -2 "$baseName"_trimmed_R2.fastq.gz \
-S "$baseName".sam

# convert to bam, sort and index alignment
samtools view -bSu "$baseName".sam > "$baseName".bam
samtools sort -T "$baseName"_sorted "$baseName".bam > "$baseName".sorted.bam
samtools index "$baseName".sorted.bam

# Create md5 sum of alignment file
md5sum "$baseName".sorted.bam > "$baseName".sorted.bam.md5
