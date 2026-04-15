#!/bin/bash

# Usage: ./script_name.sh <fastq_file_1> <fastq_file_2> <output_name> <genome> <threads>
# Example: ./script_name.sh sample_R1.fastq sample_R2.fastq output_sample hg19 12
#
# Required environment variables:
# HISAT2_INDEX_PREFIX=/path/to/hisat2/index/prefix
# GTF_FILE=/path/to/annotation.gtf
# SAMTOOLS_BIN=/path/to/samtools
# HISAT2_BIN=/path/to/hisat2
# FEATURECOUNTS_BIN=/path/to/featureCounts
# BAMCOVERAGE_BIN=/path/to/bamCoverage

# Input fastq files and parameters
fq1=$1
fq2=$2
name=$3
genome=$4
threads=$5

# Function to handle errors
handle_error() {
  echo "Error: $1"
  exit 1
}

# Log start time
echo "Starting processing: $(date)"

# Create a new directory with $name and navigate to it
mkdir -p "$name" || handle_error "Directory creation failed"
cd "$name" || handle_error "Could not navigate to directory $name"
echo "Working directory: $(pwd)"

# Set up logging
exec > >(tee -a "${name}.log") 2>&1

# Check for required arguments
if [ "$#" -ne 5 ]; then
  handle_error "Usage: $0 <fastq_file_1> <fastq_file_2> <output_name> <genome> <threads>"
fi

for var_name in HISAT2_INDEX_PREFIX GTF_FILE SAMTOOLS_BIN HISAT2_BIN FEATURECOUNTS_BIN BAMCOVERAGE_BIN; do
  [ -n "${!var_name}" ] || handle_error "Environment variable $var_name is required."
done

# Paths to reference genome and GTF annotation
index="$HISAT2_INDEX_PREFIX"
gtf="$GTF_FILE"

# Define paths for each tool explicitly
SAMTOOLS="$SAMTOOLS_BIN"
HISAT2="$HISAT2_BIN"
FEATURECOUNTS="$FEATURECOUNTS_BIN"
BAMCOVERAGE="$BAMCOVERAGE_BIN"

# Version checks for tools (for logging and verification)
echo "Checking versions of required tools..."
$SAMTOOLS --version
$HISAT2 --version
$FEATURECOUNTS -v
$BAMCOVERAGE --version

# Ensure all dependencies and references are available
[ -x "$SAMTOOLS" ] || handle_error "samtools not found: $SAMTOOLS"
[ -x "$HISAT2" ] || handle_error "hisat2 not found: $HISAT2"
[ -x "$FEATURECOUNTS" ] || handle_error "featureCounts not found: $FEATURECOUNTS"
[ -x "$BAMCOVERAGE" ] || handle_error "bamCoverage not found: $BAMCOVERAGE"
[ -f "${index}.1.ht2" ] || [ -f "${index}.1.ht2l" ] || handle_error "hisat2 index not found for genome '$genome': $index"
[ -f "$gtf" ] || handle_error "GTF annotation not found: $gtf"

# Step 1: Map raw reads to reference genome
echo "Mapping raw reads to the reference genome... $(date)"
$HISAT2 -x "$index" -p "$threads" --sp 1000,1000 -1 "$fq1" -2 "$fq2" | $SAMTOOLS view -bS - > "${name}.merge.bam" || handle_error "hisat2 mapping failed"

# Step 2: Filter non-unique mapped reads and sort BAM file
echo "Filtering non-unique reads and sorting BAM file... $(date)"
$SAMTOOLS view -F 1804 -b "${name}.merge.bam" | $SAMTOOLS sort -@ "$threads" -o "${name}.sorted.bam" || handle_error "BAM sorting failed"

# Step 3: Annotate and count transcripts for each gene
echo "Counting transcripts for each gene... $(date)"
$FEATURECOUNTS --ignoreDup -p -B -a "$gtf" -t exon -g gene_id -o "${name}.counts.txt" "${name}.sorted.bam" || handle_error "featureCounts failed"

# Generate BED file with selected columns from the counts file
echo "Extracting gene counts into BED format... $(date)"
cut -f1,7 "${name}.counts.txt" > "${name}.count.bed"

# Step 4: Generate CPM normalized track
echo "Generating CPM normalized coverage track... $(date)"
$SAMTOOLS index "${name}.sorted.bam" || handle_error "Samtools indexing failed"
$BAMCOVERAGE --bam "${name}.sorted.bam" --binSize 1 --normalizeUsing CPM --outFileName "${name}.sorted.bam.CPM.bw" --numberOfProcessors "$threads" || handle_error "bamCoverage failed"

echo "Processing complete! $(date)"
exit 0
