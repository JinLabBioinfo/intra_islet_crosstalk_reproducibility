#!/bin/bash

# Usage:
# ./FinalMappingHg19.sh <sample_name> <fastq_r1> <fastq_r2> <fra_downs>
#
# Example:
# ./FinalMappingHg19.sh H6_New H6_R1.fastq.gz H6_R2.fastq.gz 1
#
# nohup example:
# nohup ./FinalMappingHg19.sh H6_New H6_R1.fastq.gz H6_R2.fastq.gz 1 &
# nohup ./FinalMappingHg19.sh H9 Islet12219_R1.fastq.gz Islet12219_R2.fastq.gz 1 &
#
# Required environment variables:
# PICARD_DIR=/path/to/picard_tools_dir
# PICARD_JAR=/path/to/picard.jar
# DROPSEQ_DIR=/path/to/dropseq_tools_dir
# DROPSEQ_LIB_DIR=/path/to/dropseq_lib_dir
# STAR_BIN_DIR=/path/to/star_bin_dir
# STAR_GENOME_DIR=/path/to/star_genome_dir
# REFERENCE_FASTA=/path/to/genome.fa
# REFFLAT_FILE=/path/to/refFlat
# SAMTOOLS_BIN=/path/to/samtools
# TMP_DIR=/path/to/tmp_dir

myname=$1
genome=hg19
fq1=$2
fq2=$3
fra_downs=$4  # set "1" if you don't need downsample
cellnumber=20000
cellnumber2=40000
cellnumber5=100000

handle_error() {
  echo "Error: $1"
  exit 1
}

if [ "$#" -ne 4 ]; then
  handle_error "Usage: $0 <sample_name> <fastq_r1> <fastq_r2> <fra_downs>"
fi

for var_name in PICARD_DIR PICARD_JAR DROPSEQ_DIR DROPSEQ_LIB_DIR STAR_BIN_DIR STAR_GENOME_DIR REFERENCE_FASTA REFFLAT_FILE SAMTOOLS_BIN TMP_DIR; do
  [ -n "${!var_name}" ] || handle_error "Environment variable $var_name is required."
done

picard="$PICARD_DIR"
picardNew="$PICARD_JAR"
dropseq="$DROPSEQ_DIR"
mylib="$DROPSEQ_LIB_DIR"
star="$STAR_BIN_DIR"
genomeDir="$STAR_GENOME_DIR"
samtools="$SAMTOOLS_BIN"
TMP_DIR="$TMP_DIR"

##########################
# Step 1: Prepare uncompressed FASTQ paths
if [[ "$fq1" == *.gz ]]; then
  fq1_unzipped="${fq1%.gz}"
  echo "Decompressing $fq1 to $fq1_unzipped"
  zcat "$fq1" > "$fq1_unzipped"
else
  fq1_unzipped="$fq1"
fi

if [[ "$fq2" == *.gz ]]; then
  fq2_unzipped="${fq2%.gz}"
  echo "Decompressing $fq2 to $fq2_unzipped"
  zcat "$fq2" > "$fq2_unzipped"
else
  fq2_unzipped="$fq2"
fi

# Step 2: Merge R1 and R2
perl $mylib/fastqmerge.pl "$fq1_unzipped" "$fq2_unzipped" 20 50 > "$myname.R1R2.fastq"
java -Xmx32g -Djava.io.tmpdir=$TMP_DIR \
  -jar $picardNew FastqToSam \
  FASTQ=$myname.R1R2.fastq \
  OUTPUT=$myname.R1R2.bam \
  SAMPLE_NAME=$myname \
  READ_GROUP_NAME=A
$samtools view -h -o $myname.bam $myname.R1R2.bam
$dropseq/TagBamWithReadSequenceExtended INPUT=$myname.bam OUTPUT=$myname.unaligned_tagged_Cell.bam SUMMARY=$myname.unaligned_tagged_Cell.bam_summary.txt BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=False TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1
$dropseq/TagBamWithReadSequenceExtended INPUT=$myname.unaligned_tagged_Cell.bam  OUTPUT=$myname.unaligned_tagged_CellMolecular.bam SUMMARY=$myname.unaligned_tagged_Molecular.bam_summary.txt BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=True TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1
$dropseq/FilterBAM TAG_REJECT=XQ INPUT=$myname.unaligned_tagged_CellMolecular.bam OUTPUT=$myname.unaligned_tagged_filtered.bam
$dropseq/TrimStartingSequence  INPUT=$myname.unaligned_tagged_filtered.bam  OUTPUT=$myname.unaligned_tagged_trimmed_smart.bam OUTPUT_SUMMARY=$myname.adapter_trimming_report.txt SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5
$dropseq/PolyATrimmer INPUT=$myname.unaligned_tagged_trimmed_smart.bam OUTPUT=$myname.unaligned_mc_tagged_polyA_filtered.bam OUTPUT_SUMMARY=$myname.polyA_trimming_report.txt MISMATCHES=0 NUM_BASES=6
$samtools view $myname.unaligned_mc_tagged_polyA_filtered.bam |$mylib/trimsam.pl | $samtools view -h -o $myname.unaligned_mc_tagged_polyA_filtered.trimed.bam -
$samtools view -h $myname.unaligned_mc_tagged_polyA_filtered.trimed.bam | cat <($samtools view -H $myname.unaligned_mc_tagged_polyA_filtered.bam | grep -v "^@PG") - > $myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam
java -Xmx32g -Djava.io.tmpdir=$TMP_DIR -jar $picard/SamToFastq.jar  I=$myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam F=$myname.unaligned_mc_tagged_polyA_filtered.trimed.fastq
$star/STAR --runThreadN 8 --genomeDir $genomeDir --readFilesIn $myname.unaligned_mc_tagged_polyA_filtered.trimed.fastq --outFileNamePrefix $myname.star
java -Xmx32g -Djava.io.tmpdir=$TMP_DIR -jar $picard/SortSam.jar I=$myname.starAligned.out.sam  O=$myname.aligned.sorted.bam SO=queryname
java -Xmx32g -Djava.io.tmpdir=$TMP_DIR -jar $picard/MergeBamAlignment.jar UNMAPPED_BAM=$myname.unaligned_mc_tagged_polyA_filtered.trimed.h.bam ALIGNED_BAM=$myname.aligned.sorted.bam OUTPUT=$myname.merged.bam  REFERENCE_SEQUENCE=$REFERENCE_FASTA INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false
$dropseq/TagReadWithGeneExon I=$myname.merged.bam O=$myname.gene_exon_tagged.bam ANNOTATIONS_FILE=$REFFLAT_FILE  TAG=GE
$dropseq/DetectBeadSynthesisErrors I=$myname.gene_exon_tagged.bam O=$myname.gene_exon_tagged.cleaned.bam OUTPUT_STATS=$myname.gene_exon_tagged.STATS SUMMARY=$myname.gene_exon_tagged.synthesis_stats.summary.txt NUM_BARCODES=$cellnumber2  PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
$dropseq/DigitalExpression I=$myname.gene_exon_tagged.cleaned.bam O=$myname.gene_exon_tagged.cleaned.dge.txt.gz SUMMARY=$myname.gene_exon_tagged.cleaned.summary.txt NUM_CORE_BARCODES=$cellnumber
# Cleanup intermediate files while keeping final results and raw FASTQs
echo "Cleaning up intermediate files..."
ls ${myname}* | grep -v -E "gene_exon_tagged.cleaned|${myname}_R[12]\.fastq(\.gz)?$" | xargs rm -v
echo "Job Done!"
