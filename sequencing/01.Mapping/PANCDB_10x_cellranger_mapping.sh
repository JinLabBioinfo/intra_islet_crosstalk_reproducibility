#!/bin/bash

# Usage:
# ./PANCDB_10x_cellranger_mapping.sh <sample_id> <fastqs_dir> <sample_name> <localcores> <localmem>
#
# Example:
# ./PANCDB_10x_cellranger_mapping.sh HPAP034 /path/to/HPAP034 HPAP034 16 32
#
# Required environment variables:
# CELLRANGER_BIN=/path/to/CellRanger/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger
# TRANSCRIPTOME_DIR=/path/to/cellranger_reference

sample_id=$1
fastqs_dir=$2
sample_name=$3
localcores=$4
localmem=$5

handle_error() {
  echo "Error: $1"
  exit 1
}

if [ "$#" -ne 5 ]; then
  handle_error "Usage: $0 <sample_id> <fastqs_dir> <sample_name> <localcores> <localmem>"
fi

for var_name in CELLRANGER_BIN TRANSCRIPTOME_DIR; do
  [ -n "${!var_name}" ] || handle_error "Environment variable $var_name is required."
done

[ -x "$CELLRANGER_BIN" ] || handle_error "cellranger not found: $CELLRANGER_BIN"
[ -d "$TRANSCRIPTOME_DIR" ] || handle_error "transcriptome directory not found: $TRANSCRIPTOME_DIR"
[ -d "$fastqs_dir" ] || handle_error "fastqs directory not found: $fastqs_dir"

mkdir -p "$sample_id" || handle_error "Directory creation failed"
cd "$sample_id" || handle_error "Could not navigate to directory $sample_id"

$CELLRANGER_BIN count \
  --id="$sample_id" \
  --transcriptome="$TRANSCRIPTOME_DIR" \
  --fastqs="$fastqs_dir" \
  --sample="$sample_name" \
  --localcores="$localcores" \
  --localmem="$localmem" || handle_error "cellranger count failed"

echo "$sample_id finished successfully !!!"
