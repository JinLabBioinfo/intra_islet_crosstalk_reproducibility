#!/bin/bash

# Usage:
# ./Mouse_10X_cellanger_arc_mapping.sh <run_id> <libraries_csv> <localcores> <localmem>
#
# Example:
# ./Mouse_10X_cellanger_arc_mapping.sh Mouse /path/to/libraries.csv 16 32
#
# Required environment variables:
# CELLRANGER_ARC_BIN=/path/to/cellranger-arc
# CELLRANGER_ARC_REF=/path/to/refdata-cellranger-arc-mm10

run_id=$1
libraries_csv=$2
localcores=$3
localmem=$4

handle_error() {
  echo "Error: $1"
  exit 1
}

if [ "$#" -ne 4 ]; then
  handle_error "Usage: $0 <run_id> <libraries_csv> <localcores> <localmem>"
fi

for var_name in CELLRANGER_ARC_BIN CELLRANGER_ARC_REF; do
  [ -n "${!var_name}" ] || handle_error "Environment variable $var_name is required."
done

[ -x "$CELLRANGER_ARC_BIN" ] || handle_error "cellranger-arc not found: $CELLRANGER_ARC_BIN"
[ -d "$CELLRANGER_ARC_REF" ] || handle_error "reference directory not found: $CELLRANGER_ARC_REF"
[ -f "$libraries_csv" ] || handle_error "libraries.csv not found: $libraries_csv"

mkdir -p "$run_id" || handle_error "Directory creation failed"
cd "$run_id" || handle_error "Could not navigate to directory $run_id"

$CELLRANGER_ARC_BIN count \
  --id="$run_id" \
  --reference="$CELLRANGER_ARC_REF" \
  --libraries="$libraries_csv" \
  --localcores="$localcores" \
  --localmem="$localmem" || handle_error "cellranger-arc count failed"

echo "$run_id finished successfully !!!"
