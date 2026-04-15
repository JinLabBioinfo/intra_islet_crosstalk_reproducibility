# Author: Jiachen Sun

#!/bin/bash
R_SCRIPT=$(mktemp --suffix=.r)

# The text between 'EOF' and 'EOF' is written to the temp R file
cat << 'EOF' > "$R_SCRIPT"
library(argparse)

source("/home/jxs2269/islet_cosmx/for_jasmine/src/pipeline_reproducibility_common/_utility_funs.r")

parser <- ArgumentParser()
parser$add_argument("--input_dir", type = "character")
parser$add_argument("--prefix", type = "character")
parser$add_argument("--output_file", type = "character")
args <- parser$parse_args()

SaveSeuratRds(load_nanostring(args$input_dir, prefix = args$prefix), args$output_file)
EOF

# ND sample 1
Rscript "$R_SCRIPT" \
    --input_dir /home/jxs2269/islet_cosmx/nanostring_outs/RQ25389-002_AK090T_TBH3/flatFiles/YaniLiu_AK090T_TBH3 \
    --prefix AK090T \
    --output_file /home/jxs2269/islet_cosmx/for_jasmine/obj_raw_AK090T.rds

# T2D sample 1
Rscript "$R_SCRIPT" \
    --input_dir /home/jxs2269/islet_cosmx/nanostring_outs/RQ25725-001_AL189T/flatFiles/RQ25725001_AL189T \
    --prefix AL189T \
    --output_file /home/jxs2269/islet_cosmx/for_jasmine/obj_raw_AL189T.rds

# T2D sample 2
Rscript "$R_SCRIPT" \
    --input_dir /home/jxs2269/islet_cosmx/nanostring_outs/RQ25725-008_AM010/flatFiles/RQ25725008_AM010 \
    --prefix AM010 \
    --output_file /home/jxs2269/islet_cosmx/for_jasmine/obj_raw_AM010.rds

rm "$R_SCRIPT"
