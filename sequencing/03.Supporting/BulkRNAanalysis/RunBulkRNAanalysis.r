#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/BulkRNAanalysis
# Rscript RunBulkRNAanalysis.r /path/to/count_files
#
# Expected files in input_dir:
# Ctrl_KeepDu_BachOne.count.bed
# Ctrl_KeepDu_NEW.count.bed
# Retinol_KeepDu_BatchOne.count.bed
# Retinol_KeepDu_New.count.bed

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(tidyverse)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript RunBulkRNAanalysis.r /path/to/count_files")
}

input_dir <- normalizePath(args[1], mustWork = TRUE)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)

Ctrl_Islets_1 <- read.table(
  file = file.path(input_dir, "Ctrl_KeepDu_BachOne.count.bed"),
  header = TRUE, as.is = TRUE, sep = "\t", stringsAsFactors = FALSE
)
Ctrl_Islets_2 <- read.table(
  file = file.path(input_dir, "Ctrl_KeepDu_NEW.count.bed"),
  header = TRUE, as.is = TRUE, sep = "\t", stringsAsFactors = FALSE
)
Retinol_Islets_1 <- read.table(
  file = file.path(input_dir, "Retinol_KeepDu_BatchOne.count.bed"),
  header = TRUE, as.is = TRUE, sep = "\t", stringsAsFactors = FALSE
)
Retinol_Islets_2 <- read.table(
  file = file.path(input_dir, "Retinol_KeepDu_New.count.bed"),
  header = TRUE, as.is = TRUE, sep = "\t", stringsAsFactors = FALSE
)

merge_counts <- merge(Ctrl_Islets_1, Ctrl_Islets_2)
merge_counts <- merge(merge_counts, Retinol_Islets_1)
merge_counts <- merge(merge_counts, Retinol_Islets_2)

colnames(merge_counts) <- c(
  "Geneid",
  "Ctrl_Islets_1",
  "Ctrl_Islets_2",
  "Retinol_Islets_1",
  "Retinol_Islets_2"
)

write.csv(merge_counts, "HumanIslets_Sorted_CtrlRetinol.csv", row.names = FALSE)

stopifnot(all(c("Ctrl_Islets_1", "Ctrl_Islets_2", "Retinol_Islets_1", "Retinol_Islets_2") %in% colnames(merge_counts)))

counts_mat <- as.matrix(merge_counts[, c("Ctrl_Islets_1", "Ctrl_Islets_2", "Retinol_Islets_1", "Retinol_Islets_2")])
rownames(counts_mat) <- merge_counts$Geneid

metadata <- data.frame(
  id = colnames(counts_mat),
  dex = c("Ctrl", "Ctrl", "Retinol", "Retinol"),
  batch = c("batch1", "batch2", "batch1", "batch2"),
  celltype = "Islets",
  geo_id = ""
)
metadata$dex <- factor(metadata$dex, levels = c("Ctrl", "Retinol"))
metadata$batch <- factor(metadata$batch)
rownames(metadata) <- metadata$id

stopifnot(identical(colnames(counts_mat), rownames(metadata)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = metadata,
  design = ~ batch + dex
)
dds <- DESeq(dds)

res <- results(dds, contrast = c("dex", "Retinol", "Ctrl"))

res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)
res_df <- res_df[, c("Geneid", setdiff(names(res_df), "Geneid"))]

write.csv(
  res_df[order(res_df$pvalue), ],
  "DESeq2_Retinol_vs_Ctrl.batchAdjusted.results.csv",
  row.names = FALSE
)

res_all_no_na <- na.omit(as.data.frame(res_df))
write.csv(
  res_all_no_na,
  "WihoutNA.DESeq2_Retinol_vs_Ctrl.batchAdjusted.results.csv",
  quote = FALSE,
  row.names = FALSE
)

res_all_no_na$significance <- with(
  res_all_no_na,
  ifelse(padj < 0.05 & abs(log2FoldChange) > 0.5, "Significant", "Not Significant")
)

res_all_no_na$color <- with(
  res_all_no_na,
  ifelse(
    significance == "Not Significant",
    "grey",
    ifelse(log2FoldChange < -0.5, "#158BB8", ifelse(log2FoldChange > 0.5, "red", "grey"))
  )
)

write.csv(
  res_all_no_na,
  "Sig.HumanIslets_Sorted_CtrlRetinol.DEseq2.all.noNA.csv",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(dds, "HumanIslets_Sorted_CtrlRetinol.batchAdjusted.DESeq2.rds")
