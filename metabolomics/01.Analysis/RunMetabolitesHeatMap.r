#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/MetabolitesAnalysis/03.MetabolitesHeatMap
# Rscript RunMetabolitesHeatMap.r \
#   /path/to/all_sample_data_new_lysate.csv \
#   PSC_ctrl_Lys,PSC_RBP4_Lys,PSC_Rol_Lys,PSC_tx_Lys \
#   Ctrl,RBP4,Rol,Tx \
#   PSC_Lysate \
#   "Retinyl ester|Retinal|11-Cis-Retinol|11-cis-Retinal|Vitamin A acid|Isotretinoin|4-Oxoretinoic acid|all-trans-5,6-Epoxyretinoic acid|Retinoyl beta-glucuronide" \
#   "Vitamin A acid" \
#   "Retinyl ester|Retinal|11-Cis-Retinol|11-cis-Retinal" \
#   60

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(pheatmap)
library(scales)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop(
    paste(
      "Usage: Rscript RunMetabolitesHeatMap.r",
      "/path/to/input.csv",
      "sample_cols_csv",
      "group_labels_csv",
      "output_prefix",
      "label_metabolites_pipe",
      "target_A",
      "target_B_pipe",
      "k_clusters"
    )
  )
}

input_file <- normalizePath(args[1], mustWork = TRUE)
sample_cols <- strsplit(args[2], ",", fixed = TRUE)[[1]]
group_labels <- strsplit(args[3], ",", fixed = TRUE)[[1]]
output_prefix <- args[4]
label_metabolites <- strsplit(args[5], "|", fixed = TRUE)[[1]]
target_A <- args[6]
target_B <- strsplit(args[7], "|", fixed = TRUE)[[1]]
k_clusters <- as.integer(args[8])

if (length(sample_cols) != length(group_labels)) {
  stop("sample_cols_csv and group_labels_csv must have the same length.")
}

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

outdir <- file.path(script_dir, output_prefix)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

Lysate_Data <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("Compounds", sample_cols)
missing_cols <- setdiff(required_cols, colnames(Lysate_Data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

heatmap_df <- Lysate_Data[, sample_cols, drop = FALSE]

Lysate_Data$Compounds[is.na(Lysate_Data$Compounds)] <- "Unknown"
rownames(heatmap_df) <- make.unique(Lysate_Data$Compounds)

heatmap_mat <- as.matrix(data.frame(lapply(heatmap_df, as.numeric)))
rownames(heatmap_mat) <- rownames(heatmap_df)

heatmap_mat <- heatmap_mat[rowSums(is.na(heatmap_mat)) < ncol(heatmap_mat), , drop = FALSE]
log_heatmap_mat <- log2(heatmap_mat + 1)

scaled_mat <- t(apply(log_heatmap_mat, 1, function(x) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  }
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - m) / s)
  }
}))

colnames(scaled_mat) <- colnames(log_heatmap_mat)
rownames(scaled_mat) <- rownames(log_heatmap_mat)
scaled_mat <- scaled_mat[complete.cases(scaled_mat), , drop = FALSE]

col_anno <- data.frame(Group = group_labels)
rownames(col_anno) <- colnames(scaled_mat)

ph_all <- pheatmap(
  scaled_mat,
  annotation_col = col_anno,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 10,
  main = paste0(output_prefix, " Heatmap (All Metabolites)"),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA
)

pdf(paste0(output_prefix, "_all.pdf"), width = 6, height = 12)
pheatmap(
  scaled_mat,
  annotation_col = col_anno,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize = 10,
  main = paste0(output_prefix, " Heatmap (All Metabolites)"),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA
)
dev.off()

present_labels <- intersect(label_metabolites, rownames(scaled_mat))
missing_labels <- setdiff(label_metabolites, rownames(scaled_mat))

write.csv(
  data.frame(Present = present_labels),
  "present_labels.csv",
  row.names = FALSE
)
write.csv(
  data.frame(Missing = missing_labels),
  "missing_labels.csv",
  row.names = FALSE
)

if (length(present_labels) == 0) {
  stop("None of the label metabolites were found in rownames(scaled_mat).")
}

row_labels <- ifelse(rownames(scaled_mat) %in% label_metabolites, rownames(scaled_mat), "")

pdf(paste0(output_prefix, "_retinoid_labels.pdf"), width = 6, height = 12)
ph_retinoid <- pheatmap(
  scaled_mat,
  annotation_col = col_anno,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  labels_row = row_labels,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  main = paste0(output_prefix, " Heatmap (Selected Labeling)"),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA
)
dev.off()

row_tree <- ph_retinoid$tree_row
col_tree <- ph_retinoid$tree_col

saveRDS(row_tree, paste0(output_prefix, "_row_tree.rds"))
saveRDS(col_tree, paste0(output_prefix, "_col_tree.rds"))
saveRDS(ph_retinoid, paste0(output_prefix, "_pheatmap_object.rds"))

pdf(paste0(output_prefix, "_row_dendrogram.pdf"), width = 46, height = 32)
plot(row_tree, main = paste0(output_prefix, " Row Hierarchical Clustering Tree"))
dev.off()

pdf(paste0(output_prefix, "_column_dendrogram.pdf"), width = 15, height = 14)
plot(col_tree, main = paste0(output_prefix, " Column Hierarchical Clustering Tree"))
dev.off()

target_B_present <- intersect(target_B, rownames(scaled_mat))
target_B_missing <- setdiff(target_B, rownames(scaled_mat))

if (!(target_A %in% rownames(scaled_mat))) {
  stop(target_A, " is not present in rownames(scaled_mat).")
}

if (length(target_B_present) == 0) {
  stop("None of the target_B metabolites are present in rownames(scaled_mat).")
}

row_clusters <- cutree(row_tree, k = k_clusters)

cluster_df <- data.frame(
  Metabolite = names(row_clusters),
  Cluster = as.integer(row_clusters),
  stringsAsFactors = FALSE
)
write.csv(cluster_df, paste0(output_prefix, "_row_clusters.csv"), row.names = FALSE)

cluster_A_id <- cluster_df$Cluster[cluster_df$Metabolite == target_A][1]
cluster_B_df <- cluster_df[cluster_df$Metabolite %in% target_B_present, , drop = FALSE]

write.csv(
  data.frame(Metabolite = target_B_present, stringsAsFactors = FALSE),
  "target_B_present.csv",
  row.names = FALSE
)
write.csv(
  data.frame(Metabolite = target_B_missing, stringsAsFactors = FALSE),
  "target_B_missing.csv",
  row.names = FALSE
)
write.csv(
  data.frame(target_A = target_A, cluster_A_id = cluster_A_id),
  "target_A_cluster.csv",
  row.names = FALSE
)
write.csv(cluster_B_df, "target_B_clusters.csv", row.names = FALSE)

cluster_A_metabolites <- cluster_df$Metabolite[cluster_df$Cluster == cluster_A_id]
cluster_B_ids <- unique(cluster_B_df$Cluster)
cluster_B_metabolites <- cluster_df$Metabolite[cluster_df$Cluster %in% cluster_B_ids]

write.csv(
  data.frame(Metabolite = cluster_A_metabolites),
  paste0(output_prefix, "_target_A_cluster_metabolites.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Metabolite = cluster_B_metabolites),
  paste0(output_prefix, "_target_B_cluster_metabolites.csv"),
  row.names = FALSE
)

cluster_A_mat <- scaled_mat[cluster_A_metabolites, , drop = FALSE]
cluster_B_mat <- scaled_mat[cluster_B_metabolites, , drop = FALSE]

write.csv(cluster_A_mat, paste0(output_prefix, "_target_A_cluster_matrix.csv"))
write.csv(cluster_B_mat, paste0(output_prefix, "_target_B_cluster_matrix.csv"))
