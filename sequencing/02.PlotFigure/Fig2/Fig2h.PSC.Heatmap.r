#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig2
# Rscript Fig2h.PSC.Heatmap.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(pheatmap)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_root <- normalizePath(file.path(script_dir, "..", "..", "SourceData"), mustWork = FALSE)
raw_obj_dir <- file.path(source_data_root, "GetStatisticsDonorScore")
gene_dir <- file.path(source_data_root, "Fig2e_h_Heatmap")

t2d_file <- file.path(raw_obj_dir, "00.RawObj", "Disease_Status_PSC_All.rds")
aging_file <- file.path(raw_obj_dir, "01.GSEA", "Age_PSC_Healthy.rds")
bmi_file <- file.path(raw_obj_dir, "01.GSEA", "BMI_PSC_Healthy.rds")
genes_up_file <- file.path(gene_dir, "Up.PSC.genes")
genes_dn_file <- file.path(gene_dir, "Dn.PSC.genes")

setwd(script_dir)

normalize_data <- function(data_subset) {
  normalize_01 <- function(vector) {
    (vector - min(vector)) / (max(vector) - min(vector))
  }

  data_subset_scaled <- t(apply(data_subset, 1, normalize_01))
  data_subset_scaled <- as.data.frame(data_subset_scaled)
  rownames(data_subset_scaled) <- rownames(data_subset)
  colnames(data_subset_scaled) <- colnames(data_subset)
  data_subset_scaled
}

T2D_Status <- readRDS(t2d_file)
Aging <- readRDS(aging_file)
BMI <- readRDS(bmi_file)

T2D_Count_Table <- T2D_Status[["beta.RNA.PCA.20bin.ob"]][["cellvsPeak.m.aggr"]]
Aging_Table <- Aging[["beta.RNA.PCA.20bin.ob"]][["cellvsPeak.m.aggr"]]
BMI_Table <- BMI[["beta.RNA.PCA.20bin.ob"]][["cellvsPeak.m.aggr"]]

T2D_Count_Table_Normalized <- normalize_data(T2D_Count_Table)
Aging_Table_Normalized <- normalize_data(Aging_Table)
BMI_Table_Normalized <- normalize_data(BMI_Table)

top_genes <- c(readLines(genes_up_file), readLines(genes_dn_file))

T2D_Count_Table_Normalized_top <- T2D_Count_Table_Normalized[top_genes, ]
Aging_Table_Normalized_top <- Aging_Table_Normalized[top_genes, ]
BMI_Table_Normalized_top <- BMI_Table_Normalized[top_genes, ]

colnames(T2D_Count_Table_Normalized_top) <- paste0("T2D_", colnames(T2D_Count_Table_Normalized_top))
colnames(Aging_Table_Normalized_top) <- paste0("Aging_", colnames(Aging_Table_Normalized_top))
colnames(BMI_Table_Normalized_top) <- paste0("BMI_", colnames(BMI_Table_Normalized_top))

T2D_Count_Table_Normalized_top <- as.data.frame(T2D_Count_Table_Normalized_top)
Aging_Table_Normalized_top <- as.data.frame(Aging_Table_Normalized_top)
BMI_Table_Normalized_top <- as.data.frame(BMI_Table_Normalized_top)

combined_data <- merge(
  T2D_Count_Table_Normalized_top,
  Aging_Table_Normalized_top,
  by = "row.names",
  all = TRUE
)
rownames(combined_data) <- combined_data$Row.names
combined_data$Row.names <- NULL

combined_data <- merge(
  combined_data,
  BMI_Table_Normalized_top,
  by = "row.names",
  all = TRUE
)
rownames(combined_data) <- combined_data$Row.names
combined_data$Row.names <- NULL

combined_data_ordered <- combined_data[top_genes, ]

row_dir <- data.frame(Direction = c(rep("Up", 15), rep("Down", nrow(combined_data_ordered) - 15)))
rownames(row_dir) <- rownames(combined_data_ordered)

ann_colors <- list(Direction = c(Up = "#d73027", Down = "#4575b4"))

pdf(file = "Fig2h.PSC.Heatmap.pdf", width = 5, height = 5)
pheatmap(
  combined_data_ordered,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  show_colnames = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  annotation_row = row_dir,
  annotation_colors = ann_colors,
  gaps_col = c(20, 40),
  gaps_row = 15,
  border_color = NA
)
dev.off()
