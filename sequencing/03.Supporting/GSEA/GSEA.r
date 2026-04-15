#!/usr/bin/env Rscript

# Demo usage:
# cd /29.Scripts/03.Supporting/GSEA_Figure2
# Rscript GSEA_Figure2.r Beta_T2D /path/to/msigdb.v2024.1.Hs.symbols.gmt


args <- commandArgs(trailingOnly = TRUE)
analysis_name <- if (length(args) >= 1) args[1] else "Beta_T2D"

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

scripts_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
source_data_root <- file.path(scripts_root, "SourceData", "GetStatisticsDonorScore")

analysis_map <- data.frame(
  Analysis = c(
    "Beta_T2D",
    "Beta_NonT2D_Age",
    "Beta_NonT2D_BMI",
    "PSC_T2D",
    "PSC_NonT2D_Age",
    "PSC_NonT2D_BMI"
  ),
  InputFile = c(
    file.path(source_data_root, "00.RawObj", "Disease_Status_Beta_All.rds"),
    file.path(source_data_root, "01.GSEA", "Age_Beta_Healthy.rds"),
    file.path(source_data_root, "01.GSEA", "BMI_Beta_Healthy.rds"),
    file.path(source_data_root, "00.RawObj", "Disease_Status_PSC_All.rds"),
    file.path(source_data_root, "01.GSEA", "Age_PSC_Healthy.rds"),
    file.path(source_data_root, "01.GSEA", "BMI_PSC_Healthy.rds")
  ),
  SampleName = c(
    "Beta_T2D_Signatures",
    "Beta_Aging_Healthy_Signatures",
    "Beta_Aging_T2D_Signatures",
    "PSC_Aging_Healthy_Signatures",
    "PSC_Aging_Healthy_Signatures",
    "PSC_Aging_T2D_Signatures"
  ),
  stringsAsFactors = FALSE
)

if (tolower(analysis_name) == "list") {
  print(analysis_map)
  quit(save = "no", status = 0)
}

if (!analysis_name %in% analysis_map$Analysis) {
  stop("Unsupported analysis: ", analysis_name)
}

analysis_row <- analysis_map[analysis_map$Analysis == analysis_name, , drop = FALSE]
input_file <- analysis_row$InputFile[1]
sample_name <- analysis_row$SampleName[1]

if (!file.exists(input_file)) {
  stop("Input RDS not found: ", input_file)
}

gmt_candidates <- c(
  if (length(args) >= 2) args[2] else NA_character_,
  Sys.getenv("MSIGDB_GMT_FILE", "")
)
gmt_candidates <- gmt_candidates[!is.na(gmt_candidates) & nzchar(gmt_candidates)]
gmt_matches <- gmt_candidates[file.exists(gmt_candidates)]
gmt_file <- if (length(gmt_matches) > 0) gmt_matches[1] else ""

if (!nzchar(gmt_file)) {
  stop("MSigDB GMT file not found. Pass it as the second argument or set MSIGDB_GMT_FILE.")
}

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(fgsea)
library(GSEABase)
library(data.table)

output_dir <- file.path(script_dir, analysis_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)

Obj <- readRDS(input_file)
Signatures <- Obj[["betaT2D.diffGene.20bin.PCA"]]$pseudoregress.all

rankings <- sign(Signatures$slope) * (-log10(Signatures$qvalue))
names(rankings) <- rownames(Signatures)
rankings <- sort(rankings, decreasing = TRUE)

Ranklist <- data.frame(rankings)
mytime <- format(Sys.time(), "%b_%d_%Y")
myfile <- file.path(paste0(sample_name, "_", mytime, ".rnk"))
write.table(
  Ranklist,
  file = myfile,
  sep = "\t",
  row.names = TRUE,
  col.names = FALSE,
  quote = FALSE,
  append = FALSE
)

hallmark_pathways <- getGmt(gmt_file)
pathways_list <- lapply(hallmark_pathways, geneIds)
names(pathways_list) <- sapply(hallmark_pathways, function(x) x@setName)

fgseaRes <- fgsea(pathways = pathways_list, stats = rankings)
fgseaRes_ordered <- fgseaRes[order(NES)]
fgseaRes_ordered[, direction := ifelse(NES > 0, "Up", "Down")]

fwrite(fgseaRes_ordered, file = "fgsea_results_ordered.csv")
saveRDS(fgseaRes, file = "fgsea_results_hallmark_pathways.rds")

top_pathways <- fgseaRes[order(fgseaRes$padj), ][1:10, ]
top_pathways$direction <- ifelse(top_pathways$NES > 0, "Up", "Down")

p <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue")) +
  labs(
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)",
    title = paste("Top Enriched Hallmark Pathways in", analysis_name)
  )

ggsave("Top10_Hallmark_Pathways.pdf", p, width = 8, height = 6)
