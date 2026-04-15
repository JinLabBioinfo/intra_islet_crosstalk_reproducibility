#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/MetabolitesAnalysis/01.VenmCode
# Rscript GetSeparateVenn.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(VennDiagram)
library(grid)

setwd("/path/to/29.Scripts/03.Supporting/MetabolitesAnalysis/01.VenmCode")

PSC_Medium_Data <- read.csv("/path/to/all_sample_data_new_medium.csv")
PSC_Lysate_Data <- read.csv("/path/to/all_sample_data_new_lysate.csv")
ENDOC_Medium_Data <- read.csv("/path/to/EndoC_all_sample_data_new_medium.csv")
ENDOC_Lysate_Data <- read.csv("/path/to/EndoC_all_sample_data_new_lysate.csv")

PSC_Medium_subset <- PSC_Medium_Data
PSC_Lysate_subset <- PSC_Lysate_Data
ENDOC_Medium_subset <- ENDOC_Medium_Data
ENDOC_Lysate_subset <- ENDOC_Lysate_Data

Medium_compare <- merge(
  PSC_Medium_subset,
  ENDOC_Medium_subset,
  by = c("Compounds", "Class.I", "PubChem.CID"),
  suffixes = c("_PSC", "_ENDOC")
)

Lysate_compare <- merge(
  PSC_Lysate_subset,
  ENDOC_Lysate_subset,
  by = c("Compounds", "Class.I", "PubChem.CID"),
  suffixes = c("_PSC", "_ENDOC")
)

Medium_common <- length(intersect(unique(PSC_Medium_subset$Compounds), unique(ENDOC_Medium_subset$Compounds)))
Medium_psc_only <- length(setdiff(unique(PSC_Medium_subset$Compounds), unique(ENDOC_Medium_subset$Compounds)))
Medium_endoc_only <- length(setdiff(unique(ENDOC_Medium_subset$Compounds), unique(PSC_Medium_subset$Compounds)))

Lysate_common <- length(intersect(unique(PSC_Lysate_subset$Compounds), unique(ENDOC_Lysate_subset$Compounds)))
Lysate_psc_only <- length(setdiff(unique(PSC_Lysate_subset$Compounds), unique(ENDOC_Lysate_subset$Compounds)))
Lysate_endoc_only <- length(setdiff(unique(ENDOC_Lysate_subset$Compounds), unique(PSC_Lysate_subset$Compounds)))

PSC_Medium_compounds <- unique(PSC_Medium_subset$Compounds)
ENDOC_Medium_compounds <- unique(ENDOC_Medium_subset$Compounds)
PSC_Lysate_compounds <- unique(PSC_Lysate_subset$Compounds)
ENDOC_Lysate_compounds <- unique(ENDOC_Lysate_subset$Compounds)

Medium_area1 <- length(PSC_Medium_compounds)
Medium_area2 <- length(ENDOC_Medium_compounds)
Medium_cross <- length(intersect(PSC_Medium_compounds, ENDOC_Medium_compounds))

Lysate_area1 <- length(PSC_Lysate_compounds)
Lysate_area2 <- length(ENDOC_Lysate_compounds)
Lysate_cross <- length(intersect(PSC_Lysate_compounds, ENDOC_Lysate_compounds))

Venn_summary <- data.frame(
  Group = c("Medium", "Lysate"),
  Common = c(Medium_common, Lysate_common),
  PSC_only = c(Medium_psc_only, Lysate_psc_only),
  ENDOC_only = c(Medium_endoc_only, Lysate_endoc_only),
  stringsAsFactors = FALSE
)

write.csv(Venn_summary, "Venn_summary.csv", row.names = FALSE)

Medium_common_compounds <- intersect(PSC_Medium_compounds, ENDOC_Medium_compounds)
Medium_psc_only_compounds <- setdiff(PSC_Medium_compounds, ENDOC_Medium_compounds)
Medium_endoc_only_compounds <- setdiff(ENDOC_Medium_compounds, PSC_Medium_compounds)

Lysate_common_compounds <- intersect(PSC_Lysate_compounds, ENDOC_Lysate_compounds)
Lysate_psc_only_compounds <- setdiff(PSC_Lysate_compounds, ENDOC_Lysate_compounds)
Lysate_endoc_only_compounds <- setdiff(ENDOC_Lysate_compounds, PSC_Lysate_compounds)

add_region_info <- function(df, region_name, region_count) {
  out <- df
  out$Region <- region_name
  out$Region_Count <- region_count
  out
}

Medium_common_PSC <- add_region_info(
  PSC_Medium_subset[PSC_Medium_subset$Compounds %in% Medium_common_compounds, , drop = FALSE],
  "Medium_common",
  length(Medium_common_compounds)
)
Medium_common_ENDOC <- add_region_info(
  ENDOC_Medium_subset[ENDOC_Medium_subset$Compounds %in% Medium_common_compounds, , drop = FALSE],
  "Medium_common",
  length(Medium_common_compounds)
)
Medium_psc_only_df <- add_region_info(
  PSC_Medium_subset[PSC_Medium_subset$Compounds %in% Medium_psc_only_compounds, , drop = FALSE],
  "Medium_PSC_only",
  length(Medium_psc_only_compounds)
)
Medium_endoc_only_df <- add_region_info(
  ENDOC_Medium_subset[ENDOC_Medium_subset$Compounds %in% Medium_endoc_only_compounds, , drop = FALSE],
  "Medium_ENDOC_only",
  length(Medium_endoc_only_compounds)
)

Lysate_common_PSC <- add_region_info(
  PSC_Lysate_subset[PSC_Lysate_subset$Compounds %in% Lysate_common_compounds, , drop = FALSE],
  "Lysate_common",
  length(Lysate_common_compounds)
)
Lysate_common_ENDOC <- add_region_info(
  ENDOC_Lysate_subset[ENDOC_Lysate_subset$Compounds %in% Lysate_common_compounds, , drop = FALSE],
  "Lysate_common",
  length(Lysate_common_compounds)
)
Lysate_psc_only_df <- add_region_info(
  PSC_Lysate_subset[PSC_Lysate_subset$Compounds %in% Lysate_psc_only_compounds, , drop = FALSE],
  "Lysate_PSC_only",
  length(Lysate_psc_only_compounds)
)
Lysate_endoc_only_df <- add_region_info(
  ENDOC_Lysate_subset[ENDOC_Lysate_subset$Compounds %in% Lysate_endoc_only_compounds, , drop = FALSE],
  "Lysate_ENDOC_only",
  length(Lysate_endoc_only_compounds)
)

write.csv(Medium_common_PSC, paste0("Medium_common_n", length(Medium_common_compounds), "_PSC.csv"), row.names = FALSE)
write.csv(Medium_common_ENDOC, paste0("Medium_common_n", length(Medium_common_compounds), "_ENDOC.csv"), row.names = FALSE)
write.csv(Medium_psc_only_df, paste0("Medium_PSC_only_n", length(Medium_psc_only_compounds), ".csv"), row.names = FALSE)
write.csv(Medium_endoc_only_df, paste0("Medium_ENDOC_only_n", length(Medium_endoc_only_compounds), ".csv"), row.names = FALSE)

write.csv(Lysate_common_PSC, paste0("Lysate_common_n", length(Lysate_common_compounds), "_PSC.csv"), row.names = FALSE)
write.csv(Lysate_common_ENDOC, paste0("Lysate_common_n", length(Lysate_common_compounds), "_ENDOC.csv"), row.names = FALSE)
write.csv(Lysate_psc_only_df, paste0("Lysate_PSC_only_n", length(Lysate_psc_only_compounds), ".csv"), row.names = FALSE)
write.csv(Lysate_endoc_only_df, paste0("Lysate_ENDOC_only_n", length(Lysate_endoc_only_compounds), ".csv"), row.names = FALSE)

pdf("Venn_Medium.pdf", width = 7, height = 7)
grid.newpage()
draw.pairwise.venn(
  area1 = Medium_area1,
  area2 = Medium_area2,
  cross.area = Medium_cross,
  category = c("PSC Medium", "ENDOC Medium"),
  fill = c("skyblue", "salmon"),
  alpha = c(0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)
dev.off()

pdf("Venn_Lysate.pdf", width = 7, height = 7)
grid.newpage()
draw.pairwise.venn(
  area1 = Lysate_area1,
  area2 = Lysate_area2,
  cross.area = Lysate_cross,
  category = c("PSC Lysate", "ENDOC Lysate"),
  fill = c("lightgreen", "orange"),
  alpha = c(0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)
dev.off()
