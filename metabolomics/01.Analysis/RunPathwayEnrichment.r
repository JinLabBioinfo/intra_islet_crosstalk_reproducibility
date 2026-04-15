#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/MetabolitesAnalysis/02.PathwayEnrichment
# Rscript RunPathwayEnrichment.r /path/to/Lysate_ENDOC_only_n627.csv 01.Lysate_EndoC_Only_analytes

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(RaMP)
library(DT)
library(dplyr)
library(magrittr)
library(tidytext)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript RunPathwayEnrichment.r /path/to/analyte_table.csv output_folder_name")
}

input_file <- normalizePath(args[1], mustWork = TRUE)
output_name <- args[2]

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)
outdir <- file.path(script_dir, output_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rampDB <- RaMP()
listAvailableRaMPDbVersions()

input_df <- read.csv(input_file, stringsAsFactors = FALSE)

input_subset <- input_df[, c("Compounds", "PubChem.CID", "HMDB", "CAS", "ChEBI"), drop = FALSE]

Map <- c(
  "Compounds" = "Compounds",
  "PubChem.CID" = "pubchem",
  "HMDB" = "hmdb",
  "CAS" = "CAS",
  "ChEBI" = "chebi"
)

colnames(input_subset) <- Map[colnames(input_subset)]

metabprefixes <- getPrefixesFromAnalytes("metabolite", db = rampDB)
geneprefixes <- getPrefixesFromAnalytes("gene", db = rampDB)
print(rbind(metabprefixes, geneprefixes))

id_priority <- c("pubchem", "hmdb", "CAS", "chebi")

input_subset[input_subset == "-"] <- NA
input_subset[input_subset == ""] <- NA

pick_first_id <- function(x, id_priority) {
  idx <- which(!is.na(x))[1]
  if (length(idx) == 0 || is.na(idx)) {
    return(NA_character_)
  }
  paste0(id_priority[idx], ":", x[idx])
}

myanalytes_all <- apply(
  input_subset[, id_priority, drop = FALSE],
  1,
  pick_first_id,
  id_priority = id_priority
)

selected_id_type <- apply(
  input_subset[, id_priority, drop = FALSE],
  1,
  function(x) {
    idx <- which(!is.na(x))[1]
    if (length(idx) == 0 || is.na(idx)) {
      return(NA_character_)
    }
    id_priority[idx]
  }
)

myanalytes <- unique(myanalytes_all[!is.na(myanalytes_all)])
dropped_compounds <- input_subset$Compounds[is.na(myanalytes_all)]

n_total <- nrow(input_subset)
n_found <- sum(!is.na(myanalytes_all))
n_dropped <- sum(is.na(myanalytes_all))

selection_report <- data.frame(
  total_input_rows = n_total,
  rows_with_selected_analyte = n_found,
  rows_dropped_no_supported_id = n_dropped,
  unique_analytes = length(myanalytes)
)

selection_by_type <- as.data.frame(table(selected_id_type, useNA = "ifany"))
colnames(selection_by_type) <- c("selected_id_type", "Count")

dropped_rows <- input_subset[is.na(myanalytes_all), ]

check_df <- data.frame(
  Compounds = input_subset$Compounds,
  pubchem = input_subset$pubchem,
  hmdb = input_subset$hmdb,
  CAS = input_subset$CAS,
  chebi = input_subset$chebi,
  selected_id_type = selected_id_type,
  selected_analyte = myanalytes_all
)

pathways.enriched <- runEnrichPathways(analytes = myanalytes, db = rampDB)
filtered.pathways.enriched <- filterEnrichResults(
  enrichResults = pathways.enriched,
  pValType = "fdr",
  pValCutoff = 0.05
)

clusters <- findCluster(
  filtered.pathways.enriched,
  percAnalyteOverlap = 0.2,
  percPathwayOverlap = 0.2,
  db = rampDB
)

objects_to_save <- list(
  myanalytes = myanalytes,
  pathways.enriched = pathways.enriched,
  filtered.pathways.enriched = filtered.pathways.enriched,
  clusters = clusters
)

for (obj_name in names(objects_to_save)) {
  saveRDS(
    objects_to_save[[obj_name]],
    file = file.path(outdir, paste0(obj_name, ".rds"))
  )
}

write.csv(
  filtered.pathways.enriched$fishresults,
  file = file.path(outdir, "filtered.pathways.enriched_fishresults.csv"),
  row.names = FALSE
)

write.csv(
  pathways.enriched$fishresults,
  file = file.path(outdir, "total.pathways.enriched_fishresults.csv"),
  row.names = FALSE
)

write.csv(
  selection_report,
  file = file.path(outdir, "selection_report.csv"),
  row.names = FALSE
)

write.csv(
  selection_by_type,
  file = file.path(outdir, "selection_by_type.csv"),
  row.names = FALSE
)

write.csv(
  dropped_rows,
  file = file.path(outdir, "dropped_rows.csv"),
  row.names = FALSE
)

write.csv(
  check_df,
  file = file.path(outdir, "analyte_selection_check.csv"),
  row.names = FALSE
)
