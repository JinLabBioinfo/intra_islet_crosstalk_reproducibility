#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/GetStatisticsDonorScore
# Rscript GetStatisticsDonorScore.r

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "GetStatisticsDonorScore"), mustWork = FALSE)
raw_obj_dir <- file.path(source_data_dir, "00.RawObj")
donor_id_file <- file.path(source_data_dir, "Combine_DonorID.csv")

setwd(script_dir)

cell_types <- c("Beta", "Alpha", "Delta", "PP", "PSC", "Endothelial", "Duct", "Acinar")

all_summary <- data.frame(
  CellType = character(),
  KS_Test_p_value = numeric(),
  D_value = numeric(),
  stringsAsFactors = FALSE
)

CombineDonorID <- read.csv(donor_id_file, stringsAsFactors = FALSE)
colnames(CombineDonorID)[1] <- "Sample"
CombineDonorID$Sample <- trimws(CombineDonorID$Sample)
CombineDonorID <- CombineDonorID[, c(1, 8)]

for (CellType in cell_types) {
  cat("\nProcessing CellType:", CellType, "\n")

  obj_path <- file.path(raw_obj_dir, paste0("Disease_Status_", CellType, "_All.rds"))

  if (!file.exists(obj_path)) {
    cat("WARNING: File not found for", CellType, "- Skipping.\n")
    next
  }

  Obj <- readRDS(obj_path)
  DonorData <- Obj[["beta.RNA.PCA.20bin.ob"]][["data.info.withbin"]]
  DonorData <- merge(DonorData, CombineDonorID, by = "Sample", all.x = TRUE)

  if (!"pseudo.index.balanced" %in% colnames(DonorData)) {
    cat("WARNING: pseudo.index.balanced missing in", CellType, "- Skipping.\n")
    next
  }

  DonorData$pseudo.index.balanced.scaled <- (DonorData$pseudo.index.balanced - min(DonorData$pseudo.index.balanced)) /
    (max(DonorData$pseudo.index.balanced) - min(DonorData$pseudo.index.balanced))

  T2D_Risk <- DonorData %>%
    group_by(Donor) %>%
    summarize(Beta = median(pseudo.index.balanced.scaled, na.rm = TRUE), .groups = "drop")

  T2D_Risk$Group <- ifelse(grepl("T2D", T2D_Risk$Donor), "T2D", "Healthy")

  if (nrow(T2D_Risk) == 0) {
    cat("WARNING: No valid data for", CellType, "- Skipping.\n")
    next
  }

  ks_result <- ks.test(
    T2D_Risk$Beta[T2D_Risk$Group == "Healthy"],
    T2D_Risk$Beta[T2D_Risk$Group == "T2D"]
  )

  all_summary <- rbind(all_summary, data.frame(
    CellType = CellType,
    KS_Test_p_value = ks_result$p.value,
    D_value = ks_result$statistic
  ))
}

write.csv(all_summary, "All_CellTypes_Summary.csv", row.names = FALSE)
