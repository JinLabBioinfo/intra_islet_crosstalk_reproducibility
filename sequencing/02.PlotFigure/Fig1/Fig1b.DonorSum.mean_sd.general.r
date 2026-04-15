# Demo usage:
# cd 29.Scripts/02.PlotFigure/Fig1
# Rscript DonorSum.mean_sd.general.r DropSeq
# Rscript DonorSum.mean_sd.general.r PANCDB

# Load environment script when available
optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)

# Setup parameters
args <- commandArgs(trailingOnly = TRUE)

if (!exists("dataset_name")) {
  if (length(args) >= 1) {
    dataset_name <- args[1]
  } else {
    dataset_name <- "DropSeq"
  }
}

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

if (!exists("workPath")) {
  workPath <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
}

if (!exists("output_dir")) {
  output_dir <- normalizePath(script_dir, mustWork = FALSE)
}

if (!exists("input_file")) {
  input_file <- normalizePath(file.path(script_dir, "..", "..", "SourceData", "Combine_DonorID.csv"), mustWork = FALSE)
}

if (!exists("summary_csv")) {
  summary_csv <- file.path(output_dir, paste0(dataset_name, ".summary_stats.mean_sd.csv"))
}

variables_to_summarize <- c("Age", "BMI", "HbA1C")

safe_as_numeric <- function(x) {
  as.numeric(trimws(as.character(x)))
}

calculate_mean_sd <- function(data, column_name) {
  values <- data[[column_name]]
  c(
    Mean = mean(values, na.rm = TRUE),
    SD = sd(values, na.rm = TRUE),
    N = sum(!is.na(values))
  )
}

# Read and prepare donor table
Total_DonorInfo <- read.csv(input_file, stringsAsFactors = FALSE)
colnames(Total_DonorInfo)[1] <- "Sample"

donor_data <- Total_DonorInfo %>%
  filter(DataSets == dataset_name) %>%
  mutate(
    Sample = trimws(Sample),
    Disease_Status = trimws(Disease_Status),
    Sex = trimws(Sex),
    BMI = safe_as_numeric(BMI),
    Age = safe_as_numeric(Age),
    HbA1C = safe_as_numeric(HbA1C)
  )

data_N <- donor_data %>% filter(Disease_Status == "N")
data_Y <- donor_data %>% filter(Disease_Status == "Y")

# Mean +/- SD summary table
stats_N <- do.call(
  rbind,
  lapply(variables_to_summarize, function(var_name) {
    stats <- calculate_mean_sd(data_N, var_name)
    data.frame(
      Variable = var_name,
      Mean_N = stats["Mean"],
      SD_N = stats["SD"],
      N_N = stats["N"],
      row.names = NULL
    )
  })
)

stats_Y <- do.call(
  rbind,
  lapply(variables_to_summarize, function(var_name) {
    stats <- calculate_mean_sd(data_Y, var_name)
    data.frame(
      Variable = var_name,
      Mean_Y = stats["Mean"],
      SD_Y = stats["SD"],
      N_Y = stats["N"],
      row.names = NULL
    )
  })
)

summary_stats <- merge(stats_N, stats_Y, by = "Variable")
print(summary_stats)
write.csv(summary_stats, summary_csv, row.names = FALSE)
