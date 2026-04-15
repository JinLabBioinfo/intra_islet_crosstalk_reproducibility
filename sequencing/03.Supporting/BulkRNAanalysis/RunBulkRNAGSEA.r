#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/BulkRNAanalysis
# Rscript RunBulkRNAGSEA.r /path/to/Sig.HumanIslets_Sorted_CtrlRetinol.DEseq2.all.noNA.csv /path/to/msigdb.v2024.1.Hs.symbols.gmt
#
# Or:
# export MSIGDB_GMT_FILE=/path/to/msigdb.v2024.1.Hs.symbols.gmt
# Rscript RunBulkRNAGSEA.r /path/to/Sig.HumanIslets_Sorted_CtrlRetinol.DEseq2.all.noNA.csv

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(ggplot2)
library(dplyr)
library(fgsea)
library(GSEABase)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript RunBulkRNAGSEA.r /path/to/Sig.HumanIslets_Sorted_CtrlRetinol.DEseq2.all.noNA.csv [/path/to/msigdb.v2024.1.Hs.symbols.gmt]")
}

input_file <- normalizePath(args[1], mustWork = TRUE)
gmt_file <- if (length(args) >= 2) args[2] else Sys.getenv("MSIGDB_GMT_FILE", "")

if (gmt_file == "") {
  stop("Provide the MSigDB GMT file as the second argument or set MSIGDB_GMT_FILE.")
}

gmt_file <- normalizePath(gmt_file, mustWork = TRUE)

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

setwd(script_dir)

T2DSignatures <- read.csv(input_file, row.names = 1, stringsAsFactors = FALSE)

rankings <- sign(T2DSignatures$log2FoldChange) * (-log10(T2DSignatures$pvalue))
names(rankings) <- rownames(T2DSignatures)
rankings <- sort(rankings, decreasing = TRUE)

hallmark_pathways <- getGmt(gmt_file)
pathways_list <- lapply(hallmark_pathways, geneIds)
names(pathways_list) <- sapply(hallmark_pathways, function(x) x@setName)

fgseaRes <- fgsea(pathways = pathways_list, stats = rankings)
fgseaRes_ordered <- fgseaRes[order(NES)]
fgseaRes_ordered[, direction := ifelse(NES > 0, "Up", "Down")]

fwrite(fgseaRes_ordered, file = "fgsea_results_ordered.csv")

Ranklist <- data.frame(rankings)
SampleName <- "BulkRNA_T2D_Signatures"
mytime <- format(Sys.time(), "%b_%d_%Y")
myfile <- file.path(paste0(SampleName, "_", mytime, ".rnk"))
write.table(
  Ranklist,
  file = myfile,
  sep = "\t",
  row.names = TRUE,
  col.names = FALSE,
  quote = FALSE,
  append = FALSE
)

saveRDS(fgseaRes, file = "fgsea_results.rds")
