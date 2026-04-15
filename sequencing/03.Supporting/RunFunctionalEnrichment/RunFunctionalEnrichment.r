#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/RunFunctionalEnrichment
# Rscript RunFunctionalEnrichment.r Up /path/to/msigdb.v2024.1.Hs.symbols.gmt
# Rscript RunFunctionalEnrichment.r Dn /path/to/msigdb.v2024.1.Hs.symbols.gmt

optional_env_script <- Sys.getenv("DONOR_SUM_ENV_FILE", "")
if (file.exists(optional_env_script)) {
  source(optional_env_script)
}

library(dplyr)
library(qvalue)
library(fgsea)
library(GSA)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)

if (!exists("enrichment_mode")) {
  enrichment_mode <- if (length(args) >= 1) args[1] else "Up"
}

if (!(enrichment_mode %in% c("Up", "Dn"))) {
  stop("enrichment_mode must be either 'Up' or 'Dn'.")
}

script_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_file_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
} else {
  script_dir <- getwd()
}

source_data_dir <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData", "RunFunctionalEnrichment", enrichment_mode),
  mustWork = FALSE
)
shared_source_data_dir <- normalizePath(
  file.path(script_dir, "..", "..", "SourceData", "RunFunctionalEnrichment"),
  mustWork = FALSE
)

helper_file <- file.path(script_dir, "Functional.Combine.Fun.r")
if (!file.exists(helper_file)) {
  stop("Missing helper file: ", helper_file)
}
source(helper_file)

all.genes.refer <- read.table(
  file.path(shared_source_data_dir, "all.genes.refer.txt"),
  stringsAsFactors = FALSE
)[, 1]

initialize_environment <- function(workPath, DataName) {
  setwd(workPath)
  new_directory <- file.path(workPath, DataName)
  if (!dir.exists(new_directory)) {
    dir.create(new_directory)
  }
  setwd(new_directory)
  cat("Current working directory:", getwd(), "\n")
  return(new_directory)
}

workPath <- script_dir
DataName <- paste0("02.GSEA.Total.", enrichment_mode)
initialize_environment(workPath, DataName)

cluster_num <- 8
genes <- file.path(source_data_dir, paste0("CLU.", 1:cluster_num, ".genes.txt"))

load_msigdb_gmt <- function(gmt_path, save_path = NULL) {
  suppressMessages(library(GSA))

  if (!file.exists(gmt_path)) {
    stop("Error: File not found at ", gmt_path)
  }

  suppressMessages({
    suppressWarnings({
      capture.output({
        new_msig.db <- GSA.read.gmt(gmt_path)
      })
    })
  })

  suppressMessages({
    suppressWarnings({
      capture.output({
        msig.db <- setNames(new_msig.db$genesets, new_msig.db$geneset.names)
      })
    })
  })

  message("Successfully loaded ", length(msig.db), " gene sets from MSigDB.")

  if (!is.null(save_path)) {
    save(msig.db, file = save_path)
    message("Saved msig.db to ", save_path)
  }

  return(invisible(msig.db))
}

gmt_file_path <- if (length(args) >= 2) args[2] else Sys.getenv("MSIGDB_GMT_FILE", "")

if (gmt_file_path == "") {
  stop("Provide the MSigDB GMT path as the second argument or set MSIGDB_GMT_FILE.")
}

msig.db <- load_msigdb_gmt(gmt_file_path, save_path = NULL)

gene.list <- list()
for (i in seq_along(genes)) {
  gene.list[[i]] <- read.table(genes[i], stringsAsFactors = FALSE)[, 1]
}
names(gene.list) <- paste(rep("CLU", cluster_num), 1:cluster_num, sep = "")

beta.GSEA <- Fun.enrich_withFC(markergenesList = gene.list, All.genes = all.genes.refer, db = msig.db)

write.csv(beta.GSEA[[1]], "Disease.gseaterm.csv")

mat_term_gene <- list()

for (i in seq_along(beta.GSEA[[2]])) {
  list1 <- c()

  for (n in names(beta.GSEA[[2]][[i]])) {
    tmp <- str_c(beta.GSEA[[2]][[i]][[n]], collapse = ",")
    list1 <- c(list1, tmp)
  }

  list1_df <- as.data.frame(t(list1), stringsAsFactors = FALSE)
  rownames(list1_df) <- names(beta.GSEA[[2]])[i]
  colnames(list1_df) <- names(beta.GSEA[[2]][[i]])
  mat_term_gene[[i]] <- list1_df
}

mat_term_gene <- do.call(rbind, mat_term_gene)
write.csv(mat_term_gene, "Disease.gsea_term_genes.csv")
