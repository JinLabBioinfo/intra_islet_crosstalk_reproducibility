Fun.enrich_withFC <- function (markergenesList, All.genes = all.genes.refer, db = msig.db, 
    qcutoff = 0.05, top = NULL) 
{
    require(qvalue)
    binomial.test.resultsList <- list()
    resultname <- c()
    for (markname in names(markergenesList)) {
        if (length((markergenesList[[markname]])) == 0) {
            next
        }
        print(markname)
        binomialtest.msig.result <- c()
        for (name in names(db)) {
            binomialtest.msig.result <- rbind(binomialtest.msig.result, 
                binomialtest.msig.enrch_deplet(as.character(markergenesList[[markname]]), 
                  All.genes, name, thedatabase = db))
        }
        fdr.enrich <- qvalue(binomialtest.msig.result[, 2])$qvalues
        fdr.deplete <- qvalue(binomialtest.msig.result[, 3])$qvalues
        FDR.combined <- c()
        FDR.combined[which(is.na(binomialtest.msig.result$ratio))] <- NA
        FDR.combined[which(binomialtest.msig.result$ratio > 1)] <- fdr.enrich[which(binomialtest.msig.result$ratio > 
            1)]
        FDR.combined[which(binomialtest.msig.result$ratio <= 
            1)] <- fdr.deplete[which(binomialtest.msig.result$ratio <= 
            1)]
        if (is.null(FDR.combined)) {
            next
        }
        binomialtest.msig.result <- cbind(binomialtest.msig.result, 
            FDR = FDR.combined)
        binomialtest.msig.result <- binomialtest.msig.result[order(binomialtest.msig.result$FDR), 
            ]
        binomial.test.resultsList <- c(binomial.test.resultsList, 
            list(binomialtest.msig.result))
        resultname <- c(resultname, markname)
    }
    names(binomial.test.resultsList) = resultname
    DRscombo <- data.frame()
    for (markname in names(binomial.test.resultsList)) {
        tmp <- subset(binomial.test.resultsList[[markname]], 
            FDR < qcutoff)
        row.names(tmp) <- tmp[, "name"]
        colnames(tmp)[c(4, 5)] <- paste(markname, colnames(tmp)[c(4, 
            5)], sep = "_")
        for (othermarkname in setdiff(names(binomial.test.resultsList), 
            markname)) {
            tmp2 <- binomial.test.resultsList[[othermarkname]]
            row.names(tmp2) <- tmp2[, "name"]
            tmp2 <- tmp2[row.names(tmp), c("ratio", "FDR")]
            colnames(tmp2) <- paste(othermarkname, colnames(tmp2), 
                sep = "_")
            tmp <- cbind(tmp, tmp2)
        }
        tmp <- tmp[, c(-1, -2, -3)]
        if (!is.null(top)) {
            tmp <- head(tmp, n = top)
        }
        DRscombo <- Tomerge.col(DRscombo, tmp)
    }
    all.list <- c()
    for (term in row.names(DRscombo)) {
        cur.genelist <- c()
        for (i in 1:length(markergenesList)) {
            cur.genes <- markergenesList[[i]][markergenesList[[i]] %in% 
                as.character(db[[term]])]
            cur.genelist <- c(cur.genelist, list(cur.genes))
        }
        names(cur.genelist) <- names(markergenesList)
        all.list <- c(all.list, list(cur.genelist))
    }
    names(all.list) <- row.names(DRscombo)
    return(list(DRscombo, all.list = all.list))
}


Tomerge.col <- function (df1, df2) 
{
    if (length(df1) == 0) {
        newdf <- df2
        return(newdf)
    }
    else {
        df2 <- df2[!row.names(df2) %in% row.names(df1), ]
        newdf <- data.frame(row.names = c(row.names(df1), row.names(df2)))
        for (name in names(df1)) {
            if (is.factor(df1[, name])) {
                assign(name, c(as.character(df1[, name]), as.character(df2[, 
                  name])))
            }
            else {
                assign(name, c(df1[, name], df2[, name]))
            }
            newdf <- data.frame(newdf, get(name))
            names(newdf)[ncol(newdf)] <- name
        }
        return(newdf)
    }
}


binomialtest.msig.enrch_deplet <- function (mylist, All = All.genes, name, thedatabase = db) 
{
    n <- length(mylist)
    x <- length(which(mylist %in% thedatabase[[name]]))
    p <- length(setdiff(thedatabase[[name]], setdiff(thedatabase[[name]], 
        All)))/length(All)
    binomitest.enrich <- binom.test(x, n, p, alternative = "greater")
    binomitest.deplete <- binom.test(x, n, p, alternative = "less")
    statistics <- data.frame(name, enrichsc = binomitest.enrich$p.value, 
        depletesc = binomitest.deplete$p.value, ratio = (x/n)/p)
    return(statistics)
}



# Function to load multiple MSigDB GMT files and combine them
load_msigdb_gmt <- function(gmt_dir, msigNames, common_suffix, save_path = NULL) {
  # Load required package silently
  if (!requireNamespace("GSA", quietly = TRUE)) {
    install.packages("GSA", quiet = TRUE)
  }
  suppressMessages(library(GSA))

  # Initialize an empty list to store gene sets from multiple GMT files
  combined_msig.db <- list()

  # Iterate over each MSigDB file name
  for (msigName in msigNames) {
    gmt_path <- file.path(gmt_dir, paste0(msigName, common_suffix))
    
    # Check if the file exists
    if (!file.exists(gmt_path)) {
      warning("Skipping missing file: ", gmt_path)
      next
    }

    # Load the GMT file, suppressing excessive output
    suppressMessages({
      suppressWarnings({
        temp_output <- capture.output({
          new_msig.db <- GSA.read.gmt(gmt_path)
        })
      })
    })

    # Convert and store each dataset separately
    suppressMessages({
      suppressWarnings({
        combined_msig.db[[msigName]] <- setNames(new_msig.db$genesets, new_msig.db$geneset.names)
      })
    })
  }

  # Print summary without flooding output
  message("Successfully loaded ", length(combined_msig.db), " MSigDB categories.")

  # Optional: Save the combined MSigDB object
  if (!is.null(save_path)) {
    save(combined_msig.db, file = save_path)
    message("Saved combined msig.db to ", save_path)
  }

  return(invisible(combined_msig.db))  # Prevents printing large output to the console
}
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



#all.genes.refer <- read.table("/path/to/all.genes.refer.txt", stringsAsFactors = FALSE)[,1]
#str(all.genes.refer)
#head(all.genes.refer)

# Define file paths
#gmt_file_path <- "/path/to/msigdb/"
#gmt_CommonNames <- ".all.v2024.1.Hs.symbols.gmt"
#msigNames <- c("c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "h")

# Load all MSigDB GMT files
#msig.db <- load_msigdb_gmt(gmt_file_path, msigNames, gmt_CommonNames, save_path = NULL)

# Check the structure (only showing the first category for clarity)
#str(msig.db[[1]])  # Adjust index as needed to inspect different sets
