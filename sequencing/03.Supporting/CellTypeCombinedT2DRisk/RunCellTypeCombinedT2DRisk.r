#!/usr/bin/env Rscript

# Demo usage:
# cd 29.Scripts/03.Supporting/CellTypeCombinedT2DRisk
# Rscript RunCellTypeCombinedT2DRisk.r Beta
# Rscript RunCellTypeCombinedT2DRisk.r Alpha
# Rscript RunCellTypeCombinedT2DRisk.r check

args <- commandArgs(trailingOnly = TRUE)
cell_type <- if (length(args) >= 1) args[1] else "Beta"

path_map <- data.frame(
  CellType = c("Beta", "Alpha", "Delta", "PP", "PSC", "Endothelial", "Duct", "Acinar"),
  CombinedObj = c(
    "/CombinedObjPath/Beta/",
    "/CombinedObjPath/Alpha/",
    "/CombinedObjPath/Delta/",
    "/CombinedObjPath/PP/",
    "/CombinedObjPath/PSC/",
    "/CombinedObjPath/Endothelial/",
    "/CombinedObjPath/Duct/",
    "/CombinedObjPath/Acinar/"
  ),
  ObjName = c(
    "Beta.All.integrated.1.rds",
    "Alpha.All.integrated.1.rds",
    "Delta.All.integrated.1.rds",
    "PP.All.integrated.1.rds",
    "PSC.All.integrated.1.rds",
    "Endothelial.All.integrated.1.rds",
    "Duct.All.integrated.1.rds",
    "Acinar.All.integrated.1.rds"
  ),
  stringsAsFactors = FALSE
)

if (tolower(cell_type) == "check") {
  path_map$FilePath <- file.path(path_map$CombinedObj, path_map$ObjName)
  path_map$Exists <- file.exists(path_map$FilePath)
  print(path_map[, c("CellType", "FilePath", "Exists")])
  quit(save = "no", status = 0)
}

if (!cell_type %in% path_map$CellType) {
  stop("Unsupported cell type: ", cell_type)
}

match_idx <- match(cell_type, path_map$CellType)
input_dir <- path_map$CombinedObj[match_idx]
obj_name <- path_map$ObjName[match_idx]
filePath <- file.path(input_dir, obj_name)
output_name <- paste0("Disease_Status_", cell_type, "_All.rds")

if (!file.exists(filePath)) {
  stop("Input object not found: ", filePath)
}

script_dir <- getwd()

library(Seurat)
library(ggplot2)
library(RePACT)
library(plyr)
library(dplyr)
library(ggrepel)
library(pscl)
library(qvalue)
library(rlist)

options(Seurat.object.assay.version = "v3")

CallT2Dpeak_qvalue_SpeedUP <- function(cellvsPeak.m.aggr = ATAConCCA.betaT2D.tjct.3nd.10bin.ob$cellvsPeak.m.aggr,
                                       depths = ATAConCCA.betaT2D.tjct.3nd.10bin.ob$depths,
                                       index = ATAConCCA.betaT2D.tjct.3nd.10bin.ob$index,
                                       qcut = 0.1, slopecut1 = 0.5, slopecut2 = -0.5, doscale = TRUE) {
  require(parallel)
  require(qvalue)
  cellvsPeak.m.aggr.norm <- t(t(cellvsPeak.m.aggr) / depths)
  cellvsPeak.m.aggr.norm.scale <- apply(cellvsPeak.m.aggr.norm, 1, function(x) {
    (x - mean(x)) / sd(x)
  }) %>% t
  index.scale <- scale(index)

  Get_stats <- function(gene, cellvsPeak.m.aggr, cellvsPeak.m.aggr.norm.scale, index.scale) {
    if (any(is.na(as.numeric(cellvsPeak.m.aggr.norm.scale[gene, ])))) {
      return(c(NA, NA, NA, NA, NA))
    } else {
      gmodel1 <- glm(as.numeric(cellvsPeak.m.aggr.norm.scale[gene, ]) ~ index.scale, family = gaussian)
      cur.slope1 <- summary(gmodel1)$coefficients[, 1][2]
      cur.pvalues1 <- summary(gmodel1)$coefficients[, 4][2]
      cur.corr <- cor(10000 * as.numeric(cellvsPeak.m.aggr.norm.scale[gene, ]), index.scale)
      cur.MaxMedian.FC <- max(cellvsPeak.m.aggr.norm[gene, ]) / median(cellvsPeak.m.aggr.norm[gene, ])
      cur.Max <- 1e+06 * max(cellvsPeak.m.aggr.norm[gene, ])
      names(cur.slope1) <- "slope"
      names(cur.pvalues1) <- "pvalue"
      names(cur.corr) <- "cor"
      names(cur.MaxMedian.FC) <- "MaxMedianFC"
      names(cur.Max) <- "MaxRPKM"
      return(c(cur.slope1, cur.pvalues1, cur.corr, cur.MaxMedian.FC, cur.Max))
    }
  }

  if (doscale == TRUE) {
    startTime <- Sys.time()
    pseudoregress.all.lis <- mclapply(
      rownames(cellvsPeak.m.aggr),
      Get_stats,
      cellvsPeak.m.aggr,
      cellvsPeak.m.aggr.norm.scale,
      index.scale,
      mc.cores = 40
    )
    endTime <- Sys.time()
    print(endTime - startTime)
  } else {
    startTime <- Sys.time()
    pseudoregress.all.lis <- mclapply(
      rownames(cellvsPeak.m.aggr),
      Get_stats,
      cellvsPeak.m.aggr,
      cellvsPeak.m.aggr.norm,
      index,
      mc.cores = 40
    )
    endTime <- Sys.time()
    print(endTime - startTime)
  }

  pseudoregress.all <- as.data.frame(list.rbind(pseudoregress.all.lis))
  row.names(pseudoregress.all) <- row.names(cellvsPeak.m.aggr)
  pseudoregress.all$pvalue <- pseudoregress.all$pvalue + 1e-10
  pseudoregress.all$qvalue <- p.adjust(pseudoregress.all$pvalue, method = "BH")
  pseudoregress.all <- pseudoregress.all[complete.cases(pseudoregress.all), ]
  UP <- subset(pseudoregress.all, qvalue < qcut & slope > slopecut1) %>%
    .[order(.$slope, decreasing = TRUE), ]
  DN <- subset(pseudoregress.all, qvalue < qcut & slope < slopecut2) %>%
    .[order(.$slope, decreasing = FALSE), ]
  UPDN.toplot <- t(cellvsPeak.m.aggr) / depths

  return(list(
    pseudoregress.all = pseudoregress.all,
    UP = UP,
    DN = DN,
    UPDN.toplot = UPDN.toplot
  ))
}

scRNA.RePACT <- function(OBJ, Sample, pheno, pheno_levels, is_continuous = FALSE,
                         if_donorWise = FALSE, binnumber = 20, PCrange = "",
                         RePACT_qvalCut = 0.005, donorWise_qvalCut = 0.01) {
  require(Seurat)
  require(plyr)
  require(dplyr)
  require(ggrepel)
  require(pscl)
  require(qvalue)
  require(rlist)

  BetaPCA <- OBJ@reductions$pca@cell.embeddings %>% Tomerge_v2(., OBJ@meta.data)
  BetaPCA <- BetaPCA[, 1:50] %>% apply(., 2, function(x) {
    x / sqrt(sum(x^2))
  }) %>% cbind(., BetaPCA[, 51:ncol(BetaPCA)])
  BetaPCA[, Sample] <- factor(BetaPCA[, Sample], levels = unique(BetaPCA[, Sample]))

  if (is_continuous == FALSE) {
    BetaPCA[, pheno] <- as.character(BetaPCA[, pheno])
    BetaPCA[, pheno] <- factor(BetaPCA[, pheno], levels = pheno_levels)
  }

  beta.rna.pca.withinfo.subs <- list()
  print("Calculating cell pseudo-index for 100 times")
  for (i in 1:100) {
    beta.rna.pca.withinfo.sub <- c()
    for (d in levels(BetaPCA[, Sample])) {
      tmp <- BetaPCA[which(BetaPCA[, Sample] == d), ]
      if (nrow(tmp) <= 200) {
        beta.rna.pca.withinfo.sub <- rbind(beta.rna.pca.withinfo.sub, tmp)
      } else {
        beta.rna.pca.withinfo.sub <- rbind(beta.rna.pca.withinfo.sub, tmp[sample(1:nrow(tmp), 200), ])
      }
    }
    beta.rna.pca.withinfo.subs <- c(beta.rna.pca.withinfo.subs, list(beta.rna.pca.withinfo.sub))
  }

  pseudo.indexes <- list()
  if (PCrange == "") {
    sigPCs <- SelectSigPCs(BetaPCA, Sample, pheno, is_continuous = FALSE)
    sigPCs_number <- as.numeric(sigPCs[, 1])
  } else {
    sigPCs_number <- PCrange
  }

  for (i in 1:100) {
    if (is_continuous == FALSE) {
      md <- GetRePACTmodel.cca(
        ccaWithinfo = beta.rna.pca.withinfo.subs[[i]],
        prefix = "PC",
        pheno = pheno,
        CCrange = c(sigPCs_number)
      )
    } else {
      md <- GetRePACTLinearmodel.cca(
        ccaWithinfo = beta.rna.pca.withinfo.subs[[i]],
        prefix = "PC",
        pheno = pheno,
        CCrange = c(sigPCs_number)
      )
    }

    trainingdata <- beta.rna.pca.withinfo.subs[[i]]
    Restdata <- BetaPCA[setdiff(row.names(BetaPCA), row.names(beta.rna.pca.withinfo.subs[[i]])), ]
    Alldata <- BetaPCA

    beta.rna.pca.withinfo.subs[[i]]$pseudo.index <- apply(trainingdata[, sigPCs_number], 1, function(x) {
      md$coefficients[[1]] + sum(x * md$coefficients[2:11])
    })
    Restdata$pseudo.index <- apply(Restdata[, sigPCs_number], 1, function(x) {
      md$coefficients[[1]] + sum(x * md$coefficients[2:11])
    })
    pseudo.indexes <- c(pseudo.indexes, list(apply(Alldata[, sigPCs_number], 1, function(x) {
      md$coefficients[[1]] + sum(x * md$coefficients[2:11])
    })))
  }

  BetaPCA$pseudo.index.balanced <- do.call(cbind, pseudo.indexes) %>% rowMeans()
  BetaPCA$rank <- rank(BetaPCA$pseudo.index.balanced)
  print("binning the cells along the trajectory")
  beta.RNA.PCA.20bin.ob <- MakeEvenBinBydepth_SpeedUP(
    OBJ = OBJ,
    data.info = BetaPCA[, 51:ncol(BetaPCA)],
    binnumber.tmp = binnumber
  )
  print("Getting statistics")
  betaT2D.diffGene.20bin.PCA <- CallT2Dpeak_qvalue_SpeedUP(
    beta.RNA.PCA.20bin.ob$cellvsPeak.m.aggr,
    beta.RNA.PCA.20bin.ob$depths,
    beta.RNA.PCA.20bin.ob$index,
    qcut = 0.2,
    slopecut1 = 0.3,
    slopecut2 = -0.3,
    doscale = TRUE
  )

  bin20.g.up <- row.names(subset(
    betaT2D.diffGene.20bin.PCA$UP,
    qvalue < RePACT_qvalCut &
      MaxRPKM > quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM, 0.25)
  ))
  bin20.g.dn <- row.names(subset(
    betaT2D.diffGene.20bin.PCA$DN,
    qvalue < RePACT_qvalCut &
      MaxRPKM > quantile(betaT2D.diffGene.20bin.PCA$pseudoregress.all$MaxRPKM, 0.25)
  ))

  RePACT_call <- list(UP = bin20.g.up, DN = bin20.g.dn)

  if (if_donorWise == TRUE) {
    print("Start donor-wise RePACT")
    PCAInfo <- BetaPCA
    RePACT.diff.genes.bydonor <- list()
    for (curDonor in levels(PCAInfo[, Sample])) {
      print(paste("processing", which(levels(PCAInfo[, Sample]) == curDonor), "/", length(levels(PCAInfo[, Sample]))))
      PCAInfo.donor <- subset(PCAInfo, Sample == curDonor)
      PCAInfo.donor$rank <- rank(PCAInfo.donor$pseudo.index.balanced)
      OBJ.donor <- subset(OBJ, cells = rownames(PCAInfo.donor))
      PCAInfo.20bin.donor.ob <- MakeEvenBinBydepth_SpeedUP(
        OBJ = OBJ.donor,
        data.info = PCAInfo.donor[, 51:ncol(PCAInfo.donor)],
        binnumber = binnumber
      )
      PCAInfo.20bin.donor.PCA <- CallT2Dpeak_pvalueOneTail(
        PCAInfo.20bin.donor.ob$cellvsPeak.m.aggr,
        PCAInfo.20bin.donor.ob$depths,
        PCAInfo.20bin.donor.ob$index,
        doscale = TRUE,
        GlobalSlopes = betaT2D.diffGene.20bin.PCA$pseudoregress.all[, 1, drop = FALSE]
      )
      RePACT.diff.genes.bydonor <- c(RePACT.diff.genes.bydonor, list(PCAInfo.20bin.donor.PCA))
    }

    print("Combining donor-wise results")
    names(RePACT.diff.genes.bydonor) <- levels(PCAInfo[, Sample])
    ps <- lapply(RePACT.diff.genes.bydonor, function(x) x$pvalue.onetail) %>% do.call(cbind, .)
    row.names(ps) <- row.names(RePACT.diff.genes.bydonor[[1]])
    ps[is.na(ps)] <- 1
    slopes <- lapply(RePACT.diff.genes.bydonor, function(x) x$slope) %>% do.call(cbind, .)
    row.names(slopes) <- row.names(RePACT.diff.genes.bydonor[[1]])
    FishersMethod.p <- apply(ps, 1, function(x) {
      -2 * sum(log(x))
    }) %>% pchisq(., 2 * 11, lower.tail = FALSE)
    FishersMethod.p <- FishersMethod.p[c(bin20.g.up, bin20.g.dn)]
    FishersMethod.q <- qvalue(FishersMethod.p)$qvalues %>% .[order(.)]
    FishersMethod.q.df <- data.frame(gene = names(FishersMethod.q), qvalueInSample = FishersMethod.q)
    Globaltag <- c()
    Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.up)] <- "UP"
    Globaltag[which(FishersMethod.q.df$gene %in% bin20.g.dn)] <- "DN"
    Globaltag[which(!FishersMethod.q.df$gene %in% c(bin20.g.up, bin20.g.dn))] <- ""
    FishersMethod.q.df$Globaltag <- Globaltag
    FishersMethod.q.df$Globaltag <- as.factor(FishersMethod.q.df$Globaltag)
    FishersMethod.q.dn.df <- subset(FishersMethod.q.df, Globaltag == "DN")
    FishersMethod.q.up.df <- subset(FishersMethod.q.df, Globaltag == "UP")
    FishersMethod.q.dn.df$rank <- rank(FishersMethod.q.dn.df$qvalueInSample)
    FishersMethod.q.up.df$rank <- rank(FishersMethod.q.up.df$qvalueInSample)
    FishersMethod.q.dn.df$InSampleTag <- ifelse(FishersMethod.q.dn.df$qvalueInSample < donorWise_qvalCut, "Intra-donor", "Inter-donor")
    FishersMethod.q.up.df$InSampleTag <- ifelse(FishersMethod.q.up.df$qvalueInSample < donorWise_qvalCut, "Intra-donor", "Inter-donor")
    DN.hetero <- names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn] <= donorWise_qvalCut]
    DN.homo <- names(FishersMethod.q[bin20.g.dn])[FishersMethod.q[bin20.g.dn] > donorWise_qvalCut]
    UP.hetero <- names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up] <= donorWise_qvalCut]
    UP.homo <- names(FishersMethod.q[bin20.g.up])[FishersMethod.q[bin20.g.up] > donorWise_qvalCut]
    RePACT_donorWise_intermediate <- list(
      FishersMethod.q.df = FishersMethod.q.df,
      FishersMethod.q.up.df = FishersMethod.q.up.df,
      FishersMethod.q.dn.df = FishersMethod.q.dn.df
    )
    RePACT_donorWise_call <- list(
      DN.interdonor = DN.hetero,
      DN.intradonor = DN.homo,
      UP.interdonor = UP.hetero,
      UP.intradonor = UP.homo
    )
    return(list(
      BetaPCA = BetaPCA,
      beta.RNA.PCA.20bin.ob = beta.RNA.PCA.20bin.ob,
      betaT2D.diffGene.20bin.PCA = betaT2D.diffGene.20bin.PCA,
      RePACT_call = RePACT_call,
      RePACT_donorWise_intermediate = RePACT_donorWise_intermediate,
      RePACT_donorWise_call = RePACT_donorWise_call
    ))
  } else {
    return(list(
      BetaPCA = BetaPCA,
      beta.RNA.PCA.20bin.ob = beta.RNA.PCA.20bin.ob,
      betaT2D.diffGene.20bin.PCA = betaT2D.diffGene.20bin.PCA,
      RePACT_call = RePACT_call
    ))
  }
}

prepareData <- function(filePath) {
  OBJ <- readRDS(filePath)
  print(OBJ)
  print("Finished loading and processing data")
  return(OBJ)
}

sessionInfo()
OBJ <- prepareData(filePath)
print(OBJ)

result <- scRNA.RePACT(
  OBJ,
  Sample = "Sample",
  pheno = "Disease_Status",
  pheno_levels = c("N", "Y"),
  is_continuous = FALSE,
  if_donorWise = FALSE,
  binnumber = 20,
  PCrange = "",
  RePACT_qvalCut = 0.05,
  donorWise_qvalCut = 0.01
)

saveRDS(result, output_name)
print(paste0(output_name, " saved successfully!"))
