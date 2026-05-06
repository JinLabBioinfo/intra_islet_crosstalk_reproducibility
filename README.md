# Islet stellate cells sustain human beta-cell function through an RBP4-dependent retinoid paracrine loop

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20059249.svg)](https://doi.org/10.5281/zenodo.20059249)

This repository contains scripts and source data supporting the reproducibility of analysis for this study.

## Layout

- `sequencing/`: sequencing-related mapping, downstream analysis, and figure-generation scripts.
- `spatial/`: spatial analysis workflows and helper utilities.
- `metabolomics/`: metabolomics analysis and plotting scripts.
- `source_data/`: source tables and intermediate inputs for further analysis and visualization.

## Notes

- For spatial analysis, if you are using a Seurat version different from the one reported in the paper, you may need to run `UpdateSeuratObject()` after reading the `.rds` files.
- `spatial/_utility_funs.r` contains helper functions for transcript coordinate operations (adding/subsetting molecules) that are not yet implemented in Seurat, which some may find useful for their own purposes.
- PANC-DB cohort reads downsampling utility is in https://github.com/JinLabBioinfo/encapsulated_downsampling.

## Data Availability

- Raw in-house human islet scRNA-seq and bulk RNA-seq data generated in this study are available in NCBI GEO under SuperSeries accession `GSE310437`.
- The human bulk RNA-seq SubSeries is `GSE310435`.
- The human scRNA-seq SubSeries is `GSE310436`.
- In-house mouse single-cell multiomics sequencing data generated in this study are available in NCBI GEO under accession `GSE327516`.
- Count matrices and processed Seurat objects of spatial data are available on Zenodo: <https://zenodo.org/records/18474974>.
