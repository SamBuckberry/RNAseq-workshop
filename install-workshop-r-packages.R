#!/bin/R

### install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c("Rfastp", "QuasR", "Rsubread", "wiggleplotr",
                   "edgeR", "gprofiler2")

BiocManager::install(bioc_packages)

### Install from CRAN using bioc to ensure version compatability
cran_package <- c("magrittr", "tidyr", "stringr", "ggplot2", "pheatmap")

BiocManager::install(pkgs = cran_package, site_repository = "CRAN")

### Install packages that are required outside of R
### Note that this requires pip3 to be installed

system(command = "pip3 install multiqc")


