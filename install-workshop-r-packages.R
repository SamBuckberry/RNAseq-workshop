#!/bin/R

### install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_packages <- c('QuasR', 'Rsubread', 'GenomicRanges', 'GenomicAlignments',
  'GenomicFeatures', 'Rsamtools', 'rtracklayer', 'Rfastp',
  'edgeR', 'limma', 'ggbio', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
  'wiggleplotr', 'gprofiler2', 'biomaRt')

BiocManager::install(bioc_packages)

### Install from CRAN using bioc to ensure version compatability
cran_package <-  c('stringr', 'magrittr', 'DT', 'markdown', 'mime', 'rmarkdown', 'tinytex', 'knitr', 'pheatmap', 'tidyr', 'ggplot2', 'ggrepel', 'upsetR', 'cowplot'. 'dplyr')

BiocManager::install(pkgs = cran_package, site_repository = "CRAN")
