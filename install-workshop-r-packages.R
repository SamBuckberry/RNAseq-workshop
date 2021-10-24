#!/bin/R

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rfastp", "QuasR", "Rsubread"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Rhisat2")
BiocManager::install("Rsubread")

library(Rhisat2)
library(QuasR)

system(command = "pip3 install multiqc")


