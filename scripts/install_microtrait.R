#!/usr/bin/Rscript

# https://github.com/ukaraoz/microtrait

############################## LOG
sink(file=file(snakemake@log$out, open="wt"), type=c("output", "message"))

############################## LIBS
# installing packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install(version = "3.19")
suppressMessages(library(BiocManager))

list_of_packages = c("R.utils", "RColorBrewer", "MASS", "Matrix", "mgcv", "ggplot2", "ape", "assertthat", "checkmate",
 "corrplot", "doParallel", "dplyr", "futile.logger", "grid", "gtools", 
 "kmed", "lazyeval", "magrittr", "parallel", "pheatmap", "readr", "stringr", "tibble", 
  "tictoc", "tidyr")
newpackages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(newpackages)) install.packages(newpackages, repos = "http://cran.us.r-project.org")

# if (!requireNamespace("BiocManager", quietly = TRUE))
# 	install.packages("BiocManager")
#BiocManager::install("curl")

install.packages("devtools")
suppressMessages(library(devtools))

BiocManager::install("Biostrings")
BiocManager::install("coRdon")
BiocManager::install("ComplexHeatmap")

devtools::install_github("jlw-ecoevo/gRodon")

BiocManager::install("seqinr")
BiocManager::install("vegan")
devtools::install_github("ukaraoz/microtrait")

# loading library
library(microtrait)

# preparing database
microtrait::prep.hmmmodels()

list.files(file.path(.libPaths(), "microtrait/extdata/hmm/hmmpress"))

# writing dummy file
write(sprintf("Done: %s", Sys.time()), file=snakemake@output$done, append=TRUE)
