# Immune Quadrant
Code for making clusters in Tumuluru and Godfrey et al.

Functions can be found in `make_clusters.R`.

How to run:
```r
library(GSVA)
library(ggfortify)
library(tidyverse)

source("make_clusters.R")

# Read in log2CPM matrix
x <- readRDS("log2cpm_matrix.rds")

# Read in gene sets
tib_genesets <- read_csv(file = "genesets.csv",  col_types = "cc")

# Run clustering pipeline
# List includes GSVA scores, PCA data, clusters, PCA plot, and PCA biplot
clusters_list <- fun_make_clusters(x, genesets = tib_genesets)
```
