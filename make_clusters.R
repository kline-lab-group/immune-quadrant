# Functions to make clusters in Tumuluru and Godfrey et al.

#######################################################
# Dependencies:
#     - GSVA
#     - tidyverse
#     - ggfortify
#######################################################


# The pipeline contains 3 functions (gsva, pca, and kmeans) which can be run all together in make_clusters

fun_gsva <- function(lcpm, genesets = tib_genesets) {
  # Run GSVA
  
  # Parameters:
  #   lcpm - a log2 counts per million (CPM) matrix
  #   genesets - a dataframe of genesets to be used for GSVA
  #            - expects "gene" column with HUGO gene IDs and "gs" column with the name of the corresponding gene set 
  
  # Returns:
  #   gsva_dge - a matrix of GSVA scores by sample
  
  vec_gs <- unique(genesets$gs)
  gsva_dge <- gsva(
    expr = lcpm,
    kcdf = "Gaussian",
    gset.idx.list = lapply(vec_gs,
                           function(gs1) genesets %>%
                             filter(gs == gs1) %>%
                             pull(gene)
    ),
    annotation = NULL,
    method = "gsva",
    verbose = TRUE
  )
  
  rownames(gsva_dge) <- vec_gs
  
  return(gsva_dge)
  
} 

fun_pca <- function(gsva_dge, genesets = tib_genesets){
  # Run PCA on GSVA scores
  # Generally, PC1 is the "immune axis," and PC2 is the "cell of origin axis"
  # Reflects the output such that samples with higher immune scores are on the right and samples with higher ABC-DLBCL scores are on the top
  
  # Parameters:
  #   gsva_dge - output from fun_gsva()
  #   genesets - a dataframe of gene sets as used in fun_gsva()
  
  # Returns:
  #   pca_dge - an object of class "prcomp" where PCA scores are located in the "x" matrix
  
  vec_gs <- unique(genesets$gs)
  
  pca_dge <- prcomp(t(gsva_dge))
  rownames(pca_dge$rotation) <- vec_gs
  
  
  if (pca_dge$rotation["ifng","PC1"] < 0) {
    pca_dge$x[,"PC1"] <- -pca_dge$x[,"PC1"]
    pca_dge$rotation[,"PC1"] <- -pca_dge$rotation[,"PC1"]
  }
  
  if (pca_dge$rotation["ABCDLBCL-1","PC2"] < 0) {
    pca_dge$x[,"PC2"] <- -pca_dge$x[,"PC2"]
    pca_dge$rotation[,"PC2"] <- -pca_dge$rotation[,"PC2"]
  }
  
  return(pca_dge)
}  

fun_kmeans <- function(pca_dge, bool_plot = FALSE){  
  # Perform k-means clustering on PCA scores where k = 4
  # Defines 4 clusters:
  #   1. ABC Cold - cluster with low PC1 scores and high PC2 scores
  #   2. ABC Hot - cluster with high PC1 scores and high PC2 scores
  #   3. GCB Cold - cluster with low PC1 scores and low PC2 scores
  #   4. GCB Hot - cluster with high PC1 scores and low PC2 scores
  
  # Parameters:
  #   pca_dge - output from fun_pca()
  #   k - number of clusters to partition data into
  #   seed - seed for reproducible results
  #   bool_plot - if TRUE, makes a PCA plot and a PCA biplot using autoplot()
  
  # Returns a list with:
  #   kmeans - a dataframe with the cluster assignments of each sample
  #   biplot - PCA biplot (optional)
  #   pcaplot - PCA plot (optional)
  
  kmeans_dge <- kmeans(pca_dge$x, centers = 4, nstart = 50)
  
  axis_1 <- c("Cold", "Hot")
  axis_2 <- c("GCB", "ABC")
  fun_label_cluster <- function(clust) {
    axis_1_mean <- mean(pca_dge$x[kmeans_dge$cluster == clust,"PC1"])
    axis_2_mean <- mean(pca_dge$x[kmeans_dge$cluster == clust,"PC2"])
    return(paste(axis_2[(axis_2_mean > 0) + 1],
                 axis_1[(axis_1_mean > 0) + 1],
                 sep = " "))
  }
  
  df_kmeans <- data.frame(
    row.names = names(kmeans_dge$cluster),
    cluster = sapply(1:4, fun_label_cluster)[kmeans_dge$cluster]
  )
  
  if (bool_plot) {
    vec_cluster_colors <- c(
      "ABC Cold" = "blue4",
      "ABC Hot" = "#DC3220",
      "GCB Cold" = "#56B4E9",
      "GCB Hot" = "#E69F00"
    )
    
    plot_biplot <- autoplot(
      pca_dge, data = df_kmeans, colour = "cluster",
      loadings = TRUE, loadings.label = TRUE, frame = FALSE,
    ) +
      scale_color_manual(values = vec_cluster_colors)
    
    plot_pca <- autoplot(
      pca_dge, data = df_kmeans, colour = "cluster", frame = FALSE,
    ) +
      scale_color_manual(values = vec_cluster_colors)
  }else{
    plot_biplot <- NULL
    plot_pca <- NULL
  }
  
  return(
    list(
      kmeans = df_kmeans,
      biplot = plot_biplot,
      pcaplot = plot_pca
    )
  )
}

fun_make_clusters <- function(lcpm, genesets = tib_genesets, seed = 900, bool_plot = TRUE) {
  # Full pipeline for making clusters
  # Computes GSVA scores for gene sets, runs PCA, and performs k-means clustering
  
  # Parameters:
  #   lcpm - a log2 counts per million (CPM) matrix
  #   genesets - a dataframe of genesets to be used for GSVA
  #            - expects "gene" column with HUGO gene IDs and "gs" column with the name of the corresponding gene set
  #   seed - seed for reproducible results
  #   bool_plot - if TRUE, makes a PCA plot and a PCA biplot
  
  set.seed(seed)
  
  # Run GSVA
  gsva_dge <- fun_gsva(lcpm = lcpm, genesets = tib_genesets)
  
  # Run GSVA
  pca_dge <- fun_pca(gsva_dge, genesets = tib_genesets)
  
  # Run PCA
  kmeans_l <- fun_kmeans(pca_dge, bool_plot = bool_plot)
  
  # Output
  clusters_l <- list(gsva = gsva_dge,
                 pca = pca_dge,
                 kmeans = kmeans_l$kmeans,
                 pcaplot = kmeans_l$pcaplot,
                 biplot = kmeans_l$biplot)
  
  return(clusters_l)
}


