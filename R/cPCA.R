# cPCA --------------------------------------------------------------------

# run contrastive principal components analysis, first removing average cluster expression, to
# estimate the average intra-cluster covariance
# if pc_dims = NULL, all cPCs are calculated. Faster cPCA can be run by setting pc_dims to a
# value >=4 and approximating just those cPCs
run_cPCA_analysis <- function(TCGA_dat, CCLE_dat, tumor_cluster_df, CL_cluster_df, pc_dims=NULL) {
  tumor_clust_avgs <- get_cluster_averages(TCGA_dat, tumor_cluster_df)
  CL_clust_avgs <- get_cluster_averages(CCLE_dat, CL_cluster_df)

  TCGA_subtype_ms <- TCGA_dat - tumor_clust_avgs[tumor_cluster_df$seurat_clusters,]
  CCLE_subtype_ms <- CCLE_dat - CL_clust_avgs[CL_cluster_df$seurat_clusters,]

  TCGA_cov <- cov(TCGA_subtype_ms)
  CCLE_cov <- cov(CCLE_subtype_ms)

  if(!is.null(pc_dims)) {
    cov_diff_eig <- irlba::prcomp_irlba(TCGA_cov - CCLE_cov, n = pc_dims)
  } else {
    cov_diff_eig <- eigen(TCGA_cov - CCLE_cov)
  }
  return(cov_diff_eig)
}

# run contrastive principal components analysis
# set pc_dims to a value >= 4 to run fast cPCA by just calculating the top contrastive principle components
run_cPCA <- function(TCGA_obj, CCLE_obj, pc_dims = NULL) {
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')),
                                    t(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data')),
                                    TCGA_obj@meta.data, CCLE_obj@meta.data, pc_dims=pc_dims)
  return(cov_diff_eig)
}

# calculate the average expression per cluster
get_cluster_averages <- function(mat, cluster_df) {
  n_clusts <- nlevels(cluster_df$seurat_clusters)
  clust_avgs <- matrix(NA, nrow = n_clusts, ncol = ncol(mat)) %>%
    magrittr::set_colnames(colnames(mat)) %>%
    magrittr::set_rownames(levels(cluster_df$seurat_clusters))
  for (ii in levels(cluster_df$seurat_clusters)) {
    clust_avgs[ii,] <- colMeans(mat[cluster_df$seurat_clusters == ii,], na.rm=T)
  }
  return(clust_avgs)
}
