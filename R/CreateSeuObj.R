
#' Create Seurat object and run dimension reduction
#'
#' @param exp_mat expression dataframe
#' @param ann annotation data
#' @param type source of the expression data: bulk/cell
#' @param ndims number of principle components
#'
#' @return seurat object
#' @export
#'
#' @examples
#'
#' seu = CreatSeuObj(exp_mat = exp_mat, ann = ann, type = type)

CreatSeuObj = function(exp_mat, ann, type=NULL, ndims = 70){
  library(dplyr)
  library(magrittr)
  seu_obj <- Seurat::CreateSeuratObject(exp_mat,
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann )
  if (!is.null(type)) {
    seu_obj@meta.data$type = type
  }

  # mean center the data, important for PCA
  print("scale data...")

  seu_obj <- Seurat::FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)

  seu_obj <- Seurat::ScaleData(seu_obj,
                               features = rownames(Seurat::GetAssayData(seu_obj)),
                               do.scale = F)

  print("dimension reduction...")
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = Seurat::VariableFeatures(object = seu_obj),
                              npcs = ndims,
                              verbose = F)

  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA',
                               dims = 1:ndims,
                               reduction = 'pca',
                               n.neighbors = 10,
                               min.dist =  0.5,
                               metric = 'euclidean', verbose=F)
  try(
  seu_obj %<>% Seurat::RunTSNE(assay = 'RNA', dims = 1:ndims,
                                 check_duplicates = FALSE
  )

  )
  return(seu_obj)
}

# Cluster expression using Seurat clustering function
cluster_data <- function(seu_obj, ndims = 70) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca',
                                   dims = 1:ndims,
                                   k.param = 20,
                                   force.recalc = TRUE,
                                   verbose = FALSE)

  seu_obj %<>% Seurat::FindClusters(reduction = 'pca',
                                    resolution = 5)

  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters

  return(seu_obj)

}
