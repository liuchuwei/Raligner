
#' Using seurat method for dimension reduction
#'
#' @param obj raligner object
#' @param assay "raw", "combat", "mnn", "cPCA", "svd", "cPCA_mnn","svd_mnn"
#'
#' @return raligner object
#' @export
#'
#' @examples
#' ralign = DimReduce(ralign, assay = "raw")

DimReduce = function(obj, assay = "raw"){

  exp_mat = cbind(obj@assay@raw@bulk, obj@assay@raw@cell)
  ann = rbind(obj@meta@bulk, obj@meta@cell)
  type = c(rep("Tumor", ncol(obj@assay@raw@bulk)),
           rep("Cell", ncol(obj@assay@raw@cell)))

  seu = CreatSeuObj(exp_mat = exp_mat, ann = ann, type = type)

  if (!assay %in% c("raw", "combat",
                    "mnn", "cPCA",
                    "svd", "cPCA_mnn",
                    "svd_mnn")) {

    stop("assay is not found")
  }

  if (assay == "raw") {
    obj@reduction@umap@raw = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@raw = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "combat") {
    obj@reduction@umap@combat = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@combat = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "mnn") {
    obj@reduction@umap@mnn = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@mnn = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "cPCA") {
    obj@reduction@umap@cPCA = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@cPCA = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "svd") {
    obj@reduction@umap@svd = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@svd = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "cPCA_mnn") {
    obj@reduction@umap@cPCA_mnn = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@cPCA_mnn = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  if (assay == "svd_mnn") {
    obj@reduction@umap@svd_mnn = seu@reductions$umap@cell.embeddings %>% data.frame()
    obj@reduction@tnse@svd_mnn = seu@reductions$tsne@cell.embeddings %>% data.frame()
  }

  return(obj)
}
