
#' Using seurat method for dimension reduction
#'
#' @param obj raligner object
#' @param assay "raw", "correct"
#'
#' @return raligner object
#' @export
#'
#' @examples
#' ralign = DimReduce(ralign, assay = "raw")

DimReduce = function(obj, assay = "raw"){

  subobj = slot(obj, "assay")
  subobj = slot(subobj, assay)

  exp_mat = cbind(subobj@bulk, subobj@cell)
  ann = obj@meta@comb

  type = c(rep("Tumor", ncol(subobj@bulk)),
           rep("Cell", ncol(subobj@cell)))

  obj@meta@comb$type = type

  seu = CreatSeuObj(exp_mat = exp_mat, ann = ann, type = type)

  if (!assay %in% c("raw", "correct")) {
    stop("assay is not found")
  }


  obj@reduction@umap = seu@reductions$umap@cell.embeddings %>% data.frame()
  obj@reduction@tnse = seu@reductions$tsne@cell.embeddings %>% data.frame()

  return(obj)
}
