
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

DimReduce = function(obj, assay = "raw", mnn = FALSE){

  if (!mnn) {
    subobj = slot(obj, "assay")
    subobj = subobj[[assay]]
  }else{
    subobj = slot(obj, "assay")
    subobj = subobj$mnn
    subobj = subobj[[assay]]
  }


  exp_mat = cbind(subobj@bulk, subobj@cell)
  ann = obj@pData@comb

  Type = c(rep("Tumor", ncol(subobj@bulk)),
           rep("Cell", ncol(subobj@cell)))

  obj@pData@comb$Type = Type

  seu = CreatSeuObj(exp_mat = exp_mat, ann = ann, type = Type)

  if (!assay %in% c("raw", "correct", "mnn")) {
    stop("assay is not found")
  }

  dimRe.res = new("DimRe")

  dimRe.res@umap = seu@reductions$umap@cell.embeddings %>% data.frame()
  dimRe.res@tnse = seu@reductions$tsne@cell.embeddings %>% data.frame()

  if (!mnn) {
    obj@dimRe[[assay]] = dimRe.res
  }else{
    obj@dimRe[[paste("mnn", assay, sep = "_")]] = dimRe.res
  }

  return(obj)
}
