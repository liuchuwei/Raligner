# combat----
combatAdj = function(dat, meta){
  require(dplyr)
  combat<- sva::ComBat(dat = as.matrix(dat),
                       batch = meta$type)
  bulk_dat = combat[,subset(meta, type == "Tumor")$id]
  cell_dat = combat[,subset(meta, type == "Cell")$id]
  res = new("SingleAssay")
  res@cell = cell_dat %>% data.frame()
  res@bulk = bulk_dat %>% data.frame()
  return(res)
}

# cPCA----
cPCAadj = function(obj, npc = 4){
  require(dplyr)

  print("Create Seurat obj")
  tumor_obj <- CreatSeuObj(obj@assay$raw@bulk, obj@pData@bulk, type='tumor')
  cell_obj <- CreatSeuObj(obj@assay$raw@cell, obj@pData@cell, type='CL')

  tumor_obj <- cluster_data(tumor_obj)
  cell_obj <- cluster_data(cell_obj)

  print("Run cPCA method")
  cov_diff_eig <- run_cPCA(tumor_obj, cell_obj, npc)

  print("cPCA correction")
  cur_vecs <- cov_diff_eig$rotation[,  c(1:npc), drop = FALSE]

  DE_gene_set <- obj@fData %>%
    dplyr::filter(best_rank < 1000) %>%
    .[['Gene']]

  rownames(cur_vecs) <- row.names(obj@assay$raw@bulk)
  tumor_cor <- resid(lm(as.matrix(obj@assay$raw@bulk) ~ 0 + cur_vecs)) %>% t()
  cell_cor <- resid(lm(as.matrix(obj@assay$raw@cell) ~ 0 + cur_vecs)) %>% t()

  res = new("SingleAssay")
  res@cell = cell_cor %>% t() %>% data.frame()
  res@bulk = tumor_cor %>% t() %>% data.frame()

  return(res)
}

# svd----
SvdAdj = function(obj, nv = 3, type = NULL){

  dat = obj@assay$raw

  # Weighting ranks
  combdat = cbind(dat@bulk, dat@cell) %>% t()
  Rankdat = as.data.frame(apply(combdat, 2, function(x){x = rank(x)}))
  dat = data.frame(rank = unlist(as.vector(Rankdat)), exp = unlist(as.vector(combdat)))

  # Adjusted ranking matrix
  dat$rank2 = dat$rank^2
  quadraticModel <- lm(exp ~ rank + rank2, data=dat)
  RankW = 2*quadraticModel$coefficients[3]*Rankdat + quadraticModel$coefficients[2]
  NormalRank = Rankdat * RankW

  # SVD to remove non-biological effect
  BioRank = NormalRank
  if (type == "sc") {
    BioMean = apply(NormalRank, 2, mean)
  }else{
    BioMean = apply(NormalRank, 1, mean)
  }
  nonBio = BioRank - BioMean
  nonSVD = irlba::irlba(as.matrix(nonBio), nv = nv)
  svdAdj = BioRank -  nonSVD$u %*% diag(nonSVD$d) %*% t(nonSVD$v)

  bulk_dat = svdAdj[1:ncol(obj@assay$raw@bulk),]
  cell_dat = svdAdj[(ncol(obj@assay$raw@bulk)+1):nrow(svdAdj),]

  res = new("SingleAssay")
  res@cell = cell_dat %>% t() %>% data.frame()
  res@bulk = bulk_dat %>% t() %>% data.frame()
  return(res)
}

# main method----

#' AssayCorrect
#'
#' @param obj raligner object
#' @param method correction method: combat, svd, cPCA
#' @param nv represent the correct strength when chosing the svd method
#' @param npc represent the correct strength when chosing the cPCA method
#'
#' @return raligner object
#' @export
#'
#' @examples
#' ralginer = AssayCorrect(raligner, method = "combat")


AssayCorrect = function(obj, method, nv = 3, npc = 4, Seq_type = NULL){

  dat = cbind(obj@assay$raw@bulk,
              obj@assay$raw@cell)

  type = c(rep("Tumor", ncol(obj@assay$raw@bulk)),
           rep("Cell", ncol(obj@assay$raw@cell)))
  meta = data.frame(id = obj@pData@comb$id, type = type)
  row.names(meta) = meta$id

  if (!method %in% c("combat", "cPCA", "svd")) {
    stop(paste0("Don't suppot the method:", method, sep= ""))
  }
  if (method == "combat") {
    res = combatAdj(dat = dat, meta = meta)
  }

  if (method == "cPCA") {
    res = cPCAadj(obj, npc = npc)
  }

  if (method == "svd") {
    res = SvdAdj(obj, nv = nv, type = Seq_type)
  }

  obj@assay$correct = res
  return(obj)
}
