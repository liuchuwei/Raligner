# MNN --------------------------------------------------------------------

# Modification of the scran::fastMNN (https://github.com/MarioniLab/scran)
# Allows for separate k values per dataset, and simplifies some of the IO and doesn't use PCA reduction
modified_mnnCorrect <- function(ref_mat, targ_mat, k1 = 20, k2 = 20,
                                ndist = 3, subset_genes = NULL) {
  if (is.null(subset_genes)) {
    subset_genes <- colnames(ref_mat)
  }

  sets <- batchelor::findMutualNN(ref_mat[, subset_genes],
                                  targ_mat[, subset_genes],
                                  k1 = k2, k2 = k1,
                                  BPPARAM = BiocParallel::SerialParam())
  mnn_pairs <- as.data.frame(sets) %>%
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
                  targ_ID = rownames(targ_mat)[second],
                  pair = seq(nrow(.))) %>%
    dplyr::select(-first, -second)

  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)

  #remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)

  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat, re.ave.out$averaged, re.ave.out$second, k=k2, ndist=ndist, subset_genes, BPPARAM=BiocParallel::SerialParam())

  final <- list(corrected = targ_mat,
                pairs = mnn_pairs)
  return(final)
}

# Copied from dev version of scran (2018-10-28) with slight modifications as noted
#https://github.com/MarioniLab/scran
.average_correction <- function(refdata, mnn1, curdata, mnn2)
  # Computes correction vectors for each MNN pair, and then
  # averages them for each MNN-involved cell in the second batch.
{
  corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
  corvec <- rowsum(corvec, mnn2)
  npairs <- table(mnn2)
  stopifnot(identical(names(npairs), rownames(corvec)))
  corvec <- unname(corvec)/as.vector(npairs)
  # corvec = apply(corvec, 2, function(x){x = x/npairs})
  list(averaged=corvec, second=as.integer(names(npairs)))
}


.center_along_batch_vector <- function(mat, batch.vec)
  # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
  # This removes any variation along the overall batch vector within each matrix.
{
  batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
  batch.loc <- as.vector(mat %*% batch.vec)
  central.loc <- mean(batch.loc)
  mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
  return(mat)
}

#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, subset_genes, BNPARAM=NULL, BPPARAM=BiocParallel::SerialParam())
  # Computing tricube-weighted correction vectors for individual cells,
  # using the nearest neighbouring cells _involved in MNN pairs_.
  # Modified to use FNN rather than queryKNN for nearest neighbor finding
{
  cur.uniq <- curdata[in.mnn,,drop=FALSE]
  safe.k <- min(k, nrow(cur.uniq))
  # closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest <- FNN::get.knnx(cur.uniq[, subset_genes], query=curdata[, subset_genes], k=safe.k)
  # weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
  weighted.correction <- .compute_tricube_average(correction, closest$nn.index, closest$nn.dist, ndist=ndist)
  curdata + weighted.correction
}

.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3)
  # Centralized function to compute tricube averages.
  # Bandwidth is set at 'ndist' times the median distance, if not specified.
{
  if (is.null(bandwidth)) {
    middle <- ceiling(ncol(indices)/2L)
    mid.dist <- distances[,middle]
    bandwidth <- mid.dist * ndist
  }
  bandwidth <- pmax(1e-8, bandwidth)

  rel.dist <- distances/bandwidth
  rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
  tricube <- (1 - rel.dist^3)^3
  weight <- tricube/rowSums(tricube)

  output <- 0
  for (kdx in seq_len(ncol(indices))) {
    output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
  }

  if (is.null(dim(output))) {
    matrix(0, nrow(vals), ncol(vals))
  } else {
    output
  }
}

# run MNN----

#' MNN correction
#'
#' @param obj raligner object
#' @param assay raw or correct
#' @param top_DE_genes top differential genes to use
#' @param k1 knn nubmer of tumors
#' @param k2 knn number of cell lines
#' @param ndist number of knn
#'
#' @return raligner object with MNN correction result
#' @export
#'
#' @examples
#' raligner = run_MNN(raligner)
#'
run_MNN = function(obj, assay="raw", top_DE_genes = 1000, k1 = 3, k2 = 10, ndist = 3){

  # prepare data
  subobj = slot(obj, "assay")
  subobj = slot(subobj, assay)
  ref_mat = subobj@cell %>% t()
  targ_mat = subobj@bulk %>% t()
  DE_genes = obj@meta@gene
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < top_DE_genes) %>%
    .[["Gene"]]

  print("Find mutual near neighbour...")
  sets <- batchelor::findMutualNN(ref_mat[,DE_gene_set],
                                  targ_mat[,DE_gene_set],
                                  k1 = k1, k2 = k2,
                                  BPPARAM = BiocParallel::SerialParam())

  mnn_pairs <- as.data.frame(sets) %>%
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
                  targ_ID = rownames(targ_mat)[second],
                  pair = seq(nrow(.))) %>%
    dplyr::select(-first, -second)

  print("MNN correction...")
  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)

  # remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)

  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat,
                                           re.ave.out$averaged,
                                           re.ave.out$second,
                                           k=k2, ndist=ndist,
                                           subset_genes=DE_gene_set,
                                           BPPARAM=BiocParallel::SerialParam())

  mnn_slot = slot(obj, "mnn_pair")
  res = new("PairData")
  res@pair = mnn_pairs
  res@correction = re.ave.out
  res@bulk_mat =  ref_mat %>% t %>% data.frame()
  res@cell_mat = targ_mat %>% t %>% data.frame()

  if (!assay %in% c("raw", "correct")) {
    stop("assay don't find")
  }

  if (assay == "raw") {
    mnn_slot@raw = res
  }else{
    mnn_slot@correct = res

  }


  obj@mnn_pair = mnn_slot

  return(obj)
}
