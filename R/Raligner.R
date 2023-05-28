require(pryr)

# define single assay class----
setClass("SingleAssay", slots = list(cell = "data.frame",
                                bulk = "data.frame"))

# define multiple-assay class----
setClass("MultAssay", slots = list(raw = "SingleAssay",
                                correct = "SingleAssay"))

# define meta-information class----
setClass("metaIn", slot = list(gene = "data.frame",
                             cell = "data.frame",
                             bulk = "data.frame",
                             comb = "data.frame"))

# define dimension reduction class----
setClass("DimRe", slot = list(umap = "data.frame",
                              tnse = "data.frame"))

# define pair class----
setClass("PairData",slots=list(
  pair = "data.frame",
  correction = "list",
  bulk_mat = "data.frame",
  cell_mat = "data.frame"
))

# define MNN class----
setClass("MNN",slots=list(
  raw = "PairData",
  correct = "PairData"
))

# define Raligner class----
setClass("Raligner",slots=list(project="character",
                               assay="MultAssay",
                               meta ="metaIn",
                               reduction = "DimRe",
                               mnn_pair = "MNN"
                               ))



