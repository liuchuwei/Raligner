require(pryr)

# define single assay class----
setClass("SingleAssay", slots = list(cell = "data.frame",
                                bulk = "data.frame"))
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

# define multiple-assay class----
setClass("MultAssay", slots = list(raw = "SingleAssay",
                                correct = "SingleAssay",
                                mnn = "MNN"))

# define meta-information class----
setClass("metaIn", slot = list(
                             cell = "data.frame",
                             bulk = "data.frame",
                             comb = "data.frame"))

# define dimension reduction class----
setClass("DimRe", slot = list(umap = "data.frame",
                              tnse = "data.frame"))


# define Raligner class----
setClass("Raligner",slots=list(project="character",
                               assay="MultAssay",
                               pData ="metaIn",
                               fData = "data.frame",
                               dimRe = "DimRe"
))


