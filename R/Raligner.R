require(pryr)

# define single assay class----
setClass("assay", slots = list(cell = "data.frame",
                                bulk = "data.frame"))

# define multiple-assay class----
setClass("Assays", slots = list(
                                default = "character",
                                raw = "assay",
                                combat = "assay",
                                cPCA = "assay",
                                svd = "assay"))

# define meta-information class----
setClass("metaIn", slot = list(gene = "data.frame",
                             cell = "data.frame",
                             bulk = "data.frame"))

# define dimension reduction class----
setClass("UMAP", slot = list(
         raw = "data.frame",
         combat = "data.frame",
         mnn = "data.frame",
         cPCA = "data.frame",
         svd = "data.frame",
         cPCA_mnn = "data.frame",
         svd_mnn = "data.frame"))

setClass("TNSE", slot = list(
         raw = "data.frame",
         combat = "data.frame",
         mnn = "data.frame",
         cPCA = "data.frame",
         svd = "data.frame",
         cPCA_mnn = "data.frame",
         svd_mnn = "data.frame"))

setClass("DimRe", slot = list(umap = "UMAP",
                              tnse = "TNSE"))

# define pair class----
setClass("PairData",slots=list(
  pair = "data.frame",
  correction = "list"
))

# define mnn class----
setClass("MNN",slots=list(
  raw = "PairData",
  cPCA = "PairData",
  svd = "PairData"
))

# define Raligner class----
setClass("Raligner",slots=list(project="character",
                               assay="Assays",
                               meta ="metaIn",
                               reduction = "DimRe",
                               mnn_pair = "MNN"
                               ))
