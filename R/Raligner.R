require(pryr)

# define single assay class----
setClass("assay", slots = list(cell = "data.frame",
                                bulk = "data.frame"))

# define multiple-assay class----
setClass("Assays", slots = list(raw = "assay",
                                correct = "assay"))

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
  correction = "list"
))


# define Raligner class----
setClass("Raligner",slots=list(project="character",
                               assay="Assays",
                               meta ="metaIn",
                               reduction = "DimRe",
                               mnn_pair = "PairData"
                               ))



