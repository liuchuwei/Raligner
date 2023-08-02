library(devtools)
library(dplyr)
load_all()

# load data----
bulk = readRDS("data/array_dat.rds")
ccle = readRDS("data/ccle_dat.rds")

# create raligner project----
bulk_ann = bulk$phe %>% data.frame()
row.names(bulk_ann) = bulk_ann$id
bulk_dat = bulk$exp[, row.names(bulk_ann)]

cell_ann = ccle$cell_ann
cell_ann = subset(cell_ann, lineage %in% bulk_ann$lineage)
row.names(cell_ann) = cell_ann$sampleID
cell_dat = ccle$cell_dat[, cell_ann$sampleID]

combann = data.frame(id = c(row.names(bulk_ann), row.names(cell_ann)),
                     Lineage = c(bulk_ann$lineage, cell_ann$lineage),
                     Type = c(rep("Tumor", nrow(bulk_ann)), rep("Cell", nrow(cell_ann))))

name = "array_cell_line"

ralign = createObj(cell_dat = cell_dat, cell_ann = cell_ann,
                   bulk_dat = bulk_dat, bulk_ann = bulk_ann,
                   comb_ann = comb_ann, project_name = name)

# ---------------------------raligner correct------------------------
# get gene station----
ralign = CalGeneStation(ralign)

# raw assay----
## density
# DenPlot(obj = ralign, assay = "raw")

## dimension reduction
# ralign = DimReduce(ralign, assay = "raw")
# DimPlot(obj = ralign, assay = "raw", method = "umap")

# correct assay----
ralign = AssayCorrect(obj = ralign, method = "svd")
ralign = AssayCorrect(obj = ralign, method = "cPCA")
ralign = AssayCorrect(obj = ralign, method = "combat")

## density
# DenPlot(obj = ralign, assay = "correct")
ralign = DimReduce(ralign, assay = "raw")
ralign = DimReduce(ralign, assay = "svd_correct")
ralign = DimReduce(ralign, assay = "cPCA_correct")
ralign = DimReduce(ralign, assay = "combat_correct")
# DimPlot(obj = ralign, assay = "correct", method = "umap")

# mnn correct----
ralign = run_MNN(obj = ralign, assay = "raw")
ralign = run_MNN(obj = ralign, assay = "svd_correct")
ralign = run_MNN(obj = ralign, assay = "cPCA_correct")

## density
# DenPlot(obj = ralign, assay = "correct", mnn = T)
ralign = DimReduce(obj = ralign, assay = "svd_correct", mnn = TRUE)
ralign = DimReduce(obj = ralign, assay = "cPCA_correct", mnn = TRUE)
ralign = DimReduce(obj = ralign, assay = "raw", mnn = TRUE)
# DimPlot(obj = ralign, assay = "mnn_correct", method = "umap")

saveRDS(ralign, file = "log/array_ccle.rds")
