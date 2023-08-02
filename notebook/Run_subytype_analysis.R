library(dplyr)
library(devtools)
load_all()

# load data
tcga_ccle = readRDS("log/tcga_ccle.rds")

# breast
ann = subset(tcga_ccle@pData@comb, Lineage == "breast")
row.names(ann) = ann$id
dat = cbind(tcga_ccle@assay$mnn$svd_correct@bulk, tcga_ccle@assay$mnn$svd_correct@cell)
colnames(dat) = stringr::str_replace_all(colnames(dat), "[.]", "-")
dat = dat[,ann$id]

seu = CreatSeuObj(exp_mat = dat, ann = ann)

seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                             dims = 1:70,
                             k.param = 20,
                             force.recalc = TRUE,
                             verbose = FALSE)

seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                             resolution = 0.8)

saveRDS(seu, file = "log/tcga_ccle_breast_immune_subtype.rds")

# skin
ann = subset(tcga_ccle@pData@comb, Lineage == "skin")
row.names(ann) = ann$id
dat = cbind(tcga_ccle@assay$mnn$svd_correct@bulk, tcga_ccle@assay$mnn$svd_correct@cell)
colnames(dat) = stringr::str_replace_all(colnames(dat), "[.]", "-")
dat = dat[,ann$id]

seu = CreatSeuObj(exp_mat = dat, ann = ann)

seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                             dims = 1:70,
                             k.param = 20,
                             force.recalc = TRUE,
                             verbose = FALSE)

seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                             resolution = 1.5)

saveRDS(seu, file = "log/tcga_ccle_skin_immune_subtype.rds")

# lung
ann = subset(tcga_ccle@pData@comb, Lineage == "lung")
row.names(ann) = ann$id
dat = cbind(tcga_ccle@assay$mnn$svd_correct@bulk, tcga_ccle@assay$mnn$svd_correct@cell)
colnames(dat) = stringr::str_replace_all(colnames(dat), "[.]", "-")
dat = dat[,ann$id]

seu = CreatSeuObj(exp_mat = dat, ann = ann)

seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                             dims = 1:70,
                             k.param = 20,
                             force.recalc = TRUE,
                             verbose = FALSE)

seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                             resolution = 1.5)

saveRDS(seu, file = "log/tcga_ccle_lung_immune_subtype.rds")
