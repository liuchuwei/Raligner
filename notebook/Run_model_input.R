# load data
library(Seurat)
library(SeuratDisk)

# -----------------------tcga ccle------------------
# read data
seu = readRDS("log/tcga_ccle_svd_seu.rds")

# seurat to loop
seu.loom <- as.loom(x = seu, filename = "model/input/tcga_cell.loom", verbose = FALSE)
seu.loom$close_all()


# -----------------------array ccle------------------
# read data
seu = readRDS("log/array_ccle_svd_seu.rds")

# seurat to loop
seu.loom <- as.loom(x = seu, filename = "model/input/array_cell.loom", verbose = FALSE)
seu.loom$close_all()
