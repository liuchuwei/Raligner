library(devtools)
load_all()

# devtools::install_local("tmp/CytoTRACE_0.3.3.tar.gz")

# read data
tcga_ccle = readRDS("log/tcga_ccle.rds")

# CytoTRACE
library(CytoTRACE)
results <- CytoTRACE(tcga_ccle@assay$raw@cell)
saveRDS(results, file = "log/cell_cytotrace.rds")

out.table = data.frame(id = names(results$CytoTRACE), results$CytoTRACE)
write.table(out.table, file = "log/cytotrace.xls", sep = "\t", col.names = T, row.names = F)
