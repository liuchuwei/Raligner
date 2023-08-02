library(devtools)
library(glmnet)
load_all()

# load project
tcga_cell = readRDS("log/tcga_ccle.rds")

# MET drug: Crizotinib;PHA665752;Foretinib
drugs = c("PHA-665752.rds", "Crizotinib.rds", "Foretinib.rds")
for (subdrug in drugs) {
  res = readRDS(paste("model/result/tcga/", subdrug, sep = ""))

  # choose best model
  unlist(lapply(res, function(u)as.numeric(u$model_summary[5]))) -> cv_R2_avg
  which.max(cv_R2_avg) -> idx
  res.list = res[[idx]]

  fit <- res.list$model

  # predict drug response
  path = paste("data/vega/tcga/latent_", (idx-1),".csv", sep = "")
  bulk.test.data = data.table::fread(path)[tcga_cell@pData@comb$Type == "Tumor",]
  bulk.probabilities = predict(fit, as.matrix(bulk.test.data), s = 'lambda.min')
  df = data.frame(id = tcga_cell@pData@bulk$sampleID, drug = bulk.probabilities,
                  lineage = tcga_cell@pData@bulk$lineage)

  # cnv analysis
  cnv = data.table::fread("data/GDC-PANCAN.gistic.tsv.gz")
  gene_ann = data.table::fread("data/gencode.v22.annotation.gene.probeMap")
  gene_ann = subset(gene_ann, gene_ann$gene == "MET")
  cnv = subset(cnv, cnv$V1 == gene_ann$id)

  cnv_dat = data.frame(id = stringr::str_sub(colnames(cnv), 1, 15), cnv_num = t(cnv[1,]))
  cnv_dat = dplyr::left_join(df, cnv_dat, by = "id")
  cnv_dat = na.omit(cnv_dat)
  # cnv_dat$cnv_num = ifelse(cnv_dat$cnv_num>0, "Gain", "Other")

  # boxplot
  require(ggplot2)
  require(ggpubr)
  ggplot(cnv_dat, aes(x = cnv_num, y = lambda.min)) +
    geom_boxplot() + stat_compare_means()

  saveRDS(cnv_dat, file = paste("model/TCGA_", subdrug, sep = "") )
}

