library(devtools)
load_all()

# gsea analysis
ralign = readRDS("log/tcga_ccle.rds")
fgsea = DegGsea(ralign, genesets = "KEGG")
saveRDS(fgsea, file = "log/tcga_ccle_raw_svd_deg_gsea_kegg.rds")
write.table(fgsea[,1:5], file = "log/tcga_ccle_fgsea.xls", sep = "\t", col.names = T, row.names = F)

# stromal immune score
stroImm.res = stroImm(ralign)
saveRDS(stroImm.res, file = "log/tcga_ccle_raw_svd_stroImm.rds")
out.table = do.call(rbind, stroImm.res)
out.table = t(out.table)

out.table = do.call(rbind, stroImm.res)
out.table = t(out.table)
out.table = data.frame(id = row.names(out.table), out.table)
colnames(out.table) = c("id", "raw_stromal_score", "raw_immune_score",
                        "raw_purity", "correct_stromal_score", "correct_immune_score",
                        "correct_purity")
write.table(out.table, file = "log/tcga_cell_stroImm.xls",sep = "\t", row.names = F, col.names = T)

# cluster analysis
expmat = cbind(ralign@assay$mnn$svd_correct@bulk, ralign@assay$mnn$svd_correct@cell)
ann = ralign@pData@comb
row.names(ann) = ann$id

seu = CreatSeuObj(exp_mat = expmat, ann = ann)

seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                             dims = 1:70,
                             k.param = 20,
                             force.recalc = TRUE,
                             verbose = FALSE)

seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                             resolution = 5)


seu_cluster = seu@meta.data
tmp = data.frame(Lineage = ralign@pData@comb$Lineage, clusters = seu@meta.data$seurat_clusters,
                 type = ralign@pData@comb$Type)
# tmp = subset(tmp, type == "Tumor")
tmp = table(tmp$Lineage, tmp$clusters) %>% data.frame()
tmp = reshape2::dcast(tmp, Var1 ~ Var2)
tmp = tibble::column_to_rownames(tmp, var = "Var1")

get_ann = function(item){
  item_name = row.names(tmp)[which(item == max(item))]
  item_propotion = round(item[which(item == max(item))]/sum(item), 2)
  item_name = stringr::str_replace_all(item_name, "_", " ")
  out = paste(item_name, item_propotion, sep = "_")
  return(out)
}

ann = apply(tmp, 2, get_ann)
ann = stringr::str_split(ann, "_", simplify = T)
ann = as.data.frame(ann)
colnames(ann) = c("lineage" , "max_proportion")
ann$cluster = c(1:nrow(ann))
# ann$lineage[which(ann$max_proportion<0.75)]= "mixture"
ann = ann[,c(3,1,2)]
write.table(ann, file = "log/tcga_ccle_cluster_proportion.xls", sep = "\t", row.names = F, col.names = T)

new.cluster.ids <- ann$lineage
names(new.cluster.ids) <- levels(seu)
seu <- Seurat::RenameIdents(seu, new.cluster.ids)

Seurat::DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.5) + Seurat::NoLegend()

saveRDS(seu, file = "log/tcga_ccle_svd_seu.rds")

# write meta information
meta = data.frame(id = ralign@pData@comb$id,
                  type = ralign@pData@comb$Type,
                  lineage = c(ralign@pData@bulk$lineage,
                              ralign@pData@cell$lineage),
                  subtype = c(ralign@pData@bulk$subtype,
                             ralign@pData@cell$subtype),
                  seu_cluster = seu@meta.data$seurat_clusters,
                  seu_ann = seu@active.ident)
meta = cbind(meta, ralign@dimRe$mnn_svd_correct@umap)
saveRDS(meta, file = "log/tcga_ccle_cluster_meta.rds")
write.table(meta, file = "log/tcga_ccle_cluster_meta.xls", sep = "\t", row.names = F, col.names = T)
