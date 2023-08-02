# load latent
latent = list.files("data/vega/tcga/", full.names = T)
latent = latent[stringr::str_detect(latent, "latent")]

# REACTOME_SIGNALING_BY_ERBB2
getERBB2 = function(sub_latent){
  sub_latent = data.table::fread(sub_latent)
  sub_latent = sub_latent$REACTOME_SIGNALING_BY_ERBB2
  out = sub_latent
  return(out)
}

latent = as.list(latent)
ERBB2_latent = lapply(latent, getERBB2)
ERBB2_latent = do.call(cbind, ERBB2_latent)

saveRDS(ERBB2_latent, file = "model/ERBB2_latent.rds")

cor_data <- cor(ERBB2_latent, method="pearson")
corrplot::corrplot(cor_data)

# ERBB2 Gene weight
weight = list.files("data/vega/tcga/", full.names = T)
weight = weight[stringr::str_detect(weight, "weight")]

# Ann weight
ann_weight = readRDS("model/tcga_ann_weight.rds")
getWeight = function(sub_weight){
  sub_weight = data.table::fread(sub_weight)[-1,]
  sub_weight = as.data.frame(sub_weight)
  sub_weight = sub_weight[, which(ann_weight$path == "REACTOME_SIGNALING_BY_ERBB2")]
  sub_weight = as.data.frame(sub_weight)
  row.names(sub_weight) = ann_weight$gene
  out = sub_weight
}
weight = as.list(weight)
ERBB2_weight = lapply(weight, getWeight)
ERBB2_weight = do.call(cbind, ERBB2_weight)
saveRDS(ERBB2_weight, file = "model/ERBB2_weight.rds")

mean_weight = apply(ERBB2_weight, 1, sum)
ERBB2_weight = ERBB2_weight[order(mean_weight, decreasing = T),]
