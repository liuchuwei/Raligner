library(devtools)
load_all()
source("model/nest_elastic_net.R")

# load data
tcga_cell = readRDS("log/tcga_ccle.rds")

# drug information
drug = data.table::fread("data/PANCANCER_IC_Mon Jun 19 04_50_40 2023_gdsc1.csv")
drug$uniqueID = c(1:nrow(drug))
cell_model = data.table::fread("data/Model.csv")
drug1 = subset(drug, `Cell Line Name` %in% cell_model$CellLineName)
drug2 = subset(drug, `Cell Line Name` %in% cell_model$StrippedCellLineName)

drug1$id = cell_model$ModelID[match(drug1$`Cell Line Name`, cell_model$CellLineName)]
drug2$id = cell_model$ModelID[match(drug2$`Cell Line Name`, cell_model$StrippedCellLineName)]
drug = rbind(drug1, drug2)
drug = drug[!duplicated(drug$uniqueID),]
colnames(drug)[1] = "DRUG_NAME"

bigMatrix = matrix(nrow = length(unique(drug$id)), ncol = length(unique(drug$DRUG_NAME)))
row.names(bigMatrix) = unique(drug$id)
colnames(bigMatrix) = unique(drug$DRUG_NAME)

for (i in 1:nrow(drug)) {
  bigMatrix[drug$id[i], drug$DRUG_NAME[i]] = drug$IC50[i]
}

cell_drug = bigMatrix[match(tcga_cell@pData@cell$sampleID, row.names(bigMatrix)),]

drug_in = as.list(as.data.frame(cell_drug))

out.res = list()

# load latent
latents = paste("data/vega/tcga/latent_", c(0:99), ".csv", sep = "")
latents = as.list(latents)
names(latents) = paste("L", c(0:99), sep = "_")

# pool elastic net
sub_drug = drug_in$Lapatinib
judge_cell = tcga_cell@pData@comb
judge_cell = c(judge_cell$Type == "Cell")
judge_blood =  c( tcga_cell@pData@cell$lineage != "blood")
drug_name = "Lapatinib"
pool_elastic_net = function(latent_path, sub_drug){
  X = data.table::fread(latent_path)
  X = X[judge_cell,]
  Y = sub_drug

  X = X[judge_blood, ]
  Y = Y[judge_blood]
  X = X[!is.na(Y),]
  Y = na.omit(Y)

  res = tryCatch({
    nest_lasso(X=X, Y=Y, n_folds=10, n_train_test_folds=5,
               seed=666, alpha=0.5, null_testing=FALSE, drug=drug_name)
  }, warning = function(w){
    res.list = list()
    res.list$model_summary[5] = 0
  }, error = function(e){
    res.list = list()
    res.list$model_summary[5] = 0
  })
  return(res)
}
pool.res = lapply(latents, pool_elastic_net, sub_drug=sub_drug)

saveRDS(pool.res, file = "model/lapatinib.rds")
# test result
res = readRDS("model/lapatinib.rds")

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
                lineage = tcga_cell@pData@bulk$lineage, subtype = tcga_cell@pData@bulk$subtype)

# drug response subtype
phe = data.table::fread("model/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
id = stringr::str_sub(phe$sampleID, 1, 15)
phe = as.data.frame(phe)
phe = phe[, stringr::str_detect(colnames(phe), "her2")]
phe$id = id

dat = dplyr::left_join(phe, df, by = "id")
dat = data.frame(id = dat$id, lambda.min = dat$lambda.min, subtype = dat$her2_immunohistochemistry_level_result)
dat = na.omit(dat)
dat = subset(dat, subtype != "")
## boxplot
ggplot(dat, aes(x = subtype, y = lambda.min)) +
  geom_boxplot() + stat_compare_means()

## explain
dat_path = bulk.test.data[match(dat$id, tcga_cell@pData@bulk$sampleID),]
dat$erbb2 = dat_path$REACTOME_SIGNALING_BY_ERBB2

## boxplot
ggplot(dat, aes(x = subtype, y = erbb2)) +
  geom_boxplot() + stat_compare_means()

saveRDS(dat, file = "model/TCGA_Lapatinib.rds")

