library(devtools)
load_all()
source("model/nest_elastic_net.R")

# load data
tcga_cell = readRDS("log/tcga_cell.rds")

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
sub_drug = drug_in$`5-Fluorouracil`
judge_cell = tcga_cell@pData@comb
judge_cell = c(judge_cell$Type == "Cell")
judge_blood =  c( tcga_cell@pData@cell$lineage != "blood")
drug_name = "5FU"
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

saveRDS(pool.res, file = "model/Fluorouracil.rds")
# test result
res = readRDS("model/Fluorouracil.rds")

# choose best model
unlist(lapply(res, function(u)as.numeric(u$model_summary[5]))) -> cv_R2_avg
which.max(cv_R2_avg) -> idx
res.list = res[[idx]]

fit <- res.list$model

# predict drug response
library(glmnet)
path = paste("data/vega/tcga/latent_", (idx-1),".csv", sep = "")
bulk.test.data = data.table::fread(path)[tcga_cell@pData@comb$Type == "Tumor",]
bulk.probabilities = predict(fit, as.matrix(bulk.test.data), s = 'lambda.min')
df = data.frame(id = tcga_cell@pData@bulk$sampleID, drug = bulk.probabilities)

# survival analysis
sur = data.table::fread("data/GDC-PANCAN.survival.tsv")
sur$id = stringr::str_sub(sur$sample, 1, 15)
sur_dat = dplyr::left_join(df, sur, by="id")
sur_dat = na.omit(sur_dat)
sur_dat$id = stringr::str_sub(sur_dat$id, 1, 12)

## drug response information
tcga_drug = readRDS("data/tcga_drug.rds")
sur_dat = dplyr::left_join(sur_dat, tcga_drug, by="id")

dat = subset(sur_dat, stringr::str_detect(sur_dat$Drug_list,"Fluorouracil"))
dat = subset(dat, dat$lineage %in% c("colorectal", "esophagus", "gastric"))
dat$type = ifelse(dat$lambda.min > median(dat$lambda.min), "resistant", "sensistive")

library(survival)
library(survminer)
fit <- survfit(Surv(OS.time, OS) ~ type, data = dat)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)

saveRDS(dat, file = "model/TCGA_5FU.rds")

