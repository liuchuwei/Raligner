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
saveRDS(drug_in, file = "tmp/drug_in.rds")
drug_in = readRDS("tmp/drug_in.rds")

out.res = list()

# load latent
latents = paste("data/vega/tcga/latent_", c(0:99), ".csv", sep = "")
latents = as.list(latents)
names(latents) = paste("L", c(0:99), sep = "_")

# judge_cell = tcga_cell@pData@comb$Type == "Cell"
# saveRDS(judge_cell, file = "tmp/judge_cell.rds")
judge_cell = readRDS("tmp/judge_cell.rds")
# pool elastic net
for (i in 1:length(drug_in)) {
  drug_name = names(drug_in)[i]
  sub_drug = drug_in[[i]]
  pool_elastic_net = function(latent_path, sub_drug){
    X = data.table::fread(latent_path)
    X = X[judge_cell,]
    Y = sub_drug
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
  save_path = paste("model/result/tcga/", drug_name, ".rds", sep = "")
  saveRDS(pool.res, file = save_path)
}


