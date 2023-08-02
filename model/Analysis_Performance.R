# load file
fls = list.files("model/result/tcga/", full.names = T)
fls = as.list(fls)

# get in-sample pearsion correlation
pcc.res = list()
for (i in 1:length(fls)) {
  print(i)
  tryCatch({res = readRDS(fls[[i]])},
           error = function(e) {
             load(fls[[i]])
             res = pool.res
             saveRDS(res, file = fls[[i]])
             return(res)})

  tmp_list = list()
  tmp_list$model_summary[[5]] = 0
  for (j in 1:length(res)) {
    if (class(res[[j]]) == "numeric") {
      res[[j]] = tmp_list

    }
  }

  unlist(lapply(res, function(u)as.numeric(u$model_summary[5]))) -> cv_R2_avg
  which.max(cv_R2_avg) -> idx

  res = res[[idx]]
  PCC = res$model_summary[[11]]
  drug = res$model_summary[[1]]
  out = data.frame(drug = drug, pcc = as.numeric(PCC))
  pcc.res[[i]] = out
}
pcc.res = do.call(rbind, pcc.res)
saveRDS(pcc.res, file = "model/pcc_performance.rds")
