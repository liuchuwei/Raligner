# Find genes which are differentially expressed between
# clusters in the data
find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
                                       seu_obj@meta.data$cluster,
                                       limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }

  return(cur_DE_genes)

}

# Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
run_lm_stats_limma <- function (mat, vec, covars = NULL, weights = NULL, target_type = "Gene",
                                limma_trend = FALSE)
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2)
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]],
                                , drop = F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]],
                               , drop = F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  }
  else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata, ]))
  }
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata, , drop = FALSE]
    combined[["pred"]] <- pred
    form <- as.formula(paste("~", paste0(colnames(combined),
                                         collapse = " + ")))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  }
  else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep("pred", colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  if (colnames(results)[1] == "ID") {
    colnames(results)[1] <- target_type
  }
  else {
    results %<>% tibble::rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    one_sided_p <- two_sided_p/2
    if (test_dir == "right") {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat <
                                                 0]
    }
    else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat >
                                                 0]
    }
    return(one_sided_p)
  }
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(logFC = "EffectSize",
                                                             AveExpr = "Avg", t = "t_stat", B = "log_odds", P.Value = "p.value",
                                                             adj.P.Val = "q.value", min_samples = "min_samples"))) %>%
    na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value,
                                                       EffectSize, "left"), p.right = two_to_one_sided(p.value,
                                                                                                       EffectSize, "right"), q.left = p.adjust(p.left, method = "BH"),
                             q.right = p.adjust(p.right, method = "BH"))
  return(results)
}


# Differentially expressed genes --------------------------------------------------------------------

# Estimate linear-model stats for a matrix of data with respect to a group of phenotype variables
# using limma with empirical Bayes moderated F-stats for p-values
run_lm_stats_limma_group <- function (mat, phenos, covars = NULL, weights = NULL, target_type = "Gene",
                                      limma_trend = FALSE)
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- rownames(mat) %>% intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  form <- as.formula(paste("~", paste0(colnames(phenos), collapse = " + ")))
  design <- model.matrix(form, data = phenos[udata, , drop = F])
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste("~", paste0(colnames(covars),
                                         collapse = " + ")))
    Cdesign <- model.matrix(form, data = covars[udata, ,
                                                drop = F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), "(Intercept)"),
                       drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) ==
                0)
    design %<>% cbind(Cdesign)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  design <- design[, colSums(design) > 2, drop = FALSE]
  targ_coefs <- setdiff(colnames(design), "(Intercept)")
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf,
                             sort.by = "F", genelist = colnames(mat))
  results %<>% tibble::rownames_to_column(var = target_type)
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(AveExpr = "Avg",
                                                             F = "F_stat", P.Value = "p.value", adj.P.Val = "q.value"))) %>%
    na.omit() %>% dplyr::select(-ProbeID)
  return(results)
}
