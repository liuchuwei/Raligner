# calculate gene station----
#' Calculate gene station
#'
#' @param obj raligner object
#'
#' @return gene station information
#' @export
#'
#' @examples raligner = CalGeneStation(raligner)
#'
CalGeneStation = function(obj){
  library(dplyr)
  print("get gene mean expression and standard variance...")
  gene_stats <- data.frame(
    Tumor_SD = apply(obj@assay@raw@bulk, 1, sd, na.rm=T),
    CCLE_SD = apply(obj@assay@raw@cell, 1, sd, na.rm=T),
    Tumor_mean = rowMeans(obj@assay@raw@bulk, na.rm=T),
    CCLE_mean = rowMeans(obj@assay@raw@cell, na.rm=T),
    Gene = row.names(obj@assay@raw@bulk),
    stringsAsFactors = F) %>%
    dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene


  print("get differential expression gene among clusters...")
  Tumor_obj = CreatSeuObj(obj@assay@raw@bulk, obj@pData@bulk)
  Cell_obj = CreatSeuObj(obj@assay@raw@cell, obj@pData@cell)

  Tumor_obj <- cluster_data(Tumor_obj)
  Cell_obj <- cluster_data(Cell_obj)

  tumor_DE_genes <- find_differentially_expressed_genes(Tumor_obj)
  CL_DE_genes <- find_differentially_expressed_genes(Cell_obj)

  DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(
      tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats, by = 'Gene')

  obj@fData$Gene = obj@fData$symbol
  obj@fData = dplyr::left_join(DE_genes, obj@fData, by = "Gene")
  return(obj)
}
