#' Title
#'
#' @param obj raligner object
#' @param assay raw or corrected assay
#' @param top_DE_genes top genes to use
#'
#' @return density plot
#' @export
#'
#' @examples
#'
#' DenPlot(obj = ralign, assay = "raw")

DenPlot <- function(obj,
                    assay = "raw",
                    top_DE_genes = 1000, mnn = FALSE) {


  DE_genes = obj@fData
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < top_DE_genes) %>%
    .[["Gene"]]

  if (mnn) {
    obj = slot(obj, "assay")$mnn
  }else{
    obj = slot(obj, "assay")
  }

  cell = reshape2::melt(obj[[assay]]@cell[DE_gene_set,])
  cell$type = "Cell"
  bulk = reshape2::melt(obj[[assay]]@bulk[DE_gene_set,])
  bulk$type = "Tumor"

  dat = rbind(cell, bulk)
  require(ggplot2)
  p = ggplot(dat, aes(x = value, colour = type, fill = type)) +
    geom_density(alpha = 0.7) +
    theme_classic() +
    xlab("Gene Expression Intensity") +
    ylab("Density")
  print(p)
  return(p)
}
