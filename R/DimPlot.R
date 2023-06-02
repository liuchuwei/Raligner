
#' Visualization of dimension reduction result
#'
#' @param obj raligner object
#' @param method umap, tnse
#'
#' @return ggplot result
#' @export
#'
#' @examples
#' DimPlot(ralign, method = "umap")

DimPlot =  function(obj, assay = "raw", method = "umap", type = NULL) {
  subobj = slot(obj, "dimRe")
  subobj = subobj[[assay]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, method))
  colnames(dat)[c(4, 5)] = c("x", "y")

  if (is.null(type)) {
    dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))
  }

  require(ggplot2)
  p = ggplot(dat, aes(x, y, fill=Lineage, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Tumor`=0.75)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(fill=FALSE, color=FALSE) +
    scale_fill_manual(values=tissue_colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2")
  print(p)
  return(p)
}

# source of code: https://github.com/broadinstitute/celligner/blob/master/R/global_params.R
get_tissue_colors = function(){
  tissue_colors <- c(`central_nervous_system`= "#f5899e",`engineered_central_nervous_system` = "#f5899e",
                     `teratoma` = "#f5899e",
                     `bone` = "#9f55bb",
                     `pancreas` = "#b644dc",
                     `soft_tissue` = "#5fdb69",
                     `skin` = "#6c55e2",
                     `liver` = "#9c5e2b",
                     `blood` = "#da45bb",
                     `lymphocyte`=  "#abd23f",
                     `peripheral_nervous_system` = "#73e03d",
                     `ovary` = "#56e79d",`engineered_ovary` = "#56e79d",
                     `adrenal` = "#e13978",  `adrenal_cortex` = "#e13978",
                     `upper_aerodigestive` = "#5da134",
                     `kidney` = "#1f8fff",`engineered_kidney` = "#1f8fff",
                     `gastric` = "#dfbc3a",
                     `eye` = "#349077",
                     `nasopharynx` = "#a9e082",
                     `nerve` = "#c44c90",
                     `unknown` = "#999999",
                     `cervix` = "#5ab172",
                     `thyroid` = "#d74829",
                     `lung` = "#51d5e0",`engineered_lung` = "#51d5e0",
                     `rhabdoid` = "#d04850",
                     `germ_cell` = "#75dfbb",   `embryo` = "#75dfbb",
                     `colorectal` = "#96568e",
                     `endocrine` = "#d1d684",
                     `bile_duct` = "#c091e3",
                     `pineal` = "#949031",
                     `thymus` = "#659fd9",
                     `mesothelioma` = "#dc882d",
                     `prostate` = "#3870c9", `engineered_prostate` = "#3870c9",
                     `uterus` = "#e491c1",
                     `breast` = "#45a132",`engineered_breast` = "#45a132",
                     `urinary_tract` = "#e08571",
                     `esophagus` = "#6a6c2c",
                     `fibroblast` = "#d8ab6a",
                     `plasma_cell` = "#e6c241")
  return(tissue_colors)
}
