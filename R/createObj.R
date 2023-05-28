
#' create raligner object
#'
#' @param raw_dir directory of raw data
#' @param project_name name of project. eg: tcga_cell
#' @param default default assays: raw、combat、cPCA、svd
#'
#' @return raligner object
#' @export
#'
#' @examples
#' ralign = createObj(raw_dir = "data/Demo1_tcga_cell/rawdata.rds",
#'                   project_name = "tcga_ccle",
#'                   default = "raw")


createObj = function(raw_dir,
                     project_name = "TumorCell"){

  print("create project...")
  ralign = new("Raligner")
  ralign@project = project_name

  print("load data...")
  rawdata = readRDS(raw_dir)

  # judge if meta information and exp information mapping
  if(!identical(colnames(rawdata$bulk_dat), rawdata$bulk_ann$Id) &
     identical(colnames(rawdata$cell_dat), rawdata$cell_ann$Id)){
    stop("meta information and exp information don't map")
  }

  ralign@assay@raw@cell = rawdata$cell_dat
  ralign@assay@raw@bulk = rawdata$bulk_dat
  ralign@meta@gene = rawdata$gene_ann
  ralign@meta@cell = rawdata$cell_ann
  ralign@meta@bulk = rawdata$bulk_ann
  ralign@meta@comb = rbind(rawdata$bulk_ann, rawdata$cell_ann)


  return(ralign)
}
