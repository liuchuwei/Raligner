
#' create raligner object
#'
#' @param cell_dat cell line expression information
#' @param cell_ann cell line meta information
#' @param bulk_dat tumor expression information
#' @param bulk_ann tumor meta information
#' @param comb_ann combined annotaion information
#' @param project_name name of project. eg: tcga_cell
#'
#' @return raligner object
#' @export
#'
#' @examples
#' ralign = createObj(raw_dir = "data/Demo1_tcga_cell/rawdata.rds",
#'                   project_name = "tcga_ccle",
#'                   default = "raw")



createObj = function(cell_dat, cell_ann, bulk_dat, bulk_ann, comb_ann,
                     project_name = "TumorCell"){

  print("create project...")

  require(dplyr)
  obj = new("Raligner", project = project_name)

  obj@assay = list()
  obj@dimRe = list()

  # judge if meta information and exp information mapping
  if(!identical(colnames(bulk_dat), row.names(bulk_ann)) &
     identical(colnames(cell_dat), row.names(cell_ann))){
    stop("meta information and exp information don't map")
  }

  samGene = intersect(row.names(bulk_dat), row.names(cell_dat))
  rawAssay = new("SingleAssay")
  obj@assay$raw = rawAssay

  obj@assay$raw@cell = cell_dat[samGene,]
  obj@assay$raw@bulk = bulk_dat[samGene,]

  hugo = readRDS("data/hugo_2023.rds")
  obj@fData = hugo[match(samGene, hugo$symbol),] %>% data.frame()
  obj@pData@cell = cell_ann %>% data.frame()
  obj@pData@bulk = bulk_ann %>% data.frame()
  obj@pData@comb = combann %>% data.frame()

  return(obj)
}
