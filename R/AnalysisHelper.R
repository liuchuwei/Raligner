
#' Analysis gene expression change after correct method of Raligner
#'
#' @param obj raligner object after correct analysis
#' @param gensets "KEGG", "GO", "Reactome", "Hallmark"
#'
#' @return gsea result
#' @export
#'
#' @examples
DegGsea = function(obj, genesets = "KEGG"){

  require(limma)
  require(dplyr)
  require(org.Hs.eg.db)

  df = cbind(obj@assay$raw@bulk, obj@assay$correct@bulk)
  list <- c(rep("raw", ncol(df)/2), rep("correct", ncol(df)/2)) %>% factor(., levels = c("raw", "correct"), ordered = F)

  list <- model.matrix(~factor(list)+0)
  colnames(list) <- c("raw", "correct")
  df.fit <- lmFit(df, list)

  df.matrix <- makeContrasts(correct - raw, levels = list)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,n = Inf, adjust = "fdr")

  genename = row.names(tempOutput)
  gene_map <- select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
  colnames(gene_map)[1]<-"Gene"

  genelist_input = data.frame(Gene = row.names(tempOutput), logFC = tempOutput$logFC)
  aaa<-inner_join(gene_map,genelist_input,by = "Gene")
  aaa<-aaa[,-1]
  aaa<-na.omit(aaa)
  aaa$logFC<-sort(aaa$logFC,decreasing = T)

  geneList = aaa[,2]
  names(geneList) = as.character(aaa[,1])

  if (!genesets %in% c("KEGG", "GO", "Reactome", "Hallmark")) {
    stop("Genesets not support!")
  }

  if (genesets == "KEGG") {
    genset = fgsea::gmtPathways("data/path/c2.cp.kegg.v2023.1.Hs.entrez.gmt")
  }
  if (genesets == "GO") {
    genset = fgsea::gmtPathways("data/path/c5.go.v2023.1.Hs.entrez.gmt")
  }
  if (genesets == "Rectome") {
    genset = fgsea::gmtPathways("data/path/c2.cp.reactome.v2023.1.Hs.entrez.gmt")
  }
  if (genesets == "Hallmark") {
    genset = fgsea::gmtPathways("data/path/h.all.v2023.1.Hs.entrez.gmt")
  }
  fgRes <- fgsea::fgsea(pathways = genset,
                        stats = geneList,
                        minSize=3, ## minimum gene set size
                        maxSize=300, ## maximum gene set size
                        nperm=10000) %>%
    as.data.frame()

  return(fgRes)
}


#'  Analysis stromal and immune score change after correct method of Raligner
#'
#' @param obj raligner objcet
#'
#' @return stromal and immune score of raw data and correct data
#' @export
#'
#' @examples
#'
stroImm = function(obj){
  genesets = data.table::fread("data/path/stromal_immune.txt")
  genesets = as.list(genesets)

  raw = GSVA::gsva(as.matrix(ralign@assay$raw@bulk), genesets, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'
  correct = GSVA::gsva(as.matrix(ralign@assay$correct@bulk), genesets, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)   # signature 'matrix,list'

  res = list(raw = raw, correct = correct)
  return(res)
}
