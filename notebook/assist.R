##------------------------assist for Figure1 visualization------------------
library(ggplot2)
Fig1_DimPlot = function(dat){

  p = ggplot(dat, aes(x, y, fill=type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.75, `Tumor`=1.5)) +
    theme_classic() +
    # theme(legend.position = 'bottom',
    #       text=element_text(size=12),
    #       legend.margin =margin(0,0,0,0)) +
    theme(legend.position = c(0.9, 0.9)) +
    # theme(legend.position = "right") +
    guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=c(`Cell` = "#ff7f00", `Tumor` = "#984ea3")) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme(legend.background = element_blank())

  return(p)
}


Figs1_DenPlot = function(dat){

  p = ggplot(dat, aes(x = value, colour = type, fill = type)) +
    geom_density(alpha = 0.7) +
    theme_classic() +
    xlab("Gene Expression Intensity") +
    ylab("Density") +
    scale_fill_manual(values=c(`Cell` = "#ff7f00", `Tumor` = "#984ea3"))+
    scale_color_manual(values=c(`Cell` = "#ff7f00", `Tumor` = "#984ea3"))+
    theme(
      legend.position=c(0.8,0.9),
      # legend.position= "top",
      legend.key.size = unit(10, "pt")
    )
  return(p)
}

##------------------------assist for FigureS1 visualization------------------
FigS1_Den = function(obj, assays = "correct", mnn = FALSE, title = "title"){

  DE_genes = obj@fData
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < 1000) %>%
    .[["Gene"]]


  if (mnn) {
    obj = slot(obj, "assay")$mnn
  }else{
    obj = slot(obj, "assay")
  }

  cell = reshape2::melt(obj[[assays]]@cell[DE_gene_set,])
  cell$type = "Cell"
  bulk = reshape2::melt(obj[[assays]]@bulk[DE_gene_set,])
  bulk$type = "Tumor"

  dat = rbind(cell, bulk)

  p = ggplot(dat, aes(x = value, colour = type, fill = type)) +
    geom_density(alpha = 0.7) +
    theme_classic() +
    xlab("Gene Expression Intensity") +
    ylab("Density") +
    scale_fill_manual(values=c(`Cell` = "#ff7f00", `Tumor` = "#984ea3"))+
    scale_color_manual(values=c(`Cell` = "#ff7f00", `Tumor` = "#984ea3"))+
    theme(
      legend.position=c(0.8,0.9),
      # legend.position= "top",
      legend.key.size = unit(10, "pt")
    ) +
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}



FigS1_Dim =  function(obj, assay = "raw", method = "umap", type = NULL, title = "title") {

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
    scale_size_manual(values=c(`Cell`=1.75, `Tumor`= 1.5)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(fill=FALSE, color=FALSE) +
    # guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=tissue_colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2")+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))

  print(p)
  return(p)
}


##------------------------assist for Figure2 visualization------------------
Fig2_tcga_ccle = function(obj){

  subobj = slot(obj, "dimRe")
  subobj = subobj[['mnn_svd_correct']]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, 'umap'))
  colnames(dat)[c(4, 5)] = c("x", "y")


  dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))


  require(ggplot2)
  # subdat = subset(dat, cluster == 31)
  p = ggplot(dat, aes(x, y, fill=Lineage, size=type, color = type)) +
  # ggplot(subdat, aes(x, y, fill=Lineage, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.75, `Tumor`=1.5)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(fill=FALSE, color=FALSE) +
    # guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=tissue_colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    annotate("text", x=-6.5, y=-10.5, label= "breast",size = 3)+
    annotate("text", x=-10.5, y=0, label= "blood",size = 3)+
    annotate("text", x=1, y=-12.5, label= "thymus",size = 3) +
    annotate("text", x= 1.5, y=-11, label= "thyroid",size = 3) +
    annotate("text", x=0, y=15, label= "germ cell",size = 3) +
    annotate("text", x=9, y=7, label= "kidney",size = 3)+
    annotate("text", x=5, y=-10, label= "skin",size = 3)+
    annotate("text", x=4, y=-14, label= "eye",size = 3)+
    annotate("text", x=9, y=0.5, label= "fibroblast",size = 3)+
    annotate("text", x=-4.5, y=2, label= "uterus",size = 3)+
    annotate("text", x=2.5, y=4, label= "pancreas",size = 3)+
    annotate("text", x=-6.8, y=11, label= "liver",size = 3)+
    annotate("text", x=7.5, y=11, label= "ovary",size = 3)+
    annotate("text", x=-4.5, y=4, label= "urinary",size = 3)+
    annotate("text", x=-9, y=6.5, label= "lymphocyte",size = 3)+
    annotate("text", x=-5, y=5.5, label= "plasma cell",size = 3)+
    annotate("text", x=12, y=-11, label= "bone",size = 3)+
    annotate("text", x=-8, y=-6, label= "adrenal",size = 3)+
    annotate("text", x=11, y=4, label= "soft tissue",size = 3)+
    annotate("text", x=11, y=-2, label= "peripheral nervous system",size = 3)+
    annotate("text", x=3, y=10, label= "colorectal",size = 3)+
    annotate("text", x=2.5, y=7.5, label= "gastric",size = 3)+
    annotate("text", x=0, y=-3, label= "upper aerodigestive",size = 3)+
    annotate("text", x=8.5, y=-6.5, label= "central nervous system",size = 3)+
    annotate("text", x=-1.2, y=3, label= "lung",size = 3)+
    annotate("text", x=-2.5, y=-1, label= "esophagus",size = 3)+
    annotate("text", x=-12.5, y=5, label= "prostate",size = 3)
    return(p)
}

Figure2_lineage_breast = function(obj){


  subtype = data.frame(rawid = c(obj@pData@bulk$sampleID,  obj@pData@cell$sampleID),
                       subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype))

  subdat = data.frame(rawid = c(obj@pData@bulk$sampleID,  obj@pData@cell$sampleID),
                      subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype),
                      lineage = obj@pData@comb$Lineage,
                      obj@dimRe$mnn_svd_correct@umap,
                      type = obj@pData@comb$Type)

  subdat= subset(subdat,  lineage == "breast")

  subdat$subtype[subdat$subtype == "HER2 amp"] = "HER2-enriched"
  subdat$subtype[subdat$subtype == "luminal HER2 amp"] = "HER2-enriched"
  subdat$subtype[subdat$subtype == "luminal A"] = "luminal"
  subdat$subtype[subdat$subtype == "luminal B"] = "luminal"


  p = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type
                          ,shape = subtype)) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`basal` = 21,
                                  `basal A` = 22,
                                  `basal B` = 24,
                                  `HER2-enriched` = 24,
                                  `luminal` = 21,
                                  `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    scale_fill_manual(values = c(`basal` = "#ff7f00",
                                 `basal A` = "#ff7f00",
                                 `basal B` = "#ff7f00",
                                 `HER2-enriched` = "#33a02c",
                                 `luminal` = "#6a3d9a",
                                 `normal` = "#a6cee3")) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.3),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE)+
    theme(legend.background = element_blank())+
    guides(shape = guide_legend(ncol = 2))

  return(p)
}
Fig2_array_cell = function(obj){
  # obj = array_cell
  subobj = slot(obj, "dimRe")
  subobj = subobj[["mnn_svd_correct"]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))
  dat$subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype)

  dat$Lineage[which(dat$subtype == "Acute myeloid leukemia")] = "blood"
  require(ggplot2)
  # subdat = subset(dat, Lineage == "plasma_cell")
  p = ggplot(dat, aes(x, y, fill=Lineage, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.75, `Tumor`=1.5)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(fill=FALSE, color=FALSE) +
    # guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=tissue_colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2")+
    annotate("text", x=-5, y=-5, label= "blood",size = 4)+
    annotate("text", x=-2, y=10, label= "blood",size = 4)+
    annotate("text", x=-10, y=12, label= "blood",size = 4)+
    annotate("text", x=-5.5, y=11, label= "bone",size = 4)+
    annotate("text", x=-2, y=-12, label= "bone",size = 4)+
    annotate("text", x= 2, y=-9, label= "colorectal",size = 4)+
    annotate("text", x=-2.5, y=-6.5, label= "ovary",size = 4)+
    annotate("text", x=-4, y=3.5, label= "plasma cell",size = 4)+
    annotate("text", x=6, y=3, label= "skin",size = 4)+
    theme(legend.position = c(0.15, 0.3),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0))+
    theme(legend.background = element_blank())

  return(p)
}

Fig2_array_cell_subtype = function(obj){
  # obj = array_cell
  subobj = slot(obj, "dimRe")
  subobj = subobj[["mnn_svd_correct"]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))
  dat$subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype)
  dat$Lineage[which(dat$subtype == "Acute myeloid leukemia")] = "blood"

  require(ggplot2)
  # subdat = subset(dat, Lineage == "blood")
  subdat = dat
  subdat$subtype[subdat$Lineage != "blood"] = "other"
  subdat$subtype = tolower(subdat$subtype)
  table(subdat$subtype)
  # subdat1 = subset(subdat, subtype %in% c("acute myeloid leukemia"))
  subdat1 = subdat
  subdat1$subtype = stringr::str_split(subdat1$subtype, ",", simplify = T)[,1]
  table(subdat1$subtype)
  # subdat1 = subset(subdat1, subtype %in% c("acute leukemia","acute myeloid leukemia",
  #                                          "acute lymphoblastic leukemia",
  #                                          "chronic myelogenous leukemia",
  #                                          "chronic lymphoblastic leukemia"))
  p = ggplot(subdat1, aes(x, y, fill=subtype, size=type, color = type)) +
  # ggplot(dat, aes(x, y, fill=subtype, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.75, `Tumor`=1.5)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    # guides(fill=FALSE, color=FALSE) +
    guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=c(`acute leukemia` = "#fb9a99",
                               `acute lymphoblastic leukemia` = "#fdbf6f",
                               `chronic myelogenous leukemia` = "#984ea3",
                               `other` = "#969696",
                               `acute myeloid leukemia` = "#fb9a99",
                               `b-cell` = "#a6cee3",
                               `chronic lymphoblastic leukemia` = "#1f78b4")) +
    xlab("UMAP 1") +
    ylab("UMAP 2")+
    annotate("text", x=-5, y=-5, label= "acute myeloid leukemia",size = 4)+
    annotate("text", x=0, y=10, label= "acute lymphoblastic leukemia",size = 4)+
    annotate("text", x=1, y=5, label= "chronic lymphoblastic leukemia",size = 4)+
    annotate("text", x=-8, y=7, label= "chronic myelogenous leukemia",size = 4)+
    # annotate("text", x=-5.5, y=13, label= "acute myeloid leukemia",size = 4)+
    guides(fill = guide_legend(ncol = 1, title = NULL))+
    theme(legend.position = c(0.15, 0.3),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0))+
    theme(legend.background = element_blank())

  return(p)
}

Fig2_scProstate = function(obj){
  subobj = slot(obj, "dimRe")
  subobj = subobj[["correct"]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  # dat$type = c(rep("Tumor", nrow(tcga_scProstate@pData@bulk)), rep("Cell", nrow(tcga_scProstate@pData@cell)))

  # dat = dat[order(dat$type),]
  dat$type = dat$Type
  dat$Type = obj@pData@cell$Subtype[match(dat$id,
                                          obj@pData@cell$id)]
  dat$Type[is.na(dat$Type)] = "Tumor"
  dat$Type = stringr::str_replace_all(dat$Type, " ", "")

  require(ggplot2)
  p = ggplot(dat, aes(x, y, fill = Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Tumor`='white',`Cell`='black')) +
    scale_size_manual(values=c(`Tumor`=1, `Cell`=1.25)) +
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(fill=FALSE, color=FALSE) +
    scale_fill_manual(values=c(`DU145_Resistant` = "#6a3d9a",
                               `DU145_Sensitive` = "#fb9a99",
                               `Tumor` = "#1f78b4",
                               `PC3_Resistant` = "#33a02c",
                               `PC3_Sensitive` = "#ff7f00")) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    annotate("text", x=3, y=-1, label= "DU145 resistant",size = 4)+
    annotate("text", x=2.5, y=4, label= "DU145 sensitive",size = 4)+
    annotate("text", x=1, y=-3, label= "prostate",size = 4)+
    annotate("text", x=-3, y=4, label= "PC3 resistant",size = 4)+
    annotate("text", x=-3, y=-3, label= "PC3 sensistive",size = 4)
  return(p)
}


Fig2_scLung = function(obj){

  subobj = slot(obj, "dimRe")
  subobj = subobj[["correct"]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")
  #
  # tcga_scLung@pData@cell$drug = paste(tcga_scLung@pData@cell$drug,
  #                                     tcga_scLung@pData@cell$time, sep = "_")
  dat$type = obj@pData@cell$time[match(dat$id,tcga_scLung@pData@cell$id)]
  dat$type[is.na(dat$type)] = "Tumor"
  dat$type = factor(dat$type, levels = c('Tumor','untreated','one day','two days','four days','nine days','eleven days'))
  dat = dat[order(dat$type),]
  colors = c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#4575b4','#91bfdb')
  colors = rev(colors)
  names(colors) = c('Tumor','untreated','one day','two days','four days','nine days','eleven days')
  require(ggplot2)
  p4 = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
    scale_size_manual(values=c(`Tumor`=1, `Cell`=0.75)) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.15),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values=colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    annotate("text", x=-6, y=2.5, label= "lung",size = 4)+
    guides(fill = guide_legend(ncol = 2, title = NULL))+
    theme(legend.background = element_blank())
  return(p4)
}

Fig2_scSkin = function(obj){

  subobj = slot(obj, "dimRe")
  subobj = subobj[["correct"]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  # dat$type = c(rep("Tumor", nrow(tcga_scSkin@pData@bulk)), rep("Cell", nrow(tcga_scSkin@pData@cell)))

  dat$type = obj@pData@cell$drug[match(dat$id,
                                               tcga_scSkin@pData@cell$id)]
  dat$type[is.na(dat$type)] = "Tumor"

  dat = dat[order(dat$Type),]

  colors = c('#6a3d9a','#1b9e77','#d95f02')
  names(colors) = c('Tumor','untreated','vemurafenib')

  require(ggplot2)
  p = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
    scale_size_manual(values=c(`Tumor`=1.25, `Cell`=0.5)) +
    theme_classic() +
    theme(legend.position = c(0.1, 0.15),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(size=FALSE, color=FALSE) +
    scale_fill_manual(values = colors) +
    xlab("UMAP 1") +
    ylab("UMAP 2")+
    annotate("text", x=-3, y=3, label= "skin",size = 4)+
    guides(fill = guide_legend(title = NULL)) +
    theme(legend.background = element_blank())
  return(p)
}

##------------------------assist for FigureS2 visualization------------------

FigS2_scProstate = function(obj, assays = "correct", title = "text"){
  subobj = slot(obj, "dimRe")
  subobj = subobj[[assays]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  # dat$type = c(rep("Tumor", nrow(tcga_scProstate@pData@bulk)), rep("Cell", nrow(tcga_scProstate@pData@cell)))

  # dat = dat[order(dat$type),]
  dat$type = dat$Type
  dat$Type = obj@pData@cell$Subtype[match(dat$id,
                                          obj@pData@cell$id)]
  dat$Type[is.na(dat$Type)] = "Tumor"
  dat$Type = stringr::str_replace_all(dat$Type, " ", "")

  require(ggplot2)
  if (assays == "correct") {
    p = ggplot(dat, aes(x, y, fill = Type, size=type, color = type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='white',`Cell`='black')) +
      scale_size_manual(values=c(`Tumor`=1, `Cell`=1.25)) +
      theme_classic() +
      theme(legend.position = 'bottom',
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(fill=FALSE, color=FALSE) +
      scale_fill_manual(values=c(`DU145_Resistant` = "#6a3d9a",
                                 `DU145_Sensitive` = "#fb9a99",
                                 `Tumor` = "#1f78b4",
                                 `PC3_Resistant` = "#33a02c",
                                 `PC3_Sensitive` = "#ff7f00")) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      annotate("text", x=-3, y=11, label= "DU145 resistant",size = 4)+
      annotate("text", x=-13, y=7, label= "DU145 sensitive",size = 4)+
      annotate("text", x=4, y=6, label= "prostate",size = 4)+
      annotate("text", x=0, y=-13, label= "PC3 resistant",size = 4)+
      annotate("text", x=-10, y=-8, label= "PC3 sensistive",size = 4)+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }else{
    p = ggplot(dat, aes(x, y, fill = Type, size=type, color = type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='white',`Cell`='black')) +
      scale_size_manual(values=c(`Tumor`=1, `Cell`=1.25)) +
      theme_classic() +
      theme(legend.position = 'bottom',
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(fill=FALSE, color=FALSE) +
      scale_fill_manual(values=c(`DU145_Resistant` = "#6a3d9a",
                                 `DU145_Sensitive` = "#fb9a99",
                                 `Tumor` = "#1f78b4",
                                 `PC3_Resistant` = "#33a02c",
                                 `PC3_Sensitive` = "#ff7f00")) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      annotate("text", x=-5, y=-8, label= "DU145 resistant",size = 4)+
      annotate("text", x=-8, y=-15, label= "DU145 sensitive",size = 4)+
      annotate("text", x=6, y=5, label= "prostate",size = 4)+
      annotate("text", x=-4, y=13, label= "PC3 resistant",size = 4)+
      annotate("text", x=-10, y=6, label= "PC3 sensistive",size = 4)+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)
}


FigS2_scLung = function(obj, assays = "correct", title = "text"){

  subobj = slot(obj, "dimRe")
  subobj = subobj[[assays]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")
  #
  # tcga_scLung@pData@cell$drug = paste(tcga_scLung@pData@cell$drug,
  #                                     tcga_scLung@pData@cell$time, sep = "_")
  dat$type = obj@pData@cell$time[match(dat$id,obj@pData@cell$id)]
  dat$type[is.na(dat$type)] = "Tumor"
  dat$type = factor(dat$type, levels = c('Tumor','untreated','one day','two days','four days','nine days','eleven days'))
  dat = dat[order(dat$type),]
  colors = c('#d73027','#fc8d59','#fee090','#ffffbf','#e0f3f8','#4575b4','#91bfdb')
  colors = rev(colors)
  names(colors) = c('Tumor','untreated','one day','two days','four days','nine days','eleven days')
  require(ggplot2)
  if (assays == "correct") {
    p4 = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
      scale_size_manual(values=c(`Tumor`=1, `Cell`=0.75)) +
      theme_classic() +
      theme(legend.position = c(0.8, 0.8),
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(size=FALSE, color=FALSE) +
      scale_fill_manual(values=colors) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      annotate("text", x=-3, y=5, label= "lung",size = 4)+
      guides(fill = guide_legend(ncol = 2, title = NULL))+
      theme(legend.background = element_blank())+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }else{
    p4 = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
      scale_size_manual(values=c(`Tumor`=1, `Cell`=0.75)) +
      theme_classic() +
      theme(legend.position = c(0.2, 0.15),
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(size=FALSE, color=FALSE) +
      scale_fill_manual(values=colors) +
      xlab("UMAP 1") +
      ylab("UMAP 2") +
      annotate("text", x=-10, y=3, label= "lung",size = 4)+
      guides(fill = guide_legend(ncol = 2, title = NULL))+
      theme(legend.background = element_blank())+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }

  return(p4)
}

FigS2_scSkin = function(obj, assays = "correct", title = "text"){

  subobj = slot(obj, "dimRe")
  subobj = subobj[[assays]]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, "umap"))
  colnames(dat)[c(4, 5)] = c("x", "y")

  # dat$type = c(rep("Tumor", nrow(tcga_scSkin@pData@bulk)), rep("Cell", nrow(tcga_scSkin@pData@cell)))

  dat$type = obj@pData@cell$drug[match(dat$id,
                                       obj@pData@cell$id)]
  dat$type[is.na(dat$type)] = "Tumor"

  dat = dat[order(dat$Type),]

  colors = c('#6a3d9a','#1b9e77','#d95f02')
  names(colors) = c('Tumor','untreated','vemurafenib')

  require(ggplot2)
  if (assays == "raw") {
    p = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
      scale_size_manual(values=c(`Tumor`=1.25, `Cell`=0.5)) +
      theme_classic() +
      theme(legend.position = c(0.25, 0.15),
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(size=FALSE, color=FALSE) +
      scale_fill_manual(values = colors) +
      xlab("UMAP 1") +
      ylab("UMAP 2")+
      annotate("text", x=10, y=-4, label= "skin",size = 4)+
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.background = element_blank())+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))

  }else{
    p = ggplot(dat, aes(x, y, fill = type, size=Type, color = Type)) +
      geom_point(pch=21, alpha=0.7)  +
      scale_color_manual(values=c(`Tumor`='black',`Cell`='white')) +
      scale_size_manual(values=c(`Tumor`=1.25, `Cell`=0.5)) +
      theme_classic() +
      theme(legend.position = c(0.1, 0.15),
            text=element_text(size=8),
            legend.margin =margin(0,0,0,0)) +
      guides(size=FALSE, color=FALSE) +
      scale_fill_manual(values = colors) +
      xlab("UMAP 1") +
      ylab("UMAP 2")+
      annotate("text", x=0, y=3, label= "skin",size = 4)+
      guides(fill = guide_legend(title = NULL)) +
      theme(legend.background = element_blank())+
      ggtitle(title)+
      theme(plot.title = element_text(hjust = 0.5))
  }

  return(p)
}

# -----------------------assist for figure3 visualization---------------------------
Figure3_StromImmBox = function(stromImm, obj){
  require(ggplot2)
  require(ggpubr)
  purity = data.frame(raw = stromImm$raw[1,],
                       correct = stromImm$correct[1,])

  immune = data.frame(raw = stromImm$raw[2,],
                      correct = stromImm$correct[2,])

  bulk_ann = obj@pData@bulk
  colnames(bulk_ann)[1] = "sampleID"
  blood = subset(bulk_ann, lineage == "blood")
  blood_strom = purity[blood$sampleID,]
  blood_immune = immune[blood$sampleID,]

  tissue = subset(bulk_ann, lineage != "blood")
  tissue_strom = purity[tissue$sampleID,]
  tissue_immune = immune[tissue$sampleID,]

  table(bulk_ann$lineage)

  dat = reshape2::melt(blood_strom)
  p1 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Blood Purity Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(blood_immune)
  p2 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Blood Immune Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(tissue_strom)
  p3 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Tissue Stromal Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(tissue_immune)
  p4 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Tissue Immune Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  return(list(blood_stromal = p1, blood_immune = p2,
              tissue_stromal = p3, tissue_immune = p4))
}



Figure3_DrawGseaTable = function(DeGsea){
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)

  dat = DeGsea %>% arrange(padj) %>% select("pathway", "padj", "ES")
  dat = dat[1:5,]
  colnames(dat) = c("pathway", "adjusted p-val", "NES")
  dat$`adjusted p-val` = round(dat$`adjusted p-val`,6)
  dat$NES = round(dat$NES,2)
  dat = as.data.frame(dat)

  p = ggplot(dat) +
    annotate(geom = "table",
             label = list(dat),
             x = 0,
             y = 0)+
    theme_classic() +
    theme(panel.grid=element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank())

  return(p)
}


Figure3_PlotCellClass = function(cell_ann, type = "seq"){
  # dat = subset(cell_ann,  cell_ann$lineage %in% cell_ann$seu_cluster)
  dat = cell_ann
  dat = table(dat$lineage,  dat$seu_cluster) %>% data.frame()
  dat = reshape2::dcast(dat,  Var1~Var2)
  dat = tibble::column_to_rownames(dat,  var = "Var1")

  row.names(dat) = stringr::str_replace_all(row.names(dat), "_", " ")
  common = intersect(row.names(dat), colnames(dat))
  dat = dat[common,common]
  row.names(dat) = colnames(dat)
  # dat[is.na(dat)] = 0

  sum_num = apply(dat,  1,  sum)
  proportion_dat = apply(dat,  2,  function(x){x/sum_num})

  dat = reshape2::melt(proportion_dat)
  colnames(dat) = c("x",  "y",  "propotion")
  dat$propotion[is.nan(dat$propotion)] = 0
  gradCols = c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7')

  if(type == "seq"){
    p = ggplot(dat,  aes(x = x,  y = y,  fill = propotion)) +
      geom_tile(color = "white") +
      theme(axis.title = element_blank()) +
      theme(axis.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45,  hjust = 1)) +
      scale_fill_gradientn(colours = rev(gradCols))}else{
    p = ggplot(dat,  aes(x = x,  y = y,  fill = propotion)) +
      geom_tile(color = "white") +
      theme(axis.title = element_blank()) +
      theme(axis.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45,  hjust = 1)) +
      scale_fill_gradientn(colours = rev(gradCols)) +
      theme(legend.position = "none")

  }


  return(p)
}


Fig3_CellTrace = function(obj, cell_ann){

  cell_trace = readRDS("../log/cell_cytotrace.rds")
  cell_res = obj@dimRe$mnn_svd_correct@umap
  row.names(cell_res) = stringr::str_replace_all(row.names(cell_res), "[.]", "-")
  cell_res = cell_res[obj@pData@cell$sampleID,]
  cell_res$class = cell_ann$seu_cluster
  cell_res$Differentiation = cell_trace$CytoTRACE


  dat = cell_res
  gradCols = c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')
  p3 = ggplot(dat, aes(x = UMAP_1, y = UMAP_2, color = Differentiation)) +
    geom_point() +
    theme_classic()+
    scale_color_gradientn(colours = gradCols)

  diff_grad = cell_trace$CytoTRACE %>% data.frame()

  diff_grad$mixture = ifelse(cell_ann$mixture == "mixture", "mixture", "non mixture")
  # diff_grad$mixture = cell_ann$mixture
  cell_ann$lineage = stringr::str_replace_all(cell_ann$lineage, "_", " ")
  diff_grad$type = ifelse(cell_ann$seu_cluster == cell_ann$lineage, "matching", "mismatching")
  colnames(diff_grad)[1] = "Differentiation"

  dat = diff_grad

  dat$mixture = factor(dat$mixture, levels = c("non mixture", "mixture"))
  p4 = ggplot(dat, aes(x = mixture, y = Differentiation, color = mixture)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Differentiation Index") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`non mixture` = "#33a02c", `mixture` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)

  p5 = ggplot(dat, aes(x = type, y = Differentiation, color = type)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Differentiation Index") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`matching` = "#33a02c", `mismatching` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)

  return(list(DimDiff = p3, mixture = p4, matching = p5))
}


Fig3_Stemness = function(obj){

  tcga_stem = data.table::fread("../data/stem_ness.txt")

  bulk_ann = obj@pData@bulk
  bulk_ann$seu_cluster = seu@active.ident[match(bulk_ann$sampleID,
                                                           stringr::str_replace_all(colnames(seu), "[.]", "-"))]
  # bulk_ann$seu_cluster = seu_cluster$seu_cluster_ann[match(bulk_ann$sampleID,
  #                                                          seu_cluster$sampleID)]
  bulk_ann$stemness = tcga_stem$RNAss[match(bulk_ann$sampleID, tcga_stem$sample)]

  dat = data.frame(lineage = bulk_ann$lineage, class = bulk_ann$seu_cluster, stemness = bulk_ann$stemness)
  dat = na.omit(dat)
  # dat$type = ifelse(dat$class=="mixture", "matching", "mismatching")
  dat$type = ifelse(dat$class == dat$lineage, "matching", "mismatching")
  table(dat$type)

  p = ggplot(dat, aes(x = type, y = stemness, color = type)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Stemness Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`matching` = "#33a02c", `mismatching` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)
  return(p)
}

# -----------------------assist for figureS3 visualization---------------------------
FigureS3_DrawGseaTable = function(DeGsea, n = 5){
  require(ggplot2)
  require(dplyr)
  require(ggpmisc)

  dat = DeGsea %>% arrange(padj) %>% select("pathway", "padj", "ES")
  dat = dat[1:n,]
  colnames(dat) = c("pathway", "adjusted p-val", "NES")
  dat$`adjusted p-val` = round(dat$`adjusted p-val`,6)
  dat$NES = round(dat$NES,2)
  dat = as.data.frame(dat)

  p = ggplot(dat) +
    annotate(geom = "table",
             label = list(dat),
             x = 0,
             y = 0)+
    theme_classic() +
    theme(panel.grid=element_blank())+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank())

  return(p)
}


FigureS3_PlotCellClass = function(cell_ann){
  dat = subset(cell_ann,  cell_ann$lineage %in% cell_ann$seu_cluster)
  dat = table(dat$lineage,  dat$seu_cluster) %>% data.frame()
  dat = reshape2::dcast(dat,  Var1~Var2)
  dat = tibble::column_to_rownames(dat,  var = "Var1")

  sum_num = apply(dat,  1,  sum)
  proportion_dat = apply(dat,  2,  function(x){x/sum_num})

  dat = reshape2::melt(proportion_dat)
  colnames(dat) = c("x",  "y",  "propotion")

  gradCols = c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7')

  p = ggplot(dat,  aes(x = x,  y = y,  fill = propotion)) +
    geom_tile(color = "white") +
    theme(axis.title = element_blank()) +
    theme(axis.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45,  hjust = 1)) +
    scale_fill_gradientn(colours = rev(gradCols))

  return(p)
}

FigureS3_StromImmBox = function(stromImm, obj){

  require(ggplot2)
  require(ggpubr)
  stromal = data.frame(raw = stromImm$raw[1,],
                       correct = stromImm$correct[1,])

  immune = data.frame(raw = stromImm$raw[2,],
                      correct = stromImm$correct[2,])

  bulk_ann = obj@pData@bulk

  blood = subset(bulk_ann, lineage == "blood")
  blood_strom = stromal[blood$id,]
  blood_immune = immune[blood$id,]

  tissue = subset(bulk_ann, lineage != "blood")
  tissue_strom = stromal[tissue$id,]
  tissue_immune = immune[tissue$id,]

  table(bulk_ann$lineage)

  dat = reshape2::melt(blood_strom)
  p1 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Blood Stromal Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(blood_immune)
  p2 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Blood Immune Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(tissue_strom)
  p3 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Tissue Stromal Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(tissue_immune)
  p4 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Tissue Immune Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)

  return(list(blood_stromal = p1, blood_immune = p2,
              tissue_stromal = p3, tissue_immune = p4))
}


FigureS3_StromImmBoxWhole = function(stromImm, obj){

  require(ggplot2)
  require(ggpubr)
  stromal = data.frame(raw = stromImm$raw[1,],
                       correct = stromImm$correct[1,])

  immune = data.frame(raw = stromImm$raw[2,],
                      correct = stromImm$correct[2,])


  dat = reshape2::melt(stromal)
  p1 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    # geom_jitter(width = 0.2, alpha = 0.5) +
    ylab("Blood Stromal Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a")) +
    stat_compare_means(label = "p.signif", label.x = 1.4)

  dat = reshape2::melt(immune)
  p2 = ggplot(dat, aes(x = variable, y = value, color = variable)) +
    geom_violin() +
    geom_boxplot(width = 0.25, coef = 1e30) +
    ylab("Blood Immune Score") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
          legend.position = "none") +
    scale_color_manual(values = c(`raw` = "#33a02c", `correct` = "#6a3d9a"))+
    stat_compare_means(label = "p.signif", label.x = 1.4)


  return(list(stromal = p1, immune = p2))
}


# ------------------------------assit for Fig5 visualization-----------------------
FigS2_landscape = function(obj){

  subobj = slot(obj, "dimRe")
  subobj = subobj[['mnn_svd_correct']]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, 'umap'))
  colnames(dat)[c(4, 5)] = c("x", "y")

  dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))


  require(ggplot2)

  tissue_colors = get_tissue_colors()
  tissue_colors[46] = c(`other` = "#969696")

  dat$Lineage[!dat$Lineage %in% c("skin", "lung", "breast")] = "other"

  p = ggplot(dat, aes(x, y, fill=Lineage, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.25, `Miss`=0.5,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.80),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE) +
    scale_fill_manual(values=tissue_colors)+
    annotate("text", x=-6.5, y=-10.5, label= "breast",size = 3)+
    annotate("text", x=5, y=-10, label= "skin",size = 3)+
    annotate("text", x=-2.5, y=5, label= "lung",size = 3)+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())

  return(p)

}
#
Figure5_subtype_lineage_immune_breast = function(obj, cancer="breast", resolution= 0.8){
  ann = subset(tcga_ccle@pData@comb, Lineage == "breast")
  row.names(ann) = ann$id
  breast_dat = cbind(tcga_ccle@assay$mnn$correct@bulk, tcga_ccle@assay$mnn$correct@cell)
  colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
  breast_dat = breast_dat[,ann$id]

  seu = CreatSeuObj(exp_mat = breast_dat, ann = ann)

  seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                               dims = 1:70,
                               k.param = 20,
                               force.recalc = TRUE,
                               verbose = FALSE)

  seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                               resolution = resolution)

  id = colnames(seu)
  subclass = seu$seurat_clusters

  subdat = data.frame(rawid = id, subclass = subclass)
  subdat$id = stringr::str_sub(subdat$rawid, 1, 12)
  subdat = dplyr::left_join(subdat, immune, by = "id")

  dim = tcga_ccle@dimRe$mnn_correct@umap
  dim$rawid = stringr::str_replace_all(row.names(dim), "[.]", "-")

  subdat = dplyr::left_join(subdat, dim, by = "rawid")

  subtype = data.frame(rawid = c(tcga_ccle@pData@bulk$sampleID,  tcga_ccle@pData@cell$sampleID),
                       subtype = c(tcga_ccle@pData@bulk$subtype, tcga_ccle@pData@cell$subtype))

  subdat = dplyr::left_join(subdat, subtype, by = "rawid")

  if (cancer == "breast") {
    subdat$subtype[subdat$subtype == "HER2 amp"] = "HER2-enriched"
    subdat$subtype[subdat$subtype == "luminal HER2 amp"] = "HER2-enriched"
    subdat$subtype[subdat$subtype == "luminal A"] = "luminal"
    subdat$subtype[subdat$subtype == "luminal B"] = "luminal"
  }
  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  subdat$MFP[is.na(subdat$MFP) & subdat$type != "Cell"] = "miss"
  subdat$MFP[subdat$type == "Cell"] = "Cell"
  subdat = subset(subdat, (!is.na(MFP)))

  require(ggplot2)
  require(ggforce)

  # subdat = subset(subdat, seuAnn == "Mixture")
  subdat$MFP = as.character(subdat$MFP)
  subdat$MFP[subdat$MFP=="D"] = "Depleted"
  subdat$MFP[subdat$MFP=="F"] = "Fibrotic"
  subdat$MFP[subdat$MFP=="IE"] = "Immune Enriched"
  subdat$MFP[subdat$MFP=="IE/F"] = "Immune Enriched Fibrotic"
  subdat$MFP = factor(subdat$MFP, levels = c("miss", "Depleted","Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell"))
  p3 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subclass, size=type, color = type
                          ,shape = subtype
  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`basal` = 21,
                                  `basal A` = 22,
                                  `basal B` = 24,
                                  `HER2-enriched` = 24,
                                  `luminal` = 21,
                                  `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE, fill = FALSE, shape = FALSE) +
    annotate("text", x = -8.8, y = -7, label = "subclass5", size = 3)+
    annotate("text", x = -6, y = -11, label = "subclass1", size = 3)+
    annotate("text", x = -5.5, y = -7, label = "subclass2", size = 3)+
    annotate("text", x = -3, y = -10, label = "subclass0", size = 3)+
    annotate("text", x = -1, y = -5.5, label = "subclass3", size = 3)+
    annotate("text", x = -4.5, y = -3, label = "subclass4", size = 3)

  p4 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type
                          ,shape = subtype)) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`basal` = 21,
                                  `basal A` = 22,
                                  `basal B` = 24,
                                  `HER2-enriched` = 24,
                                  `luminal` = 21,
                                  `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    scale_fill_manual(values = c(`basal` = "#ff7f00",
                                 `basal A` = "#ff7f00",
                                 `basal B` = "#ff7f00",
                                 `HER2-enriched` = "#33a02c",
                                 `luminal` = "#6a3d9a",
                                 `normal` = "#a6cee3")) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.3),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE)+
    theme(legend.background = element_blank())+
    guides(shape = guide_legend(ncol = 2))
  # +
  #   annotate("text", x = -5.5, y = -7, label = "luminal", size = 3)+
  #   annotate("text", x = -8, y = -10, label = "luminal HER2-enriched", size = 3)+
  #   annotate("text", x = -4, y = -4.5, label = "basal", size = 3)+
  #   annotate("text", x = 5, y = 5, label = "basal B", size = 3)+
  #   annotate("text", x = -1, y = -5.5, label = "basal A", size = 3)

  dat = table(subdat$subclass, subdat$MFP) %>% data.frame()
  dat$Var1 = paste("subclass", dat$Var1, sep = "")
  dat = subset(dat, !Var2 %in% c("miss", "Cell"))
  dat = split(dat, dat$Var1)

  GetPercent = function(x){
    sumFreq = sum(x$Freq)
    x$Freq = x$Freq/sumFreq*100
    return(x)
  }

  dat = lapply(dat, GetPercent)
  dat = do.call(rbind, dat)
  colnames(dat)[2] = "Type"

  p5 = ggplot(dat, aes(x = Var1, y = Freq, fill = Type)) +
    geom_bar(position="stack",
             stat="identity")+
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))+
    ylab("Precentage (%)") +
    scale_fill_manual(values=c(`Depleted` = "#984ea3", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#377eb8"))+
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0),
          legend.key.size = unit(0.25, "cm"))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(fill = guide_legend(nrow = 2))

  return(list(lineage = p4, subclass_dim = p3, immune_bar = p5))
}


Fig5_plotPathHeatmap_breast = function(enrich.res){
  require(dplyr)
  dat = enrich.res[,3:ncol(enrich.res)] %>% t()
  colnames(dat) = enrich.res$id

  scale_data = function(x){
    # x = (x-mean(x))/sd(x)
    # x = (x-min(x))/(max(x) - min(x))
    x = x - mean(x)
    x = x / max(abs(x))
    return(x)
  }
  dat = apply(dat, 1, scale_data)
  dat = t(dat)

  Cluster = enrich.res$subclass

  sets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")

  columnAnn = data.frame(id = colnames(dat),
                         clusters = as.character(Cluster))

  rowAnn = data.frame(Sig = row.names(dat),
                      type = sets$Function[match(row.names(dat), sets$Signature)])

  library(ComplexHeatmap)
  library(RColorBrewer)

  # define colors
  # col_fun = c("#756bb1", "#fd8d3c", "#de2d26")
  col_fun = brewer.pal(5, "PuOr")
  col_fun = colorRampPalette(rev(col_fun))(100)
  # col_fun = colorRampPalette(col_fun)(100)

  datasets_col = colorRampPalette(brewer.pal(10, "Paired"))(12)
  names(datasets_col) = unique(columnAnn$datasets)

  # annotation
  column_ha = HeatmapAnnotation(Clusters = columnAnn$clusters,
                                col = list(Clusters = c("0" = "#ec7014",
                                                        "1" = "#88419d",
                                                        "2" = "#238b45",
                                                        "3" = "#e31a1c",
                                                        "4" = "#1f78b4",
                                                        "5" = "#fb9a99"),
                                           Datasets = datasets_col),
                                show_legend = FALSE)

  row_ha = rowAnnotation(Type = rowAnn$type,
                         col = list(
                           Type = c("Angiogenesis_Fibrosis" = "#ec7014",
                                    "Anti_Tumor_Micorenvironment" = "#e7298a",
                                    "Malignant_Cell_Properties" = "#238b45",
                                    "Pro_Tumor_Microenvironment" = "#88419d")
                         ),
                         show_annotation_name = FALSE)


  dat = dat[rowAnn$Sig,columnAnn$id]

  rowPath = split(rowAnn$Sig, rowAnn$type)

  rowType = rowAnn$type

  p = Heatmap(dat,
              name = "Activity",
              col = col_fun,
              cluster_columns = FALSE,
              clustering_distance_rows  = "pearson",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names  = FALSE,
              show_row_names = TRUE,
              top_annotation = column_ha,
              column_split = columnAnn$clusters,
              row_split = rowType,
              left_annotation = row_ha,
              row_title = NULL)
  return(p)
}


Figure5_subtype_lineage_immune_breast_enrich = function(obj, cancer="breast", resolution= 0.8){

  if (!file.exists("../log/breast_subtype.rds")) {
    ann = subset(tcga_ccle@pData@comb, Lineage == "breast")
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$mnn$correct@bulk, tcga_ccle@assay$mnn$correct@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    seu = CreatSeuObj(exp_mat = breast_dat, ann = ann)

    seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                                 dims = 1:70,
                                 k.param = 20,
                                 force.recalc = TRUE,
                                 verbose = FALSE)

    seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                                 resolution = resolution)

    id = colnames(seu)
    subclass = seu$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(tcga_ccle@pData@comb, Lineage == "breast")
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$raw@bulk, tcga_ccle@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/breast_subtype.rds")
  }else{
    subdat = readRDS("../log/breast_subtype.rds")
  }


  return(subdat)
}

# ------------------------------assit for Figure4 visualization-----------------------
Fig4_pairImm = function(obj, immune){
  require(dplyr)
  require(ggplot2)
  pair = obj@assay$mnn$correct@pair
  pair = unique(pair$targ_ID)
  pair = stringr::str_replace_all(pair, "[.]", "-")
  lineage = obj@pData@bulk$lineage[match(pair, obj@pData@bulk$sampleID)]
  pair = stringr::str_sub(pair, 1, 12)
  immType = immune$MFP[match(pair, immune$id)]
  immType = data.frame(id = pair, lineage = lineage, Immtype = immType)

  # tmp = table(immType$lineage, immType$Immtype)
  # tmp = as.data.frame(tmp)
  # tmp = reshape2::dcast(tmp, Var1 ~ Var2)

  dat = table(immType$Immtype) %>% data.frame()
  colnames(dat)[1] = "immType"
  dat$immType = as.character(dat$immType)
  dat$immType[dat$immType=="D"] = "Depleted"
  dat$immType[dat$immType=="F"] = "Fibrotic"
  dat$immType[dat$immType=="IE"] = "Immune Enriched"
  dat$immType[dat$immType=="IE/F"] = "Immune Enriched Fibrotic"

  p = ggplot(data=dat, aes(x=immType, y=Freq, fill = immType)) +
    geom_bar(stat="identity")+
    geom_text(aes(label=Freq), vjust=-0.2, size=3)+
    theme_classic() +
    ylab("Mutual Near Neighbour Pairs") +
    theme(axis.title.x = element_blank(),
          legend.position = "none",
          text=element_text(size=8)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold")) +
    scale_fill_manual(values = c(`Depleted` = "#6a3d9a",`Fibrotic`="#33a02c",
                                 `Immune Enriched` = "#ff7f00",
                                 `Immune Enriched Fibrotic` = "#fb9a99"))+
    theme(plot.margin = unit(c(0.25, 1.5, 0, 1), "cm"))
  return(p)
}

Fig4_ImmDim = function(obj, immune){

  subobj = slot(obj, "dimRe")
  subobj = subobj[['mnn_svd_correct']]
  tissue_colors = get_tissue_colors()

  dat = cbind(obj@pData@comb, slot(subobj, 'umap'))
  colnames(dat)[c(4, 5)] = c("x", "y")

  dat$type = c(rep("Tumor", nrow(obj@pData@bulk)), rep("Cell", nrow(obj@pData@cell)))
  dat$id = stringr::str_sub(dat$id, 1, 12)
  dat = dplyr::left_join(dat, immune, by = "id")

  # subdat = subset(dat, (type == "Cell" | !is.na(Response)))
  dat$type[is.na(dat$MFP) & dat$type != "Cell"] = "Miss"
  dat$MFP[dat$type == "Cell"] = "Cell"
  # subdat = subset(dat, (!is.na(MFP)))
  subdat = dat

  require(ggplot2)
  # subdat = subset(subdat, seuAnn == "Mixture")
  subdat$MFP = as.character(subdat$MFP)
  subdat$MFP[subdat$MFP=="D"] = "Depleted"
  subdat$MFP[subdat$MFP=="F"] = "Fibrotic"
  subdat$MFP[subdat$MFP=="IE"] = "Immune Enriched"
  subdat$MFP[subdat$MFP=="IE/F"] = "Immune Enriched Fibrotic"
  subdat$MFP = factor(subdat$MFP, levels = c("Miss", "Depleted","Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell"))

  subdat = subdat[order(subdat$MFP), ]
  subdat1 = data.frame(subdat)
  # subdat1 = subset(subdat, MFP != "Cell")
  # subdat1 = subset(subdat1, Lineage == "skin")
  subdat1$Type = as.character(subdat1$MFP)
  subdat1$Type[is.na(subdat1$Type)] = "other"
  # subdat1$Type[!subdat1$Lineage %in% c("skin", "lung", "breast")] = "other"
  subdat1 = subset(subdat1, Type != "Cell")
  subdat1$Type[subdat1$Type == "other"] = "Miss"
  p1 = ggplot(subdat1, aes(x, y, fill=Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Miss`=0.1,
                               `Tumor` = 0.75)) +
    theme_classic() +
    theme(legend.position = c(0.3, 0.85),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE) +
    scale_fill_manual(values=c(`Cell`="#41b6c4", `Miss`="#cccccc",
                               `Depleted` = "#984ea3", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#377eb8",
                               `color` = "grey"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())+
    annotate("text", x=-6.5, y=-10.5, label= "breast",size = 3)+
    annotate("text", x=-10.5, y=0, label= "blood",size = 3)+
    annotate("text", x=1, y=-12.5, label= "thymus",size = 3) +
    annotate("text", x= 1.5, y=-11, label= "thyroid",size = 3) +
    annotate("text", x=0, y=15, label= "germ cell",size = 3) +
    annotate("text", x=9, y=7, label= "kidney",size = 3)+
    annotate("text", x=5, y=-10, label= "skin",size = 3)+
    annotate("text", x=4, y=-14, label= "eye",size = 3)+
    annotate("text", x=9, y=0.5, label= "fibroblast",size = 3)+
    annotate("text", x=-4.5, y=2, label= "uterus",size = 3)+
    annotate("text", x=2.5, y=4, label= "pancreas",size = 3)+
    annotate("text", x=-4.5, y=12, label= "liver",size = 3)+
    annotate("text", x=7.5, y=11, label= "ovary",size = 3)+
    annotate("text", x=-4.5, y=4, label= "urinary",size = 3)+
    annotate("text", x=-9, y=6.5, label= "lymphocyte",size = 3)+
    annotate("text", x=-5, y=5.5, label= "plasma cell",size = 3)+
    annotate("text", x=12, y=-11, label= "bone",size = 3)+
    annotate("text", x=-8, y=-6, label= "adrenal",size = 3)+
    annotate("text", x=11, y=4, label= "soft tissue",size = 3)+
    annotate("text", x=11, y=-2, label= "peripheral nervous system",size = 3)+
    annotate("text", x=3, y=10, label= "colorectal",size = 3)+
    annotate("text", x=2.5, y=7.5, label= "gastric",size = 3)+
    annotate("text", x=0, y=-3, label= "upper aerodigestive",size = 3)+
    annotate("text", x=8.5, y=-6.5, label= "central nervous system",size = 3)+
    annotate("text", x=-1.2, y=3, label= "lung",size = 3)+
    annotate("text", x=-2.5, y=-1, label= "esophagus",size = 3)+
    annotate("text", x=-12.5, y=5, label= "prostate",size = 3)

  subdat2 = subset(subdat1, Type != "Cell")
  subdat2 = subset(subdat2, Lineage == "breast")
  subdat2$Type = factor(subdat2$Type, levels = c("Miss", "Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell","Depleted"))

  p2 = ggplot(subdat2, aes(x, y, fill=Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Miss`=0.1,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.4, 0.85),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE) +
    scale_fill_manual(values=c(`Cell`="#41b6c4", `Miss`="#cccccc",
                               `Depleted` = "#984ea3", `Fibrotic` = "#cccccc",
                               `Immune Enriched` = "#cccccc",
                               `Immune Enriched Fibrotic` = "#cccccc"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())

  subdat2$Type = factor(subdat2$Type, levels = c(
                                               "Fibrotic","Miss", "Immune Enriched",
                                               "Immune Enriched Fibrotic", "Cell","Depleted"))

  p3 = ggplot(subdat2, aes(x, y, fill=Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Miss`=0.1,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.4, 0.85),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE)  +
    scale_fill_manual(values=c(`Cell`="#41b6c4", `Miss`="#cccccc",
                               `Depleted` = "#cccccc", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#cccccc",
                               `Immune Enriched Fibrotic` = "#cccccc"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())

  subdat2$Type = factor(subdat2$Type, levels = c("Immune Enriched", "Miss",
                                               "Immune Enriched Fibrotic", "Cell","Depleted",
                                               "Fibrotic"))

  p4 = ggplot(subdat2, aes(x, y, fill=Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Miss`=0.1,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.4, 0.85),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE)  +
    scale_fill_manual(values=c(`Cell`="#41b6c4", `Miss`="#cccccc",
                               `Depleted` = "#cccccc", `Fibrotic` = "#cccccc",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#cccccc"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())

  subdat2$Type = factor(subdat2$Type, levels = c(
    "Immune Enriched Fibrotic", "Miss", "Cell","Depleted",
                                               "Fibrotic", "Immune Enriched"))

  p5 = ggplot(subdat2, aes(x, y, fill=Type, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1, `Miss`=0.1,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.4, 0.85),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE)  +
    scale_fill_manual(values=c(`Cell`="#41b6c4", `Miss`="#cccccc",
                               `Depleted` = "#cccccc", `Fibrotic` = "#cccccc",
                               `Immune Enriched` = "#cccccc",
                               `Immune Enriched Fibrotic` = "#377eb8"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())

  return(list(landscape = p1, sub1 = p2, sub2 = p3,
              sub3 = p4, sub4 = p5))
}

Figure4_subtype_lineage_immune_breast = function(ralign_obj, seu_obj){
  # ann = subset(tcga_ccle@pData@comb, Lineage == "breast")
  # row.names(ann) = ann$id
  # breast_dat = cbind(tcga_ccle@assay$mnn$correct@bulk, tcga_ccle@assay$mnn$correct@cell)
  # colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
  # breast_dat = breast_dat[,ann$id]
  #
  # seu = CreatSeuObj(exp_mat = breast_dat, ann = ann)
  #
  # seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
  #                              dims = 1:70,
  #                              k.param = 20,
  #                              force.recalc = TRUE,
  #                              verbose = FALSE)
  #
  # seu  <- Seurat::FindClusters(seu, reduction = 'pca',
  #                              resolution = resolution)

  id = colnames(seu)
  subclass = seu$seurat_clusters

  subdat = data.frame(rawid = id, subclass = subclass)
  subdat$id = stringr::str_sub(subdat$rawid, 1, 12)
  subdat = dplyr::left_join(subdat, immune, by = "id")

  dim = ralign_obj@dimRe$mnn_svd_correct@umap
  dim$rawid = stringr::str_replace_all(row.names(dim), "[.]", "-")

  subdat = dplyr::left_join(subdat, dim, by = "rawid")


  write.table(subdat, file = "log/immune_breast_subtype.xls", sep = "\t",
              row.names = F, col.names = T)

  subtype = data.frame(rawid = c(tcga_ccle@pData@bulk$sampleID,  tcga_ccle@pData@cell$sampleID),
                       subtype = c(tcga_ccle@pData@bulk$subtype, tcga_ccle@pData@cell$subtype))

  subdat = dplyr::left_join(subdat, subtype, by = "rawid")

  subdat$subtype[subdat$subtype == "HER2 amp"] = "HER2-enriched"
  subdat$subtype[subdat$subtype == "luminal HER2 amp"] = "HER2-enriched"
  subdat$subtype[subdat$subtype == "luminal A"] = "luminal"
  subdat$subtype[subdat$subtype == "luminal B"] = "luminal"

  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  subdat$MFP[is.na(subdat$MFP) & subdat$type != "Cell"] = "miss"
  subdat$MFP[subdat$type == "Cell"] = "Cell"
  subdat = subset(subdat, (!is.na(MFP)))

  require(ggplot2)
  require(ggforce)

  # subdat = subset(subdat, seuAnn == "Mixture")
  subdat$MFP = as.character(subdat$MFP)
  subdat$MFP[subdat$MFP=="D"] = "Depleted"
  subdat$MFP[subdat$MFP=="F"] = "Fibrotic"
  subdat$MFP[subdat$MFP=="IE"] = "Immune Enriched"
  subdat$MFP[subdat$MFP=="IE/F"] = "Immune Enriched Fibrotic"
  subdat$MFP = factor(subdat$MFP, levels = c("miss", "Depleted","Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell"))
  p3 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subclass, size=type, color = type
                          ,shape = subtype
  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`basal` = 21,
                                  `basal A` = 22,
                                  `basal B` = 24,
                                  `HER2-enriched` = 24,
                                  `luminal` = 21,
                                  `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    theme_classic() +
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE, fill = FALSE, shape = FALSE) +
    annotate("text", x = -8.8, y = -7, label = "subclass5", size = 3)+
    annotate("text", x = -6, y = -11, label = "subclass1", size = 3)+
    annotate("text", x = -5.5, y = -7, label = "subclass2", size = 3)+
    annotate("text", x = -3, y = -10, label = "subclass0", size = 3)+
    annotate("text", x = -1, y = -5.5, label = "subclass3", size = 3)+
    annotate("text", x = -4.5, y = -3, label = "subclass4", size = 3)

  p4 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type
                          ,shape = subtype)) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`basal` = 21,
                                  `basal A` = 22,
                                  `basal B` = 24,
                                  `HER2-enriched` = 24,
                                  `luminal` = 21,
                                  `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    scale_fill_manual(values = c(`basal` = "#ff7f00",
                                 `basal A` = "#ff7f00",
                                 `basal B` = "#ff7f00",
                                 `HER2-enriched` = "#33a02c",
                                 `luminal` = "#6a3d9a",
                                 `normal` = "#a6cee3")) +
    theme_classic() +
    theme(legend.position = c(0.8, 0.3),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE)+
    theme(legend.background = element_blank())+
    guides(shape = guide_legend(ncol = 2))
  # +
  #   annotate("text", x = -5.5, y = -7, label = "luminal", size = 3)+
  #   annotate("text", x = -8, y = -10, label = "luminal HER2-enriched", size = 3)+
  #   annotate("text", x = -4, y = -4.5, label = "basal", size = 3)+
  #   annotate("text", x = 5, y = 5, label = "basal B", size = 3)+
  #   annotate("text", x = -1, y = -5.5, label = "basal A", size = 3)

  dat = table(subdat$subclass, subdat$MFP) %>% data.frame()
  dat$Var1 = paste("subclass", dat$Var1, sep = "")
  dat = subset(dat, !Var2 %in% c("miss", "Cell"))
  dat = split(dat, dat$Var1)

  GetPercent = function(x){
    sumFreq = sum(x$Freq)
    x$Freq = x$Freq/sumFreq*100
    return(x)
  }

  dat = lapply(dat, GetPercent)
  dat = do.call(rbind, dat)
  colnames(dat)[2] = "Type"

  p5 = ggplot(dat, aes(x = Var1, y = Freq, fill = Type)) +
    geom_bar(position="stack",
             stat="identity")+
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))+
    ylab("Precentage (%)") +
    scale_fill_manual(values=c(`Depleted` = "#984ea3", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#377eb8"))+
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0),
          legend.key.size = unit(0.25, "cm"))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(fill = guide_legend(nrow = 2))

  return(list(lineage = p4, subclass_dim = p3, immune_bar = p5))
}


Fig4_plotPathHeatmap_breast = function(enrich.res){
  require(dplyr)
  dat = enrich.res[,3:ncol(enrich.res)] %>% t()
  colnames(dat) = enrich.res$id

  scale_data = function(x){
    # x = (x-mean(x))/sd(x)
    # x = (x-min(x))/(max(x) - min(x))
    x = x - mean(x)
    x = x / max(abs(x))
    return(x)
  }
  dat = apply(dat, 1, scale_data)
  dat = t(dat)

  Cluster = enrich.res$subclass

  sets = data.table::fread("../data/ImmGeneSets.txt")

  columnAnn = data.frame(id = colnames(dat),
                         clusters = as.character(Cluster))

  rowAnn = data.frame(Sig = row.names(dat),
                      type = sets$Function[match(row.names(dat), sets$Signature)])

  library(ComplexHeatmap)
  library(RColorBrewer)

  # define colors
  # col_fun = c("#756bb1", "#fd8d3c", "#de2d26")
  col_fun = brewer.pal(5, "PuOr")
  col_fun = colorRampPalette(rev(col_fun))(100)
  # col_fun = colorRampPalette(col_fun)(100)

  datasets_col = colorRampPalette(brewer.pal(10, "Paired"))(12)
  names(datasets_col) = unique(columnAnn$datasets)

  # annotation
  column_ha = HeatmapAnnotation(Clusters = columnAnn$clusters,
                                col = list(Clusters = c("0" = "#ec7014",
                                                        "1" = "#88419d",
                                                        "2" = "#238b45",
                                                        "3" = "#e31a1c",
                                                        "4" = "#1f78b4",
                                                        "5" = "#fb9a99"),
                                           Datasets = datasets_col),
                                show_legend = FALSE)

  row_ha = rowAnnotation(Type = rowAnn$type,
                         col = list(
                           Type = c("Angiogenesis_Fibrosis" = "#ec7014",
                                    "Anti_Tumor_Micorenvironment" = "#e7298a",
                                    "Malignant_Cell_Properties" = "#238b45",
                                    "Pro_Tumor_Microenvironment" = "#88419d")
                         ),
                         show_annotation_name = FALSE)


  dat = dat[rowAnn$Sig,columnAnn$id]

  rowPath = split(rowAnn$Sig, rowAnn$type)

  rowType = rowAnn$type

  p = Heatmap(dat,
              name = "Activity",
              col = col_fun,
              cluster_columns = FALSE,
              clustering_distance_rows  = "pearson",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names  = FALSE,
              show_row_names = TRUE,
              top_annotation = column_ha,
              column_split = columnAnn$clusters,
              row_split = rowType,
              left_annotation = row_ha,
              row_title = NULL)
  return(p)
}


Figure4_subtype_lineage_immune_breast_enrich = function(ralign_obj, seu_obj){

  if (!file.exists("../log/breast_subtype.rds")) {

    id = colnames(seu_obj)
    subclass = seu_obj$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(ralign_obj@pData@comb, Lineage == "breast")
    row.names(ann) = ann$id
    breast_dat = cbind(ralign_obj@assay$raw@bulk, ralign_obj@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/breast_subtype.rds")
  }else{
    subdat = readRDS("../log/breast_subtype.rds")
  }


  return(subdat)
}


# ------------------------------assit for Figure6 visualization-----------------------

Figure6_subtype_lineage_immune_skin = function(obj, cancer="skin", resolution= 0.8){

  if (!file.exists("../log/skin_subtype.rds")) {
    ann = subset(tcga_ccle@pData@comb, Lineage == "skin")
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$mnn$correct@bulk, tcga_ccle@assay$mnn$correct@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    seu = CreatSeuObj(exp_mat = breast_dat, ann = ann)

    seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                                 dims = 1:70,
                                 k.param = 20,
                                 force.recalc = TRUE,
                                 verbose = FALSE)

    seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                                 resolution = resolution)

    id = colnames(seu)
    subclass = seu$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(tcga_ccle@pData@comb, Lineage == cancer)
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$raw@bulk, tcga_ccle@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/skin_subtype.rds")
  }else{
    subdat = readRDS("../log/skin_subtype.rds")
  }


  return(subdat)
}


Figure6_subtype_lineage_immune_lung = function(obj, cancer="lung", resolution= 0.8){

  if (!file.exists("../log/lung_subtype.rds")) {
    ann = subset(tcga_ccle@pData@comb, Lineage == "lung")
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$mnn$correct@bulk, tcga_ccle@assay$mnn$correct@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    seu = CreatSeuObj(exp_mat = breast_dat, ann = ann)

    seu <- Seurat::FindNeighbors(seu, reduction = 'pca',
                                 dims = 1:70,
                                 k.param = 20,
                                 force.recalc = TRUE,
                                 verbose = FALSE)

    seu  <- Seurat::FindClusters(seu, reduction = 'pca',
                                 resolution = resolution)

    id = colnames(seu)
    subclass = seu$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(tcga_ccle@pData@comb, Lineage == cancer)
    row.names(ann) = ann$id
    breast_dat = cbind(tcga_ccle@assay$raw@bulk, tcga_ccle@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/lung_subtype.rds")
  }else{
    subdat = readRDS("../log/lung_subtype.rds")
  }


  return(subdat)
}


Fig6_plotPathHeatmap_skin = function(enrich.res){

  dat = enrich.res[,3:ncol(enrich.res)] %>% t()
  colnames(dat) = enrich.res$id

  scale_data = function(x){
    # x = (x-mean(x))/sd(x)
    # x = (x-min(x))/(max(x) - min(x))
    x = x - mean(x)
    x = x / max(abs(x))
    return(x)
  }
  dat = apply(dat, 1, scale_data)
  dat = t(dat)

  Cluster = enrich.res$subclass

  sets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")

  columnAnn = data.frame(id = colnames(dat),
                         clusters = as.character(Cluster))

  rowAnn = data.frame(Sig = row.names(dat),
                      type = sets$Function[match(row.names(dat), sets$Signature)])

  library(ComplexHeatmap)
  library(RColorBrewer)

  # define colors
  # col_fun = c("#756bb1", "#fd8d3c", "#de2d26")
  col_fun = brewer.pal(5, "PuOr")
  col_fun = colorRampPalette(rev(col_fun))(100)
  # col_fun = colorRampPalette(col_fun)(100)

  datasets_col = colorRampPalette(brewer.pal(10, "Paired"))(12)
  names(datasets_col) = unique(columnAnn$datasets)

  # annotation
  column_ha = HeatmapAnnotation(Clusters = columnAnn$clusters,
                                col = list(Clusters = c("0" = "#ec7014",
                                                        "1" = "#88419d",
                                                        "2" = "#238b45",
                                                        "3" = "#e31a1c",
                                                        "4" = "#1f78b4",
                                                        "5" = "#fb9a99",
                                                        "6" = "#fdbf6f",
                                                        "7" = "red",
                                                        "8" = "blue",
                                                        "9" = "yellow",
                                                        "10" = "orange",
                                                        "11" = "brown"),
                                           Datasets = datasets_col),
                                show_legend = FALSE)

  row_ha = rowAnnotation(Type = rowAnn$type,
                         col = list(
                           Type = c("Angiogenesis_Fibrosis" = "#ec7014",
                                    "Anti_Tumor_Micorenvironment" = "#e7298a",
                                    "Malignant_Cell_Properties" = "#238b45",
                                    "Pro_Tumor_Microenvironment" = "#88419d")
                         ),
                         show_annotation_name = FALSE)


  dat = dat[rowAnn$Sig,columnAnn$id]

  rowPath = split(rowAnn$Sig, rowAnn$type)

  rowType = rowAnn$type

  p = Heatmap(dat,
              name = "Activity",
              col = col_fun,
              cluster_columns = FALSE,
              clustering_distance_rows  = "pearson",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names  = FALSE,
              show_row_names = TRUE,
              top_annotation = column_ha,
              column_split = columnAnn$clusters,
              row_split = rowType,
              left_annotation = row_ha,
              row_title = NULL)
  return(p)
}

Fig6_plotPathHeatmap_lung = function(enrich.res){

  dat = enrich.res[,3:ncol(enrich.res)] %>% t()
  colnames(dat) = enrich.res$id

  scale_data = function(x){
    # x = (x-mean(x))/sd(x)
    # x = (x-min(x))/(max(x) - min(x))
    x = x - mean(x)
    x = x / max(abs(x))
    return(x)
  }
  dat = apply(dat, 1, scale_data)
  dat = t(dat)

  Cluster = enrich.res$subclass

  sets = data.table::fread("../data/tcga_cell/ImmGeneSets.txt")

  columnAnn = data.frame(id = colnames(dat),
                         clusters = as.character(Cluster))

  rowAnn = data.frame(Sig = row.names(dat),
                      type = sets$Function[match(row.names(dat), sets$Signature)])

  library(ComplexHeatmap)
  library(RColorBrewer)

  # define colors
  # col_fun = c("#756bb1", "#fd8d3c", "#de2d26")
  col_fun = brewer.pal(5, "PuOr")
  col_fun = colorRampPalette(rev(col_fun))(100)
  # col_fun = colorRampPalette(col_fun)(100)

  datasets_col = colorRampPalette(brewer.pal(10, "Paired"))(12)
  names(datasets_col) = unique(columnAnn$datasets)

  # annotation
  column_ha = HeatmapAnnotation(Clusters = columnAnn$clusters,
                                col = list(Clusters = c("0" = "#ec7014",
                                                        "1" = "#88419d",
                                                        "2" = "#238b45",
                                                        "3" = "#e31a1c",
                                                        "4" = "#1f78b4",
                                                        "5" = "#fb9a99",
                                                        "6" = "#fdbf6f",
                                                        "7" = "#b15928",
                                                        "8" = "#ffff99",
                                                        "9" = "#cab2d6",
                                                        "10" = "#b2df8a",
                                                        "11" = "#a6cee3"),
                                           Datasets = datasets_col),
                                show_legend = FALSE)

  row_ha = rowAnnotation(Type = rowAnn$type,
                         col = list(
                           Type = c("Angiogenesis_Fibrosis" = "#ec7014",
                                    "Anti_Tumor_Micorenvironment" = "#e7298a",
                                    "Malignant_Cell_Properties" = "#238b45",
                                    "Pro_Tumor_Microenvironment" = "#88419d")
                         ),
                         show_annotation_name = FALSE)


  dat = dat[rowAnn$Sig,columnAnn$id]

  rowPath = split(rowAnn$Sig, rowAnn$type)

  rowType = rowAnn$type

  p = Heatmap(dat,
              name = "Activity",
              col = col_fun,
              cluster_columns = FALSE,
              clustering_distance_rows  = "pearson",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names  = FALSE,
              show_row_names = TRUE,
              top_annotation = column_ha,
              column_split = columnAnn$clusters,
              row_split = rowType,
              left_annotation = row_ha,
              row_title = NULL)
  return(p)
}


#
#
FigS2_lineage_lung = function(obj){

  subdat = obj@dimRe$mnn_svd_correct@umap
  subdat$subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype)
  subdat$lineage = obj@pData@comb$Lineage
  subdat$rawid = obj@pData@comb$id
  subdat = subset(subdat, lineage == "lung")
  subdat$subtype[subdat$subtype == "carcinoid"] = "other"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, adenocarcinoma"] = "lung adenocarcinoma"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, adenosquamous carcinoma"] = "lung adenosquamous"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, large cell carcinoma"] = "NSCLC, large cell carcinoma"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, mucoepidermoid carcinoma"] = "NSCLC, mucoepidermoid"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, squamous cell carcinoma"] = "lung squamous cell carcinoma"
  subdat$subtype[subdat$subtype == "pleuropulmonary blastoma"] = "other"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, unspecified"] = "NSCLC, unspecified"

  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  require(ggplot2)
  require(ggforce)

  p = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type,
                          shape = subtype

  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`lung adenocarcinoma` = 21,
                                  `lung squamous cell carcinoma` = 21,
                                  `mesothelioma` = 21,
                                  `NSCLC, large cell carcinoma` = 21,
                                  `small cell lung cancer` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    scale_fill_manual(values = c(`lung adenocarcinoma` = "#fb9a99",
                                 `lung squamous cell carcinoma` = "#ff7f00",
                                 `mesothelioma` = "#33a02c",
                                 `NSCLC, large cell carcinoma` = "#1f78b4",
                                 `small cell lung cancer` = "#fdbf6f")) +
    theme_classic() +
    theme(legend.position = c(0.3, 0.3),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE) +
    guides(shape = guide_legend(ncol = 1))+
    theme(legend.background = element_blank())


  return(p)

}


FigS2_lineage_skin = function(obj){

  subdat = obj@dimRe$mnn_svd_correct@umap
  subdat$subtype = c(obj@pData@bulk$subtype, obj@pData@cell$subtype)
  subdat$lineage = obj@pData@comb$Lineage
  subdat$rawid = obj@pData@comb$id
  subdat = subset(subdat, lineage == "skin")


  subdat$subtype[subdat$subtype == "melanocytic"] = "melanoma"
  subdat$subtype[subdat$subtype == "melanoma, amelanotic"] = "amelanotic"
  subdat$subtype[subdat$subtype == "merkel cell carcinoma"] = "merkel"
  subdat$subtype[subdat$subtype == "squamous cell carcinoma"] = "squamous"


  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  p = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type
                          ,shape = subtype
  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`amelanotic` = 21,
                                  `melanoma` = 21,
                                  `merkel` = 21,
                                  `neural crest-like` = 21,
                                  `squamous` = 21,
                                  `transitory` = 21,
                                  `undifferentiated` = 24))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`= 1.5,
                               `Tumor` = 1.25))+
    scale_fill_manual(values = c(`amelanotic` = "#fdbf6f",
                                 `melanoma` = "#6a3d9a",
                                 `merkel` = "#cab2d6",
                                 `neural crest-like` = "#1f78b4",
                                 `squamous` = "#fb9a99",
                                 `transitory` = "#33a02c",
                                 `undifferentiated` = "#ff7f00")) +
    theme_classic() +
    theme(legend.position = c(0.40, 0.35),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE) +
    guides(shape = guide_legend(ncol = 2)) +
    theme(legend.background = element_blank())


  return(p)

}

#------------------assit for figureS3 visualization-------------------
FigS3_tidy_meta = function(ralign_obj, seu_obj){
  meta = ralign_obj@dimRe$mnn_svd_correct@umap
  meta$seurat_cluster = seu_obj@meta.data$seurat_clusters
  meta$seurat_ann = seu_obj@active.ident
  meta = cbind(ralign_obj@pData@comb,meta)
  meta$subtype = c(ralign_obj@pData@bulk$subtype, ralign_obj@pData@cell$subtype)
  return(meta)
}

FigureS3_mixture_mismatching = function(tcga_meta){
  dat = tcga_meta
  dat$x = dat$UMAP_1
  dat$y = dat$UMAP_2
  dat$seurat_ann = as.character(dat$seu_ann)
  dat$Lineage = stringr::str_replace_all(dat$lineage, "_", " ")
  dat$match = ifelse(dat$Lineage == dat$seurat_ann, "matching", "mismatching")
  dat$match = ifelse(dat$type == "Tumor", "tumor", dat$match)
  dat$type = dat$type
  p3 = ggplot(dat, aes(x, y, fill=match, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.25, `Miss`=0.5,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.80),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE) +
    scale_fill_manual(values= c(`mismatching`="#45a132"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank()) + theme(legend.position = "none")+
    theme(legend.background = element_blank())+
    annotate("text", x=-6.5, y=-10.5, label= "breast",size = 3)+
    annotate("text", x=-10.5, y=0, label= "blood",size = 3)+
    annotate("text", x=1, y=-12.5, label= "thymus",size = 3) +
    annotate("text", x= 1.5, y=-11, label= "thyroid",size = 3) +
    annotate("text", x=0, y=15, label= "germ cell",size = 3) +
    annotate("text", x=9, y=7, label= "kidney",size = 3)+
    annotate("text", x=5, y=-10, label= "skin",size = 3)+
    annotate("text", x=4, y=-14, label= "eye",size = 3)+
    annotate("text", x=9, y=0.5, label= "fibroblast",size = 3)+
    annotate("text", x=-4.5, y=2, label= "uterus",size = 3)+
    annotate("text", x=2.5, y=4, label= "pancreas",size = 3)+
    annotate("text", x=-4.5, y=12, label= "liver",size = 3)+
    annotate("text", x=7.5, y=11, label= "ovary",size = 3)+
    annotate("text", x=-4.5, y=4, label= "urinary",size = 3)+
    annotate("text", x=-9, y=6.5, label= "lymphocyte",size = 3)+
    annotate("text", x=-5, y=5.5, label= "plasma cell",size = 3)+
    annotate("text", x=12, y=-11, label= "bone",size = 3)+
    annotate("text", x=-8, y=-6, label= "adrenal",size = 3)+
    annotate("text", x=11, y=4, label= "soft tissue",size = 3)+
    annotate("text", x=11, y=-2, label= "peripheral nervous system",size = 3)+
    annotate("text", x=3, y=10, label= "colorectal",size = 3)+
    annotate("text", x=2.5, y=7.5, label= "gastric",size = 3)+
    annotate("text", x=0, y=-3, label= "upper aerodigestive",size = 3)+
    annotate("text", x=8.5, y=-6.5, label= "central nervous system",size = 3)+
    annotate("text", x=-1.2, y=3, label= "lung",size = 3)+
    annotate("text", x=-2.5, y=-1, label= "esophagus",size = 3)+
    annotate("text", x=-12.5, y=5, label= "prostate",size = 3)

  # 17,11,20,16,31,54,
  mixture = data.frame(table(dat$seu_cluster, dat$match))
  mixture = reshape2::dcast(mixture, Var1~Var2)
  mixture$ratio = round(mixture$mismatching/mixture$matching, 2)

  dat = tcga_meta
  dat$x = dat$UMAP_1
  dat$y = dat$UMAP_2
  dat$seu_cluster = as.character(dat$seu_cluster)
  dat$mixture = ifelse(dat$seu_cluster %in% c("17", "11", "20", "16", "31", "54"), "mixture", "other")
  dat$type = dat$type
  p4 = ggplot(dat, aes(x, y, fill=mixture, size=type, color = type)) +
    geom_point(pch=21, alpha=0.7)  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.25, `Miss`=0.5,
                               `Tumor` = 1)) +
    theme_classic() +
    theme(legend.position = c(0.2, 0.80),
          text=element_text(size=10),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, type = FALSE, size = FALSE) +
    scale_fill_manual(values= c(`mixture`="#6c55e2"))+
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank()) + theme(legend.position = "none")+
    annotate("text", x=-6.5, y=-10.5, label= "breast",size = 3)+
    annotate("text", x=-10.5, y=0, label= "blood",size = 3)+
    annotate("text", x=1, y=-12.5, label= "thymus",size = 3) +
    annotate("text", x= 1.5, y=-11, label= "thyroid",size = 3) +
    annotate("text", x=0, y=15, label= "germ cell",size = 3) +
    annotate("text", x=9, y=7, label= "kidney",size = 3)+
    annotate("text", x=5, y=-10, label= "skin",size = 3)+
    annotate("text", x=4, y=-14, label= "eye",size = 3)+
    annotate("text", x=9, y=0.5, label= "fibroblast",size = 3)+
    annotate("text", x=-4.5, y=2, label= "uterus",size = 3)+
    annotate("text", x=2.5, y=4, label= "pancreas",size = 3)+
    annotate("text", x=-4.5, y=12, label= "liver",size = 3)+
    annotate("text", x=7.5, y=11, label= "ovary",size = 3)+
    annotate("text", x=-4.5, y=4, label= "urinary",size = 3)+
    annotate("text", x=-9, y=6.5, label= "lymphocyte",size = 3)+
    annotate("text", x=-5, y=5.5, label= "plasma cell",size = 3)+
    annotate("text", x=12, y=-11, label= "bone",size = 3)+
    annotate("text", x=-8, y=-6, label= "adrenal",size = 3)+
    annotate("text", x=11, y=4, label= "soft tissue",size = 3)+
    annotate("text", x=11, y=-2, label= "peripheral nervous system",size = 3)+
    annotate("text", x=3, y=10, label= "colorectal",size = 3)+
    annotate("text", x=2.5, y=7.5, label= "gastric",size = 3)+
    annotate("text", x=0, y=-3, label= "upper aerodigestive",size = 3)+
    annotate("text", x=8.5, y=-6.5, label= "central nervous system",size = 3)+
    annotate("text", x=-1.2, y=3, label= "lung",size = 3)+
    annotate("text", x=-2.5, y=-1, label= "esophagus",size = 3)+
    annotate("text", x=-12.5, y=5, label= "prostate",size = 3)
  return(list(mismatch=p3, mixture=p4))
}

# ------------assist for figureS4 visualization---------------------
FigureS4_subtype_lineage_immune_skin = function(ralign_obj, seu_obj){

  id = colnames(seu_obj)
  subclass = seu_obj$seurat_clusters

  subdat = data.frame(rawid = id, subclass = subclass)
  subdat$id = stringr::str_sub(subdat$rawid, 1, 12)
  subdat = dplyr::left_join(subdat, immune, by = "id")

  dim = ralign_obj@dimRe$mnn_svd_correct@umap
  dim$rawid = stringr::str_replace_all(row.names(dim), "[.]", "-")

  subdat = dplyr::left_join(subdat, dim, by = "rawid")
  write.table(subdat, file = "log/immune_skin_subtype.xls", sep = "\t",
              row.names = F, col.names = T)

  subtype = data.frame(rawid = c(tcga_ccle@pData@bulk$sampleID,  tcga_ccle@pData@cell$sampleID),
                       subtype = c(tcga_ccle@pData@bulk$subtype, tcga_ccle@pData@cell$subtype))

  subdat = dplyr::left_join(subdat, subtype, by = "rawid")


  subdat$subtype[subdat$subtype == "melanocytic"] = "melanoma"
  subdat$subtype[subdat$subtype == "melanoma, amelanotic"] = "amelanotic"
  subdat$subtype[subdat$subtype == "merkel cell carcinoma"] = "merkel"
  subdat$subtype[subdat$subtype == "squamous cell carcinoma"] = "squamous"


  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  subdat$MFP[is.na(subdat$MFP) & subdat$type != "Cell"] = "miss"
  subdat$MFP[subdat$type == "Cell"] = "Cell"
  subdat = subset(subdat, (!is.na(MFP)))

  require(ggplot2)
  require(ggforce)

  # subdat = subset(subdat, seuAnn == "Mixture")
  subdat$MFP = as.character(subdat$MFP)
  subdat$MFP[subdat$MFP=="D"] = "Depleted"
  subdat$MFP[subdat$MFP=="F"] = "Fibrotic"
  subdat$MFP[subdat$MFP=="IE"] = "Immune Enriched"
  subdat$MFP[subdat$MFP=="IE/F"] = "Immune Enriched Fibrotic"
  subdat$MFP = factor(subdat$MFP, levels = c("miss", "Depleted","Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell"))

  p3 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subclass, size=type, color = type)) +
    geom_point(alpha=0.7, pch = 21)+
    # scale_shape_manual(values = c(`basal` = 21,
    #                               `basal A` = 22,
    #                               `basal B` = 24,
    #                               `HER2-enriched` = 24,
    #                               `luminal` = 21,
    #                               `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    theme_classic() +
    theme(legend.position = c(0.2, 0.3),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE) +
    guides(fill = guide_legend(ncol = 2))
  # +
  # annotate("text", x = 7.5, y = -8, label = "subclass2", size = 2)+
  # annotate("text", x = 4, y = -11, label = "subclass1", size = 2)+
  # annotate("text", x = 7.5, y = -12.5, label = "subclass0", size = 2)
  #

  p4 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type
                          ,shape = subtype
  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`amelanotic` = 21,
                                  `melanoma` = 21,
                                  `merkel` = 21,
                                  `neural crest-like` = 21,
                                  `squamous` = 21,
                                  `transitory` = 21,
                                  `undifferentiated` = 24))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`= 1,
                               `Tumor` = 0.75))+
    scale_fill_manual(values = c(`amelanotic` = "#fdbf6f",
                                 `melanoma` = "#6a3d9a",
                                 `merkel` = "#cab2d6",
                                 `neural crest-like` = "#1f78b4",
                                 `squamous` = "#fb9a99",
                                 `transitory` = "#33a02c",
                                 `undifferentiated` = "#ff7f00")) +
    theme_classic() +
    theme(legend.position = c(0.40, 0.35),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE) +
    guides(shape = guide_legend(ncol = 2)) +
    theme(legend.background = element_blank())


  dat = table(subdat$subclass, subdat$MFP) %>% data.frame()
  dat$Var1 = paste("subclass", dat$Var1, sep = "")
  dat = subset(dat, !Var2 %in% c("miss", "Cell"))
  dat = split(dat, dat$Var1)

  GetPercent = function(x){
    sumFreq = sum(x$Freq)
    x$Freq = x$Freq/sumFreq*100
    return(x)
  }

  dat = lapply(dat, GetPercent)
  dat = do.call(rbind, dat)
  colnames(dat)[2] = "Type"

  p5 = ggplot(dat, aes(x = Var1, y = Freq, fill = Type)) +
    geom_bar(position="stack",
             stat="identity")+
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))+
    ylab("Precentage (%)") +
    scale_fill_manual(values=c(`Depleted` = "#984ea3", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#377eb8"))+
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0),
          legend.key.size = unit(0.25, "cm"))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(fill = guide_legend(nrow = 2))
  return(list(lineage = p4, subclass_dim = p3, immune_bar = p5))

}


FigureS4_subtype_lineage_immune_lung = function(ralign_obj, seu_obj){

  id = colnames(seu_obj)
  subclass = seu_obj$seurat_clusters

  subdat = data.frame(rawid = id, subclass = subclass)
  subdat$id = stringr::str_sub(subdat$rawid, 1, 12)
  subdat = dplyr::left_join(subdat, immune, by = "id")

  dim = ralign_obj@dimRe$mnn_svd_correct@umap
  dim$rawid = stringr::str_replace_all(row.names(dim), "[.]", "-")

  subdat = dplyr::left_join(subdat, dim, by = "rawid")

  write.table(subdat, file = "log/immune_lung_subtype.xls", sep = "\t",
              row.names = F, col.names = T)
  subtype = data.frame(rawid = c(tcga_ccle@pData@bulk$sampleID,  tcga_ccle@pData@cell$sampleID),
                       subtype = c(tcga_ccle@pData@bulk$subtype, tcga_ccle@pData@cell$subtype))

  subdat = dplyr::left_join(subdat, subtype, by = "rawid")



  subdat$subtype[subdat$subtype == "carcinoid"] = "other"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, adenocarcinoma"] = "lung adenocarcinoma"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, adenosquamous carcinoma"] = "lung adenosquamous"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, large cell carcinoma"] = "NSCLC, large cell carcinoma"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, mucoepidermoid carcinoma"] = "NSCLC, mucoepidermoid"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, squamous cell carcinoma"] = "lung squamous cell carcinoma"
  subdat$subtype[subdat$subtype == "pleuropulmonary blastoma"] = "other"
  subdat$subtype[subdat$subtype == "non-small cell lung cancer, unspecified"] = "NSCLC, unspecified"


  subdat$type = ifelse(stringr::str_starts(subdat$rawid, "ACH"), "Cell", "Tumor")

  subdat$MFP[is.na(subdat$MFP) & subdat$type != "Cell"] = "miss"
  subdat$MFP[subdat$type == "Cell"] = "Cell"
  subdat = subset(subdat, (!is.na(MFP)))

  require(ggplot2)
  require(ggforce)

  # subdat = subset(subdat, seuAnn == "Mixture")
  subdat$MFP = as.character(subdat$MFP)
  subdat$MFP[subdat$MFP=="D"] = "Depleted"
  subdat$MFP[subdat$MFP=="F"] = "Fibrotic"
  subdat$MFP[subdat$MFP=="IE"] = "Immune Enriched"
  subdat$MFP[subdat$MFP=="IE/F"] = "Immune Enriched Fibrotic"
  subdat$MFP = factor(subdat$MFP, levels = c("miss", "Depleted","Fibrotic", "Immune Enriched",
                                             "Immune Enriched Fibrotic", "Cell"))

  p3 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subclass, size=type, color = type)) +
    geom_point(alpha=0.7, pch = 21)+
    # scale_shape_manual(values = c(`basal` = 21,
    #                               `basal A` = 22,
    #                               `basal B` = 24,
    #                               `HER2-enriched` = 24,
    #                               `luminal` = 21,
    #                               `normal` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1.5,
                               `Tumor` = 1.25))+
    theme_classic() +
    theme(legend.position = c(0.2, 0.35),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE)+
    guides(fill = guide_legend(ncol = 2))+
    theme(legend.background = element_blank())
  # +
  # annotate("text", x = 2.5, y = 5.5, label = "subclass1", size = 2)+
  # annotate("text", x = 2.5, y = -1, label = "subclass2", size = 2)+
  # annotate("text", x = -2.5, y = 5, label = "subclass0", size = 2)+
  # annotate("text", x = -1, y = -5, label = "subclass3", size = 2)+
  # annotate("text", x = -2, y = 1, label = "subclass4", size = 2)+
  # annotate("text", x = 5, y = 1, label = "subclass5", size = 2)+
  # annotate("text", x = 2, y = 1, label = "subclass6", size = 2)


  p4 = ggplot(subdat, aes(UMAP_1, UMAP_2, fill=subtype, size=type, color = type,
                          shape = subtype

  )) +
    geom_point(alpha=0.7)+
    scale_shape_manual(values = c(`lung adenocarcinoma` = 21,
                                  `lung squamous cell carcinoma` = 21,
                                  `mesothelioma` = 21,
                                  `NSCLC, large cell carcinoma` = 21,
                                  `small cell lung cancer` = 21))  +
    scale_color_manual(values=c(`Cell`='black', `Tumor`='white')) +
    scale_size_manual(values=c(`Cell`=1,
                               `Tumor` = 0.75))+
    scale_fill_manual(values = c(`lung adenocarcinoma` = "#fb9a99",
                                 `lung squamous cell carcinoma` = "#ff7f00",
                                 `mesothelioma` = "#33a02c",
                                 `NSCLC, large cell carcinoma` = "#1f78b4",
                                 `small cell lung cancer` = "#fdbf6f")) +
    theme_classic() +
    theme(legend.position = c(0.3, 0.3),
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0)) +
    guides(color=FALSE, size = FALSE) +
    guides(shape = guide_legend(ncol = 1))+
    theme(legend.background = element_blank())


  dat = table(subdat$subclass, subdat$MFP) %>% data.frame()
  dat$Var1 = paste("subclass", dat$Var1, sep = "")
  dat = subset(dat, !Var2 %in% c("miss", "Cell"))
  dat = split(dat, dat$Var1)

  GetPercent = function(x){
    sumFreq = sum(x$Freq)
    x$Freq = x$Freq/sumFreq*100
    return(x)
  }

  dat = lapply(dat, GetPercent)
  dat = do.call(rbind, dat)
  colnames(dat)[2] = "Type"

  p5 = ggplot(dat, aes(x = Var1, y = Freq, fill = Type)) +
    geom_bar(position="stack",
             stat="identity")+
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"))+
    ylab("Precentage (%)") +
    scale_fill_manual(values=c(`Depleted` = "#984ea3", `Fibrotic` = "#4daf4a",
                               `Immune Enriched` = "#ff7f00",
                               `Immune Enriched Fibrotic` = "#377eb8"))+
    theme(legend.position = 'bottom',
          text=element_text(size=8),
          legend.margin =margin(0,0,0,0),
          legend.key.size = unit(0.25, "cm"))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(fill = guide_legend(nrow = 2))

  return(list(lineage = p4, subclass_dim = p3, immune_bar = p5))

}

FigureS4_subtype_lineage_immune_skin_enrich = function(ralign_obj, seu_obj){

  if (!file.exists("../log/skin_subtype.rds")) {

    id = colnames(seu_obj)
    subclass = seu_obj$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(ralign_obj@pData@comb, Lineage == "skin")
    row.names(ann) = ann$id
    breast_dat = cbind(ralign_obj@assay$raw@bulk, ralign_obj@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/skin_subtype.rds")
  }else{
    subdat = readRDS("../log/skin_subtype.rds")
  }


  return(subdat)
}

FigureS4_subtype_lineage_immune_lung_enrich = function(ralign_obj, seu_obj){

  if (!file.exists("../log/lung_subtype.rds")) {

    id = colnames(seu_obj)
    subclass = seu_obj$seurat_clusters

    subdat = data.frame(id = id, subclass = subclass)

    Immsets = data.table::fread("../data/ImmGeneSets.txt")
    Immsets = split(Immsets$Symbol, Immsets$Signature)
    ann = subset(ralign_obj@pData@comb, Lineage == "lung")
    row.names(ann) = ann$id
    breast_dat = cbind(ralign_obj@assay$raw@bulk, ralign_obj@assay$raw@cell)
    colnames(breast_dat) = stringr::str_replace_all(colnames(breast_dat), "[.]", "-")
    breast_dat = breast_dat[,ann$id]

    Immscore = GSVA::gsva(as.matrix(breast_dat), Immsets, method = "ssgsea")
    Immscore = t(Immscore) %>% as.data.frame()
    Immscore$id = row.names(Immscore)
    subdat = dplyr::left_join(subdat, Immscore, by = "id")

    saveRDS(subdat, file = "../log/lung_subtype.rds")
  }else{
    subdat = readRDS("../log/lung_subtype.rds")
  }


  return(subdat)
}

FigS4_plotPathHeatmap_lung = function(enrich.res){

  dat = enrich.res[,3:ncol(enrich.res)] %>% t()
  colnames(dat) = enrich.res$id

  scale_data = function(x){
    # x = (x-mean(x))/sd(x)
    # x = (x-min(x))/(max(x) - min(x))
    x = x - mean(x)
    x = x / max(abs(x))
    return(x)
  }
  dat = apply(dat, 1, scale_data)
  dat = t(dat)

  Cluster = enrich.res$subclass

  sets = data.table::fread("../data/ImmGeneSets.txt")

  columnAnn = data.frame(id = colnames(dat),
                         clusters = as.character(Cluster))

  rowAnn = data.frame(Sig = row.names(dat),
                      type = sets$Function[match(row.names(dat), sets$Signature)])

  library(ComplexHeatmap)
  library(RColorBrewer)

  # define colors
  # col_fun = c("#756bb1", "#fd8d3c", "#de2d26")
  col_fun = brewer.pal(5, "PuOr")
  col_fun = colorRampPalette(rev(col_fun))(100)
  # col_fun = colorRampPalette(col_fun)(100)

  datasets_col = colorRampPalette(brewer.pal(10, "Paired"))(12)
  names(datasets_col) = unique(columnAnn$datasets)

  # annotation
  column_ha = HeatmapAnnotation(Clusters = columnAnn$clusters,
                                col = list(Clusters = c("0" = "#ec7014",
                                                        "1" = "#88419d",
                                                        "2" = "#238b45",
                                                        "3" = "#e31a1c",
                                                        "4" = "#1f78b4",
                                                        "5" = "#fb9a99",
                                                        "6" = "#fdbf6f",
                                                        "7" = "#b15928",
                                                        "8" = "#ffff99",
                                                        "9" = "#cab2d6",
                                                        "10" = "#b2df8a",
                                                        "11" = "#a6cee3"),
                                           Datasets = datasets_col),
                                show_legend = FALSE)

  row_ha = rowAnnotation(Type = rowAnn$type,
                         col = list(
                           Type = c("Angiogenesis_Fibrosis" = "#ec7014",
                                    "Anti_Tumor_Micorenvironment" = "#e7298a",
                                    "Malignant_Cell_Properties" = "#238b45",
                                    "Pro_Tumor_Microenvironment" = "#88419d")
                         ),
                         show_annotation_name = FALSE)


  dat = dat[rowAnn$Sig,columnAnn$id]

  rowPath = split(rowAnn$Sig, rowAnn$type)

  rowType = rowAnn$type

  p = Heatmap(dat,
              name = "Activity",
              col = col_fun,
              cluster_columns = FALSE,
              clustering_distance_rows  = "pearson",
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              show_column_names  = FALSE,
              show_row_names = TRUE,
              top_annotation = column_ha,
              column_split = columnAnn$clusters,
              row_split = rowType,
              left_annotation = row_ha,
              row_title = NULL)
  return(p)
  return(p)
}
