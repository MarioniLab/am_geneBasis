celltype_atlas_E8.5_colors = c("Epiblast" = "#635547", "Primitive Streak" = "#DABE99", "Caudal epiblast" = "#9e6762", "PGC" = "#FACB12",
                           "Anterior Primitive Streak" = "#c19f70", "Notochord" = "#0F4A9C", "Def. endoderm" = "#F397C0", "Gut" = "#EF5A9D",
                           "Nascent mesoderm" = "#C594BF", "Mixed mesoderm" = "#DFCDE4", "Intermediate mesoderm" = "#139992", "Caudal Mesoderm" = "#3F84AA",
                           "Paraxial mesoderm" = "#8DB5CE", "Somitic mesoderm" = "#005579", "Pharyngeal mesoderm" = "#C9EBFB", "Cardiomyocytes" = "#B51D8D",
                           "Allantois" = "#532C8A", "ExE mesoderm" = "#8870ad", "Mesenchyme" = "#cc7818", "Haematoendothelial progenitors" = "#FBBE92",
                           "Endothelium" = "#ff891c", "Blood progenitors 1" = "#f9decf", "Blood progenitors 2" = "#c9a997", "Erythroid1" = "#C72228",
                           "Erythroid2" = "#f79083", "Erythroid3" = "#EF4E22", "NMP" = "#8EC792", "Rostral neurectoderm" = "#65A83E",
                           "Caudal neurectoderm" = "#354E23", "Neural crest" = "#C3C388", "Forebrain/Midbrain/Hindbrain" = "#647a4f", "Spinal cord" = "#CDE088",
                           "Surface ectoderm" = "#f7f79e", "Visceral endoderm" = "#F6BFCB", "ExE endoderm" = "#7F6874", "ExE ectoderm" = "#989898",
                           "Parietal endoderm" = "#1A1A1A", "Doublet" = "black", "Stripped" = "black"
)


celltype_spleen_colors = c("Monocyte" = "#635547", "Plasmablast" = "#DABE99","Follicular B" = "#532C8A", 
                           "MZ B" = "#FACB12","DC2" = "#0F4A9C","Plasma cell" = "#005579","HSC" = "#EF5A9D",
                           "DC1" = "#8870ad", "Ery" = "#EF4E22","memory CD4+ ab T" = "#CDE088", "naive CD4+ ab T" = "gray84",
                           "CD4+ T" = "#CDE088","memory CD8+ ab T" = "plum1",
                           "CD8+ T" = "plum1", "FCGR3A+ NK" = "#1A1A1A", "FCGR3A- NK" = "chartreuse4", "NK" = "chartreuse4",
                           "ILC" = "aquamarine", "gd T" = "#ff891c", "EC" = "#B51D8D", "Mac" = "#FBBE92",
                           "cytotoxic CD8+ ab T" = "#F397C0","Fibroblast" = "#139992", "Unknown" = "#1A1A1A",  "Dividing T" = "azure4"
)


celltype_melanoma_colors = c("other" = "#635547", "Endo." = "#ff891c", "Macro." = "#005579", "B" = "#647a4f", "T" = "aquamarine", 
                             "CAF" = "#EF4E22" , "NK" = "#C9EBFB")

celltype_colors = list("atlas_E8.5" = celltype_atlas_E8.5_colors ,
                       "spleen" = celltype_spleen_colors , 
                       "melanoma" = celltype_melanoma_colors)


celltype_kidney_colors = c("PT" = "#635547",
                           "tIC-CNT" = "#DABE99",
                           "AVR" = "#9e6762",
                           "PC" = "#FACB12",
                           "DVR" = "#c19f70",
                           "mDC" = "#0F4A9C",
                           "Cycling" = "#F397C0",
                           "mTAL" = "#EF5A9D",
                           "aIC" = "#C594BF",
                           "CD8 T" = "#DFCDE4",
                           "NKT" = "#139992",
                           "CD4 T" = "#3F84AA",
                           "PT_VCAM1" = "#8DB5CE",
                           "TAL_unk"= "#005579",
                           "Podo"= "#C9EBFB",
                           "NK"= "#B51D8D",
                           "gEC"= "#532C8A",
                           "Mac/Mono"= "#8870ad",
                           "CNT" = "#cc7818",
                           "DTL" = "#FBBE92",
                           "B cell" = "#ff891c",
                           "cTAL"= "#f9decf",
                           "Fib"= "#c9a997",
                           "bIC" = "#C72228",
                           "DCT"= "#f79083",
                           "ATL"= "#EF4E22",
                           "Myofib"= "#8EC792",
                           "Mmrn1 EC" = "#65A83E",
                           "Msg" = "#354E23"
                           
)




plot_expr_distribution = function(sce , gene , assay = "logcounts" , title = gene){
  if (!gene %in% rownames(sce)){
    stop("Can not find gene in the counts matrix. Ensure that given entry exists.")
  }
  if (!assay %in% c("counts" , "logcounts")){
    stop("Option 'assay' have to be either 'counts' or 'logcounts'.")
  }
  if (!assay %in% names(assays(sce))){
    stop("Chosen assay option does not exist in counts matrix.")
  }
  counts = data.frame(cell = sce$cell ,
                      celltype = sce$celltype ,
                      counts = as.numeric( assay(sce[gene, ], assay)) )
  p <- ggplot(data=counts , aes(x = celltype , y = counts , fill = celltype)) +
    geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}

plot_expression_heatmap = function(sce , genes , assay = "logcounts" , title = length(genes)){
  if (!assay %in% names(assays(sce))){
    stop("Chosen assay option does not exist in counts matrix.")
  }
  sce = sce[rownames(sce) %in% genes , ]
  stat = lapply(unique(sce$celltype) , function(celltype){
    current.sce = sce[, sce$celltype == celltype]
    current.counts = as.matrix( assay(current.sce, assay))
    current.stat = data.frame(gene = rownames(sce) , mean.counts = apply(current.counts , 1 , mean))
    current.stat$celltype = celltype
    return(current.stat)
  })
  stat = do.call(rbind , stat)

  p <- ggplot(data=stat , aes(x = celltype , y = gene , fill = mean.counts)) +
    geom_tile() +
    scale_fill_viridis(discrete = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}


plot_mapping_heatmap = function(mapping , title = NULL){
  if (!.valid_mapping_df(mapping)) {
    stop()
  } else {
    if (!is.null(title) & !is(title, "character")){
      stop("Option 'title' should be either NULL or a string.")
    } else {
      mapping$celltype = as.character(mapping$celltype)
      mapping$celltype_mapped = as.character(mapping$celltype_mapped)
      tab = table(mapping$celltype , mapping$celltype_mapped)
      tab = sweep(tab, 1, rowSums(tab), "/")
      tab = as.data.frame( tab )
      colnames(tab) = c("celltype", "celltype_mapped", "n")
      tab$celltype = factor(tab$celltype , levels = unique(mapping$celltype))
      tab$celltype_mapped = factor(tab$celltype_mapped , levels = c(unique(mapping$celltype)))
      tab = tab[!is.na(tab$celltype) , ]
      p <- ggplot(tab, aes(x = celltype , y = celltype_mapped, fill = n)) +
        geom_tile() + viridis::scale_fill_viridis(discrete = F) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        ggtitle(title)
      return(p)
    }
  }
}




plot_mapping_sensitivity_trend = function(mappings , title = NULL){
  if (!is(mappings , "list")){
    stop("Input should be a list")
  } else {
    stat = lapply(1:length(mappings), function(i){
      mapping = mappings[[i]]
      if (!.valid_mapping_df(mapping)) {
        stop()
      } else {
        current_stat = get_sensitivity_mapping(mapping)
        if (!is.null(names(mappings))){
          current_stat$id = names(mappings)[i]
        } else {
          current_stat$id = i
        }
        return(current_stat)
      }
    })
    stat = do.call(rbind, stat)
    if (!is.null(names(mappings))){
      stat$id = factor(stat$id , levels = names(mappings))
    }
    pal = wes_palette("Zissou1" , length(mappings) , type = "continuous")
    p = ggplot(stat, aes( x = id, y = frac_mapped_correctly , col = id)) +
      geom_point(size=1) +
      facet_wrap(~celltype) +
      scale_color_manual(values = pal) +
      geom_hline(yintercept = 1, linetype = "dashed") +
      theme_classic() + theme( axis.text.x = element_blank() ) +
      ggtitle(title)
    return(p)
  }
}




get_umaps = function(sce , genes){
  plots = lapply(genes , function(gene){
    counts = data.frame(cell = colnames(sce) , counts = as.numeric(logcounts(sce)[gene , ]))
    current.meta = merge(meta , counts)
    current.meta = current.meta[order(current.meta$counts) , ]
    p <- ggplot(current.meta , aes(x = x , y = y , col = counts)) +
      geom_point() +
      scale_color_gradient(low = "azure3" , high = "darkgreen") + 
      theme(legend.position="none") +
      ggtitle(gene)
    return(p)
  })
  p <- ggarrange(plotlist = plots)
  return(p)
}


load_embryo_8.5 = function(dir = "local"){
  
  require(SingleCellExperiment)
  require(scater)
  
  if (dir == "cluster") {
    data.dir = "/nfs/research1/marioni/alsu/spatial/mouse_embryo/data/8_5/source/"
  } else {
    data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/"
  }
  
  sce = readRDS( paste0( data.dir , "E8.5_sce_filt_unlabelled.Rds"))
  # add normalization by libsize 
  assay(sce, "cpm") <- logcounts(scater::logNormCounts(sce, size_factors = sce$total))
  assay(sce, "cpm_wo_xist") <- logcounts(scater::logNormCounts(sce, size_factors = as.numeric( sce$total - counts( sce["Xist"] )) ))
  
  
  meta = colData(sce)
  meta = data.frame(meta)
  
  # rename Cavin3 --> Prkcdbp
  rownames.sce = rownames(sce)
  rownames.sce[rownames.sce == "Cavin3"] = "Prkcdbp"
  rownames(sce) = rownames.sce
  
  sce = sce[, sce$embryo == "embryo1" & sce$z == 2]
  return(sce)
  #assign("sce", sce, envir = .GlobalEnv)
  #assign("meta", meta, envir = .GlobalEnv)
  #invisible(0)
}

getSegmentationVerticesDF = function(DF,
                                     xname = "segmentation_vertices_x_global",
                                     yname = "segmentation_vertices_y_global",
                                     othercols = c("uniqueID","z")) {
  # DF is a DataFrame object
  # othercols is the others to keep
  
  long_x = unlist(DF[,xname])
  long_y = unlist(DF[,yname])
  
  if (length(long_x) != length(long_y)) stop("x and y need to be same length")
  
  long_xy = data.frame(
    long_x,
    long_y
  )
  colnames(long_xy) <- c(xname, yname)
  
  long_DF = cbind(
    rep(DF[,othercols], times = unlist(lapply(DF[,xname], length))),
    long_xy
  )
  
  return(as.data.frame(long_DF))
}

