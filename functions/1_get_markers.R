# This script contains all functions to get all potential markers for sce


filter_veryRareCelltypes = function(sce , nCells.thresh = 5){
  stat = table(sce$celltype)
  celltypes = names(stat)[stat >= nCells.thresh]
  sce_filtered = sce[, sce$celltype %in% celltypes]
  print(paste0( ncol(sce_filtered) , " cells are retained for the downstream analysis (out of ", ncol(sce) , ")") )
  return(sce_filtered)
}

getMarkers = function(sce , pval.type = "some", test = "t", FDR.thresh = 0.01){
  # requires there is a column called celltype and assay called logcounts
  require(scran)
  # get potential relevant genes
  markers <- findMarkers(sce , groups=sce$celltype, direction = "up", pval.type=pval.type, test = test, assay.type = "logcounts")
  # put together
  celltypes = names(markers)
  markers = lapply(1:length(celltypes), function(i){
    current.markers = as.data.frame(markers[[i]])
    current.markers = current.markers[!is.na(current.markers$FDR) & current.markers$FDR < FDR.thresh , ]
    if (nrow(current.markers) > 0){
      out = data.frame( celltype = celltypes[i], gene = rownames(current.markers))
      return(out)
    }
  })
  markers = do.call(rbind,markers)
  return(markers)
}

# Function to get rid of highly expressed genes
suitable4seq_genes = function(sce , mean.thresh = 7, max.thresh = 350){
  counts = as.matrix(counts(sce))
  stat = lapply(rownames(counts), function(gene){
    out = data.frame(gene = gene,
                     mean.counts = mean( counts[gene,] ) ,
                     max.counts = max( counts[gene,] ) )
    return(out)
  })
  stat = do.call(rbind, stat)
  suitable4seq_genes = stat$gene[stat$mean.counts < mean.thresh & stat$max.counts < max.thresh]
  return(suitable4seq_genes)
}

# combine all together
suitable4seq_markers = function(sce , nCells.thresh = 5 , mean.thresh = 7, max.thresh = 350, pval.type = "some", test = "t", FDR.thresh = 0.01){
  #sce = filter_veryRareCelltypes(sce , nCells.thresh)
  stat = getMarkers(sce , pval.type = "some", test , FDR.thresh = 0.01)
  markers = unique(stat$gene)
  markers = intersect(markers , suitable4seq_genes(sce , mean.thresh , max.thresh))
  print( paste0( length(markers) , " genes are selected as relevant for the downstream analysis (out of ", nrow(sce) , ")") )
  stat = stat[stat$gene %in% markers , ]
  out = list(stat = stat , markers = markers)
  return(out)
}


add_markers_posthoc = function(sce , genes , celltype.1 , celltype.2 = NULL ){
  if (!is.null(celltype.2)){
    current.sce = sce[, sce$celltype %in% union(celltype.1 , celltype.2)]
  } 
  else {
    current.sce = sce
  }
  current.markers = findMarkers(current.sce, current.sce$celltype == celltype.1 , assay.type = "logcounts", direction = "up" , pval.type = "all" )
  current.markers = as.data.frame(current.markers[[2]])
  current.markers = current.markers[current.markers$logFC.FALSE > 0 & !rownames(current.markers) %in% genes , ]
  current.markers = current.markers[order(current.markers$logFC.FALSE , decreasing = T),]
  return(current.markers)
}

  
  
  
  
  