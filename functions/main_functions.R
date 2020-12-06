# This script contains all functions needed for gene ranking project


get_rid_of_highly_expressed_genes = function(sce , method = "absolute" , mean.thresh = 5, max.thresh = 500){
  counts = counts(sce)
  stat = data.frame(gene = rownames(counts) , 
                    mean.counts = apply(counts , 1 , mean) , 
                    max.counts = apply(counts , 1 , max) , 
                    sd.counts = apply(counts , 1 , sd))
  if (method == "absolute"){
    genes.2keep = stat$gene[stat$mean.counts < mean.thresh & stat$max.counts < max.thresh]
    
    sce.filtered = sce[rownames(sce) %in% genes.2keep , ]
    genes.2discard = setdiff(rownames(sce) , rownames(sce.filtered))
    print(paste0( "Discarded ", length(genes.2discard) , " genes as too highly expressed.") )
    return(sce.filtered)
  }
}

detect_rare_celltypes = function(sce , nCells.thresh = 10 , FDR.thresh = 0.01){
  stat.celltype = table(sce$celltype)
  rare_celltypes = stat.celltype[stat.celltype < nCells.thresh]
  
  markers.rare_celltypes = lapply(names( rare_celltypes ) , function(celltype){
    current.sce = sce
    current.sce$celltype.bool = current.sce$celltype == celltype
    if (stat.celltype[celltype] > 1){
      markers <- scran::findMarkers(current.sce , groups=current.sce$celltype.bool, direction = "up", test = "binom", assay.type = "logcounts")
      markers = as.data.frame( markers$`TRUE` )
      markers = markers[markers$FDR < FDR.thresh , ]
      if (nrow(markers) > 0){
        markers = markers[order(markers$logFC.FALSE , decreasing = T) ,]
        markers$celltype = celltype
        return(markers)
      }
    }
  })
  markers.rare_celltypes = do.call(rbind , markers.rare_celltypes)
  out = list(celltypes = rare_celltypes , markers = markers.rare_celltypes)
  return(out)
}

get_rid_of_rare_celltypes = function(sce , nCells.thresh = 10 , FDR.thresh = 0.01){
  rare_celltypes.stat = detect_rare_celltypes(sce , nCells.thresh , FDR.thresh)
  sce.filtered = sce[, !sce$celltype %in% names( rare_celltypes.stat$celltypes )]
  return(sce.filtered)
}

getDistrPlot = function(sce , assay.type , gene, title = NULL){
  if (assay.type == "counts"){
    current.counts = as.matrix(counts(sce[gene , ]))
  }
  else if (assay.type == "logcounts"){
    current.counts = as.matrix(logcounts(sce[gene , ]))
  }
  current.counts = data.frame(cell = sce$cell , celltype = sce$celltype , counts = as.numeric( current.counts) )
  p <- ggplot(current.counts , aes(x = celltype , y = counts , fill = celltype)) + 
    geom_violin() + geom_boxplot() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    theme(legend.position = "none") + 
    ggtitle(paste0( gene , title))
  return(p)
}


getMarkers = function(sce , test = "binom", FDR.thresh = 0.01){
  # requires there is a column called celltype and assay called logcounts
  require(scran)
  # get potential relevant genes
  markers <- scran::findMarkers(sce , groups=sce$celltype, direction = "up", pval.type="some", test = test, assay.type = "logcounts")
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



get_mapping = function( genes , sce ){
  require(BiocNeighbors)
  current.sce = sce[rownames(sce) %in% genes , ]
  meta = as.data.frame(colData(current.sce))
  batchFactor = factor(current.sce$sample )
  counts = as.matrix ( cosineNorm( logcounts(current.sce) ) )
  mbpca = multiBatchPCA(counts, batch = batchFactor, d = length(unique(sce$celltype)) )
  out = do.call(reducedMNN, mbpca)
  joint.pca = out$corrected
  
  mapping = lapply(unique(sce$sample) , function(sample){
    reference.cells = colnames(sce[,!sce$sample == sample])
    query.cells = colnames(sce[, sce$sample == sample])
    knns = queryKNN( joint.pca[reference.cells ,], joint.pca[query.cells ,], k = 5, get.index = TRUE, get.distance = FALSE)
    cells.mapped = t( apply(knns$index, 1, function(x) reference.cells[x]) )
    cells.mapped = as.data.frame( cells.mapped) 
    for (i in 1:ncol(cells.mapped)){
      cells.mapped[, i] = as.character(cells.mapped[, i])
    }
    celltypes = t(apply(cells.mapped, 1, function(x) meta$celltype[match(x, meta$cell)]))
    celltype.mapped = apply(celltypes, 1, function(x) getmode(x, 1:length(x)))
    current.mapping = data.frame(cell = query.cells , celltype.mapped = celltype.mapped , n.genes = length(genes))
    current.mapping$celltype.mapped = as.character(current.mapping$celltype.mapped)
    current.mapping = merge(current.mapping , meta[, c("cell" , "celltype")] , all.x = T , all.y = F)
    current.mapping = cbind(current.mapping , cells.mapped )
    return(current.mapping)
  })
  mapping = do.call(rbind , mapping)
  return(mapping)
}

get_loadings = function(sce , genes , type){
  set.seed(32)
  require(batchelor)
  require(stats)
  current.sce = sce[rownames(sce) %in% genes , ]
  counts = as.matrix(logcounts(current.sce))
  if (type == "celltype"){
    mbpca = multiBatchPCA(counts , batch = factor( current.sce$celltype ) , get.variance = T , d = nrow(counts))
    loadings = metadata(mbpca)  
    mbpca = do.call(rbind , mbpca)
    out = list(pc = mbpca , loadings = loadings$rotation , var.explained = loadings$var.explained)
  } 
  else if (type == "none") {
    pca = prcomp(t(counts) , rank. = nrow(counts))
    out = list(pc = pca$x , loadings = pca$rotation , var.explained = pca$sdev)
  }
  return(out)
}

get_stat_one_pc = function(loadings , pcs , nPC){
  pnorm.left = pnorm(loadings , mean = median(loadings), sd = mad(loadings), lower.tail = T)
  pnorm.left = p.adjust(pnorm.left, method = "fdr")
  
  pnorm.right = pnorm(loadings , mean = median(loadings), sd = mad(loadings), lower.tail = F)
  pnorm.right = p.adjust(pnorm.right, method = "fdr")
  
  out.max = data.frame(nPC = nPC , sign = "pos" , gene = names(which.max(loadings)) , loading = max(loadings) , 
                       FDR = pnorm.right[names(pnorm.right) == names(which.max(loadings))])
  out.min = data.frame(nPC = nPC , sign = "neg" , gene = names(which.min(loadings)) , loading = min(loadings), 
                       FDR =  pnorm.left[names(pnorm.left) == names(which.min(pnorm.left))])
  out = rbind(out.max , out.min)
  
  # add anova
  current.anova = aov(pcs ~ pcs$sample + pcs$celltype)
  current.summary = summary(current.anova)[[1]]
  out = data.frame(pc = current.component , sum_sq.sample = current.summary$`Sum Sq`[1] , sum_sq.celltype = current.summary$`Sum Sq`[2] , 
                   p.sample = current.summary$`Pr(>F)`[1] , p.celltype = current.summary$`Pr(>F)`[2])
  
  return(out)
}


get_fraction_mapped_correctly = function(mapping){
  tab = table(mapping$celltype , mapping$celltype.mapped) 
  tab = sweep(tab, 1, rowSums(tab), "/" )
  tab = as.data.frame(tab)
  colnames(tab) = c("celltype" , "celltype.mapped" , "frac")
  
  celltypes = as.character(unique(tab$celltype))
  stat = lapply(celltypes, function(celltype){
    current.tab = tab[tab$celltype == celltype , ]
    if (celltype %in% current.tab$celltype.mapped){
      out = data.frame(celltype = celltype , frac_correctly_mapped = current.tab$frac[current.tab$celltype.mapped == celltype])
    }
    else {
      out = data.frame(celltype = celltype , frac_correctly_mapped = 0)
    }
    return(out)
  })
  stat = do.call(rbind, stat)
  return(stat)
}


get_graph = function(genes , sce){
  require(igraph)
  current.sce = sce[rownames(sce) %in% genes , ]
  current.counts = cosineNorm(logcounts(current.sce))
  batchFactor = factor(c(as.character(current.sce$sample)))
  mbpca = multiBatchPCA(current.counts, batch = batchFactor, d = 50)
  mbpca = do.call(rbind , mbpca)
  current.graph = buildSNNGraph( mbpca , k = 10, transposed = T)
  current.graph <- set.vertex.attribute(current.graph, "name", value = rownames(mbpca))
  return(current.graph)
}

getCorrTranscriptome = function(mapping , genes){
  cols.cells = paste0("V" , c(1:5))
  current.sce = sce[rownames(sce) %in% genes , ]
  logcounts = as.matrix(logcounts(current.sce))
  logcounts.estimated = lapply(1:nrow(mapping), function(i){
    cells = mapping[i , cols.cells]
    current.logcounts = logcounts[, colnames(logcounts) %in% cells]
    return(apply(current.logcounts , 1 , mean))
  })
  logcounts.estimated = do.call(cbind , logcounts.estimated)
  colnames(logcounts.estimated) = mapping$cell
  logcounts.real = logcounts[, colnames(logcounts) %in% mapping$cell]
  logcounts.estimated = logcounts.estimated[, order(colnames(logcounts.estimated))]
  logcounts.real = logcounts.real[, order(colnames(logcounts.real))]
  stat = lapply(rownames(logcounts.real), function(gene){
    out = data.frame(gene = gene , corr = cor(logcounts.real[gene,] , logcounts.estimated[gene,] , method = "pearson"))
  })
  stat = do.call(rbind, stat)
  stat$n.genes = mapping$n.genes[1]
  return(stat)
}
