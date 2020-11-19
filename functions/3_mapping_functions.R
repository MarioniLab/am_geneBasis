# This script contains all functions needed for mapping (and plotting functions)

getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

get_mapping = function(genes , sce){
  set.seed(32)
  require(BiocNeighbors)
  require(scran)
  require(batchelor)
  require(SingleCellExperiment)
  sce = sce[rownames(sce) %in% genes , ]
  meta = as.data.frame(colData(sce))
  batchFactor = factor(sce$sample )
  counts = as.matrix ( cosineNorm( logcounts(sce)) )
  
  mbpca = multiBatchPCA(counts, batch = batchFactor, d = 50)
  out = do.call(reducedMNN, mbpca)
  joint.pca = out$corrected
  
  mapping = lapply(unique(sce$sample) , function(sample){
    reference.cells = colnames(sce[,!sce$sample == sample])
    query.cells = colnames(sce[, sce$sample == sample])
    knns = queryKNN( joint.pca[reference.cells ,], joint.pca[query.cells ,], k = 10, get.index = TRUE, get.distance = FALSE)
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

get_fraction_mapped_correctly = function(mapping){
  tab = table(mapping$celltype , mapping$celltype.mapped) 
  tab = sweep(tab, 1, rowSums(tab), "/" )
  stat = lapply(rownames(stat) , function(celltype){
    out = data.frame(celltype = celltype , frac_correctly_mapped = tab[celltype , celltype])
  })
  stat = do.call(rbin, stat)
  return(stat)
}

get_mapping_robustness = function(genes, sce){
  mapping_all = get_mapping(genes , sce)
  stat_all = get_fraction_mapped_correctly(mapping_all)
  colnames(stat_all) = c("celltype" , "frac_correctly_mapped.all" )
  stat = lapply(genes , function(gene){
    genes.2keep = setdiff(genes, gene)
    current.mapping = get_mapping(genes.2keep , sce)
    current.stat = get_fraction_mapped_correctly(current.mapping)
    current.stat = merge(current.stat , stat_all , by = "celltype")
    current.stat$frac_correctly_mapped.ratio = current.stat$frac_correctly_mapped / current.stat$frac_correctly_mapped.all
    current.stat = current.stat[, c("celltype" , "frac_correctly_mapped.ratio")]
    current.stat$gene = gene
    return(current.stat)
  })
  stat = do.call(rbind , stat)
  return(stat)
}

getCorrTranscriptome = function(mapping , genes){
  cols.cells = paste0("V" , c(1:10))
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
