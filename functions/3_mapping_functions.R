# This script contains all functions needed for mapping (and plotting functions)


get_mapping = function(genes , sce){
  set.seed(32)
  require(BiocNeighbors)
  require(scran)
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



