# This script contains all functions needed for gene ranking project


# get_loadings = function(sce , nPC = nrow(sce) ){
#   set.seed(32)
#   require(batchelor)
#   counts = as.matrix(logcounts(sce))
#   pca = prcomp(t(counts) , rank. = nPC)
#   loadings = pca$rotation
#   return(loadings)
# }

get_loadings = function(sce , nPC = nrow(sce) ){
  set.seed(32)
  require(batchelor)
  counts = as.matrix(logcounts(sce))
  mbpca = multiBatchPCA(counts , batch = factor( sce$sample ) , d = nPC)
  loadings = metadata(mbpca)  
  loadings = loadings$rotation
  return(loadings)
}

getBrokenStickDistribution = function(n){
  out = sapply(1:n ,function(i){
    current.out = sapply(i : n , function(k){
      return(1/k)
    })
    return(sum(current.out)) 
  })
  out = out / n
  return(out)
}

get_stat_highest_loadings = function(loadings , nPC){
  pnorm.left = pnorm(loadings , mean = median(loadings), sd = mad(loadings), lower.tail = T)
  pnorm.left = p.adjust(pnorm.left, method = "fdr")
  
  pnorm.right = pnorm(loadings , mean = median(loadings), sd = mad(loadings), lower.tail = F)
  pnorm.right = p.adjust(pnorm.right, method = "fdr")
  
  out.max = data.frame(nPC = nPC , sign = "pos" , gene = names(which.max(loadings)) , loading = max(loadings) , 
                       FDR = pnorm.right[names(pnorm.right) == names(which.max(loadings))])
  out.min = data.frame(nPC = nPC , sign = "neg" , gene = names(which.min(loadings)) , loading = min(loadings), 
                       FDR =  pnorm.left[names(pnorm.left) == names(which.min(pnorm.left))])
  out = rbind(out.max , out.min)
  return(out)
}

rank_genes = function(loadings , FDR.thresh = 0.01){
  stat = lapply(1:ncol(loadings), function(i){
    current.loading = loadings[,i]
    return(get_stat_highest_loadings(current.loading , i)) 
  })
  stat = do.call(rbind , stat)
  stat = stat[stat$FDR < FDR.thresh , ]
  return(stat)
}

getBrokenStickDistribution = function(n){
  out = sapply(1:n ,function(i){
    current.out = sapply(i : n , function(k){
      return(1/k)
    })
    return(sum(current.out)) 
  })
  out = out / n
  return(out)
}

get_PCs = function(genes , sce , nPC){
  set.seed(32)
  require(batchelor)
  require(SingleCellExperiment)
  require(BiocNeighbors)
  require(scran)
  require(tibble)
  current.sce = sce[rownames(sce) %in% genes, ]
  current.counts = cosineNorm(logcounts(current.sce))
  batchFactor = factor(c(as.character(current.sce$sample)))
  mbpca = multiBatchPCA(current.counts, batch = batchFactor, d = nPC)
  mbpca = do.call(rbind , mbpca)
  #mbpca = as.data.frame(mbpca)
  #mbpca = rownames_to_column(mbpca, var = "cell")
  return(mbpca)
}
