# This script contains all functions needed for gene ranking project

get_loadings = function(sce , nPC = 100){
  set.seed(32)
  require(batchelor)
  counts = as.matrix(logcounts(sce))
  batch = factor( sce$celltype )
  mbpca = multiBatchPCA(counts , batch = batch , d = nPC)
  loadings = metadata(mbpca)  
  loadings = loadings$rotation
  return(loadings)
}

rank_genes = function(loadings){
  stat = lapply(1:ncol(loadings), function(i){
    current.loading = loadings[,i]
    out.max = data.frame(nPC = i , sign = "pos" , gene = names(which.max(current.loading)) , loading = max(current.loading))
    out.min = data.frame(nPC = i , sign = "neg" , gene = names(which.min(current.loading)) , loading = min(current.loading))
    out = rbind(out.max , out.min)
    return(out)
  })
  stat = do.call(rbind , stat)
  return(stat)
}

