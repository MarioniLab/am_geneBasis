# This script contains all functions to run search for genes that explain variation


add_logcounts = function(sce , batch){
  require(SingleCellExperiment)
  require(scuttle)
  require(batchelor)
  meta = as.data.frame(colData(sce))
  counts = counts(sce)
  umaps = reducedDim(sce , "UMAP")
  
  # add logcounts
  if (!is.null(batch)){
    unq.batches = unique(meta[, batch])
    sce = bplapply(unq.batches , function(current.batch){
      idx = which(meta[ , batch ] == current.batch)
      current.sce = SingleCellExperiment(assays = list("counts" = counts[,idx]), colData = meta[idx,])
      reducedDims(current.sce) = list( UMAP = umaps[idx,] )
      clusters <- quickCluster(current.sce, method="igraph", use.ranks=TRUE, d=50, min.mean=0.1)
      current.sce <- computeSumFactors(current.sce, clusters=clusters)
      return(current.sce)
    } , BPPARAM = mcparam)
    sce <- do.call(multiBatchNorm , sce )
  }
  else {
    sce = SingleCellExperiment(assays = list("counts" = counts), colData = meta)
    reducedDims(current.sce) = list( UMAP = umaps )
    clusters = quickCluster(sce, method="igraph", use.ranks=TRUE, d=50, min.mean=0.1)
    sce = computeSumFactors(sce, clusters=clusters)
    sce = LogNormCounts(sce)
  }
  return(sce)
}



retain_only_hvgs = function(sce, n = NULL, var.thresh = 0){
  require(scran)
  dec.sce <- modelGeneVar(sce)
  if (is.null(n)){
    hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
  } else {
    hvg.genes = getTopHVGs(dec.sce, n = n)
  }
  print(paste0(length(hvg.genes) , " genes retained (out of " , nrow(sce) , ")"))
  sce = sce[hvg.genes,]
  return(sce)
}


get_mapping = function(sce , assay = "logcounts" , genes = rownames(sce), batch = NULL, n.neigh = 3, nPC = 50 , get.dist = F){
  require(batchelor)
  require(irlba)
  current.sce = sce[genes , ]
  counts = cosineNorm( as.matrix( assay(current.sce , assay) ) )
  meta = as.data.frame(colData(current.sce))
  if (!is.null(batch)) {
    batchFactor = factor(meta[, colnames(meta) == batch])
    pcs = suppressWarnings( multiBatchPCA(counts, batch = batchFactor, d = nPC ) )
    pcs = suppressWarnings( do.call(reducedMNN, pcs) )
    pcs = pcs$corrected
  } else {
    pcs = prcomp_irlba(t(counts) , n = min(nPC, (nrow(counts)-1) ))
    pcs = pcs$x
    rownames(pcs) = colnames(current.sce)
  }
  reference_cells = colnames(current.sce)
  query_cells = colnames(current.sce)
  if (get.dist){
    knns = queryKNN( pcs[reference_cells ,], pcs[query_cells ,], k = (n.neigh+1), get.index = TRUE, get.distance = T)
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    distances = knns$distance[, 2:(n.neigh+1)]
    rownames(distances) = query_cells
    out = list(cells_mapped = cells_mapped , distances = distances)
  } else {
    knns = queryKNN( pcs[reference_cells ,], pcs[query_cells ,], k = (n.neigh+1), get.index = TRUE)
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    out = cells_mapped
  }
  return(out)
}


denoise_logcounts = function(sce, batch = NULL, n.neigh = 3, nPC = 100){
  neighs = get_mapping(sce , assay = "logcounts", genes = rownames(sce), batch = batch, n.neigh = n.neigh, nPC = nPC)
  logcounts_real = as.matrix(logcounts(sce))
  logcounts_denoised = bplapply(1:nrow(neighs) , function(i){
    cells = as.character(neighs[i ,])
    current.counts = logcounts_real[, cells]
    current.counts = apply(current.counts , 1 , median)
    return(current.counts)
  } , BPPARAM = mcparam)
  logcounts_denoised = do.call(cbind , logcounts_denoised)
  return(logcounts_denoised)
}


get_corr_transcriptome_per_gene = function(sce , assay = "logcounts", genes , batch = NULL , 
                                           n.neigh = 10 , nPC = 50 , genes.predict = rownames(sce) ,
                                           method = "pearson"){
  neighs = get_mapping(sce , assay = assay , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
  counts_predict = as.matrix(assay(sce[genes.predict , ] , assay))
  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]
  corr = lapply(1:nrow(counts_predict) , function(i){
    out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = method))
    return(out)
  }) 
  corr = do.call(rbind , corr)
  return(corr)
}


get_loadings_squared = function(sce){
  require(irlba)
  pcs = prcomp_irlba(t(as.matrix(logcounts(sce))) , n = 50)
  loadings = pcs$rotation
  rownames(loadings) = rownames(sce)
  loadings.squared = data.frame( gene = rownames(loadings) , loading.squared = apply(loadings , 1 , function(x) sum(x^2)))
  return(loadings.squared)
}





get_mapping_2_external_dataset = function(sce_reference , sce_query , cluster_id , genes, nPC = 100, n.neigh = 5){
  genes = intersect(genes,rownames(sce_reference))
  genes = intersect(genes,rownames(sce_query))
  
  sce_reference = sce_reference[genes , ]
  sce_reference = sce_reference[order(rownames(sce_reference)), ]
  sce_query = sce_query[genes , ]
  sce_query = sce_query[order(rownames(sce_query)), ]
  
  assay(sce_reference , "cosineNorm") = cosineNorm(logcounts(sce_reference))
  assay(sce_query , "cosineNorm") = cosineNorm(logcounts(sce_query))
  
  sce_joint = cbind(assay(sce_query, "cosineNorm"), assay(sce_reference, "cosineNorm"))
  batchFactor = factor(c(as.character(sce_query$sample), as.character(sce_reference$Sample)))
  
  mbpca = multiBatchPCA(sce_joint, batch = batchFactor, d = nPC)
  out = do.call(reducedMNN, mbpca)
  joint_pca = out$corrected
  current.knns = queryKNN( joint_pca[colnames(sce_reference),], joint_pca[colnames(sce_query),], k = n.neigh, 
                           get.index = TRUE, get.distance = FALSE)
  cells.mapped = t( apply(current.knns$index, 1, function(x) colnames(sce_reference)[x]) )
  
  meta_reference = as.data.frame(colData(sce_reference))
  mapped_ct = apply(cells.mapped , 1 , function(x) return(getmode(meta_reference[match(x, rownames(meta_reference)) , cluster_id] , c(1:n.neigh)) ))
  meta_query = as.data.frame(colData(sce_query))
  meta_query$mapped_ct = mapped_ct
  return(meta_query)
}


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


get_wallace_dist_between_clusters = function(clusters_1 , clusters_2){
  df = as.data.frame( cbind(clusters_1 , clusters_2) )
  colnames(df) = c("clust_1" , "clust_2")
  
  cells = rownames(df)
  N_11 = sapply(1:(nrow(df) - 1) , function(i1){
    N_11.per_cell = sapply((i1+1):nrow(df) , function(i2){
      return(df$clust_1[i1] == df$clust_1[i2] & df$clust_2[i1] == df$clust_2[i2])
    })
    return(sum(N_11.per_cell))
  })
  N_11 = sum(N_11)
  
  tab.clust_1 = as.numeric( table(df$clust_1) )
  denominator.clust_1 = sapply(tab.clust_1 , function(x) return(x*(x-1)/2))
  denominator.clust_1 = sum(denominator.clust_1)
    
  tab.clust_2 = as.numeric( table(df$clust_2) )
  denominator.clust_2 = sapply(tab.clust_2 , function(x) return(x*(x-1)/2))
  denominator.clust_2 = sum(denominator.clust_2)
  
  out = data.frame(wallace_first_over_second = N_11/denominator.clust_1 , wallace_second_over_first = N_11/denominator.clust_2)
  return(out)
}


get_cluster_comparison_stats = function(clusters_1 , clusters_2){
  out = data.frame(vi = compare(clusters_1,clusters_2,"vi"),
             nmi = compare(clusters_1,clusters_2,"nmi"),
             split.join = compare(clusters_1,clusters_2,"split.join"),
             rand = compare(clusters_1,clusters_2,"rand"),
             adjusted.rand = compare(clusters_1,clusters_2,"adjusted.rand")
  )
  out = cbind(out , get_wallace_dist_between_clusters(clusters_1 , clusters_2))
  return(out)
}



# get function that estimates whether the gene is ct marker or no
get_ct_marker_genes = function(sce , genes){
  current.logcounts = as.matrix(logcounts(sce[genes,]))
  current.logcounts = apply(current.logcounts , 1 , function(x) quantile(x,0.5))
  current.logcounts = current.logcounts > 0
  return(names(current.logcounts)[current.logcounts == F])
}
