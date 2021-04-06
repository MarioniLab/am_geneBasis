# This script contains all functions needed for gene ranking project


create_sce = function(counts , meta , batch = NULL){
  counts = counts[, order(colnames(counts))]
  meta = meta[order(meta$cell) , ]
  if (is.null(batch)){
    sce = SingleCellExperiment(assays = list("counts" = counts), colData = meta)
    clusters <- quickCluster(sce, method="igraph", use.ranks=TRUE, d=50, min.mean=0.1)
    sce <- computeSumFactors(sce, clusters=clusters)
    sce <- logNormCounts(sce)
    return(sce)
  }
  else {
    if (!(batch %in% colnames(meta))){
      stop("batch specified variable should exist in SingleCellExperiment object - are you sure there is no typo?")
    }
    else {
      batchFactor = factor(meta[, colnames(meta) == batch])
      sce = lapply(unique(batchFactor) , function(current.batch){
        idx = which(batchFactor == current.batch)
        current.sce = SingleCellExperiment(assays = list("counts" = counts[,idx]), colData = meta[idx,])
        clusters <- quickCluster(current.sce, method="igraph", use.ranks=TRUE, d=50, min.mean=0.1)
        current.sce <- computeSumFactors(current.sce, clusters=clusters)
        return(current.sce)
      })
      sce = do.call(multiBatchNorm , sce )
      sce = do.call(cbind, sce)
      return(sce)
    }
  }
}


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
    sce = do.call(cbind , sce)
  }
  else {
    sce = SingleCellExperiment(assays = list("counts" = counts), colData = meta)
    reducedDims(sce) = list( UMAP = umaps )
    clusters = quickCluster(sce, method="igraph", use.ranks=TRUE, d=50, min.mean=0.1)
    sce = computeSumFactors(sce, clusters=clusters)
    sce = logNormCounts(sce)
  }
  return(sce)
}

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

detect_rare_celltypes = function(sce , nCells.thresh = 10 , FDR.thresh = 0.05){
  stat.celltype = table(sce$celltype)
  rare_celltypes = stat.celltype[stat.celltype < nCells.thresh]
  rare_celltypes = as.data.frame(rare_celltypes)

  if (nrow(rare_celltypes) > 0){
    colnames(rare_celltypes) = c("celltype" , "n.cells")
    markers.rare_celltypes = lapply(as.character( rare_celltypes$celltype ), function(celltype){
      current.sce = sce
      current.sce$celltype.bool = current.sce$celltype == celltype
      if (stat.celltype[celltype] > 1){
        markers <- scran::findMarkers(current.sce , groups=current.sce$celltype.bool, direction = "up", test = "t", assay.type = "logcounts")
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
  }
  else {
    out = list(celltypes = NULL , markers = NULL)
  }
  return(out)
}

get_rid_of_rare_celltypes = function(sce , nCells.thresh = 10 , FDR.thresh = 0.01){
  rare_celltypes.stat = detect_rare_celltypes(sce , nCells.thresh , FDR.thresh)
  if (!is.null(rare_celltypes.stat$celltypes)){
    sce.filtered = sce[, !sce$celltype %in% unique( rare_celltypes.stat$celltypes$celltype )]
  }
  else {
    sce.filtered = sce
  }
  return(sce.filtered)
}



get_markers = function(sce , test = "binom", FDR.thresh = 0.01){
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


retain_only_hvgs = function(sce, n = NULL, var.thresh = 0){
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


get_hvgs = function(sce, n = NULL, var.thresh = 0){
  dec.sce <- modelGeneVar(sce)
  if (is.null(n)){
    hvg.genes = getTopHVGs(dec.sce, var.threshold = var.thresh)
  } else {
    hvg.genes = getTopHVGs(dec.sce, n = n)
  }
  return(hvg.genes)
}



get_mapping_single_batch = function(sce , assay = "logcounts" , genes = rownames(sce), n.neigh = 3, nPC = 50 , get.dist = F){
  require(irlba)
  require(BiocNeighbors)
  current.sce = sce[genes , ]
  counts = as.matrix( assay(current.sce , assay) )
  meta = as.data.frame(colData(current.sce))

  if (!is.null(nPC)){
    pcs = prcomp_irlba(t(counts) , n = min(nPC, (nrow(counts)-1) ))
    counts = pcs$x
    rownames(counts) = colnames(current.sce)
  }
  else {
    counts = t(counts)
  }
  
  reference_cells = colnames(current.sce)
  query_cells = colnames(current.sce)
  
  if (get.dist){
    knns = queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.index = TRUE, get.distance = T)
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    distances = knns$distance[, 2:(n.neigh+1)]
    rownames(distances) = query_cells
    out = list(cells_mapped = cells_mapped , distances = distances)
  } else {
    knns = queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.index = TRUE )
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    out = cells_mapped
  }
  return(out)
}


get_mapping = function(sce , assay = "logcounts" , genes = rownames(sce), batch = "sample", n.neigh = 3, nPC = 50 , get.dist = F){
  if (is.null(batch)){
    out = get_mapping_single_batch(sce , assay = assay , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      current.neighs = get_mapping_single_batch(current.sce , assay = assay , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist)
      return(current.neighs)
    })
    if (!get.dist){
      cells_mapped = do.call(rbind , neighs)
      cells_mapped = cells_mapped[ order(match(rownames(cells_mapped), colnames(sce))), ]
      out = cells_mapped
    }
    else {
      distances = lapply(neighs , function(current.neighs) return(current.neighs$distances))
      distances = do.call(rbind , distances)
      distances = distances[ order(match(rownames(distances), colnames(sce))), ]
      
      cells_mapped = lapply(neighs , function(current.neighs) return(current.neighs$cells_mapped))
      cells_mapped = do.call(rbind , neighs$cells_mapped)
      cells_mapped = cells_mapped[ order(match(rownames(cells_mapped), colnames(sce))), ]
      
      out = list(cells_mapped = cells_mapped , distances = distances)
    }
  }
  return(out)
}



denoise_logcounts = function(sce, batch = "sample", n.neigh = 3, nPC = 50){
  neighs = get_mapping(sce , assay = "logcounts", genes = rownames(sce), batch = batch, n.neigh = n.neigh, nPC = nPC)
  logcounts_real = as.matrix(logcounts(sce))
  logcounts_denoised = lapply(1:nrow(neighs) , function(i){
    cells = as.character(neighs[i ,])
    current.counts = logcounts_real[, cells]
    current.counts = apply(current.counts , 1 , median)
    return(current.counts)
  })
  logcounts_denoised = do.call(cbind , logcounts_denoised)
  colnames(logcounts_denoised) = colnames(sce)
  return(logcounts_denoised)
}


get_corr_transcriptome_per_gene = function(sce , genes , assay = "logcounts" , batch = "sample" , 
                                           n.neigh = 10 , nPC = 50 , genes.predict = rownames(sce) ,
                                           method = "pearson"){
  neighs = get_mapping(sce , assay = "logcounts" , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
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


get_distr_dist = function(sce , genes , assay = "logcounts" , batch = "sample" , 
                                           n.neigh = 10 , nPC = 50 , genes.predict = rownames(sce) , type ){
  eps = 0.00001
  neighs = get_mapping(sce , assay = "logcounts" , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
  counts_predict = as.matrix(assay(sce[genes.predict , ] , assay))
  
  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]
  
  stat_real = stat_real
  stat_predict = stat_predict
  
  stat = lapply(1:nrow(counts_predict) , function(i){
    if (type == "KS"){
      out = data.frame(gene = rownames(counts_predict)[i] , ks.dist = ks.dist(stat_real[i,] , stat_predict[i,]))
      out = out[ , c("gene" , "ks.dist.D")]
      colnames(out) = c("gene" , "ks")
      out$ks[is.na(out$ks)] = eps
      out$ks[out$ks < eps] = eps
    }
    else if (type == "JS"){
      vec = rbind(stat_real[i,] , stat_predict[i,])
      out = data.frame(gene = rownames(counts_predict)[i] , JS = suppressMessages(JSD(vec)))
      out$JS[out$JS < eps] = eps
    }
    else if (type == "Wasserstein"){
      out = data.frame(gene = rownames(counts_predict)[i] , wasserstein = wasserstein1d(stat_real[i,] , stat_predict[i,]))
      out$wasserstein[out$wasserstein < eps] = eps
    }  
    else if (type == "corr"){
      out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = "spearman"))
      out$corr[is.na(out$corr)] = eps
      out$corr[out$corr < eps] = eps
    }
    else if (type == "dist"){
      out = data.frame(gene = rownames(counts_predict)[i] , dist = as.numeric(dist(rbind(stat_real[i,] , stat_predict[i,]) , method = "manhattan")))
      out$dist[out$dist < eps] = eps
    }
    return(out)
  }) 
  stat = do.call(rbind , stat)
  return(stat)
}


get_cells_not_mapped_well = function(sce , stat_predict.all , genes.selection , assay = "logcounts" , batch = "sample" , 
                                     n.neigh = 10 , nPC.selection = 50,  genes.predict = rownames(sce) ){
  neighs.selection = get_mapping(sce , assay = "logcounts" , genes = genes.selection, batch = batch , n.neigh = n.neigh , nPC = nPC.selection)
  counts_predict = as.matrix(assay(sce[genes.predict , ] , assay))
  
  stat_predict.selection = lapply(1:ncol(neighs.selection) , function(j){
    cells = neighs.selection[,j]
    current.stat_predict = counts_predict[, cells]
    current.stat_predict = current.stat_predict
    return(current.stat_predict)
  })
  
  stat = lapply(1:2000 , function(i){
    stat.per_gene = sapply(1:ncol(counts_predict) , function(j) {
      vec_1 = sapply(1:length(stat_predict.all) , function(k) return(stat_predict.all[[k]][i,j]))
      vec_2 = sapply(1:length(stat_predict.selection) , function(k) return(stat_predict.selection[[k]][i,j]))
      my.t = t.test(vec_1 , vec_2)
      p = my.t$p.value
      if (is.na(p) | p > 0.05){
        return(0)
      }
      else {
        return(1)
      }
    })
    out = data.frame(gene = rownames(counts_predict)[i] , n_cells_mismapped = sum(stat.per_gene))
    return(out)
  })
  stat = do.call(rbind , stat)
  return(stat)
  
}





get_moran_score = function(sce , genes , nPC = 50 , assay = "logcounts" , genes.predict = rownames(sce) ){
  
  snn.graph = as.matrix( as_adjacency_matrix( buildSNNGraph(assay(sce[genes , ], "logcounts") , d = nPC)) ) 
  snn.graph = snn.graph*0.9999
  
  counts = as.matrix(assay(sce[genes.predict , ] , assay))
  moran.score = aquila::MoransI(values = counts , weights = snn.graph, temp_dir = "~/Develop/test")
  colnames(moran.score) = c("gene" , "moran")
  return(moran.score)
  
}


get_preservation_score = function(sce , assay = "logcounts" , genes.all = rownames(sce) , 
                                  genes.compare , batch = "sample", n.neigh = 5 , nPC.all = 50 , nPC.compare = 200 , 
                                  bandwidth.type = "per_cell"){
  neighs.compare = get_mapping(sce , assay = assay , genes = genes.compare, batch = batch, n.neigh = n.neigh, nPC = nPC.compare)
  neighs.all = get_mapping(sce , assay = assay, genes = genes.all, batch = batch, n.neigh = (ncol(sce)-1), nPC = nPC.all , get.dist = T)
  
  # add random neighbors from the celltype
  neighs.rand = lapply(colnames(sce) , function(cell){
    cells = colnames(sce)[sce$celltype == sce$celltype[colnames(sce) == cell]]
    cells = setdiff(cells , cell)
    return(sample(cells , n.neigh))
  })
  neighs.rand = do.call(rbind , neighs.rand)
  
  
  neighs.all.cells_mapped = neighs.all$cells_mapped
  neighs.all.distances = neighs.all$distances
  n.cells = ncol(sce)
  
  score = lapply(1:nrow(neighs.compare) , function(i){
    cells = neighs.all.cells_mapped[i,]
    idx_all = c(1:n.neigh)
    idx_compare = which(cells %in% neighs.compare[i,] )
    idx_rand = which(cells %in% neighs.rand[i,] )
    
    dist_all = neighs.all.distances[i, idx_all]
    dist_compare = neighs.all.distances[i, idx_compare]
    dist_rand = neighs.all.distances[i, idx_rand]
    
    if (bandwidth.type == "per_cell") {
      current.bandwidth = max(dist_all)
      dist_all = exp(-1*dist_all/current.bandwidth)
      dist_compare = exp(-1*dist_compare/current.bandwidth)
      dist_rand = exp(-1*dist_rand/current.bandwidth)
      current.score = mean(dist_compare)/mean(dist_all)
      current.score.rand = mean(dist_rand)/mean(dist_all)
    } else if (bandwidth.type == "all") {
      current.bandwidth = mean(neighs.all.distances[i, ])
      dist_all = exp(-1*dist_all/current.bandwidth)
      dist_compare = exp(-1*dist_compare/current.bandwidth)
      dist_rand = exp(-1*dist_rand/current.bandwidth)
      current.score = mean(dist_compare)/mean(dist_all)
      current.score.rand = mean(dist_rand)/mean(dist_all)
    } else if (bandwidth.type == "quantile_mean"){
      current.score = mean(idx_compare/n.cells)
      current.score.rand =  mean(idx_rand/n.cells)
    } else if (bandwidth.type == "quantile_max"){
      current.score = max(idx_compare/n.cells)
      current.score.rand =  max(idx_rand/n.cells)
    }
    out = data.frame(score = current.score , score_rand = current.score.rand)
    return(out)
  })
  score = do.call(rbind , score)
  score$cell = rownames(neighs.compare)
  return(score)
}


get_preservation_score_simple = function(sce , assay = "logcounts" , genes.all = rownames(sce) , 
                                  genes.compare , batch = "sample", n.neigh = 5 , nPC.all = 50 , nPC.compare = 200 ){
  neighs.compare = get_mapping(sce , assay = assay , genes = genes.compare, batch = batch, n.neigh = n.neigh, nPC = nPC.compare)
  neighs.all = get_mapping(sce , assay = assay, genes = genes.all, batch = batch, n.neigh = (ncol(sce)-1), nPC = nPC.all , get.dist = T)
  
  neighs.all.cells_mapped = neighs.all$cells_mapped
  neighs.all.distances = neighs.all$distances
  n.cells = ncol(sce)
  
  score = lapply(1:nrow(neighs.compare) , function(i){
    cells = neighs.all.cells_mapped[i,]
    idx_all = c(1:n.neigh)
    idx_compare = which(cells %in% neighs.compare[i,] )
    
    dist_all = neighs.all.distances[i, idx_all]
    dist_compare = neighs.all.distances[i, idx_compare]
    
    current.bandwidth = mean(neighs.all.distances[i, ])
    dist_all = exp(-1*dist_all/current.bandwidth)
    dist_compare = exp(-1*dist_compare/current.bandwidth)
    current.score = mean(dist_compare)/mean(dist_all)
    out = data.frame(score = current.score )
    return(out)
  })
  score = do.call(rbind , score)
  score$cell = rownames(neighs.compare)
  return(score)
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


get_loadings_squared = function(sce, nPC = 50 , type ){
  require(irlba)
  pcs = prcomp_irlba(t(as.matrix(logcounts(sce))) , n = nPC)
  loadings = pcs$rotation
  rownames(loadings) = rownames(sce)
  if (type == "sum"){
    loadings.squared = data.frame( gene = rownames(loadings) , loading.squared = apply(loadings , 1 , function(x) sum(x^2)))
  } 
  else if (type == "max") {
    loadings.squared = data.frame( gene = rownames(loadings) , loading.squared = apply(loadings , 1 , function(x) max(x^2)))
  }
  return(loadings.squared)
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
