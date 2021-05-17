# This script contains all functions needed for gene ranking project

get_markers = function(sce , test = "binom", type = "all", FDR.thresh = 0.01){
  # requires there is a column called celltype and assay called logcounts
  require(scran)
  # get potential relevant genes
  markers = scran::findMarkers(sce , groups=sce$celltype, direction = "up", pval.type=type, test = test, assay.type = "logcounts")
  # put together
  celltypes = names(markers)
  markers = lapply(1:length(celltypes), function(i){
    current.markers = as.data.frame(markers[[i]])
    current.markers = current.markers[!is.na(current.markers$FDR) & current.markers$FDR < FDR.thresh , ]
    if (nrow(current.markers) > 0){
      out = data.frame( celltype = celltypes[i], gene = rownames(current.markers) )
      return(out)
    }
  })
  markers = do.call(rbind,markers)
  return(markers)
}

get_fraction_mapped_correctly = function(mapping, cluster.id = "celltype" , cluster_mapped.id = "mapped_celltype"){
  tab = table(mapping[, cluster.id] , mapping[, cluster_mapped.id])
  tab = sweep(tab, 1, rowSums(tab), "/" )
  tab = as.data.frame(tab)
  colnames(tab) = c("cluster" , "mapped_cluster" , "frac")

  clusters = as.character(unique(tab$cluster))
  stat = lapply(clusters, function(cluster){
    current.tab = tab[tab$cluster == cluster , ]
    if (cluster %in% current.tab$mapped_cluster){
      out = data.frame(cluster = cluster , frac_correctly_mapped = current.tab$frac[current.tab$mapped_cluster == cluster])
    }
    else {
      out = data.frame(cluster = cluster , frac_correctly_mapped = 0)
    }
    return(out)
  })
  stat = do.call(rbind, stat)
  colnames(stat) = c(cluster.id , "frac_correctly_mapped")
  return(stat)
}


get_hvgs = function(sce, n = NULL, var.thresh = 0){
  require(scran)
  dec.sce = scran::modelGeneVar(sce)
  hvg.genes = scran::getTopHVGs(dec.sce, var.threshold = var.thresh)
  if (!is.null(n) & n < length(hvg.genes)){
    hvg.genes = scran::getTopHVGs(dec.sce, n = n)
  }
  return(hvg.genes)
}

retain_only_hvgs = function(sce, n = NULL, var.thresh = 0){
  hvgs = get_hvgs(sce, n = n, var.thresh = var.thresh)
  print(paste0(length(hvgs) , " genes retained (out of " , nrow(sce) , ")"))
  sce = sce[hvgs,]
  return(sce)
}


assign_neighbors = function(counts , reference_cells , query_cells, n.neigh = 5, get.dist = F){
  require(BiocNeighbors)
  
  knns = BiocNeighbors::queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.distance = get.dist)
  cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
  rownames(cells_mapped) = query_cells
  if (!get.dist){
    out = cells_mapped
  }
  if (get.dist){
    distances = knns$distance[, 2:(n.neigh+1)]
    rownames(distances) = query_cells
    out = list(cells_mapped = cells_mapped , distances = distances)
  }
  return(out)
}


get_mapping_single_batch = function(sce , genes = rownames(sce), n.neigh = 5, nPC = 50 , get.dist = F , cosine = F){
  require(irlba)
  require(BiocNeighbors)
  set.seed(32)
  
  current.sce = sce[genes , ]
  if (cosine){
    logcounts(current.sce) = cosineNorm(logcounts(current.sce))
  }
  counts = as.matrix( logcounts(current.sce) )
  meta = as.data.frame(colData(current.sce))
  res = tryCatch(
    {
      if (!is.null(nPC)){
        pcs = irlba::prcomp_irlba(t(counts) , n = min(nPC, (nrow(counts)-1) , (ncol(counts) - 1)))
        counts = pcs$x
        rownames(counts) = colnames(current.sce)
      } else {
        counts = t(counts)
      }
      reference_cells = colnames(current.sce)
      query_cells = colnames(current.sce)
      
      if (n.neigh == "all"){
        n.neigh = length(reference_cells) - 1
      }
      out = assign_neighbors(counts , reference_cells , query_cells, n.neigh = n.neigh, get.dist = get.dist)
      return(out)
    },
    error = function(dump){
      message("Features you selected can not be used for pca")
      return(NULL)
    }
  )
  return(res)
}


get_mapping_many_batches = function(sce , genes = rownames(sce), batch = "sample", n.neigh = 5, nPC = 50 , get.dist = F, cosine = F){
  require(batchelor)
  require(BiocNeighbors)
  set.seed(32)
  
  current.sce = sce[genes , ]
  if (cosine){
    logcounts(current.sce) = cosineNorm(logcounts(current.sce))
  }
  counts = as.matrix( logcounts(current.sce) )
  meta = as.data.frame(colData(current.sce))
  batchFactor = factor(meta[, colnames(meta) == batch])
  res = tryCatch(
    {
      if (!is.null(nPC)){
        counts = batchelor::multiBatchPCA(counts , batch = batchFactor , d = nPC)
        counts = do.call(batchelor::reducedMNN , counts)
        counts = counts$corrected
      } else {
        counts = batchelor::fastMNN(counts , batch = batchFactor , d = NA)
        counts = reducedDim(counts , "corrected")
      }
      reference_cells = colnames(current.sce)
      query_cells = colnames(current.sce)
      
      if (n.neigh == "all"){
        n.neigh = length(reference_cells) - 1
      }
      out = assign_neighbors(counts , reference_cells , query_cells, n.neigh = n.neigh, get.dist = get.dist)
      return(out)
    },
    error = function(dump){
      message("Features you use are insufficient to perform batch correction via fastMNN")
      return(NULL)
    }
  )
  return(res)
}


get_mapping = function(sce , genes = rownames(sce), batch = "sample", n.neigh = 5, nPC = 50 , get.dist = F, type = "per batch", cosine = F){
  #require(BiocSingular)
  #require(BiocParallel)
  require(BiocNeighbors)
  
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)
  #set.seed(32)
  
  if (is.null(batch)){
    out = get_mapping_single_batch(sce , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
  }
  else if (type == "together") {
    out = get_mapping_many_batches(sce , genes = genes, batch = batch , n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
  }
  else if (type == "per batch") {
    
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    
    neighs = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.neighs = get_mapping_single_batch(sce[, idx] , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist, cosine = cosine)
      return(current.neighs)
    })
    if (!get.dist){
      cells_mapped = do.call(rbind , neighs)
      cells_mapped = cells_mapped[ match(colnames(sce), rownames(cells_mapped)), ]
      out = cells_mapped
    }
    else {
      distances = lapply(neighs , function(current.neighs) return(current.neighs$distances))
      distances = do.call(rbind , distances)
      distances = distances[ match(colnames(sce), rownames(distances)), ]
      
      cells_mapped = lapply(neighs , function(current.neighs) return(current.neighs$cells_mapped))
      cells_mapped = do.call(rbind , neighs$cells_mapped)
      cells_mapped = cells_mapped[ match(colnames(sce), rownames(cells_mapped)), ]
      
      out = list(cells_mapped = cells_mapped , distances = distances)
    }
  }
  return(out)
}


get_corr_per_gene = function(sce , genes , batch = "sample" , 
                                           n.neigh = 10 , nPC = 50 , genes.predict = rownames(sce) , method = "spearman"){
  #require(BiocSingular)
  #require(BiocParallel)
  
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)
  
  eps = 0.00001
  neighs = get_mapping(sce , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
  counts_predict = as.matrix(logcounts(sce[genes.predict , ]))
  
  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]
  
  
  stat = lapply(1:nrow(counts_predict) , function(i){
    out = data.frame(gene = rownames(counts_predict)[i] , corr = cor(stat_real[i,] , stat_predict[i,] , method = method))
    out$corr[is.na(out$corr)] = eps
    out$corr[out$corr < eps] = eps
    return(out)
  }) 
  stat = do.call(rbind , stat)
  return(stat)
}


get_lp_norm_dist = function(sce , genes , batch = "sample" , n.neigh = 5 , nPC = 50 , genes.predict = rownames(sce) , p){
  #require(BiocSingular)
  #require(BiocParallel)
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)

  if (!is.null(genes)){
    neighs = get_mapping(sce , genes = genes, batch = batch , n.neigh = n.neigh , nPC = nPC)
  }
  else {
    neighs = initiate_random_mapping(sce , batch = batch , n.neigh = n.neigh)
  }
  
  counts_predict = as.matrix(logcounts(sce[genes.predict , ]))

  stat_predict = lapply(1:ncol(neighs) , function(j){
    cells = neighs[,j]
    current.stat_predict = counts_predict[, cells]
    return(current.stat_predict)
  })
  stat_predict = Reduce("+", stat_predict) / length(stat_predict)
  stat_real = counts_predict[, rownames(neighs)]
  
  if (p < Inf){
    stat = lapply(1:nrow(counts_predict) , function(i){
      out = data.frame(gene = rownames(counts_predict)[i] , dist = as.numeric(dist(rbind(stat_real[i,] , stat_predict[i,]) , method = "minkowski" , p = p)))
      return(out)
    }) 
    stat = do.call(rbind , stat)
  } else if (is.infinite(p)){
    stat = lapply(1:nrow(counts_predict) , function(i){
      out = data.frame(gene = rownames(counts_predict)[i] , dist = max(abs(stat_real[i,] - stat_predict[i,])))
      return(out)
    }) 
    stat = do.call(rbind , stat)
  }
  return(stat)
}


initiate_random_mapping = function(sce , batch = "sample", n.neigh = 5){
  require(irlba)
  require(BiocNeighbors)
  #require(BiocSingular)
  #require(BiocParallel)
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)
  
  if (is.null(batch)){
    batchFactor = factor(rep(1 , ncol(sce)))
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
  }
  
  initial_random_mtrx = suppressWarnings( abs(matrix(rnorm(10),2,ncol(sce))) )
  colnames(initial_random_mtrx) = colnames(sce)
  
  neighs = lapply(unique(batchFactor) , function(current.batch){
    idx = which(batchFactor == current.batch)
    counts = t( initial_random_mtrx[, idx] )
    
    reference_cells = colnames(sce[,idx])
    query_cells = colnames(sce[,idx])
    
    knns = suppressWarnings( queryKNN( counts[reference_cells ,], counts[query_cells ,], k = (n.neigh+1), get.distance = F) )
    cells_mapped = t( apply(knns$index, 1, function(x) reference_cells[x[2:(n.neigh+1)]]) )
    rownames(cells_mapped) = query_cells
    return(cells_mapped)
  })
  neighs = do.call(rbind , neighs)
  neighs = neighs[ match(colnames(sce), rownames(neighs)), ]
  return(neighs)
}

add_gene_to_current_selection = function(sce , stat_all, genes = NULL , batch = NULL , n.neigh = 5, nPC = NULL, p = 3){
  stat_genes = suppressWarnings( get_lp_norm_dist(sce , genes = genes , batch = batch, n.neigh = n.neigh , nPC = nPC, 
                                                    genes.predict = rownames(sce) , p = p) )
  stat_genes = stat_genes[!stat_genes$gene %in% genes , ] 
  stat_genes = merge(stat_genes , stat_all)
  stat_genes$dist_diff = stat_genes$dist - stat_genes$dist_all 
  idx = which(stat_genes$dist_diff == max(stat_genes$dist_diff))
  gene = stat_genes$gene[idx[1]]
  return(gene)
}
  
add_first_gene = function(sce , stat_all, batch = NULL , n.neigh = 5, p = 3 , K = 5){
  #require(BiocSingular)
  #require(BiocParallel)
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)
  
  first_genes = lapply(1:K , function(dump){
    gene = add_gene_to_current_selection(sce , stat_all, genes = NULL , batch = batch , n.neigh = n.neigh, nPC = NULL, p = p)
  })
  first_genes = unlist(first_genes)
  out = as.character(names(sort(table(first_genes),decreasing=TRUE)[1]))
  return(out)
}
  

gene_search = function(sce , genes_base = NULL, n_genes_total , batch = NULL, n.neigh = 5, p = 3, K = 5, nPC = NULL){
  # get baseline stat-all
  stat_all = suppressWarnings( get_lp_norm_dist(sce, genes = rownames(sce), batch = batch , n.neigh = n.neigh , nPC = 50 , 
                                                genes.predict = rownames(sce) , p = p) )
  colnames(stat_all) = c("gene" , "dist_all")
  
  # add first gene if selection is empty
  if (is.null(genes_base)){
    genes_base = add_first_gene(sce , stat_all, batch = batch , n.neigh = n.neigh, p = p, K = K)
  }
  
  genes_all = genes_base
  while(length(genes_all) < n_genes_total){
    gene = add_gene_to_current_selection(sce , stat_all, genes = genes_all , batch = batch , n.neigh = n.neigh, nPC = nPC, p = p)
    #print(as.character(gene))
    genes_all = c(genes_all , as.character(gene))
  }
  out = data.frame(rank = c(1:length(genes_all)) , gene = genes_all)
  return(out)
}


get_expr_real_and_neighbors = function(sce , genes , assay = "logcounts" , batch = "sample" , 
                                       n.neigh = 10 , nPC = 50 , genes.predict = rownames(sce)){
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
  out = list(real = stat_real , predict = stat_predict)
  return(out)
}



get_z_scaled_distances_single_batch = function(sce , genes.all = rownames(sce) , n.neigh = 5 , nPC.all = 50){
  neighs.all = get_mapping(sce , genes = genes.all, batch = NULL, n.neigh = "all", nPC = nPC.all , get.dist = T)
  distances = neighs.all$distances
  distances_scaled = t( apply(distances , 1 , function(x) scale(x)) )
  rownames(distances_scaled) = rownames(distances)
  neighs.all$distances = distances_scaled
  return(neighs.all)
}

get_z_scaled_distances = function(sce , genes.all = rownames(sce) , batch = NULL, n.neigh = 5 , nPC.all = 50){
  if (is.null(batch)){
    out = get_z_scaled_distances_single_batch(sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    return(out)
  }
  else {
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    neighs.all = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      out =  get_z_scaled_distances_single_batch(current.sce , genes.all = genes.all , n.neigh = n.neigh , nPC.all = nPC.all)
    })
    names(neighs.all) = unique(batchFactor)
    return(neighs.all)
  }
}


get_preservation_score_single_batch = function(sce , neighs.all = NULL ,  genes.all = rownames(sce) , 
                                               genes.compare , n.neigh = 5 , nPC.all = 50 , nPC.compare = 200){
  if (is.null(neighs.all)){
    neighs.all = get_z_scaled_distances(sce, genes.all = genes.all, batch = NULL, n.neigh = n.neigh, nPC.all = nPC.all)
  }
  neighs.compare = get_mapping(sce , genes = genes.compare, batch = NULL, n.neigh = n.neigh, nPC = nPC.compare)
  neighs.all.cells_mapped = neighs.all$cells_mapped
  neighs.all.distances = neighs.all$distances
  n.cells = ncol(sce)
  
  score = lapply(1:nrow(neighs.compare) , function(i){
    cells = neighs.all.cells_mapped[i,]
    idx_all = c(1:n.neigh)
    idx_compare = which(cells %in% neighs.compare[i,] )
    
    dist_all = neighs.all.distances[i, idx_all]
    dist_compare = neighs.all.distances[i, idx_compare]
    current.score = median(-dist_compare)/median(-dist_all)
    out = data.frame(score = current.score )
    return(out)
  })
  score = do.call(rbind , score)
  score$cell = rownames(neighs.compare)
  return(score)
}


get_preservation_score = function(sce , neighs.all = NULL ,  genes.all = rownames(sce) , 
                                               genes.compare , batch = "sample" ,n.neigh = 5 , nPC.all = 50 , nPC.compare = 200){
  if (is.null(batch)){
    out = get_preservation_score_single_batch(sce , neighs.all = neighs.all ,  genes.all = genes.all, 
                                                         genes.compare = genes.compare, n.neigh = n.neigh , nPC.all = nPC.all , nPC.compare = nPC.compare)
    return(out)
  }
  else {
    if (is.null(neighs.all)){
      neighs.all = get_z_scaled_distances(sce , genes.all = genes.all , batch = batch, n.neigh = n.neigh, nPC.all = nPC.all)
    }
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    score = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.sce = sce[, idx]
      current.neighs.all = neighs.all[[which(names(neighs.all) == current.batch)]]
      out = get_preservation_score_single_batch(current.sce , neighs.all = current.neighs.all , genes.all = genes.all , 
                                                genes.compare = genes.compare, n.neigh = n.neigh , nPC.all = nPC.all, nPC.compare = nPC.compare)
      return(out)
    })
    score = do.call(rbind, score)
    return(score)
  }
}



get_preservation_score_simple = function(sce , neighs.all = NULL ,  genes.all = rownames(sce) , 
                                         genes.compare , batch = "sample", n.neigh = 5 , nPC.all = 50 , nPC.compare = 200){
  
  if (is.null(neighs.all)){
    if (is.null(batch)){
      neighs.all = get_mapping(sce , genes = genes.all, batch = NULL, n.neigh = "all", nPC = nPC.all , get.dist = T)
    }
    else {
      meta = as.data.frame(colData(sce))
      batchFactor = factor(meta[, colnames(meta) == batch])
      neighs.all = lapply(unique(batchFactor) , function(current.batch){
        idx = which(batchFactor == current.batch)
        current.sce = sce[, idx]
        current.neighs.all = get_mapping(current.sce , genes = genes.all, batch = NULL, n.neigh = "all", nPC = nPC.all , get.dist = T)
      })
      names(neighs.all) = unique(batchFactor)
    }
  }
  
  
  
  neighs.compare = get_mapping(sce , genes = genes.compare, batch = batch, n.neigh = n.neigh, nPC = nPC.compare , type = "together")
  if (!is.null(neighs.compare)){
    if (is.null(neighs.all)){
      neighs.all = get_mapping(sce , genes = genes.all, batch = batch, n.neigh = "all", nPC = nPC.all , get.dist = T , type = "together")
    }
    neighs.all.cells_mapped = neighs.all$cells_mapped
    neighs.all.distances = neighs.all$distances
    n.cells = ncol(sce)
    
    score = lapply(1:nrow(neighs.compare) , function(i){
      cells = neighs.all.cells_mapped[i,]
      idx_all = c(1:n.neigh)
      idx_compare = which(cells %in% neighs.compare[i,] )
      
      dist_all = neighs.all.distances[i, idx_all]
      dist_compare = neighs.all.distances[i, idx_compare]
      current.score = median(-dist_compare)/median(-dist_all)
      
      #current.bandwidth = mean(neighs.all.distances[i, ])
      #dist_all = exp(-1*dist_all/current.bandwidth)
      #dist_compare = exp(-1*dist_compare/current.bandwidth)
      #current.score = median(dist_compare)/median(dist_all)
      out = data.frame(score = current.score )
      return(out)
    })
    score = do.call(rbind , score)
    score$cell = rownames(neighs.compare)
    return(score)
  }
  else {
    return(NULL)
  }
}


get_preservation_score_simple = function(sce , neighs.all = NULL ,  genes.all = rownames(sce) , 
                                  genes.compare , batch = "sample", n.neigh = 5 , nPC.all = 50 , nPC.compare = 200){
  neighs.compare = get_mapping(sce , genes = genes.compare, batch = batch, n.neigh = n.neigh, nPC = nPC.compare , type = "together")
  if (!is.null(neighs.compare)){
    if (is.null(neighs.all)){
      neighs.all = get_mapping(sce , genes = genes.all, batch = batch, n.neigh = "all", nPC = nPC.all , get.dist = T , type = "together")
    }
    neighs.all.cells_mapped = neighs.all$cells_mapped
    neighs.all.distances = neighs.all$distances
    n.cells = ncol(sce)
    
    score = lapply(1:nrow(neighs.compare) , function(i){
      cells = neighs.all.cells_mapped[i,]
      idx_all = c(1:n.neigh)
      idx_compare = which(cells %in% neighs.compare[i,] )
      
      dist_all = neighs.all.distances[i, idx_all]
      dist_compare = neighs.all.distances[i, idx_compare]
      current.score = median(-dist_compare)/median(-dist_all)
      
      #current.bandwidth = mean(neighs.all.distances[i, ])
      #dist_all = exp(-1*dist_all/current.bandwidth)
      #dist_compare = exp(-1*dist_compare/current.bandwidth)
      #current.score = median(dist_compare)/median(dist_all)
      out = data.frame(score = current.score )
      return(out)
    })
    score = do.call(rbind , score)
    score$cell = rownames(neighs.compare)
    return(score)
  }
  else {
    return(NULL)
  }
}



get_mapping_2_external_dataset = function(sce_reference , sce_query , cluster_id , genes, nPC = 100, n.neigh = 5, skip.first = F, cosine = F){
  genes = intersect(genes,rownames(sce_reference))
  genes = intersect(genes,rownames(sce_query))
  
  sce_reference = sce_reference[genes , ]
  sce_reference = sce_reference[order(rownames(sce_reference)), ]
  sce_query = sce_query[genes , ]
  sce_query = sce_query[order(rownames(sce_query)), ]
  
  if (cosine){
    assay(sce_reference , "logcounts") = cosineNorm(logcounts(sce_reference))
    assay(sce_query , "logcounts") = cosineNorm(logcounts(sce_query))
  }
  
  sce_joint = cbind(assay(sce_query, "logcounts"), assay(sce_reference, "logcounts"))
  batchFactor = factor(c(as.character(sce_query$sample), as.character(sce_reference$sample)))
  
  mbpca = multiBatchPCA(sce_joint, batch = batchFactor, d = nPC)
  out = do.call(reducedMNN, mbpca)
  joint_pca = out$corrected
  current.knns = queryKNN( joint_pca[colnames(sce_reference),], joint_pca[colnames(sce_query),], k = n.neigh, 
                           get.index = TRUE, get.distance = FALSE)
  if (skip.first){
    cells.mapped = t( apply(current.knns$index, 1, function(x) colnames(sce_reference)[x[2:n.neigh]]) )
  }
  else {
    cells.mapped = t( apply(current.knns$index, 1, function(x) colnames(sce_reference)[x]) )
  }
  
  meta_reference = as.data.frame(colData(sce_reference))
  mapped_cluster = apply(cells.mapped , 1 , function(x) return(getmode(meta_reference[match(x, rownames(meta_reference)) , cluster_id] , c(1:n.neigh)) ))
  meta_query = as.data.frame(colData(sce_query))
  meta_query$mapped_cluster = mapped_cluster
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

generateSimilarity = function(SCE, k = 50, batchFactor = NULL, HVGs = NULL) {
  # SCE is a single cell experiment object containing the gene expression
  # in "logcounts" slot, otherwise a genes x cells matrix of logcounts
  # k is the number of nearest neighbours in the estimated KNN network
  # batchFactor is a factor matching columns of SCE specifying batches for MNN correction
  # after PCA
  # HVGs is optional set of genes to calculate similarity
  
  require(scran)
  require(SingleCellExperiment)
  require(bluster)
  require(igraph)
  require(scater)
  require(batchelor)
  # 
  # if (!"logcounts" %in% names(assays(SCE))) {
  #   require(scuttle)
  #   SCE <- logNormCounts(SCE)
  # }
  
  if (is.null(HVGs)) {
    fit = modelGeneVar(logcounts(SCE))
    HVGs = getTopHVGs(fit)
  }
  
  SCE <- runPCA(SCE, subset_row = HVGs)
  
  if (!is.null(batchFactor)) {
    SCE_corrected <- fastMNN(SCE, batch = batchFactor)
    PCs = reducedDim(SCE_corrected, "corrected")
  } else {
    PCs = reducedDim(SCE, "PCA")
  }
  
  graph = bluster::makeKNNGraph(PCs, k = k)
  V(graph)$name <- colnames(SCE)
  
  graph_sim = igraph::similarity(graph, method = "jaccard")
  rownames(graph_sim) <- V(graph)$name
  colnames(graph_sim) <- V(graph)$name
  
  return(graph_sim)
}

getSubsetUncertainty = function(SCE, batchFactor.sim = NULL,
                                querySCE = NULL,
                                subsetGenes = NULL,
                                k = 50, 
                                full_sim = NULL,
                                plot = FALSE,
                                plotAdditional = NULL,
                                verbose = FALSE,
                                jointBatchFactor = NULL,
                                returnAdditional = NULL,
                                ...) {
  
  # output is a named numeric vector of uncertainty values for each cell
  # if querySCE is provided, uncertainty values will include these
  # cells too
  
  # SCE is a SingleCellExperiment object of the reference dataset
  # querySCE is a SingleCellExperiment object of the query dataset,
  # if subsetGenes is NULL then the rownames of these are given as the 
  # subset
  # subsetGenes is a character vector of genes to subset with, this can
  # be NULL if querySCE is provided
  # k integer is the number of nearest neighbours
  # full_sim is a square matrix assumed to be the similarity of the 
  # reference data given as SCE, which can be generated a priori 
  # using generateSimilarity(SCE)
  # jointBatchFactor is a named factor that should have values for SCE and
  # querySCE if provided, which will be included as a batch via interaction
  # with the Reference and Query batch
  # returnAdditional is a character vector of any additional objects to be 
  # returned along with the uncertainty score, e.g. to also extract the 
  # joint PCs set returnAdditional = "jointPCs", or "g" for the plot
  
  require(igraph)
  require(BiocNeighbors)
  
  if (is.null(full_sim)) {
    full_sim = generateSimilarity(SCE, ...)
  }
  
  # combining SCE objects is nontrivial in general
  if (is.null(querySCE)) {
    if (is.null(subsetGenes)) stop("Either querySCE or subsetGenes needs to be provided")
    jointSCE = SCE[subsetGenes,]
    batchFactor = rep(c("Reference"), times = c(ncol(SCE)))
  } else {
    if (is.null(subsetGenes)) {
      jointSCE = cbind(SCE, querySCE)[rownames(querySCE),]
      batchFactor = rep(c("Reference", "Query"), times = c(ncol(SCE), ncol(querySCE)))
    } else {
      jointSCE = cbind(SCE, querySCE)[subsetGenes,]
      batchFactor = rep(c("Reference"), times = c(ncol(SCE)))
    }
  }
  
  # add the additional batch factor if given
  if (!is.null(jointBatchFactor)) {
    batchFactor <- interaction(batchFactor, jointBatchFactor[colnames(jointSCE)])
  }
  
  # extract similarity of the subsetted genes
  subset_sim = generateSimilarity(SCE, batchFactor = batchFactor.sim, HVGs = rownames(jointSCE))
  
  # concatenate and batch correct the reference and query datasets (if applicable)
  jointSCE <- logNormCounts(jointSCE)
  jointSCE <- runPCA(jointSCE)
  if (length(unique(batchFactor)) != 1) {
    jointSCE_corrected <- fastMNN(jointSCE, batch = batchFactor)
    jointPCs = reducedDim(jointSCE_corrected, "corrected")
  } else {
    jointPCs = reducedDim(jointSCE, "PCA")
  }
  
  # identify nearest neighbours
  tmp_r = jointPCs[colnames(SCE),]
  
  ref_knn = queryKNN(tmp_r,
                     query = jointPCs,
                     k = k)$index
  ref_knn_name = apply(ref_knn, 2, function(x) rownames(tmp_r)[x])
  rownames(ref_knn_name) <- rownames(jointPCs)
  
  # extract cell-specific uncertainty score
  uncertainty_scores = sapply(rownames(ref_knn_name), function(i) {
    if (verbose) print(i)
    ref_sim_nn = full_sim[ref_knn_name[i,], ref_knn_name[i,]]
    ref_sim_sub_nn = subset_sim[ref_knn_name[i,], ref_knn_name[i,]]
    
    stat = suppressWarnings({ks.test(c(ref_sim_nn[lower.tri(ref_sim_nn)]),
                                     c(ref_sim_sub_nn[lower.tri(ref_sim_sub_nn)]))$stat})
    names(stat) <- NULL
    return(stat)
  })
  
  # generate UMAP for plotting
  if (plot) {
    jointSCE$uncertainty = uncertainty_scores
    reducedDim(jointSCE, "UMAP") <- calculateUMAP(t(jointPCs))
    g = plotUMAP(jointSCE, colour_by = "uncertainty")
    if (!is.null(plotAdditional)) {
      # e.g. plotAdditional = list("celltype", scale_colour_manual(values = celltype_colours))
      require(patchwork)
      gAdditional = plotUMAP(jointSCE, colour_by = plotAdditional[[1]])
      gAll = gAdditional + plotAdditional[[2]] + labs(colour = plotAdditional[[1]]) + g
      print(gAll)
    } else {
      print(g)
    }
  }
  
  if (!is.null(returnAdditional)) {
    out = mget(c("uncertainty_scores", intersect(ls(), returnAdditional)))
    return(out)
  } else {
    return(uncertainty_scores)
  }
}


compare_gene_selections = function(sce , genes , genes_to_compare , corr.thresh = 1, rank.thresh = 1000){
  genes = genes[genes$rank <= rank.thresh , ]
  genes_to_compare = genes_to_compare[genes_to_compare$rank <= rank.thresh,]
  sce = sce[rownames(sce) %in% unique( c(as.character(genes$gene) , as.character(genes_to_compare$gene)) ) ,]
  counts = as.matrix(logcounts(sce))
  corr_stat = as.data.frame( cor(t(counts), method = "pearson") )
  
  stat_genes = lapply(1:nrow(genes) , function(i){
    current.corr_stat = corr_stat[, colnames(corr_stat) %in% genes_to_compare$gene]
    current.corr = current.corr_stat[rownames(current.corr_stat) == genes$gene[i], ]
    out = data.frame(corr = max(current.corr) , gene.which = colnames(current.corr_stat)[which.max(current.corr)])
    return(out)
  })
  stat_genes = do.call(rbind , stat_genes)
  genes = cbind(genes , stat_genes)
  genes = genes[genes$corr <= corr.thresh , ]
  
  return(genes)
}



