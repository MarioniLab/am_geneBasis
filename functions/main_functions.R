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
  colnames(stat) = c(cluster.id , frac_correctly_mapped)
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


get_mapping_single_batch = function(sce , genes = rownames(sce), n.neigh = 5, nPC = 50 , get.dist = F){
  require(irlba)
  require(BiocNeighbors)
  set.seed(32)
  
  current.sce = sce[genes , ]
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


get_mapping_many_batches = function(sce , genes = rownames(sce), batch = "sample", n.neigh = 5, nPC = 50 , get.dist = F){
  require(batchelor)
  require(BiocNeighbors)
  set.seed(32)
  
  current.sce = sce[genes , ]
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


get_mapping = function(sce , genes = rownames(sce), batch = "sample", n.neigh = 5, nPC = 50 , get.dist = F, type = "per batch"){
  #require(BiocSingular)
  #require(BiocParallel)
  require(BiocNeighbors)
  
  #mcparam = MulticoreParam(workers = ncores)
  #register(mcparam)
  #set.seed(32)
  
  if (is.null(batch)){
    out = get_mapping_single_batch(sce , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist)
  }
  else if (type == "together") {
    out = get_mapping_many_batches(sce , genes = genes, batch = batch , n.neigh = n.neigh, nPC = nPC , get.dist = get.dist)
  }
  else if (type == "per batch") {
    
    meta = as.data.frame(colData(sce))
    batchFactor = factor(meta[, colnames(meta) == batch])
    
    neighs = lapply(unique(batchFactor) , function(current.batch){
      idx = which(batchFactor == current.batch)
      current.neighs = get_mapping_single_batch(sce[, idx] , genes = genes, n.neigh = n.neigh, nPC = nPC , get.dist = get.dist)
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



get_mapping_2_external_dataset = function(sce_reference , sce_query , cluster_id , genes, nPC = 100, n.neigh = 5, skip.first = F){
  genes = intersect(genes,rownames(sce_reference))
  genes = intersect(genes,rownames(sce_query))
  
  sce_reference = sce_reference[genes , ]
  sce_reference = sce_reference[order(rownames(sce_reference)), ]
  sce_query = sce_query[genes , ]
  sce_query = sce_query[order(rownames(sce_query)), ]
  
  assay(sce_reference , "cosineNorm") = cosineNorm(logcounts(sce_reference))
  assay(sce_query , "cosineNorm") = cosineNorm(logcounts(sce_query))
  
  sce_joint = cbind(assay(sce_query, "cosineNorm"), assay(sce_reference, "cosineNorm"))
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
