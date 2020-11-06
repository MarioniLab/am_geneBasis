# GENE SELECTION SCRIPTS
getHVGs = function(sce, min.mean = 1e-3, dir = "cluster", FDR = 0.05){
  require(biomaRt)
  trend = scran::trendVar(sce, use.spikes = FALSE, loess.args = list(span = 0.05))
  decomp = scran::decomposeVar(sce, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  if (dir == "cluster"){
    ychr = read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
  } else if (dir == "local"){
    ychr = read.table("/Users/alsu/Develop/FetalAlcoholSyndrome/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
  }
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$FDR < FDR])
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


getMarkers = function(sce , pval.type = "some", test , FDR.thresh){
  # assumes there is a column called celltype and assay called logcounts
  require(scran)
  # get potential relevant genes
  markers <- findMarkers(sce , groups=sce$celltype, direction = "up", pval.type=pval.type, test = test, assay.type = "logcounts")
  
  # put together
  CTs = names(markers)
  markers = lapply(1:length(CTs), function(i){
    current.markers = as.data.frame(markers[[i]])
    current.markers = current.markers[!is.na(current.markers$FDR) & current.markers$FDR < FDR.thresh , ]
    if (nrow(current.markers) > 0){
      out = data.frame( celltype = CTs[i], gene = rownames(current.markers))
      return(out)
    }
  })
  markers = do.call(rbind,markers)
  return(markers)
}




# get genes that are highly correlated with the chosen one

getCorrelatedGenes = function(corr.thresh , counts , gene){
  counts.rest = counts[, !colnames(counts) == gene]
  counts.gene = counts[, colnames(counts) == gene]
  corr.genes = lapply(1:ncol(counts.rest) , function(i){
    out = data.frame(gene = colnames(counts.rest)[i] , corr = cor(counts.gene , counts.rest[,i] , method = "pearson"))
    return(out)
  })
  corr.genes = do.call(rbind , corr.genes)
  corr.genes = corr.genes[corr.genes$corr > corr.thresh , ]
  return( corr.genes$gene) 
}



#MAPPING FUNCTIONS

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


# load embryo 8.5 data
load_embryo_8.5 = function(dir = "local"){
  
  require(SingleCellExperiment)
  require(scater)
  
  if (dir == "cluster") {
    data.dir = "/nfs/research1/marioni/alsu/spatial/mouse_embryo/data/8_5/source/"
  } else {
    data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/"
  }
  
  sce = readRDS( paste0( data.dir , "E8.5_sce_filt_unlabelled.Rds"))
  # add normalization by libsize 
  assay(sce, "cpm") <- logcounts(scater::logNormCounts(sce, size_factors = sce$total))
  assay(sce, "cpm_wo_xist") <- logcounts(scater::logNormCounts(sce, size_factors = as.numeric( sce$total - counts( sce["Xist"] )) ))
  
  
  meta = colData(sce)
  meta = data.frame(meta)
  
  # rename Cavin3 --> Prkcdbp
  rownames.sce = rownames(sce)
  rownames.sce[rownames.sce == "Cavin3"] = "Prkcdbp"
  rownames(sce) = rownames.sce
  
  assign("sce", sce, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  invisible(0)
}


getSegmentationVerticesDF = function(DF,
                                     xname = "segmentation_vertices_x_global",
                                     yname = "segmentation_vertices_y_global",
                                     othercols = c("uniqueID","z")) {
  # DF is a DataFrame object
  # othercols is the others to keep
  
  long_x = unlist(DF[,xname])
  long_y = unlist(DF[,yname])
  
  if (length(long_x) != length(long_y)) stop("x and y need to be same length")
  
  long_xy = data.frame(
    long_x,
    long_y
  )
  colnames(long_xy) <- c(xname, yname)
  
  long_DF = cbind(
    rep(DF[,othercols], times = unlist(lapply(DF[,xname], length))),
    long_xy
  )
  
  return(as.data.frame(long_DF))
}

add_scalebar = function(dist_um = 250, x = 2.75, y = -3.20, ...) {
  # this is a quantity to add to an existing ggplot
  
  # dist_um is the distance for the scalebar in um, default 250um
  # x and y are the coordinates to place the scalebar
  # usage: to add to a ggplot object like "g + add_scalebar()"
  
  # useful optional argument is box.col = "white"
  
  # need to make sure that ggplot of interest doesn't have group aesthetic
  # defined in the ggplot, but in the geom_polygon itself
  
  require(ggsn)
  add = ggsn::scalebar(location = "bottomright",
                       dist = dist_um/227.74, dist_unit = "units", 
                       transform = FALSE,
                       x.min = x, x.max = x,
                       y.min = y, y.max = y,
                       height = 0.2, 
                       box.fill = "black",
                       st.size = 0,
                       inherit.aes = FALSE,
                       ...)
  return(add)
}

# function to rotate the embryos
rotateDF = function(DF, 
                    xname = "segmentation_vertices_x_global_affine", 
                    yname = "segmentation_vertices_y_global_affine", 
                    ang = 0) {
  # ang is a numeric vector named three values corresponding to embryo1, embryo2, and embryo3
  
  ang_long = ang[as.character(DF$embryo)]
  ang_rad = ang_long/180
  
  x = DF[,xname]
  y = DF[,yname]
  
  x_turn = x*cos(ang_rad) - y*sin(ang_rad)
  y_turn = x*sin(ang_rad) + y*cos(ang_rad)
  
  # reset the columns and then return the DF
  
  DF[,xname] <- x_turn
  DF[,yname] <- y_turn
  
  return(DF)
  
}


get_suitable4seq_genes = function(dir = "local", system , mean.thresh , max.thresh){
  if (dir == "cluster"){
    root.dir = "/nfs/research1/marioni/alsu/geneBasis/"
  } 
  else if (dir == "local"){
    root.dir = "/Users/alsu/Develop/geneBasis/"
  }
  
  stat = read.table(paste0(root.dir , "data/expression_stat_perGene/stat__" , system , ".tab"), header = T, sep = "\t")
    
  idx = stat$mean.counts < mean.thresh & stat$max.counts < max.thresh
  print(paste0(sum(idx) , " genes are kept for the analysis"))
  return(stat$gene[idx])
}