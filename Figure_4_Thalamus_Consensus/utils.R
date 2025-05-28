calc_entropy <- function(u) {
  p <- u[u>0]
  p <- p/sum(p)
  -sum(p*log(p))
}


# given a matrix of labels, calculate all pairwise ARIs
calc_aris <- function(m, flavour="ARI") {
  a <- diag(ncol(m))
  for(i in 1:(ncol(m)-1))
    for(j in 2:ncol(m)) {
      if(flavour=="ARI") {
        require(mclust)
        a[i,j] <- a[j,i] <- mclust::adjustedRandIndex(m[,i], m[,j])
      } else if(flavour=="sARI") {
        
      }
    }
  rownames(a) <- colnames(a) <- colnames(m)
  a
}

num_methods <- function(u) {
  s <- sapply(strsplit(u, "[_\\.]"), .subset, 1)
  g <- grep("label",s)
  if(length(g)) 
    s <- s[-g]
  s %>% table %>% length
}


calc_summary <- function(tree, arism, k = 3, gtcol = "label") {
  ct <- cutree(tree, k=k)
  tct <- table(ct)
  data.frame(group = names(tct), 
             number = as.integer(tct),
             ari_to_truth = sapply(names(tct),
                                   function(u) median(arism[ct==u & names(ct)!=gtcol,gtcol])),
             ari_block = sapply(names(tct),
                                function(u) {z <- arism[ct==u,ct==u]; median(z[upper.tri(z)]) }),
             k = k,
             gt_in_group = (names(tct) == ct[gtcol]),
             number_distinct_methods = sapply(names(tct),
                                              function(u) num_methods(names(ct)[ct==u])))
}


save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Add an option here to allow alignment of only part of the columns 
align_classes <- function(d, ref, columns=NULL) {
  
  bcs <- d %>% as.data.frame

  # Adding a columns filtering here
  if (!is.null(columns)){
    all_cols <- colnames(bcs)
    bcs <- bcs[,union(columns, ref)]
  }

  for(j in 1:ncol(bcs)){
    bcs[,j] <- factor(bcs[,j], levels = unique(bcs[,j])) %>% as.numeric %>% as.factor()
  }
  levs <- lapply(bcs, levels)
  n_levs <- sapply(levs, length)
  
  cols_to_change <- setdiff(colnames(bcs), ref)
  
  for(i in cols_to_change) {
    if(n_levs[i] != n_levs[ref]) {
      message(paste0(i," has ",n_levs[i],
                     " levels; reference (", ref, 
                     ") has ",n_levs[ref], " (not modifying)"))
      next
    }
    hung <- clue::solve_LSAP(table(bcs[,ref], bcs[,i]), 
                             maximum = TRUE)
    lookup <- cbind(seq_along(hung), levs[[i]][hung])
    levels(bcs[,i]) <- lookup[order(as.integer(lookup[,2])),1]
    bcs[,i] <- as.factor(as.character(bcs[,i]))
  }

  # Adding back the remaining columns 
  if (!is.null(columns)){
    bcs <- cbind(bcs, as.data.frame(d)[, setdiff(all_cols, colnames(bcs))])
  }
  
  return(bcs)
}

spot_entropy <- function(spatial_coords, label, k=6) {
  suppressPackageStartupMessages(require(dbscan))
  knns <- kNN(spatial_coords, k=k)
  label <- as.factor(label)
  neighb_labels <- apply(knns$id, 2, function(u) label[u])
  apply(neighb_labels, 
        1, function(u) calc_entropy(table(factor(u,
                                                 levels=levels(label)))))
}


get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    reducedDim_file = NA,
    assay_name = "counts",
    reducedDim_name = "reducedDim") {
  # Make sure row name columns are treated as characters 
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1, 
                        colClasses = c("character", rep(NA, ncol(read.delim(feature_file, header=TRUE))-1)))
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1,
                        colClasses = c("character", rep(NA, ncol(read.delim(observation_file, header=TRUE))-1)))
  
  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1,
                            colClasses = c("character", rep(NA, ncol(read.delim(coord_file, header=TRUE))-1)))
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])
  
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowData, colData = colData, spatialCoords = coordinates
  )
  
  if (!is.na(matrix_file)) {
    assay(spe, assay_name, withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
    #assay(spe, "logcounts", withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
  }
  
  # Filter features and samples
  if ("selected" %in% colnames(rowData(spe))) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if ("selected" %in% colnames(colData(spe))) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }
  
  if (!is.na(reducedDim_file)) {
    dimRed <- read.delim(reducedDim_file, stringsAsFactors = FALSE, row.names = 1)
    reducedDim(spe, reducedDim_name) <- as.matrix(dimRed[colnames(spe), ])
  }
  return(spe)
}



plotter_fun <- function(u,v,z, point_size = 3) {
  
  df <- data.frame(colData(spe)[,c("row","col")], 
                   gene=logcounts(spe)[u,],
                   ground_truth = gt)
  df$highlight <- df$ground_truth %in% v
  
  
  ggplot(df, aes(x=row, y=col, fill=gene, colour=highlight)) + 
    geom_point(shape=21, size = point_size) +
    theme_classic() +
    # theme(legend.position = "none") +
    xlab("") + ylab("") +
    scale_colour_manual(values = gt_col) +
    scale_fill_gradient(low="gray95", high="deeppink4") +
    scale_x_reverse(expand = expansion(0,1)) +
    scale_y_continuous(expand = expansion(0,1)) +
    ggtitle(paste0(u,": ",z," (",paste0(v,collapse = ","),")"))
}


## Adopted from https://github.com/keyalone/EnSDD/blob/main/R/utils.R, modified for JSD instead of L2 norm


########### ensemble strategy ###############
#' The adaptive weighted ensemble-based learning method to integrate the multiple binary spots similarity matrix
#'
#'
#' @importFrom parallel makeCluster stopCluster parApply
#' @importFrom abind abind
#'
#' @TODO: Make Reuslts.clustering a sparse matrix object. Does it need to change anything?
#' @param Results.clustering a list contains all the results of individual similarity matrix. The elements of list is a matrix, spots * spots.
#' @param lambda hyper-parameter constrain the weight of individual methods for ensemble. If the parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param prob.quantile numeric of probabilities with values in [0,1]. Default setting is 0.5.
#' @param niter a positive integer represents the maximum number of updating algorithm. Default setting is 100.
#' @param epsilon a parameter represents the stop criterion.
#'
#' @return a list contains a matrix of the ensemble similarity of spots and a vector of the weight assigned to base results.
#'
#'@export

solve_ensemble <- function(Results.clustering, 
                          lambda = NULL, 
                          prob.quantile = 0.5,
                          niter = 100, 
                          epsilon = 1e-5,
                          verbose = FALSE){
  suppressPackageStartupMessages(require(future.apply))
  
  plan(multisession, workers = parallelly::availableCores() - 1)

  options(digits = 7)
  # Results.clustering <- Results.clustering.all[[1]]
  num.methods <- length(Results.clustering)
  num.spots <- nrow(Results.clustering[[1]])
  num.cell.type <- ncol(Results.clustering[[1]])

  ## initialization V by the mean of individual values
  w <- c(rep(1/num.methods, num.methods))
  H <-  Reduce("+", Map("*", Results.clustering, w))

  if(is.null(lambda)){
    cat("We will adpote a value for lambda in our algorithm...", "\n")
  }

  k <- 1

  while (k <= niter) {
    st <- proc.time()
    if(k == 1){
      loss_all_temp <- 0
      # Generate the first loss value 
      temp2 <-  future_sapply(Results.clustering, JSD_Matrix, Y = H)
      # Empricial estimation of lambda in the paper
      if(is.null(lambda)){
        lambda <- quantile(temp2, probs = prob.quantile)
      }
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  future_sapply(Results.clustering, JSD_Matrix, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    ##### update H
    H <-  Reduce("+", Map("*", Results.clustering, w))

    # Objective function loss
    loss_main <- sum(sapply(Results.clustering, JSD_Matrix, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy

    if(k == niter){
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")}

    # Stopping criteria
    diff_iter <- abs(loss_all - loss_all_temp)
    delta <- proc.time() - st 
    if (verbose){
      cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
          "loss_all: ", loss_all, "lambda: ", lambda, "diff:", 
          diff_iter, "epoch_time:", delta[3], "\n")
    }

    if(diff_iter < epsilon | k >= niter){
      break
      }
    k <- k + 1

  }
  colnames(H) <- colnames(Results.clustering[[1]])
  return(list(H = H, w = w))

}

JSD_Matrix <- function(X, Y, epi=1e-10){
  # Flatten the matrix, get probability, add a pseudo-count
  x <- pmax(as.vector(X)/sum(X), epi)
  y <- pmax(as.vector(Y)/sum(Y), epi)
  m <- (x + y)/2

  # Make the JSD bounded by 1 by using log2
  JSD <- sqrt(0.5 * (sum(x * log2(x / m)) +
                     sum(y * log2(y / m))))
  return(JSD)
}
########### ensemble strategy ###########

#' Get binary similarity matrix from a cluster vector,
get_binary_matrix <- function(cluster_vector){
  suppressPackageStartupMessages(require(Matrix))
  N <- length(cluster_vector)

  pairs <- which(outer(cluster_vector, cluster_vector, FUN = "=="), arr.ind = TRUE)
  # Remove diagonal entries
  pairs <- pairs[pairs[,1] != pairs[,2], ]

  # Create a sparse matrix using the pairs of indices
  bm <- Matrix::sparseMatrix(
    i = pairs[,1],
    j = pairs[,2],
    x = rep(1, nrow(pairs)),
    dims = c(N, N)
  )

  return(bm)
}

#' Get cluster label from Leiden clustering of the consensus binary matrix.
#' Using also binary_search function to ensure the proper number of clusters
get_cluster_label <- function(binary_matrix,
                              n_clust_target,
                              resolution_update = 2,
                              resolution_init = 0.5,
                              num_rs = 100,
                              tolerance = 1e-5,
                              verbose=FALSE,
                              min_target=NULL){
  suppressPackageStartupMessages(require(igraph))
  # require(leidenalg)

  graph <- igraph::graph_from_adjacency_matrix(binary_matrix,
                                       mode = "upper",
                                       weighted = TRUE,
                                       diag = FALSE)
  write.table(as_data_frame(graph, what="edges") %>% sort("from"), 
              "edgeList.txt", sep="\t", 
              col.names=FALSE, row.names=FALSE)
  if (verbose){cat("Weighted neighborhood graph created... \n")}

  # Initialize boundaries
  lb <- rb <- NULL
  n_clust <- -1
  if (is.null(min_target)){
    min_target <- nrow(binary_matrix)/(10*n_clust_target)
  }
  print(min_target)

  get_clusters <- function(graph, resolution, mt=min_target){

    if (!file.exists("networkanalysis-1.3.0.jar")){
      system("wget https://repo1.maven.org/maven2/nl/cwts/networkanalysis/1.3.0/networkanalysis-1.3.0.jar")
    }
    system(sprintf("java -cp networkanalysis-1.3.0.jar nl.cwts.networkanalysis.run.RunNetworkClustering -w -m %1.0f -r %s -o cluster.txt edgeList.txt", mt, resolution))

    results <- read.table("cluster.txt", header=TRUE, row.names=NULL)[,2]

    return(results)
  }

  res <-  resolution_init
  result <- get_clusters(graph, resolution = res)
  # Adjust cluster_ids extraction per method
  n_clust <- length(unique(result))
  if (verbose){cat(sprintf("Boundary search starts..res = %s, n_clust=%s \n", res, n_clust))}
  if (n_clust > n_clust_target) {
    while (n_clust > n_clust_target && res > 1e-5) {
      rb <- res
      res <- res / resolution_update
      result <- get_clusters(graph, resolution = res)
      n_clust <- length(unique(result))
      if (verbose){cat(sprintf("Boundary search..lb = %s, rb = %s, n_clust=%s \n", res, rb, n_clust))}
    }
    lb <- res
  } else if (n_clust < n_clust_target) {
    while (n_clust < n_clust_target) {
      lb <- res
      res <- res * resolution_update
      result <- get_clusters(graph, resolution = res)
      n_clust <- length(unique(result))
      if (verbose){cat(sprintf("Boundary search..lb = %s, rb = %s, n_clust=%s \n", lb, res, n_clust))}
    }
    rb <- res
  }
  if (n_clust == n_clust_target) {lb = rb = res}

  i <- 0
  if (verbose){cat(sprintf("Boundary search done. lb = %s, rb = %s, res = %s, n_clust=%s \n", lb, rb, res, n_clust))}

  while ((rb - lb > tolerance || lb == rb) && i < num_rs) {
    mid <- sqrt(lb * rb)
    # message("Resolution: ", mid)
    result <- get_clusters(graph, resolution = mid)
    n_clust <- length(unique(result))
    min_clust_size <- min(table(result))
    if (verbose){cat(sprintf("iter: %s, res: %s, n_clust: %s, min_clust_size: %s \n", i,  mid, n_clust, min_clust_size))} # nolint
    if (n_clust == n_clust_target && min_clust_size >= min_target) break
    if (n_clust > n_clust_target) {
      rb <- mid
    } else if (n_clust < n_clust_target){
      lb <- mid
    } else if (n_clust == n_clust_target && min_clust_size < min_target){
      rb <- mid
      if (rb == lb){lb <- lb*0.9}
    }
    i <- i + 1
  }

  # Warning if target not met
  if (n_clust != n_clust_target) {
    warning(sprintf("Warning: n_clust = %d not found in binary search, return best approximation with res = %f and n_clust = %d. (rb = %f, lb = %f, i = %d)", n_clust_target, mid, n_clust, rb, lb, i))
  }

  cluster_vector <- result
  return(cluster_vector)
}

