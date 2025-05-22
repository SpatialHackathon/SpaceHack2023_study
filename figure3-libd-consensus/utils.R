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


align_classes <- function(d, ref) {
  
  bcs <- d %>% as.data.frame
  for(j in 1:ncol(bcs))
    bcs[,j] <- as.factor(bcs[,j])
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
  
  bcs
  
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
    if(k == 1){
      loss_all_temp <- 0
      # Generate the first loss value 
      temp2 <-  sapply(Results.clustering, L2_norm, Y = H)
      # Empricial estimation of lambda in the paper
      if(is.null(lambda)){
        lambda <- quantile(temp2, probs = prob.quantile)
      }
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  sapply(Results.clustering, L2_norm, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    ##### update H
    H <-  Reduce("+", Map("*", Results.clustering, w))

    # Objective function loss
    loss_main <- sum(sapply(Results.clustering, L2_norm, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy

    if(k == niter){
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")}

    # Stopping criteria
    diff_iter <- abs(loss_all - loss_all_temp)
    if (verbose){
      cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
          "loss_all: ", loss_all, "lambda: ", lambda, "diff",
          diff_iter, "\n")
    }

    if(diff_iter < epsilon | k >= niter){
      break
      }
    k <- k + 1

  }
  colnames(H) <- colnames(Results.clustering[[1]])
  return(list(H = H, w = w))

}

L2_norm <- function(X, Y){
  return(sqrt(sum((X-Y)^2)))
}
########### ensemble strategy ###########
