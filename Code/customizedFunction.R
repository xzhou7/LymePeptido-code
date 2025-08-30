library("future")
library("pbapply")
library("Matrix")
library("sctransform")
library("Seurat")
library("SeuratObject")

PrepSCTFindMarkers_force <- function(object, assay = "SCT", verbose = TRUE) {
  if (verbose && nbrOfWorkers() == 1) {
    my.lapply <- pblapply
  } else {
    my.lapply <- future_lapply
  }
  if (length(x = levels(x = object[[assay]])) == 1) {
    if (verbose) {
      message("Only one SCT model is stored - skipping recalculating corrected counts")
    }
    return(object)
  }
  observed_median_umis <- lapply(
    X = SCTResults(object = object[[assay]], slot = "cell.attributes"),
    FUN = function(x) median(x[, "umi"])
  )
  model.list <- slot(object = object[[assay]], name = "SCTModel.list")
  median_umi.status <- lapply(X = model.list,
                              FUN = function(x) { return(tryCatch(
                                expr = slot(object = x, name = 'median_umi'),
                                error = function(...) {return(NULL)})
                              )})
  if (any(is.null(x = unlist(x = median_umi.status)))){
    # For old SCT objects  median_umi is set to median umi as calculated from obserbed UMIs
    slot(object = object[[assay]], name = "SCTModel.list") <- lapply(X = model.list,
                                                                     FUN = UpdateSlots)
    SCTResults(object = object[[assay]], slot = "median_umi") <- observed_median_umis
    
  }
  model_median_umis <- SCTResults(object = object[[assay]], slot = "median_umi")
  min_median_umi <- min(unlist(x = observed_median_umis), na.rm = TRUE)
  if (verbose) {
    message(paste0("Found ",
                   length(x = levels(x = object[[assay]])),
                   " SCT models.",
                   " Recorrecting SCT counts using minimum median counts: ",
                   min_median_umi))
  }
  umi.assay <- unique(
    x = unlist(
      x = SCTResults(object = object[[assay]], slot = "umi.assay")
    )
  )
  if (length(x = umi.assay) > 1) {
    stop("Multiple UMI assays are used for SCTransform: ",
         paste(umi.assay, collapse = ", ")
    )
  }
  umi.layers <- Layers(object = object, assay = umi.assay, search = 'counts')
  if (length(x = umi.layers) > 1) {
    object[[umi.assay]] <- JoinLayers(
      object = object[[umi.assay]],
      layers = "counts", new = "counts")
  }
  raw_umi <- GetAssayData(object = object, assay = umi.assay, slot = "counts")
  corrected_counts <- Matrix(
    nrow = nrow(x = raw_umi),
    ncol = ncol(x = raw_umi),
    data = 0,
    dimnames = dimnames(x = raw_umi),
    sparse = TRUE
  )
  cell_attr <- SCTResults(object = object[[assay]], slot = "cell.attributes")
  model_pars_fit <- lapply(
    X = SCTResults(object = object[[assay]], slot = "feature.attributes"),
    FUN = function(x) x[, c("theta", "(Intercept)", "log_umi")]
  )
  arguments <- SCTResults(object = object[[assay]], slot = "arguments")
  model_str <- SCTResults(object = object[[assay]], slot = "model")
  set_median_umi <- rep(min_median_umi, length(levels(x = object[[assay]])))
  names(set_median_umi) <- levels(x = object[[assay]])
  set_median_umi <- as.list(set_median_umi)
  all_genes <- rownames(x = object[[assay]])
  # correct counts
  my.correct_counts <- function(model_name){
    model_genes <- rownames(x = model_pars_fit[[model_name]])
    x <- list(
      model_str = model_str[[model_name]],
      arguments = arguments[[model_name]],
      model_pars_fit = as.matrix(x = model_pars_fit[[model_name]]),
      cell_attr = cell_attr[[model_name]]
    )
    cells <- rownames(x = cell_attr[[model_name]])
    umi <- raw_umi[all_genes, cells]
    
    umi_corrected <- correct_counts(
      x = x,
      umi = umi,
      verbosity = 0,
      scale_factor = min_median_umi
    )
    missing_features <- setdiff(x = all_genes, y = rownames(x = umi_corrected))
    corrected_counts.list <- NULL
    gc(verbose = FALSE)
    empty <- SparseEmptyMatrix(nrow = length(x = missing_features), ncol = ncol(x = umi_corrected))
    rownames(x = empty) <- missing_features
    colnames(x = umi_corrected) <- colnames(x = umi_corrected)
    
    umi_corrected <- rbind(umi_corrected, empty)[all_genes,]
    
    return(umi_corrected)
  }
  corrected_counts.list <- my.lapply(X = levels(x = object[[assay]]),
                                     FUN = my.correct_counts)
  names(x = corrected_counts.list) <- levels(x = object[[assay]])
  
  corrected_counts <- do.call(what = MergeSparseMatrices, args = corrected_counts.list)
  corrected_counts <- as.sparse(x = corrected_counts)
  corrected_data <- log1p(x = corrected_counts)
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "counts",
                                           new.data = corrected_counts)})
  suppressWarnings({object <- SetAssayData(object = object,
                                           assay = assay,
                                           slot = "data",
                                           new.data = corrected_data)})
  SCTResults(object = object[[assay]], slot = "median_umi") <- set_median_umi
  return(object)
}

PrepSCTFindMarkers.V5 <- function(object, assay = "SCT", umi.assay = "RNA", layer = "counts", verbose = TRUE) {
  layers <- Layers(object = object[[umi.assay]], search = layer)
  dataset.names <- gsub(pattern = paste0(layer, "."), replacement = "", x = layers)
  for (i in seq_along(along.with = layers)) {
    l <- layers[i]
    counts <- LayerData(
      object = object[[umi.assay]],
      layer = l
    )
  }
  cells.grid <- DelayedArray::colAutoGrid(x = counts, ncol = min(length(Cells(object)), ncol(counts)))
}

MergeSparseMatrices <- function(...) {
  
  colname.new <- character()
  rowname.new <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()
  
  for (mat in list(...)) {
    colname.old <- colnames(x = mat)
    rowname.old <- rownames(x = mat)
    
    # does not check if there are overlapping cells
    colname.new <- union(x = colname.new, y = colname.old)
    rowname.new <- union(x = rowname.new, y = rowname.old)
    
    colindex.new <- match(x = colname.old, table = colname.new)
    rowindex.new <- match(x = rowname.old, table = rowname.new)
    
    ind <- summary(object = mat)
    # Expand the list of indices and x
    i <- c(i, rowindex.new[ind[,1]])
    j <- c(j, colindex.new[ind[,2]])
    x <- c(x, ind[,3])
  }
  
  merged.mat <- sparseMatrix(i=i,
                             j=j,
                             x=x,
                             dims=c(length(rowname.new), length(colname.new)),
                             dimnames=list(rowname.new, colname.new))
  return (merged.mat)
}





