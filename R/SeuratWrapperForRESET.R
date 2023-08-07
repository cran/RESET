#
# SeuratWrapperForRESET.R
#
# A Seurat-specific wrapper around RESET.R
#
# @author rob.frost@dartmouth.edu
#

#
# Computes RESET scores using the normalized gene expression values for the cells in a Seurat object.
# 
# Inputs:
#
# -seurat.data: The Seurat object that holds the scRNA-seq data. 
#  Assumes PCs have already computed on centered and scaled normalized counts.
# -num.pcs: Number of PCs to use in the REST method. If not specified, will use all computed PCs. 
#  The projection of the scRNA-seq data onto these PCs will be used as the X.test matrix for the reset() call.
# -gene.set.collection: List of m gene sets for which scores are computed.
#  Each element in the list corresponds to a gene set and the list element is a vector
#  of indices for the genes in the set. The index value is defined relative to the
#  order of genes in the gene.exprs matrix. Gene set names should be specified as list names.
#  See createVarSetCollection() function in REST.R for utility function that can be used to 
#  help generate this list of indices.
# 
#  See reset() docs for other arguments
#
# Output: Updated Seurat object with a new Assay object called "RESET" that contains:
#
#   The matrix of computed cell-level gene set scores in slot "data"
#   The vector of overall variable set scores in the feature metadata column "RESET"
#
resetForSeurat = function(seurat.data, num.pcs, gene.set.collection,
    k=10, random.threshold, k.buff=0, q=0, test.dist="normal", norm.type="2",
    per.var=FALSE) {
  
  if (!requireNamespace("Seurat", quietly=TRUE)) {
    stop("Seurat package not available!")
  }
  
  if (missing(seurat.data)) {
    stop("Missing Seurat object!")
  }
  
  if (missing(gene.set.collection)) {
    stop("Missing gene set collection list!")
  }
  
  if (missing(random.threshold)) {
    random.threshold = k
    message("Setting random.threhold to specified k of ", k)
  } else if (random.threshold < k) {
    stop("random.threshold cannot be than k!")
  }
  
  # Use the active.assay slot to determine whether SCTransform or the standard normalization logic was used
  if (seurat.data@active.assay == "RNA") {

    # Use the normalized counts
    message("Using log-normalized RNA counts...")
    normalized.counts = seurat.data@assays$RNA@data
    
  } else if (seurat.data@active.assay == "SCT") {
    
    # Use the corrected counts. The SCT correction process reverses the regression model
    # to generate counts that approximate what would be found if all cells were sequenced to the same depth.
    message("Using SCTransform normalized RNA counts...")
    normalized.counts = seurat.data@assays$SCT@counts
        
  } else {    
    stop("Unsupported active assay: ", seurat.data@active.assay)
  }    

  # Get the projection onto the to PCs  
  X.pca = seurat.data@reductions$pca@cell.embeddings
  if (!missing(num.pcs)) {
    if (num.pcs < ncol(X.pca)) {
      message("Testing reconstruction against top ", num.pcs, " PCs")
      X.pca = X.pca[,1:num.pcs]      
    } else {
      warning("num.pcs larger than number of PCs computed on Seurat object")
    }
  } else {
    message("Testing reconstruction against all ", ncol(X.pca), " computed PCs")
  }
  
  # Execute RESET 
  reset.results = reset(X=Matrix::t(normalized.counts), 
                      X.test=X.pca, 
                      center.X=TRUE,
                      scale.X=TRUE,
                      center.X.test=FALSE, # assume PCs were computed on centered data
                      scale.X.test=FALSE, # assume PCs were computed on scaled data
                      k=k, random.threshold=random.threshold, k.buff=k.buff, q=q,
                      var.sets=gene.set.collection,
                      test.dist=test.dist,
                      norm.type=norm.type,
                      per.var=per.var)

  # Create Assay object to store the sample level scores
  reset.assay = Seurat::CreateAssayObject(counts = t(reset.results$S))
  # Add the overall scores as metadata using name "RESET
  reset.assay = Seurat::AddMetaData(reset.assay, reset.results$v, "RESET")
  seurat.data[["RESET"]] = reset.assay    
  
  return (seurat.data)
}  

#
# Convert the overall and cell-level RESET scores to/from per-variable scores.
#
# Inputs:
#
#   seurat.data: Seurat object returned from resetForSeurat()
#   var.sets: List containing variable set indices (this must be the same list used to call reset()).
#   to.per.var: If true, converts to per-variable scores, i.e., divides scores by variable set size.
#             If false, converts from per-variable scores, i.e., multiplies scores by variable set size.
#
# Returns a modified Seurat object
#
convertToPerVarScoresForSeurat = function(seurat.data, gene.set.collection, to.per.var=TRUE) {
  if (missing(seurat.data)) {
    stop("seurat.data must be specified!")
  }
  if (missing(gene.set.collection)) {
    stop("gene.set.collection must be specified!")
  }
  
  # Compute the size of each variable set
  gene.set.size = unlist(lapply(gene.set.collection, length))
  
  # Scale the sizes by the mean size
  # Thiw will prevent the scores from dramatically changing in magnitude.
  gene.set.size = gene.set.size/mean(gene.set.size)
  
  # Extract the overall and cell-level scores
  overall.scores = seurat.data@assays$RESET@meta.features$RESET
  cell.scores = seurat.data@assays$RESET@data
  if (to.per.var) {
    overall.scores = overall.scores/gene.set.size  
    # scale each column (which correspond to cells)
    cell.scores = apply(cell.scores, 2, function(x) {
      return (x/gene.set.size)        
    })
  } else {
    overall.scores = overall.scores*gene.set.size  
    cell.scores = apply(cell.scores, 2, function(x) {
      return (x*gene.set.size)        
    })
  }
  # Update the Seurat object
  seurat.data@assays$RESET@meta.features$RESET = overall.scores
  seurat.data@assays$RESET@data = cell.scores
  
  return (seurat.data)
}
