#
# Reconstruction Set Test (RESET)
#
# @author rob.frost@dartmouth.edu
#

#
# Implementation of the RESET method.
#
# X: n-by-p target matrix; columns represent variables and rows represent samples
# X.test: matrix (subset or transformation of X, e.g., projection on top PCs) that will be combined with the var.set variables 
#         to compute the reduced rank reconstruction via randomized SVD. 
#         Reconstruction error will be measured on the variables in X.test. 
#         If not specified, the X matrix will be used.
# center.X: Flag which controls whether the values in X are mean centered during execution of the algorithm.
#         If only X is specified and center.X=TRUE, then all columns in X will be centered. If both X and X.test are specified,
#         then centering performed on just the columns of X contained in the specified variable sets.
# scale.X: Flag which controls whether the values in X are scaled to have variance 1 during execution of the algorithm.
#         If only X is specified and scale.X=TRUE, then all columns in X will be scaled. If both X and X.test are specified,
#         then scaling is performed on just the columns of X contained in the specified variable sets.
# center.X.test: Flag which controls whether the values in X.test are mean centered during execution of the algorithm.
# scale.X.test: Flag which controls whether the values in X.test are scaled to have variance 1 during execution of the algorithm.
# var.sets: list of m variable sets, each element is a vector of indices of variables in the set
# k: Rank of reconstruction. Default to 2. Must be greater than or equal to the minimum variable set size.
# random.threshold: If specified, indicates the variable set size above which a randomized reduced-rank reconstruction is used. If
#  the variable set size is less or equal to random.threshold, then a non-random reconstruction is computed. Defaults to k 
#  and cannot be less than k.
# k.buff: Additional dimensions used in randomized reduced-rank construction algorithm. Defaults to 0. 
# q: Number of power iterations for randomized SVD. Defaults to 0.
# test.dist: Distribution for elements of random test matrix used in randomized algorithm. Must be one of "normal" or "uniform".
# norm.type: The type of norm to use for computing reconstruction error. Defaults to "2" for Euclidean/Frobenius norm. Other supported option
#    is "1" for L1 norm.
# weight.type: If specified, then the generated sample-level scores are weighted by a statistic computed as a function of set values for each sample.
#    Supported types include "mean" for mean of the set values, "mean.abs" for mean of the absolute set values, "stand.mean" for the
#    standardized mean of the set values and "stand.mean.abs" for the standardized mean of the absolute set values.
#    All of the supported weights will result in larger magnitude scores for samples with large values.
#
# Returns a list with two elements:
# 
# S: n-by-m matrix of sample-level variable set scores computed by RESET
# v: length m vector of the overall variable set scores computed by RESET
#
reset = function(X, X.test, center.X=TRUE, scale.X=FALSE, center.X.test=TRUE, scale.X.test=FALSE, var.sets, k=2, random.threshold, k.buff=0, q=0, test.dist="normal",
                 norm.type="2") {
  
  if (missing(X)) {
    stop("Matrix X not specified!")
  }
  if (missing(var.sets)) {
    stop("List of variable set indices var.sets must be specified!")
  }
  num.sets = length(var.sets)
  if (num.sets == 0) {
    stop("Must specify at least one variable set!")
  }
  min.set.size = min(sapply(var.sets, length))
  if (min.set.size == 0) {
    stop("All variable sets must contain at least one member!")
  }
  if (k < 1) {
    stop("k must be greater than 1!")
  }
  if ((k+k.buff) > min.set.size) {
    stop("k+k.buff is ", k+k.buff, ", which is larger than the minimum set size of ", min.set.size, "!")
  }
  if (missing(random.threshold)) {
    random.threshold = k
    message("Setting random.threhold to specified k of ", k)
  } else if (random.threshold < k) {
    stop("random.threshold cannot be than k!")
  }
  if (k.buff < 0) {
    stop("k.buff cannot be negative!")
  }
  if (q < 0) {
    stop("q cannot be negative!")
  }
  if (norm.type != "2" && norm.type != "1") {
    stop("unsupported norm.type ", norm.type)
  }
  # if (!missing(weight.type)) {
  #   if (!(weight.type %in% c("mean", "stand.mean", "mean.abs", "stand.mean.abs"))) {
  #     stop("weight.type must be one of 'mean', 'stand.mean', 'mean.abs' or 'stand.mean.abs'!")
  #   }
  # }
  
  p = ncol(X)
  n = nrow(X)
  m = length(var.sets)
  
  # If requested, centering X or X.test
  if (!missing(X.test)) {
    if (center.X.test || scale.X.test) {
      X.test = scale(X.test, scale=scale.X.test, center=center.X.test)
    }
  } else {
    if (center.X || scale.X) {
      X = scale(X, scale=scale.X, center=center.X)
    }
    X.test = X
  }

  # Compute the L1 or L2 norm of each row and the test matrix
  X.row.norms = apply(X.test, 1, function(x){vectorNorm(x,norm.type)})
  # Get correct matrix norm type
  m.norm = norm.type
  if (norm.type == "2") {
    m.norm = "F"
  }
  X.norm = norm(X.test, type=m.norm)
    
  # matrix that holds sample level scores
  S = matrix(0,nrow=n,ncol=m)
  
  # vector that holds overall set scores
  v = rep(0, m)
  
  message("Computing scores for collection of ", m, " sets")
  
  for (i in 1:m) {
    if (i %% 50 == 0) {
      message("Computing scores for set ", i)
    }
    var.set = var.sets[[i]]
    var.set.size = length(var.set)
    
    # Subset X for the variable set
    X.var.set = X[,var.set]
    
    # If X.test was specified and center.X is true, then
    # the subset needs to be centered
    if (!missing(X.test) && (center.X || scale.X)) {
      X.var.set = scale(X.var.set, center=center.X, scale=scale.X)
    }
    
    # Compute sample-level weights if needed
    # if (!missing(weight.type)) {
    #   var.vals = X.var.set
    #   if (endsWith(weight.type, ".abs")) {
    #     var.vals = abs(var.vals)
    #   }
    #   var.set.weights = apply(var.vals, 1, mean)
    #   if (startsWith(weight.type, "stand.")) {
    #     var.set.weights = scale(var.set.weights, center=FALSE, scale=TRUE)        
    #   }
    # }

    # Get rank k basis for column space of X.var.set
    if (var.set.size <= random.threshold) {
      # Variable set size is <= random.threshold,
      # so compute orthonormal basis directly using column-pivoted QR
      Q = qr.Q(qr(X.var.set, LAPACK=TRUE))
    } else {
      # Variable set size is > random.threshold, so use a randomized algorithm
      # implemented by randomColumnSpace to compute the rank k orthonormal basis
      Q = randomColumnSpace(X=X.var.set,k=(k+k.buff),q=q,test.dist=test.dist)
    }
    
    # Use the first k columns of Q as the rank k basis
    Q = Q[,1:k]
    # Compute reconstruction of X.test by projecting onto Q
    X.recon = Q %*% (t(Q) %*% X.test)
    
    # Compute reconstruction error
    X.diff = X.test - X.recon
    
    # Compute the row-level scores as log2 fold-change of the L1 or L2 norm of original row relative to the equivalent norm of the error. 
    # Larger values reflect improved reconstruction performance.
    X.diff.norms = apply(X.diff, 1, function(x){vectorNorm(x,norm.type)})
    S[,i] = log2(X.row.norms/X.diff.norms)
    
    # if weight.type is specified, then multiple scores by the appropriate weights
    # if (!missing(weight.type)) {
    #   S[,i] = S[,i] * var.set.weights
    # }
    
    # Compute overall score as the log2 fold-change of the L1 or Frobineus norm of the original test matrix relative to the relative to the equivalent norm 
    # of the error
    v[i] = log2(X.norm/norm(X.diff,type=m.norm))
  }
  
  S[which(is.nan(S))] = 0
  rownames(S) = rownames(X)
  colnames(S) = names(var.sets)
  names(v) = names(var.sets)
  
  results = list()
  results$S = S
  results$v = v
  
  return (results)
}

#
# Wrapper around the reset() function that uses the top PCs as X.test.
#
resetViaPCA = function(X, center=TRUE, scale=FALSE, num.pcs=2, pca.buff=2, pca.q=1, 
                       var.sets, k=2, random.threshold, k.buff=0, q=0, test.dist="normal", norm.type="2") {
  if (missing(X)) {
    stop("Matrix X not specified!")
  }
  if (num.pcs < 1) {
    stop("num.pcs cannot be less than 1!")
  }
  if (pca.buff < 0) {
    stop("pca.buff cannot be negative!")
  }
  if (pca.q < 0) {
    stop("pca.q cannot be negative!")
  }
  if (missing(random.threshold)) {
    random.threshold = k
  }

  # center/scale as desired
  if (center || scale) {
    X = scale(X, center=center, scale=scale)
  }
  
  # use randomized SVD algorithm to compute approximate projection onto the top PCs
  X.pcs = X %*% randomSVD(X=X,k=(num.pcs+pca.buff),q=pca.q,test.dist=test.dist)$v[,1:num.pcs]
  
  # call reset()
  # if the data was not centered before calling randomSVD, then center it when 
  # needed in reset()
  reset.results = reset(X=X,X.test=X.pcs,center.X=!center,center.X.test=!center,
                        scale.X=FALSE,scale.X.test=FALSE, var.sets=var.sets,
                        k=k, random.threshold=random.threshold, k.buff=k.buff,
                        q=q,test.dist=test.dist,norm.type=norm.type)
  return (reset.results)
}

#
# Utility function that creates a variable set list in the format required
# by rest() (i.e., list of variable indices) given the variable names and a
# list of variable names.
#
# Inputs:
#
# -var.names: Vector of variable names. This should correspond to the order of variables in
#             the target matrix X.
# -var.sets: List of m variable sets where each element in the list corresponds to
#            a set and the list element is a vector variable names. List names are variable set names.
# -min.size: Minimum set size after filtering out variable not in the var.names vector.
#            Sets whose post-filtering size is below this are removed from the final
#            collection list. Default is 1 and cannot be set to less than 1.
# -max.size: Maximum variable set size after filtering out variables not in the var.names vector.
#            Sets whose post-filtering size is above this are removed from the final
#            collection list. If not specified, no filtering is performed.
#
# Output:
#
#   Version of the input variable set collection list where variable names have been replaced by position indices,
#   variables not present in the var.names vector have been removed and sets failing the min/max size
#   constraints have been removed.
#
createVarSetCollection = function(var.names, var.sets, min.size=1, max.size) {
  
  # min.size must be at least 1
  if (min.size < 1) {
    stop("Invalid min.size value! Must be 1 or greater.")
  }
  # If max size is set, make sure it is not less than min size
  if (!missing(max.size)) {
    if (max.size < min.size) {
      stop("max.size cannot be less than min.size!")
    }
  }
  
  num.vars = length(var.names)
  if (num.vars < 1) {
    stop("var.names must contain at least one variable!")
  }
  
  num.sets = length(var.sets)
  if (num.sets < 1) {
    stop("var.sets must contain at least one set!")
  }
  
  set.names = names(var.sets)
  var.set.indices = list()
  for (i in 1:num.sets) {
    set.var.names = var.sets[[i]]
    # map names to indices
    set.indices = unlist(sapply(set.var.names, function(x){which(var.names == x)}))
    set.size = length(set.indices)
    if (set.size < min.size) {
      next
    }
    if (!missing(max.size)) {
      if (set.size > max.size) {
        next
      }
    }
    current.index = length(var.set.indices)+1
    var.set.indices[[current.index]] = set.indices
    names(var.set.indices)[current.index] = set.names[i]
  }
  
  return (var.set.indices)
}

#
# Utility function that computes the L1 or L2 norm of a vector
#
vectorNorm = function(x, norm.type="2") {
  if (norm.type == "1") {
    return (sum(abs(x)))
  } else if (norm.type == "2") {
    return (sqrt(sum(x^2)))
  } else {
    stop("Invalid norm.type ", norm.type)
  }
}