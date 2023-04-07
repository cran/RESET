#
# RNLAUtils.R
#
# @author rob.frost@dartmouth.edu
#

#
# Implementation of a randomized algorithm that computes an approximate rank k
# orthonormal basis for the column space of the input method.
# Includes the option for structured sparsity.
#
# X: an n x p target matrix 
# k: target rank, defaults to 2
# q: number of power iterations, defaults to 0
# sparsity.structure: If specified, is a vector whose elements are
#    are the indices (in column-oriented format) of the non-sparse elements in the p x k random test matrix.
# test.dist: The distribution for elements of the test matrix. Default is "normal". Supported options:
#    normal: Elements are generated as N(0,1) RVs
#    constant: Elements are set to 1
#    uniform: Elements are set to U(0,1) RVs
#
randomColumnSpace = function(X, k=2, q=0, sparsity.structure=NULL, test.dist="normal") {
  
  if (missing(X)) {
    stop("Target matrix X must be specified!")
  }
  
  n = nrow(X)
  p = ncol(X)
  
  if (k > min(n,p)) {
    stop("k cannot be larger than the minimum dimension of X!")
  }
  
  if (k < 1) {
    stop("k must be >= 1!")
  }
  
  if (q < 0) {
    stop("q cannot be less than 0!?")
  }
  
  if (is.null(sparsity.structure)) {
    # Create a dense random test matrix
    if (test.dist == "normal") {
      O = matrix(rnorm(k*p), nrow=p) 
    } else if (test.dist == "uniform") {
      O = matrix(runif(k*p), nrow=p)
    } else {
      stop("Only normal or uniform are supported for dense test matrix!")
    }
  } else {
    # Create a sparse random test matrix
    O = Matrix::Matrix(0, nrow=p, ncol=k, sparse=TRUE)
    if (test.dist == "normal") {
      O[sparsity.structure] = rnorm(length(sparsity.structure))
    } else if (test.dist == "constant") {
      O[sparsity.structure] = 1
    } else if (test.dist == "uniform") {
      O[sparsity.structure] = runif(length(sparsity.structure))
    } else {
      stop("test.dist ", test.dist, " not supported!")
    }
  }
  
  # Create sketch matrix
  Y = X %*% O

  # Subspace power iterations via algorithm 2 from rSVD paper
  if (q > 0) {
    for (i in 1:q) {
      Y = qr.Q(qr(Y))
      Y = X %*% qr.Q(qr(t(X) %*% Y))
    }
  }
  
  # Get orthonormal basis via QR decomp of sketch
  Q = qr.Q(qr(Y))
  
  return (Q)
}

#
# Computes a rank-k randomized SVD. Uses randomColumnSpace() with identical arguments.
#
randomSVD = function(X, k=2, q=0, sparsity.structure=NULL, test.dist="normal") {
  
  Q = randomColumnSpace(X,k,q,sparsity.structure, test.dist=test.dist)
  B = t(Q) %*% X
  svd.B = svd(B)
  u = Q %*% svd.B$u
  
  results = list()
  results$d = svd.B$d
  results$u = u
  results$v = svd.B$v
  
  return (results)
}

#
# Computea a rank-k reconstruction using the specified SVD output.
#
recon = function(svd.out, k) {
  if (k == 1) {
    X.recon = (svd.out$u[,1] * svd.out$d[1]) %*% t(svd.out$v[,1])
  } else {
    X.recon = svd.out$u[,1:k] %*% diag(svd.out$d[1:k]) %*% t(svd.out$v[,1:k])
  }
  return (X.recon)
}

