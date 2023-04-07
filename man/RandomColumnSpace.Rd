\name{randomColumnSpace}
\alias{randomColumnSpace}
\title{Implementation of a sparse powered randomized algorithm for computing a basis for the column space of a matrix.}
\description{
Computes a rank \code{k} approximation of the column space of an n-by-p input matrix \code{X} using a sparse randomized embedding with optional subspace power iterations.
Specifically, a p-by-k random text matrix \code{O} is created where all elements are generated as independent N(0,1) or U(0,1) random variables except for elements designated as sparse via the specified \code{sparsity.structure}, which are set to 0. If a sparse structure is used, the non-zero elements can alternatively be set to the constant value of 1 for a non-random embedding. The test matrix is used to create an n-by-k sketch matrix \code{Y} as \code{Y=XO}. If \code{q>0}, subspace power iterations are performed on \code{Y} via algorithm 2 in the paper by Erichson, et al. associated with the \code{rsvd} R package (https://doi.org/10.18637/jss.v089.i11). The returned rank \code{k} column space approximation of \code{X} is then generated via a column-pivoted QR decomposition of \code{Y}.
}
\usage{
    randomColumnSpace(X, k=2, q=0, sparsity.structure=NULL, test.dist="normal")
}
\arguments{
\item{X}{An n-by-p target matrix.}
\item{k}{Target rank. Defaults to 2.}
\item{q}{Number of power iterations. Defaults to 0.}
\item{sparsity.structure}{Optional sparsity structure. Should be specified as a vector whose elements are 
      the indices (in column-oriented format) of the non-sparse elements in the p x k random test matrix \code{O}. If not specified, \code{O} will be dense.}
\item{test.dist}{Type of random variable used to populate non-sparse elements of random test matrix \code{O}. Must be either 'normal' for N(0,1) RVs, 'uniform' for U(0,1) RVs or 'constant' for the value of 1. Note that 'constant' should only be used if \code{sparsity.structure} is specified.}
}
\value{A n-by-k estimate of the column space of \code{X}.}
\seealso{\code{\link{randomSVD}}}
\keyword{file}
\examples{
  # Simulate a 100-by-100 matrix of random Poisson data
  X = matrix(rpois(10000, lambda=2), nrow=100)
  # Create a random sparsity structure for 100-by-5 random test matrix; half elements will be 0
  sparsity.structure = sample(1:500, 250, replace=TRUE)
  # Compute rank 5 estimate of column space of X using a sparse test matrix
  Q = randomColumnSpace(X,k=5,sparsity.structure=sparsity.structure)
  # Compute using a dense test matrix with U(0,1) RVs
  Q = randomColumnSpace(X,k=5,test.dist="uniform")
}