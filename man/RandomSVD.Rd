\name{randomSVD}
\alias{randomSVD}
\title{Implementation of a sparse powered randomized singular value decomposition.}
\description{
Computes an approximate rank \code{k} singular value decomposition (SVD) of an n-by-p input matrix \code{X} using a sparse randomized embedding with optional subspace power iterations. The \code{randomColumnSpace} method is used to generate an rank \code{k} approximation of the column space of \code{X}. This n-by-k approximation \code{Y} is then used to create a k-by-p projection \code{B} of \code{X} onto this rank \code{k} subspace via \code{B=Y^TX}. A non-random SVD is computed for \code{B} and this SVD solution is used to generate an approximate rank \code{k} SVD of \code{X}.
}
\usage{
    randomSVD(X, k=2, q=0, sparsity.structure=NULL, test.dist="normal")
}
\arguments{
\item{X}{An n-by-p target matrix.}
\item{k}{Target rank. Defaults to 2. See description in \code{\link{randomColumnSpace}}.}
\item{q}{Number of power iterations. Defaults to 0. See description in \code{\link{randomColumnSpace}}.}
\item{sparsity.structure}{Optional sparsity structure. See description in \code{\link{randomColumnSpace}}.}
\item{test.dist}{Type of random variable used to populate non-sparse elements of random test matrix.
See description in \code{\link{randomColumnSpace}}.}
}
\value{
List with the following elements:
    \itemize{
    \item \code{u} a matrix whose columns are the top k approximate left singular vectors of \code{X}.
    \item \code{d} a vector containing the top k approximate singular values of \code{X}.
    \item \code{v} a matrix whose columns are the top k approximate right singular vectors of \code{X}.
     }
}
\seealso{\code{\link{randomColumnSpace}}}
\keyword{file}
\examples{
  # Simulate a 100-by-100 matrix of random Poisson data
  X = matrix(rpois(10000, lambda=2), nrow=100)
  # Create a random sparsity structure for 100-by-5 random test matrix; half elements will be 0
  sparsity.structure = sample(1:500, 250, replace=TRUE)
  # Compute rank 5 SVD of X using a sparse test matrix
  svd.out = randomSVD(X,k=5,sparsity.structure=sparsity.structure)
  # Compute using a dense test matrix with U(0,1) RVs
  svd.out = randomSVD(X,k=5,test.dist="uniform")
}