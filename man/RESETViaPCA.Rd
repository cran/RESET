\name{resetViaPCA}
\alias{resetViaPCA}
\title{Reconstruction Set Test (RESET) via PCA}
\description{
Wrapper around the \code{\link{reset}} method that uses the projection of \code{X} onto the top \code{num.pcs} principal components as \code{X.test}.
This PC projection is computed using a randomized reduced rank SVD as implemented by \code{\link{randomSVD}}.
}
\usage{
resetViaPCA(X, center=TRUE, scale=FALSE, num.pcs=2, pca.buff=2, pca.q=1, var.sets, k=2, 
      random.threshold, k.buff=0, q=0, test.dist="normal", norm.type="2", per.var=FALSE)
}
\arguments{
\item{X}{See description in \code{\link{reset}}}
\item{center}{Flag which controls whether the values in \code{X} are mean centered.
     Note that if \code{center} is set to true, centering is performed on the entire \code{X} matrix prior to 
     calling \code{\link{randomSVD}}, which may have significant performance and memory implications of \code{X} is large and/or sparse.
     If \code{center} is false, then \code{X} will be projected onto the uncentered PCs and \code{center.X} and \code{center.X.test} will be set to
     TRUE in the call to \code{\link{reset}}.} 
\item{scale}{Flag which controls whether the values in \code{X} are are scaled to have variance 1.
     Note that if \code{scale} is set to true, scaling is performed on the entire \code{X} matrix prior to 
     calling \code{\link{randomSVD}}, which may have significant performance and memory implications of \code{X} is large and/or sparse.}
\item{num.pcs}{Number of principal components used for computing the projection of \code{X}.}
\item{pca.buff}{Number of extra dimensions used when calling \code{\link{randomSVD}} to compute the PCs. See \code{k.buff} parameter for
      \code{\link{randomSVD}} function.}
\item{pca.q}{Number of power iterations used when calling \code{\link{randomSVD}} to compute the PCs. See \code{q} parameter for
      \code{\link{randomSVD}} function.}      
\item{var.sets}{See description in \code{\link{reset}}}
\item{k}{See description in \code{\link{reset}}}
\item{random.threshold}{See description in \code{\link{reset}}}
\item{k.buff}{See description in \code{\link{reset}}}
\item{q}{See description in \code{\link{reset}}}
\item{test.dist}{See description in \code{\link{reset}}}
\item{norm.type}{See description in \code{\link{reset}}}
\item{per.var}{See description in \code{\link{reset}}}
}
\value{
A list with the following elements:
    \itemize{
      \item \code{S} an n-by-m matrix of sample-level variable set scores.
      \item \code{v} a length m vector of overall variable set scores.
     }
}
\seealso{\code{\link{reset}},\code{\link{createVarSetCollection}},\code{\link{randomColumnSpace}}} 
\keyword{file}
\examples{
  # Create a collection of 5 variable sets each of size 10
  var.sets = list(set1=1:10, 
                  set2=11:20,
                  set3=21:30,
                  set4=31:40,
                  set5=41:50)                  

  # Simulate a 100-by-100 matrix of random Poisson data
  X = matrix(rpois(10000, lambda=1), nrow=100)

  # Inflate first 10 rows for first 10 variables, i.e., the first
  # 10 samples should have elevated scores for the first variable set
  X[1:10,1:10] = rpois(100, lambda=5)

  # Execute RESET when reconstruction measured on top 10 PCs 
  # with mean centering performed before computing PCs
  resetViaPCA(X, num.pcs=10, var.sets=var.sets, k=2, random.threshold=10)
  
  # Execute RESET when reconstruction measured on top 10 
  # uncentered PCs with centering performed as needed inside reset()
  resetViaPCA(X, center=FALSE, num.pcs=10, var.sets=var.sets, k=2, random.threshold=10)
}
