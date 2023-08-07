\name{reset}
\alias{reset}
\title{Reconstruction Set Test (RESET)}
\description{
Implementation of the Reconstruction Set Test (RESET) method, which transforms an n-by-p input matrix \code{X} into an n-by-m matrix of sample-level variable set scores and a length m vector of overall variable set scores. Execution of RESET involves the following sequence of steps:
    \itemize{
      \item If \code{center.X=TRUE}, mean center the columns of \code{X}. If \code{X.test} is specified, the centering is instead 
            performed on just the columns of \code{X} corresponding to each variable set. 
            See documentation for the \code{X} and \code{center.X} parameters for more details.
      \item If \code{scale.X=TRUE}, scale the columns of \code{X} to have variance 1. If \code{X.test} is specified, the scaling is instead 
            performed on just the columns of \code{X} corresponding to each variable set. 
            See documentation for the \code{X} and \code{scale.X} parameters for more details.            
      \item If \code{center.X.test=TRUE}, mean center the columns of \code{X.test}. 
            See documentation for the \code{X.test} and \code{center.X.test} parameters for more details.
      \item If \code{scale.X.test=TRUE}, scale the columns of \code{X.test}. 
            See documentation for the \code{X.test} and \code{scale.X.test} parameters for more details.            
      \item Set the reconstruction target matrix \code{T} to \code{X} or, if \code{X.test} is specified, to \code{X.test}.
      \item Compute the norm of \code{T} and norm of each row of \code{T}. By default, these are the Frobenius and Euclidean norms respectively. 
      \item For each set in \code{var.sets}, sample-level and matrix level scores are generated as follows:
         \itemize{
          \item Create a subset of \code{X} called \code{X.var.set} that only includes the columns of \code{X} correponding to the variables
                in the set.
          \item Compute a rank \code{k} orthonormal basis \code{Q} for the column space of \code{X.var.set}. 
                If the size of the set is less then or equal to \code{random.threshold}, then this is computed as the top \code{k} columns
                of the \code{Q} matrix from a column-pivoted QR decomposition of \code{X.var.set}, otherwise, it is approximated using
                a randomized algorithm implemented by \code{\link{randomColumnSpace}}.
          \item The reduced rank reconstruction of \code{T} is then created as \code{Q Q^T T}.
          \item The original \code{T} is subtracted from the reconstruction to represent the reconstruction error and the appropriate norm
                is computed on each row and the entire error matrix.
          \item The overall score is the log2 ratio of the norm of the original \code{T} to the norm of the reconstruction error matrix.
          \item The score for each sample is the log2 ratio of the norm of the corresponding row of the original \code{T} to the norm of the same row of the reconstruction error matrix.
          \item If \code{per.var=TRUE}, then the overall and sample-level scores are divided by the variable set size.
         }
      }
}
\usage{
reset(X, X.test, center.X=TRUE, scale.X=FALSE, center.X.test=TRUE, scale.X.test=FALSE, 
      var.sets, k=2, random.threshold, k.buff=0, q=0, test.dist="normal", norm.type="2",
      per.var=FALSE)
}
\arguments{
\item{X}{The n-by-p target matrix; columns represent variables and rows represent samples.}
\item{X.test}{Matrix that will be combined with the \code{var.set} variables to compute the reduced rank reconstruction.
      This is typically a subset or transformation of \code{X}, e.g., projection on top PCs. 
      Reconstruction error will be measured on the variables in \code{X.test}. 
      If not specified, the entire \code{X} matrix will be used for calculating reconstruction error.}
\item{center.X}{Flag which controls whether the values in \code{X} are mean centered during execution of the algorithm.
      If only \code{X} is specified and \code{center.X=TRUE}, then all columns in \code{X} will be centered. If both \code{X} and \code{X.test} are specified,
      then centering is performed on just the columns of \code{X} contained in the specified variable sets.
      Mean centering is especially important for accurate performance when \code{X.test} is specified as a reduced rank representation of the \code{X}, e.g, 
      as the projection of \code{X} onto the top principal components. However, mean centering the entire matrix \code{X} can have a dramatic impact on memory 
      requirements if \code{X} is a large sparse matrix. 
      In this case, a non-centered \code{X} and appropriate \code{X.test} (e.g., project onto top PCs of \code{X}) can be provided and mean centering performed
      on just the needed variables during execution of RESET. This "just-in-time" centering is enabled by setting \code{center.X=TRUE} and 
      providing both \code{X} and \code{X.test}. If \code{X} has already been mean-centered (and \code{X.test} is a subset of this mean-centered matrix or
      computed using this mean-centered matrix), then center should be specified as FALSE.}
\item{scale.X}{Flag which controls whether the values in \code{X} are are scaled to have variance 1 during execution of the algorithm.
      Defaults to false.
      If only \code{X} is specified and \code{scale.X=TRUE}, then all columns in \code{X} will be scaled. If both \code{X} and \code{X.test} are specified,
      then scaling is performed on just the columns of \code{X} contained in the specified variable sets.}
\item{center.X.test}{Flag which controls whether the values in \code{X.test}, if specified, are mean centered during execution of the algorithm.
      Centering should be performed consistently for \code{X} and \code{X.test}, i.e., if \code{center.X} is true or \code{X} was previously centered, 
      then \code{center.X.test} should te true unless \code{X.test} previously centered or generated from a centered \code{X}.}
\item{scale.X.test}{Flag which controls whether the values in \code{X.test}, if specified, are scaled to have variance 1 during execution of the algorithm.
      Similar to centering, scaling should be performed consistently for \code{X} and \code{X.test}, i.e., if \code{scale.X} is true or \code{X}
      was previously scaled then \code{scale.X.test} should te true unless \code{X.test} previously scaled or generated from a scaled \code{X}.}
\item{var.sets}{List of m variable sets, each element is a vector of indices of variables in the set that correspond to columns in \code{X}. 
      If variable set information is instead available in terms of variable names, the appropriate format can be generated 
      using \code{\link{createVarSetCollection}}.}
\item{k}{Rank of reconstruction. Default to 2. Cannot be larger than the minimum variable set size.}
\item{random.threshold}{If specified, indicates the variable set size above which a randomized reduced-rank reconstruction is used. 
      If the variable set size is less or equal to random.threshold, then a non-random reconstruction is computed. Defaults to k 
      and cannot be less than k.}
\item{k.buff}{Additional dimensions used in randomized reduced-rank construction algorithm. Defaults to 0. 
      Values above 0 can improve the accuracy of the
      randomized reconstruction at the expense of additional computational complexity. If \code{k.buff}=0, then the reduced rank 
      reconstruction can be generated directly from the output of \code{\link{randomColumnSpace}}, otherwise, a reduced rank SVD must also be computed
      with the reconstruction based on the top \code{k} components.}
\item{q}{Number of power iterations for randomized SVD (see \code{\link{randomSVD}}). 
      Defaults to 0. Although power iterations can improve randomized SVD
      performance in general, it can decrease the sensitivity of the RESET method to detect mean or covariance differences.}
\item{test.dist}{Distribution for non-zero elements of random test matrix used in randomized SVD algorithm. See description for \code{test.dist} parameter
      of \code{\link{randomSVD}} method.}
\item{norm.type}{The type of norm to use for computing reconstruction error. Defaults to "2" for Euclidean/Frobenius norm. Other supported option
      is "1" for L1 norm.}
\item{per.var}{If true, the computed scores for each variable set are divided by the scaled variable set size to generate per-variable scores.
      Variable set size scaling is performed by dividing all sizes by the mean size (this will generate per-variable scores of approximately the
      same magnitude as the non-per-variable scores).}
}
\value{
A list with the following elements:
    \itemize{
      \item \code{S} an n-by-m matrix of sample-level variable set scores.
      \item \code{v} a length m vector of overall variable set scores.
     }
}
\seealso{\code{\link{createVarSetCollection}},\code{\link{randomColumnSpace}}} 
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

  # Execute RESET using non-randomized basis computation
  reset(X, var.sets=var.sets, k=2, random.threshold=10)

  # Execute RESET with randomized basis computation
  # (random.threshold will default to k value which is less
  # than the size of all variable sets)
  reset(X, var.sets=var.sets, k=2, k.buff=2)

  # Execute RESET with non-zero k.buff
  reset(X, var.sets=var.sets, k=2, k.buff=2)
  
  # Execute RESET with non-zero q
  reset(X, var.sets=var.sets, k=2, q=1)

  # Execute RESET with L1 vs L2 norm
  reset(X, var.sets=var.sets, k=2, norm.type="1")

  # Project the X matrix onto the first 5 PCs and use that as X.test
  # Scale X before calling prcomp() so that no centering or scaling
  # is needed within reset()
  X = scale(X)
  X.test = prcomp(X,center=FALSE,scale=FALSE,retx=TRUE)$x[,1:5]
  reset(X, X.test=X.test, center.X=FALSE, scale.X=FALSE, 
    center.X.test=FALSE, scale.X.test=FALSE, var.sets=var.sets, k=2)
}
