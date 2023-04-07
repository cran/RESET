\name{createVarSetCollection}
\alias{createVarSetCollection}
\title{Utility function to help create a variable set collection list object}
\description{
Utility function that creates a variable set list in the format required by \link{reset} (i.e., list of variable indices) given the variable names and an ordered vector of variable names.
}
\usage{
    createVarSetCollection(var.names, var.sets, min.size=1, max.size)
}
\arguments{
\item{var.names}{Vector of variable names. This should correspond to the order of variables as represented by the columns of the target matrix X.}
\item{var.sets}{List of m variable sets where each element in the list corresponds to a set and the list element is a vector variable names. List names are variable set names.}
\item{min.size}{Minimum set size after filtering out variable not in the var.names vector. Sets whose post-filtering size is below this are removed from the final collection list. Default is 1 and cannot be set to less than 1.}
\item{max.size}{Maximum variable set size after filtering out variables not in the var.names vector. Sets whose post-filtering size is above this are removed from the final collection list. If not specified, no filtering is performed.}
}
\value{
Version of the input variable set collection list where variable names have been replaced by position indices, variables not present in the var.names vector have been removed and sets failing the min/max size constraints have been removed.
}
\seealso{\code{\link{reset}}}
\keyword{file}
\examples{
    # Create a collection with two sets defined over 3 variables
    createVarSetCollection(var.names=c("A", "B", "C"),
        var.sets = list(set1=c("A", "B"), set2=c("B", "C")),
        min.size=2, max.size=3)                    
}
