\name{resetForSeurat}
\alias{resetForSeurat}
\title{RESET wrapper for scRNA-seq data processed using the Seurat framework}
\description{
Executes the Reconstruction Set Test (RESET) method (\code{\link{reset}}) on 
    normalized scRNA-seq data stored in a Seurat object.
}
\usage{
resetForSeurat(seurat.data, num.pcs, gene.set.collection, k=10, random.threshold, 
               k.buff=0, q=0, test.dist="normal", norm.type="2")
}
\arguments{
\item{seurat.data}{The Seurat object that holds the scRNA-seq data. Assumes PCA have already been performed 
        on a centered and scaled version of the normalized counts.}
\item{num.pcs}{Number of PCs to use in the RESET method. If not specified, will use all computed PCs. 
        The projection of the scRNA-seq data onto these PCs will be used as the \code{X.test} matrix for the \code{\link{reset}} call.}
\item{gene.set.collection}{List of m gene sets for which scores are computed.
    Each element in the list corresponds to a gene set and the list element is a vector
    of indices for the genes in the set. The index value is defined relative to the
    order of genes in the gene.exprs matrix. Gene set names should be specified as list names.
    See \code{\link{createVarSetCollection}} for utility function that can be used to help generate this list of indices.}
\item{k}{See description in \code{\link{reset}}}    
\item{random.threshold}{See description in \code{\link{reset}}}
\item{k.buff}{See description in \code{\link{reset}}}    
\item{q}{See description in \code{\link{reset}}}    
\item{test.dist}{See description in \code{\link{reset}}}
\item{norm.type}{See description in \code{\link{reset}}}
}
\value{
An updated Seurat object with the following modifications:
    \itemize{
    \item RESET sample-level gene set score matrix S stored in the "data" slot of a new "RESET" Assay object.
    \item RESET overall gene set scores stored in the feature metadata column "RESET".
     }
}
\seealso{\code{\link{reset}},\code{\link{createVarSetCollection}}} 
\keyword{file}
\examples{
    # Only run example code if Seurat package is available
    if (requireNamespace("Seurat", quietly=TRUE) & requireNamespace("SeuratObject", quietly=TRUE)) {
        # Define a collection with three gene sets
        collection=list(set1=1:10, set2=11:20, set3=21:30)
        # Execute on the pbmc_small scRNA-seq data set included with SeuratObject
        # See vignettes for more detailed Seurat examples
        reset.out = resetForSeurat(seurat.data=SeuratObject::pbmc_small,
            num.pcs=5,
            gene.set.collection=collection)
        # Retrieve the scores for the first 10 cells
        reset.out@assays$RESET[,1:10]
        # Retrieve the overall scores
        reset.out@assays$RESET@meta.features
    }
}