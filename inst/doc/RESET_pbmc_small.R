### R code from vignette source 'RESET_pbmc_small.Rnw'

###################################################
### code chunk number 1: RESET_pbmc_small.Rnw:16-17
###################################################
library(RESET)


###################################################
### code chunk number 2: RESET_pbmc_small.Rnw:24-32
###################################################
if (requireNamespace("Seurat", quietly=TRUE) && requireNamespace("SeuratObject", quietly=TRUE)) {
	SeuratObject::pbmc_small
	gene.names = rownames(SeuratObject::pbmc_small)
	gene.names[1:5]
	Seurat::VariableFeatures(SeuratObject::pbmc_small)[1:5]
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 3: RESET_pbmc_small.Rnw:39-54
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	gene.set.id.list = list()
	# Create set with top 5 variable genes
	gene.set.id.list[[1]] = c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3")
	names(gene.set.id.list)[1] = "VarGenes"
	# Create set with 5 random genes
	gene.set.id.list[[2]] = c("TREML1", "CD79B", "LRRC25", "GPX1", "CFD")
	names(gene.set.id.list)[2] = "RandomGenes"
	print(gene.set.id.list)
	# Create the list of gene indices required by resetForSeurat()
	gene.set.collection = createVarSetCollection(var.names=gene.names,
		var.sets=gene.set.id.list)
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 4: RESET_pbmc_small.Rnw:63-71
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	pbmc.reset = resetForSeurat(seurat.data=SeuratObject::pbmc_small,
	                          num.pcs=5,
	                          gene.set.collection=gene.set.collection,
	                          k=5)
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 5: RESET_pbmc_small.Rnw:76-84
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
  # Display RESET scores for first 10 cells
	print(pbmc.reset@assays$RESET[,1:10])
  # Display overall RESET scores
	pbmc.reset@assays$RESET@meta.features
} else {
	message("Seurat package not available! Not executing associated vignette content.")
}


###################################################
### code chunk number 6: RESET_pbmc_small.Rnw:91-103
###################################################
if (requireNamespace("Seurat", quietly=TRUE)) {
	Seurat::DefaultAssay(object = pbmc.reset) = "RESET"
	Seurat::FeaturePlot(pbmc.reset, reduction="tsne", features="VarGenes")
} else {
	message("Seurat package not available! Not executing associated vignette content.")
	oldpar = par(mar = c(0,0,0,0))
	plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
	text(x = 0.5, y = 0.5,paste("Seurat package not available!\n",
					 "FeaturePlot not generated."),
	cex = 1.6, col = "black")
  par(oldpar)
}


