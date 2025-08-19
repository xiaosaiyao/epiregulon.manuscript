library(FigR)
library(dplyr)
library(FNN)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SummarizedExperiment)
library(ArchR)
library(epiregulon.archr)

archR_project_path <- "/gstore/project/epigen/PBMC/saved_project"
proj <- loadArchRProject(path = archR_project_path, showLogo = FALSE)
proj <- proj[proj$cellNames[!is.na(proj$cell_type)]]

# load gene expression matrix
GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "GeneExpressionMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)

GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts")
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

#load peakMatrix
peakMatrix <- getMatrixFromProject(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)
peakMatrix <- as(peakMatrix, "SingleCellExperiment")
peakMatrix.sce <- peakMatrix[, colnames(GeneExpressionMatrix)]
geneExprMatrix <- assay(GeneExpressionMatrix)

# Remove genes with zero expression across all cells
geneExprMatrix <- geneExprMatrix[Matrix::rowSums(geneExprMatrix)!=0,]

# choose sequences present in the reference genome
selected_sequences <- seqnames(rowRanges(peakMatrix.sce)) %in% names(BSgenome.Hsapiens.UCSC.hg38)
peakMatrix.sce <- peakMatrix.sce[selected_sequences,]

# extract LSI
LSI_ATAC <- getReducedDims(proj, reducedDims = "LSI_ATAC")

rm(list=ls()[sapply(ls(), function(x) !x %in% c("LSI_ATAC", "peakMatrix.sce", "geneExprMatrix"))])
save.image(".RData")
