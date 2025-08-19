library(epiregulon.archr)
library(ArchR)
library(BiocParallel)
library(GenomicRanges)


archR_project_path <- "/gstore/project/epigen/PBMC/saved_project"
proj <- loadArchRProject(path = archR_project_path, showLogo = FALSE)

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


# Add reduced dimensionality and embedding

reducedDim(GeneExpressionMatrix, "LSI_Combined") <- getReducedDims(ArchRProj = proj,
                                                                   reducedDims = "LSI_Combined")

reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- getEmbedding(ArchRProj = proj,
                                                                  embedding = "UMAP_Combined")

GeneExpressionMatrix <- GeneExpressionMatrix[,!is.na(GeneExpressionMatrix$cell_type)]



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
peakMatrix <- peakMatrix[, colnames(GeneExpressionMatrix)]
names(assays(peakMatrix)) <- "counts"


# grl <- getTFMotifInfo(genome = "hg38", "atlas.tissue")
grl <- getTFMotifInfo(genome = "hg38")

rm(list=ls()[sapply(ls(), function(x) !x %in% c("grl", "peakMatrix", "GeneExpressionMatrix"))])
save.image(".RData")
