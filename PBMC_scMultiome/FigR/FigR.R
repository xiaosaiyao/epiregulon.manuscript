library(FigR)
library(dplyr)
library(FNN)
library(doParallel)
library(SummarizedExperiment)

load(".RData")
args <- commandArgs(trailingOnly = TRUE)
cores <- as.numeric(args[1])

# Cell kNN for later use in results smoothing
cellkNN <- get.knn(LSI_ATAC , k = 30)$nn.index

rownames(cellkNN) <- rownames(LSI_ATAC)

# use cells with both peak and gene expression data
geneExprMatrix <- geneExprMatrix[,colnames(geneExprMatrix) %in% colnames(peakMatrix.sce)]

# log normalization , following Seurat::NormalizeData

geneExprMatrix<-log(geneExprMatrix+1)

# Calculate peak to gene correlation
# default window around each gene: 50000
names(assays(peakMatrix.sce)) <- "counts"
cisCorr <- runGenePeakcorr(ATAC.se = peakMatrix.sce,
                           RNAmat = geneExprMatrix,
                           genome = "hg38",
                           nCores = cores,
                           p.cut = NULL,
                           n_bg = 100,
                           windowPadSize = 2.5e5)


# filter siginificant peak-gene associations

cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)

dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                       cutoff =10,
                       labelTop = 20,
                       returnGeneList = TRUE,
                       force = 2)


# calculate domain of regulatory chromatin (DORC) accessibility scores for each gene
dorcMat <- getDORCScores(ATAC.se = peakMatrix.sce,
                         dorcTab = cisCorr.filt,
                         geneList = dorcGenes,
                         nCores = cores)


# smooth dorc socres using cell KNN
dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30], mat = dorcMat, nCores = cores)

# smooth gene expression data
RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30], mat = geneExprMatrix, nCores = cores)

FigR_GRN <- runFigRGRN(ATAC.se = peakMatrix.sce,
                       dorcTab = cisCorr.filt,
                       genome = "hg38",
                       dorcMat = dorcMat.s,
                       rnaMat = RNAmat.s,
                       nCores = cores)

# saveRDS(FigR_GRN, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/FigR_GRN.rds")
