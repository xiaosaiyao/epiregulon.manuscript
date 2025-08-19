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

saveRDS(GeneExpressionMatrix, "/gstore/project/epigen/PBMC/OUTPUT/GeneExpressionMatrix.rds")

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
saveRDS(peakMatrix, "/gstore/project/epigen/PBMC/OUTPUT/peakMatrix.rds")


# detach("package:epiregulon.archr", unload= TRUE)
# detach("package:epiregulon", unload= TRUE)
# library(epiregulon)
grl <- getTFMotifInfo(genome = "hg38")
grl <- grl[unlist(lapply(grl, length))>=1000]


set.seed(1010, kind ="L'Ecuyer-CMRG")

# find peak to gene links
p2g <- calculateP2G(peakMatrix = peakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    reducedDim = reducedDim(GeneExpressionMatrix, "LSI_Combined"),
                    peak_assay = "counts",
                    exp_assay = "normalizedCounts",
                    cor_cutoff = 0.5
)


# Construct regulons
overlap <- addTFMotifInfo(grl = grl,
                          p2g = p2g,
                          peakMatrix = peakMatrix)

regulon_df_full <- getRegulon(p2g, overlap, aggregate = FALSE)

# prune network
pruned.regulon <- pruneRegulon(regulon = regulon_df_full,
                               expMatrix = GeneExpressionMatrix,
                               exp_assay = "normalizedCounts",
                               peakMatrix = peakMatrix,
                               peak_assay = "counts",
                               prune_value = "pval"
)

selected_tfs <- c("TCF7", "GATA3", "BCL11B", "SPI1", "CEBPA", "EBF1", "PAX5", "POU2F2", "POU2AF1", "EOMES", "RUNX3", "TBX21", "SPIB", "IRF8", "STAT6","ELK1", "GATA3", "JUN", "NFATC3", "NFKB1", "STAT3", "MAF",
                  "RUNX1", "TCF3", "STAT3", "BCL6", "FOXP3", "PRDM1", "KLF4", "GATA2", "CEBPB", "IKZF1", "NFIL3", "TCF4")

saveRDS(pruned.regulon, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/pruned.regulon_merged_no_clusters.rds")
saveRDS(pruned.regulon[pruned.regulon$tf %in% selected_tfs,], "/gstore/project/epigen/benchmark/PBMC/OUTPUT/pruned.regulon_merged_trimmed_no_clusters.rds")

regulon.w <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                        peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                        peak_assay="counts")

saveRDS(regulon.w[regulon.w$tf %in% selected_tfs,], "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_merged_trimmed_no_clusters.rds")
pruned.regulon.trimmed <- pruned.regulon[pruned.regulon$tf %in% selected_tfs,]

regulon.w2 <- addWeights(pruned.regulon.trimmed, expMatrix = GeneExpressionMatrix,
                         peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                         peak_assay="counts", clusters = GeneExpressionMatrix$cell_type,
                         method = "MI")
saveRDS(regulon.w2, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_MI_trimmed_merged_no_clusters.rds")

regulon.w3 <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                         peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                         peak_assay="counts", clusters = GeneExpressionMatrix$cell_type,
                         method = "corr")
saveRDS(regulon.w3[regulon.w3$tf %in% selected_tfs,], "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_merged_trimmed_no_clusters.rds")
saveRDS(regulon.w3, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_merged_no_clusters.rds")
