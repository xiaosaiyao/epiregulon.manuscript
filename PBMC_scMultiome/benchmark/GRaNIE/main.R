library(GRaNIE)
load(".RData")

args <- commandArgs(trailingOnly = TRUE)
cores <- as.numeric(args[1])

GRN <- initializeGRN(outputFolder = ".", genomeAssembly = "hg38")

# default normalization of peak data doesn't work since there is a zero value in each row making it impossible
# to calculate geometrical means
GRN <- addData(GRN, counts_peaks = peakMatrix_df, normalization_peaks = "limma_quantile",
               counts_rna = GeneExpressionMatrix_data, normalization_rna = "none", idColumn_peaks = "peakID",
               idColumn_RNA = "geneID", forceRerun = TRUE, sampleMetadata = sampleMetadata)


GRN <- addTFBS(GRN, motifFolder = "/gstore/project/epigen/PBMC/analysis/GRaNIE/FIMO_HOCOMOCOv11", TFs = "all", filesTFBSPattern = "_TFBS",
               fileEnding = ".bed.gz", forceRerun = TRUE)

GRN <- overlapPeaksAndTFBS(GRN, forceRerun = TRUE)

GRN <- addConnections_TF_peak(GRN, plotDiagnosticPlots = FALSE, connectionTypes = c("expression"),
                              corMethod = "pearson", forceRerun = TRUE, outputFolder = ".")



GRN <- addConnections_peak_gene(GRN, corMethod = "pearson",
                                TADs= NULL, plotDiagnosticPlots = FALSE,
                                plotGeneTypes = list(c("all")), forceRerun = TRUE,
                                outputFolder = ".", nCores = cores)



GRN = filterGRNAndConnectGenes(GRN, TF_peak.fdr.threshold = 0.2,
                               peak_gene.fdr.threshold = 0.2,
                               peak_gene.fdr.method = "BH",
                               gene.types = c("protein_coding", "lincRNA"),
                               allowMissingTFs = FALSE,
                               allowMissingGenes = FALSE)



GRN_connections.all <- getGRNConnections(GRN, type = "all.filtered",
                                         include_geneMetadata = TRUE,
                                         include_TF_gene_correlations = FALSE)
