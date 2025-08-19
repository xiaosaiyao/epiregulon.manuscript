library(ArchR)
library(epiregulon.archr)

archR_project_path <- "/gstore/project/epigen/PBMC/saved_project"
proj <- loadArchRProject(path = archR_project_path, showLogo = TRUE)
proj <- proj[proj$cellNames[!is.na(proj$cell_type)]]

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

# Add embeddings
reducedDim(GeneExpressionMatrix, "LSI_Combined") <- getReducedDims(ArchRProj = proj,
                                                                            reducedDims = "LSI_Combined",
                                                                            returnMatrix = TRUE)[colnames(GeneExpressionMatrix), ]


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

library(GenomicRanges)
peakMatrix_assay <- as.matrix(SummarizedExperiment::assay(peakMatrix))
rownames(peakMatrix_assay) <- paste0(as.character(seqnames(rowRanges(peakMatrix))),
                                     ":", as.character(ranges(rowRanges(peakMatrix))))
colnames(peakMatrix_assay) <- gsub(".*#", "", colnames(peakMatrix_assay))
colnames(peakMatrix_assay) <- gsub("-", "\\.", colnames(peakMatrix_assay))
paths_to_matrices <- "/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/peak_matrix.tsv"
write.table(peakMatrix_assay, file = paths_to_matrices,  col.names = TRUE, row.names = TRUE, sep = "\t")


barcode_tab <- data.frame(barcode = colnames(peakMatrix))
barcode_tab$sample_id <- gsub("(.*)(#.*)", "\\1",barcode_tab$barcode)
barcode_tab$barcode <- gsub(".*#", "",barcode_tab$barcode)
barcode_tab$barcode <- gsub("-", "\\.",barcode_tab$barcode)
barcode_tab[["cell_type"]] <- colData(peakMatrix)$cell_type
write.table(barcode_tab, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/barcode_tab_PBMC.csv", row.names = FALSE, sep = ",")

peakRegions <- rowRanges(peakMatrix)
peakRegions <- resize(peakRegions, fix="center", width = 1)
rtracklayer::export.bed(peakRegions, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/peakRegions.bed")

# download motif file from https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl

# https://resources.aertslab.org/cistarget/
