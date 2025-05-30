library(scran)
library(scater)
library(epiregulon.archr)
library(epiregulon.extra)
library(ArchR)
library(BiocParallel)
library(eulerr)
library(ggpubr)

#checkdata
archR_project_path <- "OUTPUT/ArchRProject/"
proj.all <- loadArchRProject(path = archR_project_path, showLogo = TRUE)


celllines <- c("LNCaP", "VCaP", "MDA")

###############################
# load gene expression matrix
GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "GeneExpressionMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)



GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts")
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
GeneExpressionMatrix$TEST_ARTICLE <- factor(as.character(GeneExpressionMatrix$TEST_ARTICLE),
                                            levels = c("DMSO", "Enza","ARV110", "A9690") )




# Add embeddings
reducedDim(GeneExpressionMatrix, "IterativeLSI_Combined") <- getReducedDims(ArchRProj = proj.all,
                                                                            reducedDims = "IterativeLSI_Combined",
                                                                            returnMatrix = TRUE)[colnames(GeneExpressionMatrix), ]

#load peakMatrix
peakMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = 1,
    logFile = createLogFile("getMatrixFromProject")
)

peakMatrix <- as(peakMatrix, "SingleCellExperiment")
peakMatrix <- peakMatrix[, colnames(GeneExpressionMatrix)]
names(assays(peakMatrix)) <- "counts"


# get TF binding
#ExperimentHub::setExperimentHubOption("CACHE", "/gstore/scratch/u/yaox19/.cache/ExperimentHub")

for (cell in c("LNCaP","VCaP", "MDA")){
    message(cell)
    selected <- which(proj.all$Cell == cell )
    proj <- proj.all[selected,]
    GeneExpressionMatrix.select <- GeneExpressionMatrix[, selected]
    peakMatrix.select <- peakMatrix[, selected]


    # peak2gene links

    set.seed(1010, kind = "Mersenne-Twister")


    # find cell line specific peaks
    cell_peaks <- lapply(unique(proj$hash_assignment2),
                         function(cluster) {readRDS(file.path(archR_project_path, "PeakCalls", paste0(make.names(cluster), "-reproduciblePeaks.gr.rds")))})
    cell_peaks <- GenomicRanges::reduce(do.call(c, GRangesList(cell_peaks)))
    peaks_overlaps <- findOverlaps(cell_peaks, rowRanges(peakMatrix.select ))
    peakMatrix.cell <- peakMatrix.select[unique(subjectHits(peaks_overlaps)),]

    grl <- getTFMotifInfo(genome = "hg38", mode = "motif", peaks = rowRanges(peakMatrix.cell) )


    # find peak to gene links
    p2g <- calculateP2G(peakMatrix = peakMatrix.cell,
                        expMatrix = GeneExpressionMatrix.select ,
                        reducedDim = reducedDim(GeneExpressionMatrix.select ),
                        peak_assay = "counts",
                        exp_assay = "normalizedCounts",
                        cor_cutoff = 0.5
    )


    # Construct regulons
    overlap <- addTFMotifInfo(grl = grl,
                              p2g = p2g,
                              peakMatrix = peakMatrix.cell)



    regulon_df_full <- getRegulon(p2g, overlap, aggregate = FALSE)




    # prune network
    pruned.regulon <- pruneRegulon(regulon = regulon_df_full,
                                   expMatrix = GeneExpressionMatrix.select ,
                                   exp_assay = "normalizedCounts",
                                   peakMatrix = peakMatrix.cell ,
                                   peak_assay = "counts",
                                   prune_value = "pval",
                                   clusters = GeneExpressionMatrix.select$TEST_ARTICLE
    )


    dim(pruned.regulon)


    # Add weights
    regulon.w <- addWeights(regulon = pruned.regulon,
                            expMatrix = GeneExpressionMatrix.select,
                            exp_assay = "normalizedCounts",
                            peakMatrix = peakMatrix.cell,
                            peak_assay = "counts",
                            method = "wilcoxon",
                            clusters = GeneExpressionMatrix.select$TEST_ARTICLE
    )

    # calculate activity
    score.combine <- calculateActivity(expMatrix = GeneExpressionMatrix.select,
                                       regulon = regulon.w,
                                       method = "weightedMean",
                                       exp_assay = "normalizedCounts",
                                       mode = "weight",
                                       FUN = "mean")


    saveRDS(regulon.w, paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.", cell, ".motif.rds"))
    saveRDS(score.combine, paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cell, ".motif.rds"))

    #########
    pdf(paste0("OUTPUT/ArchRProject/Epiregulon/activity.separate.", cell,".motif.pdf"), width = 9, height = 6)


    plotviolin <- plotActivityViolin(score.combine,
                                     tf = c( "AR","FOXA1", "SMARCA4"),
                                     clusters = GeneExpressionMatrix.select$TEST_ARTICLE,
                                     ncol = 3,
                                     nrow = 2,
                                     title = cell,
                                     colors = c("grey","red","blue","orange"),
                                     boxplot = TRUE)
    print(plotviolin)


    dev.off()
}






