library(scran)
library(scater)
library(epiregulon)
library(ArchR)
library(BiocParallel)
library(eulerr)
library(ggpubr)
library(epiregulon.archr)
library(epiregulon.extra)
#checkdata
archR_project_path <- "OUTPUT/ArchRProject/"
proj.all <- loadArchRProject(path = archR_project_path, showLogo = TRUE)


#import grl from chip
# ExperimentHub::setExperimentHubOption("CACHE", "/gstore/scratch/u/yaox19/.cache/ExperimentHub")
#grl.all <- getTFMotifInfo("hg38","atlas.sample")
grl.all <- readRDS("/gne/data/genomics/external_sources/chipAtlas/chipAtlas/chip_peaks/grl.chipatlas.sample.hg38.qc.rds")
grl <- list()
grl[["LNCaP"]] <- grl.all$LNCAP
grl[["VCaP"]] <- grl.all$VCaP
grl[["22Rv1"]] <- grl.all$`22Rv1`

grl[["MDA"]][["SMARCA4"]] <- rtracklayer::import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allSMARCA4.DBA_DESEQ2.bed")
grl[["MDA"]][["FOXA1"]] <- rtracklayer::import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allFOXA1.DBA_DESEQ2.bed")
grl[["MDA"]][["AR"]] <- rtracklayer::import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allAR.DBA_DESEQ2.bed")
grl[["MDA"]] <- GenomicRanges::GRangesList(grl[["MDA"]])

celllines <- names(grl)

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



for (cell in  celllines ){
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

    # find peak to gene links
    p2g <- calculateP2G(peakMatrix = peakMatrix.cell,
                        expMatrix = GeneExpressionMatrix.select ,
                        reducedDim = reducedDim(GeneExpressionMatrix.select ),
                        peak_assay = "counts",
                        exp_assay = "normalizedCounts",
                        cor_cutoff = 0.5
    )


    # Construct regulons
    overlap <- addTFMotifInfo(grl = grl[[cell]],
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


    saveRDS(regulon.w, paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.", cell, ".chip.rds"))
    saveRDS(score.combine, paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cell, ".chip.rds"))

    #########
    pdf(paste0("OUTPUT/ArchRProject/Epiregulon/activity.separate.", cell,".chip.pdf"), width = 9, height = 6)


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


plot.scatter <- list()

for (cell in celllines){
    score_chip <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cell, ".chip.rds"))
    score_public <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cell, ".rds"))

    for (tf in c("AR", "FOXA1")){
        message(cell, tf)
        activity.df <- data.frame(chip = score_chip[tf, ],
                                  public = score_public[tf, ])

        plot.scatter[[paste0(cell,"_",tf)]] <- ggscatter(activity.df,
                  x = "chip",
                  y = "public",
                  xlab = "cell line matched chip", ylab = "public chip",
                  color = "black", shape = 21, size = 0.5, # Points color, shape and size
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                  cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                  title = paste(tf, cell)
        )


    }


}


pdf("OUTPUT/ArchRProject/Epiregulon/overlap.targets.public.vs.chip.corr.pdf", width = 6, height=12)

gridExtra::grid.arrange(grobs=plot.scatter, ncol=2)
dev.off()







