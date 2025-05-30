library(epiregulon)
library(ArchR)
library(scran)
library(scater)
library(epiregulon.archr)
options(ggrastr.default.dpi=300)

number <-  commandArgs(trailingOnly=TRUE)
print(number)
part_no <- paste0('part_', number)

selected_cells <- c("Tumor", "Ductal-like1", "Ductal-like2", "Acinar", "Acinar REG+")

# import geneIntegrationMatrix
archR_project_path <- paste0("/gstore/scratch/u/yaox19/TGI/No_866_ATAC/", part_no, "/output_TSS4_Frags1000")
archr.proj <- loadArchRProject(path = archR_project_path, showLogo = FALSE)

# add information from scRNAseq
sce <- dsassembly::getExperiment('DS000016526', part_no)

archr.proj$clusterName <- colData(sce)[archr.proj$predictedCell,"clusterName"]

rm(sce)

# subset to epithelial cells
archr.proj <- archr.proj[archr.proj$clusterName %in% selected_cells ,]


GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = archr.proj,
    useMatrix = "GeneIntegrationMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = 1,
    logFile = createLogFile("getMatrixFromProject")
)


GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix,
                                        rename = "NormalizedCounts",
                                        transform = TRUE,
                                        transform_method = "log")
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

# remove duplicated rownames
GeneExpressionMatrix <- GeneExpressionMatrix[which(!duplicated(rownames(GeneExpressionMatrix))),]


# Add embeddings
reducedDim(GeneExpressionMatrix, "UMAP_ATAC") <- getEmbedding(ArchRProj = archr.proj,
                                                              embedding = "UMAP",
                                                              returnDF = TRUE)[colnames(GeneExpressionMatrix), ]

pdf("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/ATAC.umap.pdf", width=4, height=4)
plotReducedDim(GeneExpressionMatrix, dimred="UMAP_ATAC", color_by ="clusterName", text_by="clusterName", rasterise = TRUE, point_size = 0.1, point_alpha = 0.5)
dev.off()


# import peakMatrix
PeakMatrix <- getMatrixFromProject(
    ArchRProj = archr.proj,
    useMatrix = "PeakMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads=1
)
PeakMatrix <- ArchRMatrix2SCE(PeakMatrix,rename = "counts")

# subsample
set.seed(1010)
subsample.idx <- 1:ncol(GeneExpressionMatrix)

GeneExpressionMatrix <- GeneExpressionMatrix[, subsample.idx]
PeakMatrix <- PeakMatrix[, subsample.idx]


# epiregulon
grl <- getTFMotifInfo(genome = "hg38")
head(grl)

set.seed(1010)
p2g <- calculateP2G(peakMatrix = PeakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    exp_assay = "NormalizedCounts",
                    reducedDim = reducedDim(GeneExpressionMatrix, "UMAP_ATAC"),
                    cor_method = "spearman")


overlap <- addTFMotifInfo(grl = grl, p2g = p2g, peakMatrix = PeakMatrix)
head(overlap)

regulon <- getRegulon(p2g = p2g, overlap = overlap, aggregate = FALSE)
regulon

saveRDS(regulon, paste0("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/", part_no, ".regulon.full.rds"))


# prune regulon
pruned.regulon <- pruneRegulon(expMatrix = GeneExpressionMatrix,
                               exp_assay = "NormalizedCounts",
                               peakMatrix = PeakMatrix,
                               peak_assay = "counts",
                               test = "chi.sq",
                               regulon,
                               prune_value = "pval",
                               regulon_cutoff = 0.05
)

pruned.regulon

regulon.w <- addWeights(regulon = pruned.regulon,
                        expMatrix  = GeneExpressionMatrix,
                        exp_assay  = "NormalizedCounts",
                        peakMatrix = PeakMatrix,
                        peak_assay = "counts",
                        method = "wilcox")

regulon.w

saveRDS(regulon.w, paste0("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/", part_no, ".pruned.regulon.w.wilcox.rds"))

score.combine <- calculateActivity(expMatrix = GeneExpressionMatrix,
                                   regulon = regulon.w,
                                   mode = "weight",
                                   method = "weightedMean",
                                   exp_assay = "NormalizedCounts",
                                   normalize = FALSE)

saveRDS(score.combine, paste0("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/", part_no,".score.combine.rds"))

