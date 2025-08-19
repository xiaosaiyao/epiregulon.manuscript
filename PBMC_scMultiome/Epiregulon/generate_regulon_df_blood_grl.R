library(epiregulon)
library(BiocParallel)
library(GenomicRanges)


GeneExpressionMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/GeneExpressionMatrix.rds")
peakMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/peakMatrix.rds")
grl <- readRDS("/gne/data/genomics/external_sources/chipAtlas/chipAtlas/chip_peaks/grl.chipatlas.tissue.hg38.rds")

grl <- grl$Blood


set.seed(1010, kind ="L'Ecuyer-CMRG")

# find peak to gene links
p2g <- calculateP2G(peakMatrix = peakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    reducedDim = reducedDim(GeneExpressionMatrix, "LSI_Combined"),
                    peak_assay = "counts",
                    exp_assay = "normalizedCounts",
                    cor_cutoff = 0.5,
                    clusters = GeneExpressionMatrix$Clusters
)


# Construct regulons
overlap <- addTFMotifInfo(grl = grl,
                          p2g = p2g,
                          peakMatrix = peakMatrix)

regulon_df_full <- getRegulon(p2g, overlap, aggregate = FALSE)

selected_tfs = c("STAT6","ELK1", "GATA3", "JUN", "NFATC3", "NFKB1", "STAT3", "TCF7", "GATA3", "BCL11B", "SPI1", "CEBPA", "EBF1", "PAX5", "POU2F2", "POU2AF1", "EOMES", "RUNX3", "TBX21", "SPIB", "IRF8")


regulon_df_full <- regulon_df_full[regulon_df_full$tf %in% selected_tfs, ]

# prune network
pruned.regulon <- pruneRegulon(regulon = regulon_df_full,
                               expMatrix = GeneExpressionMatrix,
                               exp_assay = "normalizedCounts",
                               peakMatrix = peakMatrix,
                               peak_assay = "counts",
                               prune_value = "pval",
                               clusters = GeneExpressionMatrix$cell_type
)


regulon.w <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                        peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                        peak_assay="counts",
                        clusters = GeneExpressionMatrix$cell_type)


saveRDS(regulon.w, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_blood_trimmed.rds")


regulon.w2 <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                         peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                         peak_assay="counts", clusters = GeneExpressionMatrix$cell_type,
                         method = "MI")
saveRDS(regulon.w2, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_MI_blood_trimmed.rds")

regulon.w3 <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                         peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                         peak_assay="counts", clusters = GeneExpressionMatrix$cell_type,
                         method = "corr")
saveRDS(regulon.w3, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_blood_trimmed.rds")


