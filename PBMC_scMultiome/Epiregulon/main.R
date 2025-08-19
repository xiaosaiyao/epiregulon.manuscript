library(epiregulon)
set.seed(1010, kind ="L'Ecuyer-CMRG")

load(".RData")
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



regulon.w <- addWeights(pruned.regulon, expMatrix = GeneExpressionMatrix,
                        peakMatrix = peakMatrix, exp_assay="normalizedCounts",
                        peak_assay="counts", method = "corr", clusters=GeneExpressionMatrix$cell_type)


activity.matrix <- calculateActivity(expMatrix = GeneExpressionMatrix, exp_assay="normalizedCounts",
                                     regulon = regulon.w)

