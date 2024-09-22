rm(list=ls())
library(epiregulon)
library(scMultiome)
library(igraph)
library(ggplot2)
library(ArchR)
library(ggpubr)
library(epiregulon.extra)
RNGkind()
mae <- scMultiome::reprogramSeq()

# peak matrix
PeakMatrix <- mae[["PeakMatrix"]]

# expression matrix
GeneExpressionMatrix <- mae[["GeneExpressionMatrix"]]
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

# dimensional reduction matrix
reducedDimMatrix <- reducedDim(mae[['TileMatrix500']], "LSI_ATAC")
reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- reducedDim(mae[['TileMatrix500']], "UMAP_Combined")


# plot distribution of hash assignment
GeneExpressionMatrix$hash_assignment <- factor(GeneExpressionMatrix$hash_assignment, levels = c("HTO3_NKX2.1_v2","HTO8_NKX2.1_UTR",
                                                                                                "HTO10_GATA6_UTR","HTO2_GATA6_v2",
                                                                                                "HTO6_hFOXA1_UTR","HTO4_mFOXA1_v2",
                                                                                                "HTO1_FOXA2_v2","HTO5_NeonG_v2"))


GeneExpressionMatrix <- scuttle::logNormCounts(GeneExpressionMatrix)

pdf("reprogramseq/preprocessing/OUTPUT/hash_distribution.pdf", width=5, height=5)
ggplot(data.frame(colData(GeneExpressionMatrix)),
       aes(fill = hash_assignment, y = 1, x = Clusters)) +
    geom_bar(position = "fill", stat = "identity") + theme_classic()
ggplot(data.frame(colData(GeneExpressionMatrix)),
       aes(fill = Clusters, y = 1, x = hash_assignment)) +
    geom_bar(position = "fill", stat = "identity") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()



# plot cluster and hash assignment as UMAP
pdf("reprogramseq/preprocessing/OUTPUT/Reprogramseq.umap.pdf")
options(ggrastr.default.dpi=300)
scater::plotReducedDim(GeneExpressionMatrix,
                       dimred = "UMAP_Combined",
                       text_by = "Clusters",
                       colour_by = "Clusters",
                       rasterise = TRUE)

options(ggrastr.default.dpi=300)
scater::plotReducedDim(GeneExpressionMatrix,
                       dimred = "UMAP_Combined",
                       text_by = "hash_assignment",
                       colour_by = "hash_assignment",
                       rasterise = TRUE)
dev.off()

# perform GRN using epiregulon
grl <- getTFMotifInfo(genome = "hg38")
head(grl)

set.seed(1010, "Mersenne-Twister")
p2g <- calculateP2G(peakMatrix = PeakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    exp_assay = "normalizedCounts",
                    reducedDim = reducedDimMatrix)

p2g

saveRDS(p2g, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.p2g.rds")


overlap <- addTFMotifInfo(grl = grl, p2g = p2g, peakMatrix = PeakMatrix)
head(overlap)

saveRDS(overlap, "reprogramseq/preprocessing/OUTPUT/overlap.rds")

regulon <- getRegulon(p2g = p2g, overlap = overlap, aggregate = FALSE)
regulon

saveRDS(regulon, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.regulon.rds")

pruned.regulon <- pruneRegulon(expMatrix = GeneExpressionMatrix,
                               exp_assay = "normalizedCounts",
                               peakMatrix = PeakMatrix,
                               peak_assay = "counts",
                               test = "chi.sq",
                               regulon,
                               clusters = GeneExpressionMatrix$Clusters,
                               prune_value = "pval",
                               regulon_cutoff = 0.05
)

pruned.regulon

saveRDS(pruned.regulon, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.pruned.regulon.rds")

#wilcoxon
regulon.w <- addWeights(regulon = pruned.regulon,
                        expMatrix  = GeneExpressionMatrix,
                        exp_assay  = "normalizedCounts",
                        peakMatrix = PeakMatrix,
                        peak_assay = "counts",
                        clusters = GeneExpressionMatrix$Clusters,
                        method = "wilcox")


saveRDS(regulon.w, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.regulon.w.rds")

regulon.w.GATA6_NKX21 <- regulon.w[regulon.w$tf %in% c("GATA6","NKX2-1"),]
saveRDS(regulon.w.GATA6_NKX21, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.regulon.w.NKX2_1.GATA6.rds")


score.combine <- calculateActivity(expMatrix = GeneExpressionMatrix,
                                   regulon = regulon.w,
                                   mode = "weight",
                                   method = "weightedMean",
                                   exp_assay = "normalizedCounts")

saveRDS(score.combine, "reprogramseq/preprocessing/OUTPUT/reprogram.seq.score.combine.rds")



# perform global differential activity
markers <- findDifferentialActivity(
    score.combine,
    GeneExpressionMatrix$hash_assignment,
    pval.type = "some",
    direction = "up",
    test.type = "t"
)
markers <- getSigGenes(markers, topgenes = 8)



pdf("reprogramseq/preprocessing/OUTPUT/reprogram.seq.bubble.pdf", width=5, height = 5)
plotBubble(
    activity_matrix = score.combine[,GeneExpressionMatrix$hash_assignment!= "HTO3_NKX2.1_v2"],
    tf = markers$tf,
    clusters = GeneExpressionMatrix$hash_assignment[GeneExpressionMatrix$hash_assignment!= "HTO3_NKX2.1_v2"]
)

plotBubble(
    activity_matrix = score.combine[,GeneExpressionMatrix$hash_assignment!= "HTO3_NKX2.1_v2"],
    tf = c("GATA6","NKX2-1","FOXA2","FOXA1"),
    clusters = GeneExpressionMatrix$hash_assignment[GeneExpressionMatrix$hash_assignment!= "HTO3_NKX2.1_v2"]
)
dev.off()

hash_matrix <- matrix(0, nrow=2, ncol = ncol(score.combine))
rownames(hash_matrix) <- c("NKX2-1","GATA6")
colnames(hash_matrix) <- colnames(score.combine)
hash_matrix["NKX2-1", which(GeneExpressionMatrix$hash_assignment ==  "HTO8_NKX2.1_UTR")] <- 1
hash_matrix["GATA6", GeneExpressionMatrix$hash_assignment %in%  c("HTO10_GATA6_UTR", "HTO2_GATA6_v2")] <- 1


# plot TF activity, gex in UMAP
pdf("reprogramseq/preprocessing/OUTPUT/reprogram.seq.gex.activty.umap.pdf", width=10, height = 6)
options(ggrastr.default.dpi=300)
plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = score.combine,
                tf = c("NKX2-1","GATA6","FOXA2","FOXA1", "AR"),
                dimtype = "UMAP_Combined",
                label = "Clusters",
                point_size = 0.1,
                colors = c("darkblue","yellow"),
                limit = c(0,0.1),
                ncol = 3,
                rasterise=TRUE)

plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = assay(GeneExpressionMatrix,"logcounts"),
                tf = c("NKX2-1","GATA6","FOXA2","FOXA1", "AR"),
                dimtype = "UMAP_Combined",
                label = "Clusters",
                point_size = 0.1,
                ncol = 3,
                limit = c(0,1),
                colors = c("grey","blue"),
                legend.label = "GEX",
                rasterise=TRUE)

plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = hash_matrix,
                tf = c("NKX2-1","GATA6"),
                dimtype = "UMAP_Combined",
                label = "Clusters",
                point_size = 0.1,
                ncol = 3,
                nrow = 2,
                limit = c(0,1),
                colors = c("grey","darkgreen"),
                legend.label = "GEX",
                rasterise=TRUE)
dev.off()

# plot TF activity as violin plots
pdf("reprogramseq/preprocessing/OUTPUT/violin.plots.hash.pdf", width = 4, height = 3)

plotActivityViolin(activity_matrix = score.combine,
                   tf = c("NKX2-1","GATA6"),
                   clusters = GeneExpressionMatrix$hash_assignment)

plotActivityViolin(activity_matrix = assay(GeneExpressionMatrix,"logcounts"),
                   tf = c("NKX2-1","GATA6"),
                   clusters = GeneExpressionMatrix$hash_assignment,legend.label = "gex",
                   )

dev.off()


pdf("reprogramseq/preprocessing/OUTPUT/violin.plots.clusters.pdf", width = 4, height = 2)

plotActivityViolin(activity_matrix = score.combine,
                   tf = c("NKX2-1","GATA6"),
                   clusters = GeneExpressionMatrix$Clusters)
plotActivityViolin(activity_matrix = assay(GeneExpressionMatrix,"logcounts"),
                   tf = c("NKX2-1","GATA6"),
                   clusters = GeneExpressionMatrix$Clusters, legend.label = "gex")
dev.off()



# add logFC with respect to neonG

outdir <- "/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/"
proj.all <- loadArchRProject(path = outdir, showLogo = TRUE)

proj.all$Clusters2 <- proj.all$Clusters
proj.all$Clusters2[proj.all$hash_assignment == "HTO5_NeonG_v2"] <- "mNeonG"

set.seed(1010, "Mersenne-Twister")
differential_genes <- getMarkerFeatures(
    ArchRProj = proj.all,
    useMatrix = "GeneExpressionMatrix",
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = c("C1","C3"),
    bgdGroups = c("mNeonG"),
    logFile = "differentialgenes"
)


gene_exp <- assay(differential_genes, "Log2FC")
rownames(gene_exp) <- rowData(differential_genes)[,"name"]
gene_exp[c("NKX2-1","GATA6"),]

for (cluster in colnames(gene_exp)) {
    regulon.w[,paste0("scLogFCexp_", cluster)] <- gene_exp[match(regulon.w$target,rownames(gene_exp)), cluster]
}



# correlate regulon weights vs logFC

C3.weight.wilcox.vs.logfc <- ggscatter(data.frame(regulon.w[regulon.w$tf=="NKX2-1",]), x = "scLogFCexp_C3", y = "weight.all",
                                       color = "black", shape = 20, size = 1, # Points color, shape and size
                                       add = "reg.line",  # Add regressin line
                                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                       conf.int = TRUE, # Add confidence interval
                                       cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                       cor.coeff.args = list(method = "spearman", label.sep = "\n"),
                                       ylab = "regulon weight (wilcoxon)", xlab = "LogFC relative to uninfected",
                                       title = "logFC vs regulon weight (wilcoxon) C3")

C1.weight.wilcox.vs.logfc <- ggscatter(data.frame(regulon.w[regulon.w$tf=="GATA6",]), x = "scLogFCexp_C1", y = "weight.all",
                                       color = "black", shape = 20, size = 1, # Points color, shape and size
                                       add = "reg.line",  # Add regressin line
                                       add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                       conf.int = TRUE, # Add confidence interval
                                       cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                       cor.coeff.args = list(method = "spearman", label.sep = "\n"),
                                       ylab = "regulon weight (wilcoxon)", xlab = "LogFC relative to uninfected",
                                       title = "logFC vs regulon weight (wilcoxon) C1")



pdf("reprogramseq/preprocessing/OUTPUT/weight.vs.logFC.pdf", height = 4, width =4)
print(C3.weight.wilcox.vs.logfc)
print(C1.weight.wilcox.vs.logfc)
dev.off()


# plot target genes in regulons
regulon.w <- readRDS("reprogramseq/preprocessing/OUTPUT/reprogram.seq.regulon.w.rds")
score.combine <- readRDS("reprogramseq/preprocessing/OUTPUT/reprogram.seq.score.combine.rds")

pdf("reprogramseq/preprocessing/OUTPUT/regulon.heatmap.NKX2-1.GATA6.pdf", width=9, height = 4)

plotHeatmapRegulon(sce=GeneExpressionMatrix,
                   tfs= c("GATA6","NKX2-1"),
                   regulon=regulon.w,#[regulon.w$C1.vs.uninfected.FDR <0.1 & abs(regulon.w$C1.vs.uninfected.logFC>0.5),],
                   regulon_cutoff=0,
                   downsample=4000,
                   color_breaks = c(-0.5, 0, 0.5),
                   colors = c("darkblue", "black", "yellow"),
                   cell_attributes="Clusters",
                   col_gap="Clusters",
                   exprs_values="normalizedCounts",
                   #column_title="LNCaP AR regulon heatmap",
                   raster_quality = 10,
                   show_row_names = FALSE)

dev.off()



# chromvar chromVar of NEPC matrix
# load chromvar matrix
NEPCMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "NEPCMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)

NEPCMatrix <- as(NEPCMatrix,"SingleCellExperiment")
reducedDim(NEPCMatrix, "UMAP_Combined") <- getEmbedding(ArchRProj = proj.all,
                                                        embedding = "UMAP_Combined",
                                                        returnDF = TRUE)[colnames(NEPCMatrix), ]

chromvar_z <- assay(NEPCMatrix, "z")


# plot chromVar
pdf("reprogramseq/preprocessing/OUTPUT/reprogram.seq.chromvar.umap.pdf", width = 12, height = 8)
options(ggrastr.default.dpi=300)
plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = assay(NEPCMatrix,"z"),
                tf = rownames(NEPCMatrix),
                dimtype = "UMAP_Combined",
                label = "Clusters",
                point_size = 0.5,
                ncol = 3,
                nrow=2,
                colors = c("grey","red"),
                limit = c(0,2),
                rasterise = TRUE)

dev.off()

combo <- t(rbind(score.combine, chromvar_z))
combo <- data.frame(as.matrix(combo))
combo$Clusters <- GeneExpressionMatrix$Clusters

pdf("reprogramseq/preprocessing/OUTPUT/corr.activity.chromvar.pdf", width=5, height=5)


for (TF in c("GATA6", "NKX2.1")){

    # if (TF == "GATA6"){
    #     combo.select <- data.frame(combo)[colnames(GeneExpressionMatrix)[which(GeneExpressionMatrix$Clusters == "C1")],]
    # } else {
    #     combo.select <- data.frame(combo)[colnames(GeneExpressionMatrix)[which(GeneExpressionMatrix$Clusters == "C3")],]
    # }
    combo.select <- combo
    for (regions in c("Tang_CRPC_SCL", "Tang_CRPC_NE")) {
        scatterplot <- ggscatter(combo.select, x = TF, y = regions,
                                 color = "Clusters", size = 1, # Points color, shape and size
                                 add = "reg.line",  # Add regressin line
                                 add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                 conf.int = TRUE, # Add confidence interval
                                 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                                 cor.coeff.args = list(method = "pearson",  label.sep = "\n"),
                                 title =  paste(TF, regions),
                                 xlab = "TF activity(epiregulon)",
                                 ylab = "chromatin accessibility (chromvar)",
                                 cor.coef.size = 8 ) + theme(text = element_text(size = 20))
        print(scatterplot)
    }

}


dev.off()


## chromvar TF binding
TF_bindingMatrix <- mae[["TF_bindingMatrix"]]
reducedDim(TF_bindingMatrix, "UMAP_Combined") <- reducedDim(mae[['TileMatrix500']], "UMAP_Combined")



# plot chromVar
pdf("reprogramseq/preprocessing/OUTPUT/reprogram.seq.chromvar.TFbinding.umap.pdf", width = 12, height = 8)
options(ggrastr.default.dpi=300)
plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = assay(TF_bindingMatrix,"z"),
                tf = c("NKX2-1","GATA6","FOXA1","FOXA2", "AR"),
                dimtype = "UMAP_Combined",
                label = "Clusters",
                point_size = 0.5,
                ncol = 3,
                nrow=2,
                colors = c("grey","red"),
                limit = c(0,2),
                rasterise = TRUE)

dev.off()





# # differential network
# topTFplot <- list()
# topTFplot_ratio <- list()
# topTFplot_diff <- list()
#
# cluster <- c(GATA6="C1",NKX2.1="C3", FOXA2="C4")
#
# for (TF in c("GATA6", "NKX2-1","FOXA2")) {
#
#     # build graphs
#     TF_network <- buildGraph(regulon.w, weights = "weight", cluster=cluster[make.names(TF)])
#     uninfected_network <- buildGraph(regulon.w, weights = "weight", cluster="C5")
#
#     # build differential graphs
#     diff_graph <- buildDiffGraph(TF_network, uninfected_network)
#     diff_graph_filter <- subgraph.edges(diff_graph, E(diff_graph)[E(diff_graph)$weight>0], del=T)
#     diff_graph_filter <- addCentrality(diff_graph_filter)
#     diff_graph_filter <- normalizeCentrality(diff_graph_filter)
#     rank_table <- rankTfs(diff_graph_filter)
#     head(rank_table,20)
#
#     similarity_score <- calculateJaccardSimilarity(diff_graph_filter )
#     similarity_score_TF <- sort(similarity_score[, TF], decreasing = TRUE)
#
#     # permutation
#     set.seed(1010)
#     permute_matrix <- permuteGraph(diff_graph_filter, TF, 100, p=1)
#     permute_matrix <- permute_matrix[names(similarity_score_TF),]
#     ratio_matrix <- similarity_score_TF/rowMeans(permute_matrix)
#     diff_matrix <- similarity_score_TF-rowMeans(permute_matrix)
#     p_matrix <- rowMeans(apply(permute_matrix, 2, function(x) {x > similarity_score_TF}))
#
#
#     # similarity matrix
#     similarity_df <- data.frame(similarity = head(sort(similarity_score_TF, decreasing = TRUE),20),
#                                 TF = names(head(sort(similarity_score_TF,decreasing = TRUE),20)))
#     similarity_df$TF <- factor(similarity_df$TF, levels = rev(unique(similarity_df$TF)))
#     topTFplot[[TF]] <- ggplot(similarity_df, aes(x=TF, y=similarity)) +
#         geom_bar(stat="identity") +
#         coord_flip() +
#         ggtitle(TF) +
#         theme_classic()
#
#     # ratio similarity
#
#     ratio_df <- data.frame(similarity = head(sort(ratio_matrix, decreasing = TRUE), 20),
#                            TF= names(head(sort(ratio_matrix, decreasing = TRUE), 20)))
#     ratio_df$TF <- factor(ratio_df$TF, levels = rev(unique(ratio_df$TF)))
#     topTFplot_ratio[[TF]] <- ggplot(ratio_df, aes(x=TF, y=similarity)) +
#         geom_bar(stat="identity") +
#         coord_flip() +
#         ggtitle(TF) +
#         ylab("ratio in similarity") +
#         theme_classic()
#
#     # diff similarity
#     diff_df <- data.frame(similarity = head(sort(diff_matrix, decreasing = TRUE), 20),
#                           TF= names(head(sort(diff_matrix, decreasing = TRUE), 20)))
#     diff_df$TF <- factor(diff_df$TF, levels = rev(unique(diff_df$TF)))
#     topTFplot_diff[[TF]] <- ggplot(diff_df, aes(x=TF, y=similarity)) +
#         geom_bar(stat="identity") +
#         coord_flip() +
#         ggtitle(TF) +
#         ylab("diff in similarity") +
#         theme_classic()
#
#
# }
#
#
#
# pdf("reprogramseq/OUTPUT/GATA6.similarity.pdf", width=9, height=4)
# gridExtra::grid.arrange(grobs=topTFplot, ncol=3)
# gridExtra::grid.arrange(grobs=topTFplot_ratio, ncol=3)
# gridExtra::grid.arrange(grobs=topTFplot_diff, ncol=3)
# dev.off()
#
# # visualization
#
# pdf("reprogramseq/OUTPUT/plotdiffnetwork.pdf", width=5, height=4)
#
# plotDiffNetwork(regulon.w,
#                 cutoff = 0,
#                 tf = c("GATA6"),
#                 weight = "weight",
#                 clusters = c("C1","C5"),
#                 layout = "stress")
#
# plotDiffNetwork(regulon.w,
#                 cutoff = 0,
#                 tf = c("NKX2-1"),
#                 weight = "weight",
#                 clusters = c("C3","C5"),
#                 layout = "stress")
#
#
# plotDiffNetwork(regulon.w,
#                 cutoff = 0,
#                 tf = c("GATA6","YAP1"),
#                 weight = "weight",
#                 clusters = c("C1"),
#                 layout = "stress")
#
# dev.off()
#

