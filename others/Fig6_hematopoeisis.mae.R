## ----results = "hide", warning = FALSE, message = FALSE-------------------------------------
mae <- scMultiome::hematopoiesis()

# Load peak matrix
PeakMatrix <- mae[["PeakMatrix"]]

# Load expression matrix
GeneExpressionMatrix <- mae[["GeneIntegrationMatrix"]]

# Add gene symbols to rownames
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

# Transfer dimensionality reduction matrix to GeneExpression
reducedDim(GeneExpressionMatrix, "IterativeLSI") <-
  reducedDim(mae[['TileMatrix500']], "IterativeLSI")
reducedDim(GeneExpressionMatrix, "UMAP") <-
  reducedDim(mae[['TileMatrix500']], "UMAP")



## -------------------------------------------------------------------------------------------

scater::plotReducedDim(GeneExpressionMatrix,
                       dimred = "UMAP",
                       text_by = "Clusters2",
                       colour_by = "Clusters2",
                       point_size = 0.3,
                       point_alpha = 0.3)



## -------------------------------------------------------------------------------------------
library(epiregulon)
grl <- getTFMotifInfo(genome = "hg19")
grl



## -------------------------------------------------------------------------------------------

set.seed(1010)
p2g <- calculateP2G(peakMatrix = PeakMatrix,
                    expMatrix = GeneExpressionMatrix,
                    exp_assay = "normalizedCounts",
                    reducedDim = reducedDim(GeneExpressionMatrix, "IterativeLSI"))

p2g


## -------------------------------------------------------------------------------------------
overlap <- addTFMotifInfo(grl = grl, p2g = p2g, peakMatrix = PeakMatrix)
head(overlap)


## -------------------------------------------------------------------------------------------
regulon <- getRegulon(p2g, overlap, aggregate=FALSE)
regulon


## -------------------------------------------------------------------------------------------

pruned.regulon <- pruneRegulon(expMatrix = GeneExpressionMatrix,
                               exp_assay = "normalizedCounts",
                               peakMatrix = PeakMatrix,
                               peak_assay = "counts",
                               regulon = regulon,
                               prune_value = "pval",
                               regulon_cutoff = 0.05,
                               clusters = GeneExpressionMatrix$Clusters2)



## -------------------------------------------------------------------------------------------
set.seed(1010)
regulon.w <- addWeights(regulon = pruned.regulon,
                        expMatrix  = GeneExpressionMatrix,
                        exp_assay  = "normalizedCounts",
                        peakMatrix = PeakMatrix,
                        peak_assay = "counts",
                        clusters = GeneExpressionMatrix$Clusters2,
                        aggregateCells = TRUE,
                        method = "wilcox",
                        useDim = "IterativeLSI")


regulon.w

saveRDS(regulon, "/gstore/project/lineage/manuscript/epiregulon/OUTPUT/regulon.rds")
saveRDS(regulon.w, "/gstore/project/lineage/manuscript/epiregulon/OUTPUT/regulon.w.rds")

## -------------------------------------------------------------------------------------------
score.combine <- calculateActivity(expMatrix = GeneExpressionMatrix,
                                   regulon = regulon.w,
                                   mode = "weight",
                                   method = "weightedMean",
                                   exp_assay = "normalizedCounts")
head(score.combine[1:5,1:5])


## -------------------------------------------------------------------------------------------
library(epiregulon.extra)
markers  <- findDifferentialActivity(activity_matrix = score.combine,
                                    clusters = GeneExpressionMatrix$Clusters2,
                                    pval.type = "some",
                                    direction = "up",
                                    test.type = "t")


## -------------------------------------------------------------------------------------------
markers.sig <- getSigGenes(markers, topgenes = 3 )



## -------------------------------------------------------------------------------------------
options(ggrastr.default.dpi=300)
tfs_interest <- c("EBF1","PAX5", "GATA3","SPI1")
plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = score.combine[tfs_interest,],
                tf = tfs_interest,
                dimtype = "UMAP",
                nrow=2,
                ncol=2,
                point_size=0.1,
                rasterise = TRUE)


## -------------------------------------------------------------------------------------------

plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = assay(GeneExpressionMatrix, "normalizedCounts")[tfs_interest,],
                tf = tfs_interest,
                dimtype = "UMAP",
                nrow=2,
                ncol=2,
                legend.label = "Gex",
                colors = c("grey","blue"),
                point_size=0.1,
                rasterise = TRUE)


## -------------------------------------------------------------------------------------------
plotActivityViolin(activity_matrix = score.combine,
                   tf = tfs_interest,
                   clusters = GeneExpressionMatrix$Clusters2,
                   legend.label = "Gex",
                   nrow=2,
                   ncol=2)


## -------------------------------------------------------------------------------------------
plotActivityViolin(activity_matrix = assay(GeneExpressionMatrix, "normalizedCounts")[tfs_interest,],
                   tf = tfs_interest,
                   clusters = GeneExpressionMatrix$Clusters2,
                   nrow=2,
                   ncol=2,
                   legend.label = "gene expression")


## -------------------------------------------------------------------------------------------
plotBubble(activity_matrix = score.combine,
           tf = tfs_interest,
           GeneExpressionMatrix$Clusters2,
           bubblesize = "FDR")


## ----fig.height=9---------------------------------------------------------------------------
plotBubble(activity_matrix = score.combine,
           tf = markers.sig$tf,
           GeneExpressionMatrix$Clusters2,
           bubblesize = "FDR")


## ----enrichment_hematopoeisis, fig.width=12, fig.height=12----------------------------------
#retrieve genesets
H <- EnrichmentBrowser::getGenesets(org = "hsa",
                                    db = "msigdb",
                                    cat = "H",
                                    gene.id.type = "SYMBOL" )
C2 <- EnrichmentBrowser::getGenesets(org = "hsa",
                                     db = "msigdb",
                                     cat = "C2",
                                     gene.id.type = "SYMBOL" )


#combine genesets and convert genesets to be compatible with enricher
gs <- c(H, C2)
gs.list <- do.call(rbind,lapply(names(gs), function(x)
  {data.frame(gs=x, genes=gs[[x]])}))

enrichresults <- regulonEnrich(TF = tfs_interest,
                               regulon = regulon.w,
                               weight = "weight",
                               weight_cutoff = 0,
                               genesets = gs.list)

#plot results

pdf("OUTPUT/heme.mono.CD4.M.enrich.pdf", width=14, height=8)
enrichPlot(results = enrichresults, ncol=2)
dev.off()

## ----differential networks hematopoeisis----------------------------------------------------

pdf("OUTPUT/heme.EBF1.preB.CD4.M.pdf")
plotDiffNetwork(regulon.w,
                cutoff = 0,
                tf = c("EBF1"),
                weight = "weight",
                clusters = c("PreB","CD4.M"),
                layout = "stress")
dev.off()


## -------------------------------------------------------------------------------------------
library(ggplot2)
# construct a graph of the preB cells
preB_network <- buildGraph(regulon.w, weights = "weight", cluster="PreB")

# compute a similarity matrix of all TFs
similarity_score <- calculateJaccardSimilarity(preB_network)

# Focus on EBF1
similarity_score_EBF1 <- similarity_score[, "EBF1"]
similarity_df <- data.frame(similarity = head(sort(similarity_score_EBF1,
                                                   decreasing = TRUE),20),
                            TF = names(head(sort(similarity_score_EBF1,
                                                 decreasing = TRUE),20)))

similarity_df$TF <- factor(similarity_df$TF, levels = rev(unique(similarity_df$TF)))

# plot top TFs most similar to EBF1
topTFplot <- ggplot(similarity_df, aes(x=TF, y=similarity)) +
  geom_bar(stat="identity") +
  coord_flip() +
  ggtitle("EBF1 similarity") +
  theme_classic()

print(topTFplot)



## -------------------------------------------------------------------------------------------

# create a permuted graph by rewiring the edges 100 times
permute_matrix <- permuteGraph(preB_network, "EBF1", 100, p=1)
permute_matrix <- permute_matrix[names(similarity_score_EBF1),]
diff_matrix <- similarity_score_EBF1-rowMeans(permute_matrix)

diff_matrix_df <- data.frame(similarity = head(sort(diff_matrix,
                                                    decreasing = TRUE),20),
                            TF = names(head(sort(diff_matrix,
                                                 decreasing = TRUE),20)))

diff_matrix_df$TF <- factor(diff_matrix_df$TF, levels = rev(unique(diff_matrix_df$TF)))

# plot top TFs most similar to EBF1
topTFplot <- ggplot(diff_matrix_df, aes(x=TF, y=similarity)) +
            geom_bar(stat="identity") +
            coord_flip() +
            ggtitle("background subtracted EBF1 similarity ") +
            theme_classic()

pdf("OUTPUT/heme.topTF.EBF1.pdf")
print(topTFplot)
dev.off()

# obtain empirical p-values
p_matrix <- rowMeans(apply(permute_matrix, 2, function(x) {x > similarity_score_EBF1}))
p_matrix[names(head(sort(diff_matrix,decreasing = TRUE),20))]


## -------------------------------------------------------------------------------------------
#regulon.w.2 <- regulon.w
#regulon.w <- readRDS("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/regulon.w.rds")
# construct a graph of the CD4.M and NK cells respectively
CD4.M_network <- buildGraph(regulon.w, weights = "weight", cluster="CD4.M")
Mono_network <- buildGraph(regulon.w, weights = "weight", cluster="Mono")

# construct a difference graph
diff_graph <- buildDiffGraph(Mono_network,CD4.M_network, abs_diff = FALSE)
diff_graph <- addCentrality(diff_graph)
diff_graph <- normalizeCentrality(diff_graph)
rank_table <- rankTfs(diff_graph)

library(ggplot2)

pdf("OUTPUT/heme.mono.CD4.M.pdf", width=8, height=5)
ggplot(rank_table, aes(x = rank, y = centrality)) +
    geom_point() +
    ggrepel::geom_text_repel(data = rbind(head(rank_table, 10),
                                          tail(rank_table, 10)),
                             aes(label = tf),
                             nudge_x = 0, nudge_y = 0, box.padding = 0.5, max.overlaps = Inf) +
    theme_classic() + ggtitle ("differential TFs (Mono-CD4.M) ranked by degree centrality")



## ----results = "hide", warning = FALSE, message = FALSE-------------------------------------

dev.off()

library(igraph)

## -------------------------------------------------------------------------------------------
diff_graph_filter <- subgraph.edges(diff_graph,
                                    E(diff_graph)[E(diff_graph)$weight>0],
                                    del=TRUE)


# compute a similarity matrix of all TFs
similarity_score <- calculateJaccardSimilarity(diff_graph_filter)

# Focus on SPI1
similarity_score_SPI1 <- similarity_score[, "SPI1"]
similarity_df <- data.frame(similarity = head(sort(similarity_score_SPI1,
                                                   decreasing = TRUE),20),
                            TF = names(head(sort(similarity_score_SPI1,
                                                 decreasing = TRUE),20)))

similarity_df$TF <- factor(similarity_df$TF,
                           levels = rev(unique(similarity_df$TF)))

# plot top TFs most similar to SPI1
topTFplot <- ggplot(similarity_df, aes(x=TF, y=similarity)) +
  geom_bar(stat="identity") +
  coord_flip() +
  ggtitle("SPI1 similarity") +
  theme_classic()

print(topTFplot)

pdf("OUTPUT/heme.topTF.SPI.pdf")
print(topTFplot)
dev.off()
## -------------------------------------------------------------------------------------------
sessionInfo()

