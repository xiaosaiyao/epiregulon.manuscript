source("/gstore/project/epigen/benchmark/src/utils.R")
# load the MAE object
mae <- scMultiome::reprogramSeq()
# expression matrix
GeneExpressionMatrix <- mae[["GeneExpressionMatrix"]]
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

regulon <- readRDS("/gstore/project/epigen/benchmark/reprogram/OUTPUT/pruned.regulon.rds")
geneExpr <- as.matrix(assay(GeneExpressionMatrix, "normalizedCounts"))
colnames(geneExpr) <- colData(GeneExpressionMatrix)$Clusters

GATA6_tf <- regulon[regulon$tf == "GATA6",]
gene_metrics_all <- diff_expr_metrics(rownames(geneExpr), geneExpr, "C1")
gene_metrics_all$GATA6_target <- gene_metrics_all$gene %in% GATA6_tf$target
#gene_metrics_all <- gene_metrics_all[gene_metrics_all$gene %in% regulon$target,]
gene_metrics_all$combined_index <- rank(rank(-gene_metrics_all$mean_expr_diff)*rank(-gene_metrics_all$expr_prop))

GATA6_targets <- unique(GATA6_tf$target)
ranks <- setNames(rank(-gene_metrics_all$combined_index),gene_metrics_all$gene)

set.seed(1691, kind ="Mersenne-Twister")
fgseaRes <- fgsea::fgsea(list(GATA6 = GATA6_targets), ranks, minSize=15, maxSize=1100, scoreType = "pos")
# pathway  pval  padj log2err        ES      NES  size  leadingEdge
# <char> <num> <num>  <lgcl>     <num>    <num> <int>       <list>
#     1:   GATA6 1e-50 1e-50      NA 0.6139429 2.342736   490 NR3C2, Z....

dev.new()
pdf("/gstore/project/epigen/benchmark/reprogram/OUTPUT/plots/gsea_GATA6.pdf", width = 4, height = 3)
fgsea::plotEnrichment(GATA6_targets, ranks)
dev.off()

NKX21_tf <- regulon[regulon$tf == "NKX2-1",]
gene_metrics_C3 <- diff_expr_metrics(rownames(geneExpr), geneExpr, "C3")
gene_metrics_C3$NKX21_target <- gene_metrics_C3$gene %in% NKX21_tf$target
gene_metrics_C3$combined_index <- rank(rank(-gene_metrics_C3$mean_expr_diff)*rank(-gene_metrics_C3$expr_prop))
NKX21_targets <- unique(NKX21_tf$target)
ranks <- setNames(rank(-gene_metrics_C3$combined_index),gene_metrics_C3$gene)
set.seed(1691, kind ="Mersenne-Twister")
fgseaRes <- fgsea::fgsea(list(NKX21 = NKX21_targets), ranks, minSize=15, maxSize=3000, scoreType = "pos")
# pathway        pval        padj log2err        ES      NES  size  leadingEdge
# <char>       <num>       <num>   <num>     <num>    <num> <int>       <list>
#     1:   NKX21 1.35169e-27 1.35169e-27 1.36518 0.6249244 2.264333   146 LINC0057....


dev.new()
pdf("/gstore/project/epigen/benchmark/reprogram/OUTPUT/plots/gsea_NKX21.pdf", width = 4, height = 3)
fgsea::plotEnrichment(NKX21_targets,ranks)
dev.off()

