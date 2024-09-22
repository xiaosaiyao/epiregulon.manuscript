# archr
library(ArchR)

outdir <- "/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/"
archr.proj <- loadArchRProject(path = outdir, showLogo = TRUE)

set.seed(1010)
differential_genes <- getMarkerFeatures(
    ArchRProj = archr.proj,
    useMatrix = "GeneExpressionMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = "C1",
    bgdGroups = "C5",
    logFile = "differentialgenes"
)


archr.gene_exp <- assay(differential_genes, "Log2FC")
rownames(archr.gene_exp) <- rowData(differential_genes)[,"name"]


#scran
library(scran)
GeneExpressionMatrix <- getMatrixFromProject(archr.proj, useMatrix = "GeneExpressionMatrix")
GeneExpressionMatrix <- scuttle::logNormCounts(GeneExpressionMatrix, assay.type = "GeneExpressionMatrix" )
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
scran.results <- findMarkers(GeneExpressionMatrix, GeneExpressionMatrix$Clusters)
scran.gene_exp <- data.frame(gene = rownames(scran.results$C1),
               scran.C1.vs.C5.logFC = scran.results$C1$logFC.C5)

# combine
gene_exp <- scran.gene_exp
gene_exp$archr.C1.vs.C5 <- archr.gene_exp[scran.gene_exp$gene, "C1"]


ggplot(gene_exp, aes(scran.C1.vs.C5.logFC, archr.C1.vs.C5)) + geom_point()



########################
differential_genes <- getMarkerFeatures(
    ArchRProj = archr.proj,
    useMatrix = "GeneExpressionMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    useGroups = c("C1","C3"),
    bgdGroups = c("C5","C6"),
    logFile = "differentialgenes"
)


archr.gene_exp <- assay(differential_genes, "Log2FC")
rownames(archr.gene_exp) <- rowData(differential_genes)[,"name"]

#
GeneExpressionMatrix$Clusters2 <- GeneExpressionMatrix$Clusters
GeneExpressionMatrix$Clusters2[GeneExpressionMatrix$Clusters2 %in% c("C5","C6")] <- "uninfected"
scran.results <- findMarkers(GeneExpressionMatrix, GeneExpressionMatrix$Clusters2, sorted = FALSE)
scran.gene_exp <- data.frame(gene = rownames(scran.results$C1),
                             scran.C1.vs.uninfected.logFC = scran.results$C1$logFC.uninfected,
                             scran.C3.vs.uninfected.logFC = scran.results$C3$logFC.uninfected)



# combine
gene_exp <- scran.gene_exp
gene_exp$archr.C1.vs.uninfected <- archr.gene_exp[scran.gene_exp$gene, "C1"]
gene_exp$archr.C3.vs.uninfected <- archr.gene_exp[scran.gene_exp$gene, "C3"]


ggplot(gene_exp, aes(scran.C1.vs.uninfected.logFC, archr.C1.vs.uninfected)) + geom_point()
ggplot(gene_exp, aes(scran.C3.vs.uninfected.logFC, archr.C3.vs.uninfected)) + geom_point()


problem.genes <- gene_exp[gene_exp$scran.C3.vs.uninfected.logFC < 0.3 & gene_exp$scran.C3.vs.uninfected.logFC > -0.3 & gene_exp$archr.C3.vs.uninfected>2,]
problem.genes.C3 <- assay(GeneExpressionMatrix)[problem.genes$gene,colnames(GeneExpressionMatrix)[GeneExpressionMatrix$Clusters == "C3"]]
problem.genes.uninfected <- assay(GeneExpressionMatrix)[problem.genes$gene,colnames(GeneExpressionMatrix)[GeneExpressionMatrix$Clusters2 == "uninfected"]]

FC <- rowMeans(problem.genes.C3)/rowMeans(problem.genes.uninfected)

problem.genes$logFC <- log2(FC[problem.genes$gene])

plot(problem.genes$logFC, problem.genes$scran.C3.vs.uninfected.logFC)
plot(problem.genes$logFC, problem.genes$archr.C3.vs.uninfected)


