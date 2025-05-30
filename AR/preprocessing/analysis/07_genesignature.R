library(ArchR)
library(epiregulon)
library(scater)
library(epiregulon.extra)
library(epiregulon.archr)
#checkdata

archR_project_path <- "OUTPUT/ArchRProject/"
proj.all <- loadArchRProject(path = archR_project_path, showLogo = TRUE)
getAvailableMatrices(proj.all)



###############################
# load gene expression matrix
GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "GeneExpressionMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = 1,
    logFile = createLogFile("getMatrixFromProject")
)

GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts",
                                        transform = TRUE,
                                        transform_method = "log")
rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name

GeneExpressionMatrix$TEST_ARTICLE <- factor(as.character(GeneExpressionMatrix$TEST_ARTICLE),
                                            levels = c("DMSO", "Enza", "ARV110","A9690"))




# Add embeddings
reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- getEmbedding(ArchRProj = proj.all,
                                                                  embedding = "UMAP_Combined",
                                                                  returnDF = TRUE)[colnames(GeneExpressionMatrix), ]

reducedDim(GeneExpressionMatrix, "UMAP_ATAC") <- getEmbedding(ArchRProj = proj.all,
                                                              embedding = "UMAP_ATAC",
                                                              returnDF = TRUE)[colnames(GeneExpressionMatrix), ]


reducedDim(GeneExpressionMatrix, "UMAP_RNA") <- getEmbedding(ArchRProj = proj.all,
                                                             embedding = "UMAP_RNA",
                                                             returnDF = TRUE)[colnames(GeneExpressionMatrix), ]

# plot
GeneExpressionMatrix <- GeneExpressionMatrix[, which(!is.na(GeneExpressionMatrix$TEST_ARTICLE))]

options(ggrastr.default.dpi=300)

library(RColorBrewer)

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

pdf("OUTPUT/ArchRProject/Plots/umap.pdf", width = 5, height =4)
plotReducedDim(GeneExpressionMatrix,
               dimred="UMAP_ATAC",
               colour_by = "Cell",
               text_by = "Cell",
               point_size=0.3,
               point_alpha=0.3,
               rasterise = TRUE)

colourCount <- length(unique(GeneExpressionMatrix$TEST_ARTICLE))
plotReducedDim(GeneExpressionMatrix,
               dimred="UMAP_ATAC",
               colour_by = "TEST_ARTICLE",
               point_size=0.3,
               point_alpha=0.3,
               rasterise = TRUE,
               text_by="Cell") +
    scale_color_manual(values=c("grey","red","blue","orange","darkgreen","navyblue"))

plotReducedDim(GeneExpressionMatrix[,which(GeneExpressionMatrix$Cell %in% c("LNCaP", "VCaP","DU145","H660") &
                                        GeneExpressionMatrix$TEST_ARTICLE %in% c("DMSO", "Enza", "ARV110","A9690"))],
               dimred="UMAP_ATAC",
               colour_by = "TEST_ARTICLE",
               point_size=0.3,
               point_alpha=0.3,
               rasterise = TRUE,
               text_by="Cell") +
    scale_color_manual(values=c("grey","red","blue","orange"))

dev.off()


############################
# calculate signature scores
library(genomitory)
prostate <- getFeatureSetCollection("GMTY205:analysis/prostate.gmt.bz2@REVISION-4")
names(prostate) <- prostate@elementMetadata@listData[["name"]]


outdir <- archR_project_path

signature.scores <- calculateActivity(GeneExpressionMatrix, genesets = as.list(prostate), exp_assay = "normalizedCounts")

prostate.regulon <- epiregulon:::genesets2regulon(prostate)

# plot violin
pdf("OUTPUT/ArchRProject/Plots/signatures.AR.pdf", width = 12, height = 12)

for (cell in unique(GeneExpressionMatrix$Cell)){
    selected_cells <- which(GeneExpressionMatrix$Cell == cell)


    violinplot <- plotActivityViolin(activity_matrix=signature.scores[, selected_cells],
                                     tf=rownames(signature.scores[, selected_cells]),
                                     clusters=GeneExpressionMatrix$TEST_ARTICLE[selected_cells],
                                     #facet_grid_variable=GeneExpressionMatrix$Cell[selected_cells],
                                     legend.label="signature score",
                                     ncol = 5,
                                     colors = c("grey","red","blue","orange"),
                                     title = cell,
                                     boxplot = TRUE
    )
    print(violinplot)

    violinplot <- plotActivityViolin(activity_matrix=assay(GeneExpressionMatrix, "logcounts")[, selected_cells],
                                     tf=c("AR"),
                                     clusters=GeneExpressionMatrix$TEST_ARTICLE[selected_cells],
                                     #facet_grid_variable=GeneExpressionMatrix$Cell[selected_cells],
                                     legend.label="gene expression",
                                     ncol = 5, nrow=4,
                                     colors = c("grey","red","blue","orange"),
                                     title = cell,
                                     boxplot = TRUE
    )
    print(violinplot)

}

plotActivityDim(sce = GeneExpressionMatrix,
                activity_matrix = signature.scores,
                tf = c("AR_Hieronymus", "NE_Dong", "Tang_CRPC_AR",
                       "Tang_CRPC_NE", "Tang_CRPC_SCL", "basal","luminal"),
                dimtype = "UMAP_ATAC",
                label = "Cell",
                point_size = 0.5,
                ncol = 3, rasterise=TRUE)



dev.off()


pdf("OUTPUT/ArchRProject/Plots/signatures.AR.heatmap.pdf", height = 6)

for (cell in unique(GeneExpressionMatrix$Cell)){

    selected_cells <- which(GeneExpressionMatrix$Cell == cell & GeneExpressionMatrix$TEST_ARTICLE %in% c("DMSO","Enza","ARV110","A9690"))
    heatmapplot <- plotHeatmapRegulon(sce=GeneExpressionMatrix[, selected_cells],
                                      exprs_values = "normalizedCounts",
                                      regulon = prostate.regulon,
                                      tf = c("AR_10_Bluemn","AR_Dong", "AR_Hieronymus"),
                                      color_breaks = c(-0.5, 0, 0.5),
                                      colors = c("darkblue", "black", "yellow"),
                                      cell_attributes="TEST_ARTICLE",
                                      col_gap = "TEST_ARTICLE",
                                      raster_quality=10,
                                      downsample = 2000,
                                      column_col=list(TEST_ARTICLE = c("DMSO" = "grey",
                                                                       "Enza" = "red",
                                                                       "ARV110" = "blue",
                                                                       "A9690"= "orange")),
                                      row_col=list(tf = c("AR_Hieronymus" = "brown",
                                                          "AR_Dong" = "yellow",
                                                          "AR_10_Bluemn" = "pink")),
                                      column_title = cell,
                                      name = "z-score")
    print(heatmapplot)

}

dev.off()

#######################################
# compare activity vs signatures


df <- t(signature.scores)

# import the AR activity from chip

MDA.score.combine <- readRDS("OUTPUT/ArchRProject/Epiregulon/score.combine.MDA.chip.rds")
LNCaP.score.combine <- readRDS("OUTPUT/ArchRProject/Epiregulon/score.combine.LNCaP.chip.rds")
VCaP.score.combine <- readRDS("OUTPUT/ArchRProject/Epiregulon/score.combine.VCaP.chip.rds")

df <- as.matrix(df)
df <- data.frame(df)
df[,"epiregulon"] <- NA
df[colnames(MDA.score.combine), "epiregulon"] <- as.vector(MDA.score.combine["AR",])
df[colnames(LNCaP.score.combine), "epiregulon"] <- as.vector(LNCaP.score.combine["AR",])
df[colnames(VCaP.score.combine), "epiregulon"] <- as.vector(VCaP.score.combine["AR",])
df$TEST_ARTICLE <- colData(GeneExpressionMatrix[,rownames(df)])$TEST_ARTICLE
df$Cell <- colData(GeneExpressionMatrix[,rownames(df)])$Cell

df <- df[df$TEST_ARTICLE %in% c("DMSO", "Enza", "ARV110",  "A9690" ),]
library(ggplot2)

pdf("OUTPUT/ArchRProject/Epiregulon/signature.vs.activity.pdf", width=4, height=2.5)
ggplot(df[df$Cell == "LNCaP",], aes(x=AR_10_Bluemn, y=epiregulon, color=TEST_ARTICLE)) +
    geom_point(size=1) +
    scale_color_manual(values=c("grey","red","blue","orange","darkgreen","navyblue")) +
    theme_classic() + ggtitle("LNCaP")

ggplot(df[df$Cell == "VCaP",], aes(x=AR_10_Bluemn, y=epiregulon, color=TEST_ARTICLE)) +
    geom_point(size=1) +
    scale_color_manual(values=c("grey","red","blue","orange","darkgreen","navyblue")) +
    theme_classic() + ggtitle("VCaP")

ggplot(df[df$Cell == "MDA",], aes(x=AR_10_Bluemn, y=epiregulon, color=TEST_ARTICLE)) +
    geom_point(size=1) +
    scale_color_manual(values=c("grey","red","blue","orange","darkgreen","navyblue")) +
    theme_classic() + ggtitle("MDA")

dev.off()
