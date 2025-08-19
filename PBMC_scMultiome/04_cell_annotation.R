library(ArchR)
library(scater)
library(SingleCellExperiment)
library(epiregulon.archr)

proj <- loadArchRProject("saved_project")
GeneExpressionMatrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")

# Loading reference data with Ensembl annotations.
library(celldex)
bpe_data <- BlueprintEncodeData()

expr_assay <- assays(GeneExpressionMatrix)[[1]]
rownames(expr_assay) <- rowData(GeneExpressionMatrix)$name

# Performing predictions
library(SingleR)
predictions <- SingleR(test=expr_assay, assay.type.test=1,
                       ref=bpe_data, labels=bpe_data$label.fine) # could be switched to label.main

proj$cell_type_SingleR <- predictions$labels

saveArchRProject(ArchRProj = proj, outputDirectory = "saved_project", load = FALSE)

proj <- loadArchRProject("/gstore/project/epigen/PBMC/saved_project")
GeneExpressionMatrix <- getMatrixFromProject(proj, "GeneExpressionMatrix")
GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts")
reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- getEmbedding(proj, embedding = "UMAP_Combined")
GeneExpressionMatrix$cell_type[is.na(GeneExpressionMatrix$cell_type)] <- "unknown"
plotReducedDim(GeneExpressionMatrix, dimred="UMAP_Combined", colour_by = "Clusters", text_by = "Clusters")
plotReducedDim(GeneExpressionMatrix, dimred="UMAP_Combined", colour_by = "cell_type_SingleR", text_by = "cell_type_SingleR")

clusters <- proj$Clusters
manual_annotation <- rep(NA, length(clusters))
manual_annotation[clusters %in% "C6"] <- "Naive CD4+ T"
manual_annotation[clusters %in% c("C13")] <- "CD14+ Mono"
manual_annotation[clusters %in% "C12"] <- "Monocytes"
manual_annotation[clusters %in% c("C2", "C3")] <- "B"
manual_annotation[clusters %in% "C4"] <- "Memory CD8+ T"
manual_annotation[clusters %in% "C14"] <- "FCGR3A+ Mono"
manual_annotation[clusters %in% "C5"] <- "NK"
manual_annotation[clusters %in% c("C9", "C10")] <- "Memory CD4+ T"
manual_annotation[clusters %in% c("C7", "C8")] <- "Naive CD8+ T"
manual_annotation[clusters %in% c("C1")] <- "DC"
proj$cell_type <- manual_annotation
saveArchRProject(ArchRProj = proj, outputDirectory = "saved_project", load = FALSE)

proj <- loadArchRProject("saved_project")
GeneExpressionMatrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")

GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix, rename="normalizedCounts")
reducedDim(GeneExpressionMatrix, "UMAP_Combined") <- getEmbedding(proj, embedding = "UMAP_Combined")

plot_df <- reducedDim(GeneExpressionMatrix)
plot_df$cell_type <- GeneExpressionMatrix$cell_type
centers <- lapply(split(reducedDim(GeneExpressionMatrix), GeneExpressionMatrix$cell_type), function(x) colMeans(x))
centers <- as.data.frame(do.call(rbind, centers))
centers$cell_type <- as.factor(rownames(centers))
colnames(centers) <- gsub("LSI_Combined#","",colnames(centers))
colnames(plot_df) <- gsub("LSI_Combined#","",colnames(plot_df))


dev.new()
pdf("/gstore/project/epigen/PBMC/Plots/PBMC_UMAP_ArchR.pdf", width = 7, height = 6)
ggplot(plot_df, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, color = cell_type))+
  geom_point(alpha = 0.6)+
  scale_color_manual(values= c("CD14+ Mono" ="#33A02C","Memory CD4+ T" = "#E31A1C", "B"  =   "#1F78B4",  "Monocytes" ="#B2DF8A","Naive CD8+ T" = "#FDBF6F" , "FCGR3A+ Mono" = "olivedrab4" ,
                               "DC"= "darkmagenta","Naive CD4+ T" = "#FB9A99", "Memory CD8+ T" =  "#FF7F00", "NK"  ="#CAB2D6",
                               "unknown" = "grey"))+
  geom_text(data=centers, aes(x = UMAP_Dimension_1, y = UMAP_Dimension_2, label = cell_type), color = "black", vjust = -0.5)+
  theme_classic()
dev.off()


