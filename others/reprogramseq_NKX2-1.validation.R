library(gp.sa)
library(gp.sa.diff)
library(DataSetDB)
library(ggplot2)
library(superheat)
library(pheatmap)
library(gridExtra)
library(ggplotify)
library(ggpubr)
library(epiregulon)
library(epiregulon.extra)
source("/gstore/project/lineage/analysis/function/functions.R")


#import data from DSDB
se <- getDatasetAsSE("DS000013689")
assays(se)$RPKM <- normalizedRPKM(assays(se)$count, rowData(se)$size)

se <- se[,order(se$SAMPLE_ID)]

se$celllines <- sapply(strsplit(se$INVENTORY_SAMPLE_NAME, split = "_"), "[", 2)
se$gene <- sapply(strsplit(se$INVENTORY_SAMPLE_NAME, split = "_"), "[", 3)
se$media <- sapply(strsplit(se$INVENTORY_SAMPLE_NAME, split = "_"), "[", 4)
se$dox <- sapply(strsplit(se$INVENTORY_SAMPLE_NAME, split = "_"), "[", 5)
se$dox <- gsub("\\+dox", "+ dox", se$dox)
se$treatment <- paste0(se$celllines,  "_", se$gene, "_", se$media,  "_",  se$dox)
se$treatment <- factor(se$treatment,
                       levels = c("FGC_eGFP_normal_no dox", "FGC_eGFP_normal_+ dox",
                                  "FGC_NKX2-1_normal_no dox", "FGC_NKX2-1_normal_+ dox",
                                  "FGC_NKX2-1_CDSS_no dox", "FGC_NKX2-1_CDSS_+ dox",
                                  "p53koE2_eGFP_normal_no dox", "p53koE2_eGFP_normal_+ dox",
                                  "p53koE2_NKX2-1_normal_no dox", "p53koE2_NKX2-1_normal_+ dox",
                                  "p53koE2_NKX2-1_CDSS_no dox", "p53koE2_NKX2-1_CDSS_+ dox")
)

################# plot individual genes

signature.list <- list()


# Load precalculated regulons
regulon.w <- readRDS("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/reprogram.seq.regulon.w.rds")
regulon.w <- regulon.w[regulon.w$tf == "NKX2-1",]
signature.list[["NKX2-1.epiregulon"]] <- DataFrame(
    gene_id = regulon.w$target,
    weights = regulon.w$weight)


desired_order <- c(16:18, 13:15, 4:6, 1:3, 10:12, 7:9, 34:36, 31:33, 22:24, 19:21, 28:30,25:27)




# Load bulk data

outdir <- "/gstore/project/ar_ligands/NE/rnaseq/NGS4677_NKX2.1_P53KO/OUTPUT/"
voom.output <- readRDS(paste0(outdir,"NKX2.1.voom.output.rds"))
bulkdiff <- voom.output[['FGC']][['Difference between `high` vs `low`']]
bulkdiff$Symbol <- rowData(se)$symbol[match(rownames(bulkdiff), rowData(se)$ID)]

# Add bulk LFC to regulons
regulon.w$bulkLFC <- bulkdiff$LogFC[match(regulon.w$target, bulkdiff$Symbol)]

pdf("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/bulk.vs.weight.pdf", height=5, width = 5)


ggscatter(data.frame(regulon.w), x = "bulkLFC", y = "weight.C3",
          color = "black", shape = 20, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                label.x = 3, label.y = 0
                                ),
          ylab = "regulon weight", xlab = "exp logFC", title = "Wilcox C3",
          xlim = c(-1,5),  ylim=c(0, 0.6))

ggscatter(data.frame(regulon.w), x = "bulkLFC", y = "weight.all",
          color = "black", shape = 20, size = 1, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE, # Add confidence interval
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          cor.coeff.args = list(method = "pearson", label.sep = "\n",
                                label.x = 3, label.y =0.03
                                ),
          ylab = "regulon weight", xlab = "exp logFC", title = "Wilcox all",
          xlim=c(-1,5), ylim=c(0, 0.25))

dev.off()




############
#NKX2-1 specific plots

heatmaplist <- list()
for (i in c("NKX2-1.epiregulon") ){

    genes_to_plot <- signature.list[[i]]["gene_id"]
    gene_index <- na.omit(unlist(match(genes_to_plot, rowData(se)$symbol)))
    expr <- assays(se)$RPKM[gene_index, colnames(se)[desired_order]]
    rownames(expr) <- rowData(se)$symbol[gene_index]

    # remove duplicates
    expr <- expr[!duplicated(rownames(expr)),]
    # remove NA
    expr <- na.omit(expr)
    # remove Inf
    expr <- expr[!is.infinite(rowSums(expr)),]
    # remove 0 variance
    expr <- expr[rowVars(expr) > 0.001,]

    annotation_col <- data.frame(Dox = se$dox, Gene = se$gene)
    rownames(annotation_col) <- colnames(se)
    annotation_col <- annotation_col[colnames(se)[desired_order],]

    #annotate genes
    annotation_rows <- data.frame(weight=regulon.w$weight[match(rownames(expr),regulon.w$target)])
    rownames(annotation_rows) <- rownames(expr)

    #order expr
    expr <- expr[order(annotation_rows$weight),]

    expr.add <- data.frame(colData(se)[desired_order,],t(expr))
    expr <- expr[, expr.add$Sample[which(expr.add$media == "normal" & expr.add$celllines == "FGC")]]
    expr.add <- expr.add[colnames(expr),]
    heatmaplist[[i]] <- as.grob(plotHeatmap(expr = data.frame(expr),
                                            labels = se$INVENTORY_SAMPLE_NAME[
                                                match(colnames(expr),se$Sample)],
                                            title = i,
                                            min = -2,
                                            max = 2,
                                            cluster_col = FALSE,
                                            cluster_row = FALSE,
                                            annotation_row = annotation_rows))




    plotIndividual(
        make.names(rownames(expr)),
        data.frame(expr.add),
        file = paste0("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/",
                      i,
                      ".individual.plots.pdf"
        ),
        facet_name = "media+celllines",
        x = "treatment",
        nrow = 2,
        ncol = 5,
        manual_colors = c("black", "grey","orange","red"))

}

heatmaplist.m <- marrangeGrob(heatmaplist, nrow = 1, ncol = 1)
ggsave("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/NKX2.1.heatmap.pdf", heatmaplist.m, width = 7, height = 8)

# calculate activity
rownames(se) <- rowData(se)$symbol
activity <- calculateActivity(se, exp_assay = "RPKM", regulon=regulon.w, mode = "weight")
se$NKX2.1 <- as.vector(activity)

plotIndividual(
    make.names(rownames(activity)),
    data.frame(colData(se[se$media == "normal" & se$celllines == "FGC",])),
    file = paste0("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/",
                  "activity.plots.pdf"
    ),
    facet_name = "media+celllines",
    x = "treatment",
    nrow = 2,
    ncol = 2,
    manual_colors = c("black", "grey","orange","red", "black", "grey","orange","red","black", "grey","orange","red")
)



# Validate with Tang et al
#cell lines and organoids
model <- read.csv("/gstore/project/lineage/manuscript/epiregulon/data/Tang_science/GSE181374_totalDf_rlog.csv", row.names = 1)
rownames(model) <- sapply(strsplit(rownames(model), split="\\."), "[",1)

metadata <- read.csv("/gstore/project/lineage/manuscript/epiregulon/data/Tang_science/patient_class.csv")
rownames(metadata) <- metadata$Sample.alt
metadata <- metadata[colnames(model),]
samples <- intersect(metadata$Sample.alt, colnames(model))

model.se <- SingleCellExperiment(assays=list(counts=as.matrix(model[,samples])), colData = metadata[samples,])


library(EnsDb.Hsapiens.v86)

ensembl.genes <- rownames(model)
symbols <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

model.se.filtered <- model.se[symbols$GENEID,]
rowData(model.se.filtered) <- symbols
rownames(model.se.filtered) <- symbols$SYMBOL
model.se.filtered.nodup <- model.se.filtered[which(!duplicated(rownames(model.se.filtered))),]


heatmapplot <- plotHeatmapRegulon(sce=model.se.filtered,
                                  tfs=c("GATA6","NKX2-1"),
                                  regulon=regulon.w,
                                  regulon_cutoff=0.01,
                                  downsample=2000,
                                  cell_attributes="Subtype",
                                  col_gap="Subtype",
                                  exprs_values="counts",
                                  name="regulon heatmap"
                                  #column_title_rot = 45
                                  )

pdf("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/Tang_GATA6_heatmap.pdf")

print(heatmapplot)
dev.off()

regulon.w <- readRDS("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/reprogram.seq.regulon.w.rds")

score.combine <- calculateActivity(expMatrix = model.se.filtered,
                                   regulon = regulon.w,
                                   mode = "weight",
                                   method = "weightedMean",
                                   exp_assay = "counts",
                                   normalize = FALSE)
plotActivityViolin(score.combine, tf = c("NKX2-1","GATA6"), clusters=model.se.filtered$Subtype)
activity_mat <- cbind(colData(model.se.filtered), as.matrix(t(score.combine[c("NKX2-1","GATA6"),])))

gex_mat <- cbind(colData(model.se.filtered), t(counts(model.se.filtered)[c("NKX2-1","GATA6"),]))


pdf("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/GATA6_activity_CRPC.pdf")
ggplot(data.frame(activity_mat), aes(x=Subtype, y=GATA6, fill=Subtype)) +  geom_violin() +
    geom_point() + theme_classic() + ylab("GATA6 activity")
ggplot(data.frame(gex_mat), aes(x=Subtype, y=GATA6, fill=Subtype)) +  geom_violin() +
    geom_point() + theme_classic() + ylab("GATA6 expression")
dev.off()


library(genomitory)
prostate <- getFeatureSetCollection("GMTY205:analysis/prostate.gmt.bz2@REVISION-4")
names(prostate) <- prostate@elementMetadata@listData[["name"]]

gene.list.regulon <- epiregulon:::genesets2regulon(prostate)
score.combine.gene.list <- calculateActivity(expMatrix = GeneExpressionMatrix,
                                   regulon = gene.list.regulon,
                                   mode = "weight",
                                   method = "weightedMean",
                                   exp_assay = "counts",
                                   normalize = FALSE)

pdf("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/genesignature.geneIntegrationMatrix.heatmap.plots.pdf",
    width = 10, height = 4 )
plotActivityViolin(activity_matrix = score.combine.gene.list,
                   clusters = GeneExpressionMatrix$hash_assignment,
                   tf = c("Tang_CRPC_NE","Tang_CRPC_SCL"),boxplot = TRUE)
plotHeatmapRegulon(sce= GeneExpressionMatrix,
                   tfs=c("Tang_CRPC_NE","Tang_CRPC_SCL"),
                   regulon=gene.list.regulon,
                   regulon_cutoff=0,
                   downsample=4000,
                   cell_attributes="Clusters",
                   col_gap="Clusters",
                   exprs_values="normalizedCounts",
                   name="regulon heatmap",
                   column_title_rot = 45,
                   color_breaks = c(-2, 0, 2),

)

dev.off()

