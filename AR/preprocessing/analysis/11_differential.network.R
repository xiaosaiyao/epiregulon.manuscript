library(scran)
library(scater)
library(epiregulon)
library(epiregulon.extra)
#devtools::load_all("/gstore/project/lineage/xiaosai/epiregulon.extra")
library(epiregulon.archr)
library(ArchR)
library(BiocParallel)
library(igraph)
library(ggrepel)

#checkdata
archR_project_path <- "OUTPUT/ArchRProject"
proj.all <- loadArchRProject(path = archR_project_path, showLogo = TRUE)

celllines <- unique(proj.all$Cell)


###############################
# load gene expression matrix
GeneExpressionMatrix <- getMatrixFromProject(
    ArchRProj = proj.all,
    useMatrix = "GeneExpressionMatrix",
    useSeqnames = NULL,
    verbose = TRUE,
    binarize = FALSE,
    threads = getArchRThreads(),
    logFile = createLogFile("getMatrixFromProject")
)


GeneExpressionMatrix <- ArchRMatrix2SCE(GeneExpressionMatrix,
                                        rename="normalizedCounts",
                                        transform=TRUE,
                                        transform_method="log")

rownames(GeneExpressionMatrix) <- rowData(GeneExpressionMatrix)$name
GeneExpressionMatrix$TEST_ARTICLE <- factor(as.character(GeneExpressionMatrix$TEST_ARTICLE),
                                            levels = c("DMSO", "Enza","ARV110", "A9690") )


# Add embeddings
reducedDim(GeneExpressionMatrix, "UMAP_ATAC") <- getEmbedding(ArchRProj = proj.all,
                                                              embedding = "UMAP_ATAC",
                                                              returnDF = TRUE)[colnames(GeneExpressionMatrix), ]


##################################
# differential activity
markers <- list()
diff_act_plot <- list()
regulon.w <- list()
drugs <-  c("Enza", "ARV110", "A9690")


logFC_cutoff <- 0.3
sig_cutoff <- 0.05
topgenes <- 10
regulon.size.filter <- 30

findDifferentialActivityCustom <- function(activity_matrix,
                                     clusters,
                                     test.type = "t",
                                     pval.type = "some",
                                     direction = c("any","up","down"),
                                     groups = deprecated(),
                                     logvalues = TRUE,
                                     ...){

    if(lifecycle::is_present(groups)){
        lifecycle::deprecate_warn( "1.0.0", "findDifferentialActivity(groups)",
                                   "findDifferentialActivity(clusters)")
        clusters <- groups
    }

    direction <- match.arg(direction)
    activity_matrix <- stats::na.omit(as.matrix(activity_matrix))
    tf_markers <- scran::findMarkers(activity_matrix, clusters, test.type=test.type,
                                     pval.type=pval.type, direction=direction, ...)


    if (!isTRUE(logvalues)){

        for (cluster in unique(clusters)) {

            # replace logFC with diff
            colnames(tf_markers[[cluster]]) <- gsub("logFC", "diff", colnames(tf_markers[[cluster]]))

            # calculate logFC
            current <- rowMeans(activity_matrix[, which(clusters == cluster)])
            rest <- rowMeans(activity_matrix[, which(clusters != cluster)])
            summary.logFC <- log2(current/rest)
            tf_markers[[cluster]][,"summary.logFC"] <- summary.logFC[rownames(tf_markers[[cluster]])]

            tf_markers[[cluster]][,"summary.diff"] <- NULL

            # loop through all other comparisons
            for (othercluster in setdiff(unique(clusters), cluster)){
                othercluster.mean <- rowMeans(activity_matrix[, which(clusters == othercluster)])
                othercluster.logFC <- log2(current/othercluster.mean)
                tf_markers[[cluster]][,paste0("logFC.", make.names(othercluster))] <- othercluster.logFC[rownames(tf_markers[[cluster]])]
            }
        }
    }
    return(tf_markers)

}

# find the top differential TFs
for (cell in celllines ){

    # select the cell line
    selected <- which(GeneExpressionMatrix$Cell == cell)
    GeneExpressionMatrix.select <- GeneExpressionMatrix[, selected]
    regulon.w.cell <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.", cell, ".rds"))
    score.combine <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cell, ".rds"))
    score.combine <- score.combine[,colnames(GeneExpressionMatrix.select)]

    markers[[cell]] <- findDifferentialActivityCustom(activity_matrix = score.combine,
                                                clusters = GeneExpressionMatrix.select$TEST_ARTICLE,
                                                pval.type = "any", direction = "any", log.p=TRUE, logvalues = FALSE)

    markers.sig <- getSigGenes(markers[[cell]], topgenes = topgenes, direction = "any" )


    regulon.w[[cell]] <- addLogFC(expMatrix = GeneExpressionMatrix.select,
                                  clusters = GeneExpressionMatrix.select$TEST_ARTICLE,
                                  regulon = regulon.w.cell,
                                  pval.type = "some",
                                  assay.type = "logcounts",
                                  logFC_condition = drugs,
                                  logFC_ref = "DMSO")


    # filter regulons based on logFC
    Enza_index <- which(abs(regulon.w[[cell]]$`Enza.vs.DMSO.logFC`) > logFC_cutoff & regulon.w[[cell]]$`Enza.vs.DMSO.FDR` < sig_cutoff)
    ARV110_index <- which(abs(regulon.w[[cell]]$`ARV110.vs.DMSO.logFC`) > logFC_cutoff & regulon.w[[cell]]$`ARV110.vs.DMSO.FDR` < sig_cutoff)
    A9690_index <- which(abs(regulon.w[[cell]]$`A9690.vs.DMSO.logFC`) > logFC_cutoff & regulon.w[[cell]]$`A9690.vs.DMSO.FDR` < sig_cutoff)

    combined_index <- unique(c(Enza_index, ARV110_index, A9690_index))
    regulon.w.filtered <- regulon.w[[cell]][combined_index,]

    # calculate altered regulon size
    regulon.size <- table(regulon.w.filtered$tf)

    # filter top TFs based on altered regulon size
    markers.sig$regulon.size <-  as.vector(regulon.size[markers.sig$tf])
    markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]


    # plot
    # tryCatch(bubble <- plotBubble(activity_matrix = score.combine,
    #                      tf = unique(markers.sig$tf),
    #                      clusters = GeneExpressionMatrix.select$TEST_ARTICLE,
    #                      bubblesize = "FDR"))
    # set.seed(1010)
    # tfplot <- plotHeatmapActivity(activity=score.combine,
    #                               sce=GeneExpressionMatrix.select,
    #                               tfs=unique(markers.sig$tf),
    #                               downsample=1000,
    #                               cell_attributes="TEST_ARTICLE",
    #                               col_gap="TEST_ARTICLE",
    #                               name = "activity",
    #                               column_title_rot = 45,
    #                               color_breaks = c(-2, 0, 2),
    #                               colors = c("blue", "white", "red"),
    #                               cluster_rows = FALSE)
    #
    #
    # set.seed(1010)
    # gexplot <- plotHeatmapActivity(activity= assay(GeneExpressionMatrix.select, "logcounts"),
    #                                sce=GeneExpressionMatrix.select,
    #                                tfs=unique(markers.sig$tf),
    #                                downsample=1000,
    #                                cell_attributes="TEST_ARTICLE",
    #                                col_gap="TEST_ARTICLE",
    #                                name = "gex",
    #                                column_title_rot = 45,
    #                                color_breaks = c(-2, 0, 2),
    #                                colors = c("blue", "white", "red"),
    #                                cluster_rows = FALSE)

    #pdf(paste0("OUTPUT/ArchRProject/Epiregulon/diffactivity.", cell, ".pdf"), height = min(20, round(length(unique(markers.sig$tf))*0.15)+2))
    #print(bubble)
    #print(tfplot)
    #print(gexplot)
    #dev.off()
}



# plot regulon size vs differential activity
ylim.max <- c(LNCaP=300, VCaP=200, MDA=100, DU145=500, H660=400, `22Rv1` = 300)

regulon.w.filtered <- list()

for (cell in celllines){
    for (drug in drugs){
        message(drug, cell)
        rank_list <- markers[[cell]][[drug]]
        rank_list <- rank_list[rank_list$log.FDR < -2,]

        if (nrow(rank_list)>0) {
            rank_list <- rank_list[order(rank_list$logFC.DMSO),]
            rank_list$rank <- 1:nrow(rank_list)
            rank_list$label <- ""

            # add regulon size
            filtered_index <- as.numeric(which(abs(regulon.w[[cell]][[paste0(drug,".vs.DMSO.logFC")]]) > logFC_cutoff & regulon.w[[cell]][[paste0(drug,".vs.DMSO.FDR")]] < sig_cutoff))
            regulon.w.filtered[[cell]][[drug]] <- regulon.w[[cell]][filtered_index,]

            freq <- table(regulon.w.filtered[[cell]][[drug]]$tf)
            rank_list$count <- as.vector(freq)[match(rownames(rank_list), names(freq))]


            top_logFC_genes <- intersect(rownames(head(rank_list,10)), rownames(rank_list)[rank_list$count>regulon.size.filter])
            bottom_logFC_genes <- intersect(rownames(tail(rank_list,10)), rownames(rank_list)[rank_list$count>regulon.size.filter])
            most_freq_TF <- head(rownames(rank_list)[order(rank_list$count, decreasing = TRUE)],10)
            custom_genes <- c("AR", "TEAD1") #,"NKX3-1","FOXA1","HOXB13","ERG","TEAD1")
            genes_label <- c(custom_genes, top_logFC_genes, bottom_logFC_genes, most_freq_TF)


            idx.match <- na.omit(match(genes_label, rownames(rank_list)))
            rank_list$label[idx.match] <- rownames(rank_list)[idx.match]

            rank_list <- rank_list[order(rank_list$label),]
            diff_act_plot[[paste0(cell,drug)]] <- ggplot(rank_list, aes(logFC.DMSO, count, label = label)) +
                geom_point() +
                geom_text_repel(max.overlaps = Inf) +
                geom_point(color = ifelse(rank_list$label == "", "grey50", "red")) +
                ggtitle(paste(cell, "after", drug, "treatment total activity in", "regulon")) +
                ylab("Regulon size") + xlab("logFC in activity with respect to DMSO") +
                xlim(-1.5,1) + ylim(0,ylim.max[cell]) +
                geom_vline(xintercept = 0, linetype="dashed", color = "gray", size=1) + theme_classic()
        }
    }

}


library(gridExtra)

pdf("OUTPUT/ArchRProject/Epiregulon/diff.activity.pdf", width=4, height = 4)
marrangeGrob(grobs = diff_act_plot, ncol=1, nrow=1)
dev.off()


#### Geneset enrichment
H <- EnrichmentBrowser::getGenesets(org = "hsa",
                                    db = "msigdb",
                                    cat = "H",
                                    gene.id.type = "SYMBOL",
                                    cache = FALSE)
C6 <- EnrichmentBrowser::getGenesets(org = "hsa",
                                     db = "msigdb",
                                     cat = "C6",
                                     gene.id.type = "SYMBOL",
                                     cache = FALSE)

gs <- c(H,C6)
gs.list <- do.call(rbind,lapply(names(gs),
                                function(x) {data.frame(gs=x, genes=gs[[x]])}))



pdf("OUTPUT/ArchRProject/Epiregulon/gsea.pdf", width=12, height=4)

# LNCaP

enrichresults <- regulonEnrich(TF = c("AR","SMARCA4"),
                               regulon = regulon.w.filtered[["LNCaP"]][["A9690"]],
                               weight = "weight",
                               weight_cutoff = 0,
                               genesets = gs.list)
gsea_plot <- enrichPlot(results = enrichresults, title = "LNCaP A9690")
print(gsea_plot)

# MDA

enrichresults <- regulonEnrich(TF = c("AR","SMARCA4"),
                               regulon = regulon.w.filtered[["MDA"]][["A9690"]],
                               weight = "weight",
                               weight_cutoff = 0,
                               genesets = gs.list)
gsea_plot <- enrichPlot(results = enrichresults, title = "MDA A9690")
print(gsea_plot)


# DU145
enrichresults <- regulonEnrich(TF = c("NR3C1","SMARCA4"),
                               regulon = regulon.w.filtered[["DU145"]][["A9690"]],
                               weight = "weight",
                               weight_cutoff = 0,
                               genesets = gs.list)
gsea_plot <- enrichPlot(results = enrichresults, title = "DU145 A9690")
print(gsea_plot)

# H660

enrichresults <- regulonEnrich(TF = c("TCF12","SMARCA4"),
                               regulon = regulon.w.filtered[["H660"]][["A9690"]],
                               weight = "weight",
                               weight_cutoff = 0,
                               genesets = gs.list)
gsea_plot <- enrichPlot(results = enrichresults, title = "H660 A9690")
print(gsea_plot)
dev.off()
