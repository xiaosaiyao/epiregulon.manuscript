devtools::load_all("/gstore/project/lineage/xiaosai/epiregulon.extra")
library(epiregulon)
library(dsassembly)
library(ArchR)
library(scran)
library(ggrepel)
topgenes <- 20
regulon.size.filter <- 35
part_no_all <- c("part_1","part_2","part_5","part_8")
part_no_all <- "part_1"

topgenes <- c(part_1 = 200, part_2=20, part_5=20, part_8=20)

#zipped file /gstore/project/abs_datasets/data/No_866_ATAC/2.000/No_866_ATAC.tar.gz
for (part_no in part_no_all){
    # load activity matrix
    score.combine <- readRDS(paste0("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/", part_no,".score.combine.rds"))

    # import archr project
    archR_project_path <- paste0("/gstore/scratch/u/yaox19/TGI/No_866_ATAC1/", part_no, "/output_TSS4_Frags1000")
    archr.proj <- loadArchRProject(path = archR_project_path, showLogo = FALSE, force=TRUE)

    # import geneIntegrationMatrix
    sce <- dsassembly::getExperiment('DS000016526', part_no)

    # transfer clusterName info from sce onto archR based on predicted cells
    archr.proj$clusterName <- colData(sce)[archr.proj$predictedCell,"clusterName"]
    cluster_info <- as.vector(archr.proj@cellColData[colnames(score.combine), "clusterName"])

    # differential activity
    markers  <- findDifferentialActivity(activity_matrix = score.combine,
                                         clusters = cluster_info ,
                                         pval.type = "some",
                                         direction = "up",
                                         log.p=TRUE)

    markers.sig <- getSigGenes(markers, direction = "up", topgenes =topgenes[part_no])

    plotBubble(activity_matrix = score.combine,
               tf = markers.sig$tf,
               clusters = cluster_info ,
               bubblesize = "summary.logFC",
               direction = "up")


    # import regulon
    regulon.w <- readRDS(paste0("/gstore/project/lineage/kidney/DS000016526_ccRCC/OUTPUT/", part_no, ".pruned.regulon.w.wilcox.rds"))

    # prepare a downsampled gene expression matrix for computing logFC of target genes
    sce <- logNormCounts(sce)
    sce <- sce[, sce$clusterName %in% unique(cluster_info)]
    set.seed(1010)
    subsample.idx <- sample(1:ncol(sce), 10000, replace = FALSE)
    sce.test  <- sce [, subsample.idx]
    rownames(sce.test) <- rowData(sce.test)$symbol

    # add logFC onto regulons
    regulon.w <- addLogFC(expMatrix = sce.test,
                          clusters = sce.test$clusterName,
                          regulon = regulon.w,
                          sig_type = "FDR",
                          logFC_condition="Tumor")


    # filter regulons based on logFC
    regulon.w.filtered <- regulon.w[which(abs(regulon.w$Tumor.vs.rest.logFC) > 0.5 & regulon.w$Tumor.vs.rest.FDR < 0.05), ]

    # calculate altered regulon size
    regulon.size <- table(regulon.w.filtered$tf)

    # filter top TFs based on altered regulon size
    markers.sig$regulon.size <-  as.vector(regulon.size[markers.sig$tf])
    markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]


    # plot
    bubble <- plotBubble(activity_matrix = score.combine,
                         tf = markers.sig$tf,
                         clusters = cluster_info ,
                         bubblesize = "summary.logFC", direction = "up")

    # assemble a sce with colData info
    colData <- archr.proj@cellColData[colnames(score.combine),]
    sce.dummy <- SingleCellExperiment(colData=colData)
    tfplot <- plotHeatmapActivity(activity=score.combine,
                                  sce=sce.dummy,
                                  tfs=unique(markers.sig$tf),
                                  downsample=1000,
                                  cell_attributes="clusterName",
                                  col_gap="clusterName",
                                  name = "activity",
                                  column_title_rot = 45,
                                  color_breaks = c(-2, 0, 2),
                                  colors = c("blue", "white", "red"),
                                  cluster_rows = FALSE)



    gex <- assay(sce.test, "logcounts")
    gex <- as(gex,"dgCMatrix")
    gexplot <- plotHeatmapActivity(activity=gex,
                                   sce=sce.test,
                                   tfs=unique(markers.sig$tf),
                                   downsample=1000,
                                   cell_attributes="clusterName",
                                   col_gap="clusterName",
                                   name = "gex",
                                   column_title_rot = 45,
                                   color_breaks = c(-2, 0, 2),
                                   colors = c("blue", "white", "red"),
                                   cluster_rows = FALSE)

    pdf(paste0("OUTPUT/activity.", part_no, ".pdf"), height = min(20, round(length(unique(markers.sig$tf))*0.2)+2))
    print(bubble)
    print(tfplot)
    print(gexplot)
    dev.off()
}
