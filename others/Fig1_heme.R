library(ArchR)
library(epiregulon)

## load archR project
archR_project_path <- "/gstore/project/lineage/sam/heme_GRN/OUTPUT"
proj <- loadArchRProject(archR_project_path)

# read sce data
sce <- readRDS("/gstore/project/lineage/sam/heme_GRN/scRNA-Granja-2019.rds")

## construct regulon
grl <- getTFMotifInfo(genome = "hg19")
head(grl)

p2g <- calculateP2G(ArchR_path = archR_project_path)
head(p2g)

overlap <- addTFMotifInfo(p2g, grl, archR_project_path = archR_project_path)
head(overlap)

regulon.full <- getRegulon(p2g, overlap, aggregate=TRUE)
head(regulon.full)
saveRDS(regulon.full, "OUTPUT/heme.regulon.full.rds")


## select for TFs

TFs <- c("SPI1","PRDM1")

regulon <- regulon.full %>% dplyr::filter(tf %in% TFs)
nrow(regulon)




regulon.w <- addWeights(regulon=regulon,
                        sce=sce,
                        cluster_factor="BioClassification",
                        block_factor=NULL,
                        corr=TRUE,
                        MI=FALSE,
                        BPPARAM=BiocParallel::MulticoreParam())
head(regulon.w)

score.combine <- calculateActivity(sce, regulon.w, "weight", method="weightedMean",
                                   assay = "logcounts")
head(score.combine[,1:10])

da_list <- findDifferentialActivity(score.combine, sce$BioClassification, pval.type="some",
                                    direction="up", test.type= "t")
markers <- getSigGenes(da_list, fdr_cutoff = 0.05)
head(markers)

pdf("OUTPUT/Fig2.heme.umap.pdf", width = 12)

subsample <- sample(1:ncol(sce),size = 5000, replace = FALSE)
activityplot <- plotActivityDim(sce[, subsample], score.combine[, subsample], TFs, "TSNE", combine = T, point_size = 0.5)
ggsave(filename = "OUTPUT/activityplot.pdf",
       plot = activityplot, device = "pdf",
       units = "in", width=12, height=9, dpi=72)

plotActivityDim(sce, assay(sce), TFs, "TSNE", combine = T, limit = c(0,2),
                colors = c("gray","blue"), point_size = 0.5 )
dev.off()

plotActivityViolin(score.combine, TFs, sce$BioClassification)
plotActivityViolin(assay(sce), TFs, sce$BioClassification)
plotBubble(as.matrix(score.combine), TFs,
           sce$BioClassification, bubblesize = "FDR")
dev.off()
