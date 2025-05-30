library(maw.utils)
library(ArchR)
library(gridExtra)
library(DropletUtils)
library(SingleCellExperiment)
library(ggplot2)
library(scater)

# retrieve arcseq paths
arcseq_info <- getFireDBResourceSetInfo(frs.id = "FRS20569")$file

# retrieve HTO paths
HTO_frs <- c("FRS17161", "FRS17163", "FRS17165", "FRS17171")
HTO_info <- list()
for (id in HTO_frs){
    uri <- getFireDBResourceSetInfo(frs.id = id)$file$uri
    HTO_info[[id]] <- list.files(uri, recursive=TRUE, full.names = TRUE, pattern = ".hashing.csv")
}
HTO_info <- unname(do.call(c, HTO_info))
LIBID_SAMID <- unlist(lapply(strsplit(HTO_info, split = "/"), "[", 10))
SAMID  <- unlist(lapply(strsplit(LIBID_SAMID, split = "_"), "[", 2))

HTO_uri <- data.frame(SAMID = SAMID, HTO_uri = HTO_info)


# merge arcseq and HTO paths
merge_file_info <- merge(HTO_uri, arcseq_info, by.x = "SAMID" , by.y = "sampleName")
merge_file_info <- merge_file_info[order(merge_file_info$SAMID), ]



######## load matrices
# load RNA matrix

hashing_qc <- list()
umapplot <- list()
seRNA_final <- list()
for (i in seq_len(nrow(merge_file_info))) {
    message(merge_file_info$SAMID[i])
    # import gex
    seRNA <- import10xFeatureMatrix(
        input = file.path( merge_file_info$uri[i], "raw_feature_bc_matrix.h5"),
        names = merge_file_info$SAMID[i])
    names(assays(seRNA)) <- "counts"

    # import HTO and convert to sce
    HTO <- data.table::fread(merge_file_info$HTO_uri[i])
    rownames_HTO <- HTO$Antibody
    HTO <- HTO[,-1]
    colnames(HTO) <- paste0(merge_file_info$SAMID[i], "#", colnames(HTO),"-1")
    HTO <- as(as.matrix(HTO), "dgCMatrix")
    rownames(HTO) <- rownames_HTO
    HTO <- SingleCellExperiment(assays = list(counts=HTO))

    # merge GEX and HTO into a SCE
    common_cells <- intersect(colnames(HTO), colnames(seRNA))
    HTO <- HTO[, common_cells]
    seRNA <- seRNA[, common_cells]
    seRNA <- as(seRNA, "SingleCellExperiment")
    altExp(seRNA, "HTO") <- HTO

    # call empty droplets to define ambient droplets
    set.seed(10010)
    e.out.gene <- emptyDrops(counts(seRNA), by.rank = 30000 )
    is.cell <- e.out.gene$FDR <= 0.001
    summary(is.cell)

    # plot empty droplet assignments
    par(mfrow=c(1,2))
    r <- rank(-e.out.gene$Total)
    plot(r, e.out.gene$Total, log="xy", xlab="Rank", ylab="Total gene count", main="")
    abline(h=metadata(e.out.gene)$retain, col="darkgrey", lty=2, lwd=2)
    hist(log10(e.out.gene$Total[is.cell]), xlab="Log[10] gene count", main="")

    # Estimate HTO ambient proportions using empty droplets
    hto.mat <- assay(altExp(seRNA),"counts")[,which(is.cell)]
    ambient <- proportions(rowSums(assay(altExp(seRNA), "counts")[,is.na(e.out.gene$FDR)]))

    # plot ambient proportions
    barplot(ambient,las=2, main="ambient proportion")
    hash.stats <- hashedDrops(hto.mat,
                              ambient=ambient,
                              doublet.nmads = 1.5,
                              doublet.min = 2,
                              confident.nmads = 3,
                              confident.min = 1)
    #doublet.min = 3, confident.min=1, confident.nmads = 2
    table(hash.stats$Best[hash.stats$Confident])

    # examine hashing
    hash.stats$colors <- NA
    hash.stats$colors[hash.stats$Doublet] <- "doublet"
    hash.stats$colors[hash.stats$Confident] <- "confident"

    hashing_qc[[merge_file_info$SAMID[i]]] <- ggplot(data.frame(hash.stats),
                                                     aes(LogFC, LogFC2, color=colors)) +
        geom_point(size=0.5)+ scale_color_manual(values=c("black","red")) +
        xlab("Log fold-change from best to second HTO")+
        ylab("Log fold-change of second HTO over ambient") +
        ggtitle(merge_file_info$SAMID[i])


    # keep only non-empty cells
    seRNA <- seRNA[, which(is.cell)]
    colData(seRNA) <- cbind(colData(seRNA), hash.stats)
    colData(seRNA)$library <- sapply(strsplit(colnames(seRNA), split = "#"), "[",1)

    assay(altExp(seRNA), "logcounts") <- log10(assay(altExp(seRNA), "counts")+1)
    assay(altExp(seRNA), "clr") <- sweep(assay(altExp(seRNA), "logcounts"), 2,
                                         colMeans(assay(altExp(seRNA), "logcounts")), "-")
    seRNA <- runUMAP(seRNA, altexp = "HTO", name="UMAP_HTO", assay.type = "clr", exprs_values = "clr")
    seRNA$hash_assignment <- rownames_HTO[seRNA$Best]
    umapplot[[merge_file_info$SAMID[i]]] <- plotReducedDim(seRNA[, which(seRNA$Doublet == FALSE & seRNA$Confident == TRUE)],
                                                           dimred = "UMAP_HTO",
                                                           color_by = "hash_assignment",
                                                           point_size=0.3 ) + ggtitle(merge_file_info$SAMID[i])

    # save seRNA
    seRNA_final[[merge_file_info$SAMID[i]]] <- seRNA
}


seRNA_final <- do.call(cbind, seRNA_final)

saveRDS(seRNA_final, "OUTPUT/seRNA_final.rds")
saveRDS(umapplot, "OUTPUT/umapplot.rds")
saveRDS(hashing_qc, "OUTPUT/hashing_qc.rds")

pdf("OUTPUT/hashing.pdf", width = 12, height = 12)
grid.arrange(grobs=hashing_qc)
grid.arrange(grobs=umapplot)
dev.off()
