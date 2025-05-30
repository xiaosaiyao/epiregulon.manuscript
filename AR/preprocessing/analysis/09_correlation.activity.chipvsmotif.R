library(ggplot2)
library(gridExtra)

# Compare total activity

p <- list()

celllines <- c("LNCaP", "VCaP", "MDA")

for (cellline in celllines){
    activity_matched <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.", cellline,".chip.rds"))
    activity_merged <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.",cellline,".rds"))
    activity_motif <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/score.combine.",cellline,".motif.rds"))

    # matched chipseq vs merged chipseq
    common <- intersect(rownames(activity_matched), rownames(activity_merged))
    common <- intersect(common, rownames(activity_motif))
    common <- setdiff(common, "EP300")

    activity_merged <- activity_merged[common,]
    activity_matched <- activity_matched[common,]

    correlation <- cor(t(as.matrix(activity_matched)), t(as.matrix(activity_merged )))

    correlation.vector <- colSums(correlation*diag(1,length(common)))
    correlation.df <- data.frame(regulator=common,
                                 correlation = correlation.vector,
                                 source = "merged_chipseq")

    # matched chipseq vs merged chipseq
    activity_motif <- activity_motif[common,]

    correlation <- cor(t(as.matrix(activity_matched)), t(as.matrix(activity_motif )))

    correlation.vector <- colSums(correlation*diag(1,length(common)))
    correlation.df2 <- data.frame(regulator=common,
                                 correlation = correlation.vector,
                                 source = "motif")
    correlation.df <- rbind(correlation.df, correlation.df2)
    write.table(correlation.df, paste0("OUTPUT/ArchRProject/Epiregulon/correlation.merged.vs.motif", cellline,".txt"), sep="\t", quote=FALSE, row.names = FALSE)

    # plot results
    correlation.df$regulator <- factor(correlation.df$regulator,
                                       levels = unique(correlation.df$regulator[order(correlation.df$source,correlation.df$correlation)]))

    ggtitle <- paste(cellline,"activity correlation")
    p[[ggtitle ]] <- ggplot(data=correlation.df, aes(x=regulator, y=correlation, fill=source)) +
        geom_bar(stat="identity", position=position_dodge(preserve = "single")) + theme_classic()  + ggtitle(ggtitle ) +
        scale_fill_manual(values=c('darkred','darkblue')) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}



# Compare overlap in target genes
jaccard <- function(x, a, b) {
    a <- unique(a[a$tf==x, "target"])
    b <- unique(b[b$tf==x,"target"])
    intersection <- length(intersect(a, b))
    union <- length(a) + length(b) - intersection
    return (intersection/union)
}


for (cellline in celllines){
    regulon_matched <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.", cellline,".chip.rds"))
    regulon_merged <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.",cellline,".rds"))
    regulon_motif <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.",cellline,".motif.rds"))

    # matched chipseq vs merged chipseq
    common <- intersect(unique(regulon_matched$tf), unique(regulon_merged$tf))
    common <- intersect(common, unique(regulon_motif$tf))
    common <- setdiff(common, "EP300")

    jaccard.similarity <- sapply(common, jaccard, regulon_matched, regulon_merged)

    similarity.df <- data.frame(regulator=common,
                                similarity=jaccard.similarity,
                                source = "merged_chipseq")

    similarity.df$regulator <- factor(similarity.df$regulator,
                                      levels = similarity.df$regulator[order(similarity.df$similarity)])



    # matched chipseq vs merged chipseq
    jaccard.similarity <- sapply(common, jaccard, regulon_matched, regulon_motif)

    similarity.df2 <- data.frame(regulator=common,
                                similarity=jaccard.similarity,
                                source = "motif")
    similarity.df <- rbind(similarity.df,similarity.df2)
    write.table(similarity.df, paste0("OUTPUT/ArchRProject/Epiregulon/jaccard.similarity.merged.vs.motif", cellline,".txt"), sep="\t", quote=FALSE, row.names = FALSE)

    # plot results
    similarity.df$regulator <- factor(similarity.df$regulator,
                                      levels = unique(similarity.df$regulator[order(similarity.df$source, similarity.df$similarity)]))

    ggtitle <- paste(cellline," target similarity")
    p[[ggtitle]] <- ggplot(data=similarity.df, aes(x=regulator, y=similarity, fill=source)) +
        geom_bar(stat="identity", position=position_dodge(preserve = "single")) + theme_classic()  + ggtitle(ggtitle) +
        scale_fill_manual(values=c('darkred','darkblue')) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

}

pdf("OUTPUT/ArchRProject/Epiregulon/activity.regulon.groundtruth.LNCaP.pdf", width=6, height=6)
gridExtra::marrangeGrob(grobs=list(p$`LNCaP activity correlation`, p$`LNCaP  target similarity`), nrow=3, ncol=1)
dev.off()

pdf("OUTPUT/ArchRProject/Epiregulon/activity.regulon.groundtruth.VCaP.pdf", width=5, height=6)
gridExtra::marrangeGrob(grobs=list(p$`VCaP activity correlation`, p$`VCaP  target similarity`), nrow=3, ncol=1)
dev.off()

pdf("OUTPUT/ArchRProject/Epiregulon/activity.regulon.groundtruth.MDA.pdf", width=3, height=6)
gridExtra::marrangeGrob(grobs=list(p$`MDA activity correlation`, p$`MDA  target similarity`), nrow=3, ncol=1)
dev.off()

