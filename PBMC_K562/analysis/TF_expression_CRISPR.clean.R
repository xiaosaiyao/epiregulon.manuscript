#library(tidyverse)
library(yardstick)
library(gridExtra)
library(ggplotify)
library(org.Hs.eg.db)
library(annotate)
library(ggrepel)
library(epiregulon)
library(decoupleR)
library(OmnipathR)
library(Seurat)


######## load sce GSE133344 #########
sce <- readRDS("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/data/CRISPRa.perturb.gp.rds")
sce <- sce$sce
scater::plotReducedDim(sce, dimred="UMAP", colour_by = "cluster")

rownames(sce) <- rowData(sce)$symbol

unique_TFs <- unique(c(sce$first_gRNA, sce$second_gRNA))

#trim cells
cell_index <- c(which(sce$first_gRNA %in% unique_TFs& sce$second_gRNA =="NegCtrl0"),
             which(sce$second_gRNA %in% unique_TFs & sce$first_gRNA =="NegCtrl0"),
             which(sce$second_gRNA =="NegCtrl0" & sce$first_gRNA =="NegCtrl0"))

CRISPR.sce <- sce[,cell_index]

# add tf identity
CRISPR.sce$tf <- "NTC"
CRISPR.sce$tf[CRISPR.sce$first_gRNA != "NegCtrl0"] <- CRISPR.sce$first_gRNA[CRISPR.sce$first_gRNA != "NegCtrl0"]
CRISPR.sce$tf[CRISPR.sce$second_gRNA != "NegCtrl0"] <- CRISPR.sce$second_gRNA[CRISPR.sce$second_gRNA != "NegCtrl0"]

# add status
CRISPR.sce$status <- "Active"
CRISPR.sce$status[CRISPR.sce$tf == "NTC"] <- "NTC"
CRISPR.sce$status <- factor(CRISPR.sce$status, levels = c("Active", "NTC"))


saveRDS(CRISPR.sce, "/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/data/small.rds")

###################### load small sce and regulon #####################
CRISPR.sce <- readRDS("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/data/small.rds")

unique_TFs <- unique(c(CRISPR.sce$first_gRNA, CRISPR.sce$second_gRNA))
unique_TFs <- unique_TFs[unique_TFs %in% rownames(CRISPR.sce)]

# load different regulons
regulon <- list()

# epiregulon
#regulon[["epiregulon"]] <- readRDS("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/regulon.rds")
#regulon[["epiregulon.pruned"]] <- readRDS("/gstore/project/lineage/manuscript/epiregulon/OUTPUT/regulon.w.rds")

regulon[["epiregulon"]] <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/regulon.w.rds")

# dorothea
library(dorothea)
data(dorothea_hs, package = "dorothea")
regulon[["dorothea"]] <- dorothea_hs

# humanbase
regulon[["humanbase"]] <- read.table("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/data/blood_top.gz")
conversion_table <- getSYMBOL(as.character(unique(regulon[["humanbase"]]$V1)), data='org.Hs.eg')
regulon[["humanbase"]]$tf <- conversion_table[as.character(regulon[["humanbase"]]$V1)]
regulon[["humanbase"]]$target <- conversion_table[as.character(regulon[["humanbase"]]$V2)]
regulon[["humanbase"]] <- na.omit(regulon[["humanbase"]])
regulon[["humanbase"]] <- regulon[["humanbase"]][, c("tf","target")]

# collecTRI
regulon[["get_collectri"]] <- decoupleR::get_collectri(organism='human', split_complexes=FALSE)
colnames(regulon[["get_collectri"]]) <- c("tf", "target","mor")

# Pando
Seurat_obj <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/Seurat_obj.rds")
regulon[["pando"]] <- data.frame(Seurat_obj@grn@networks[[1]]@coefs[,c("tf","target")])# target genes in 'target' column, tf in 'tf'

# FigR
FigR_GRN <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/FigR_GRN.rds") # target genes in 'DORC' column, tf in 'Motif' column
regulon[["FigR"]] <- FigR_GRN[,c("Motif","DORC")]
colnames(regulon[["FigR"]]) <- c("tf","target")

# Scenic plus
scenic_GRN <- read.csv('/gstore/project/epigen/PBMC/analysis/scenicplus/OUTPUT/scenic_GRN.csv') # target genes in 'Gene' column, tf in 'TF' column
regulon[["scenicplus"]] <- scenic_GRN[,c("TF","Gene")]
colnames(regulon[["scenicplus"]]) <- c("tf","target")

# trim regulon
regulon.short <- lapply(names(regulon), function(x){ regulon[[x]][which(regulon[[x]]$tf %in% unique_TFs),]})
names(regulon.short) <- names(regulon)

################################
# compute weights
score.combine <- list()

for (regulon_name in names(regulon.short)){
    print(regulon_name)
    # addWeights
    regulon.short[[regulon_name ]] <- addWeights(regulon = regulon.short[[regulon_name]],
                                                 expMatrix  = CRISPR.sce,
                                                 clusters = CRISPR.sce$guide_identity,
                                                 method = "corr",
                                                 BPPARAM = BiocParallel::MulticoreParam())

    # calculate activity
    score.combine[[regulon_name ]] <- calculateActivity(expMatrix = CRISPR.sce,
                                                        regulon =  regulon.short[[regulon_name]])
}




results <- data.frame(matrix(NA, nrow=length(unique_TFs), ncol = length(names(regulon))+1))
colnames(results) <- c("gex", names(regulon))
rownames(results) <- unique_TFs

roc.plot <- list()
violin.plot <- list()

for (tf in unique_TFs){
    message(tf)
    CRISPR.sce.tf <- CRISPR.sce[, CRISPR.sce$tf %in% c(tf,"NTC")]

    # gene expression auc
    rocdata <- data.frame(truth = CRISPR.sce.tf$status,
                          estimate = counts(CRISPR.sce.tf)[tf, colnames(CRISPR.sce.tf)])
    colnames(rocdata) <- c("truth","estimate")
    auc <- roc_auc(rocdata, truth, estimate, na_rm = T)
    results[tf, "gex"] <- auc$.estimate

    # plot ROC activity
    roc <- roc_curve(rocdata,truth, estimate, na_rm = T)
    roc.plot[[tf]][["gex"]] <- autoplot(roc) +
        ggtitle(paste(tf, "expression AUC=", round(auc$.estimate,4)))

    # violin plot activity
    violin.plot[[tf]][["gex"]] <- ggplot(rocdata, aes_string(x="truth", y="estimate", fill = "truth")) +
        geom_violin(scale="width") + theme_classic() + ggtitle(paste(tf, "expression")) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none",
              axis.text.x = element_text(angle = 90, hjust=1)) +
        ylab("gene expression") + xlab ("cluster")

    for (regulon_name in  names(regulon.short)){

        if (tf  %in% rownames(score.combine[[regulon_name]])){
            rocdata <- data.frame(truth = CRISPR.sce.tf$status,
                                  estimate = score.combine[[regulon_name]][tf, colnames(CRISPR.sce.tf)])
            colnames(rocdata) <- c("truth","estimate")
            auc <- roc_auc(rocdata, truth, estimate, na_rm = T)
            results[tf, regulon_name] <- auc$.estimate

            # plot ROC activity
            roc <- roc_curve(rocdata,truth, estimate, na_rm = T)
            roc.plot[[tf]][[regulon_name]] <- autoplot(roc) +
                ggtitle(paste(tf, regulon_name, "AUC=", round(auc$.estimate,4)))

            # violin plot activity
            violin.plot[[tf]][[regulon_name]] <- ggplot(rocdata,
                                                        aes_string(x="truth", y="estimate", fill = "truth")) +
                geom_violin(scale="width") + theme_classic() + ggtitle(paste(tf, regulon_name)) +
                theme(plot.title = element_text(hjust = 0.5), legend.position="none",
                      axis.text.x = element_text(angle = 90, hjust=1)) +
                ylab("activity") + xlab ("cluster")
        }

    }

}

results$TF <- rownames(results)
saveRDS(violin.plot,"/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/violin.plot.rds")
saveRDS(roc.plot,"/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/roc.plot.rds")
saveRDS(results,"/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/results.rds")


# plot summary plot
summary_plots <- list()
for (regulon_name in  names(regulon.short)){
    summary_plots[[regulon_name]] <- ggplot(results, aes_string("gex", regulon_name)) +
        geom_point() +
        geom_text_repel(aes(label=TF), max.overlaps = Inf) +
        geom_abline(slope = 1, intercept = 0) +
        theme_classic() + ylim (0.4, 1.2) + xlim (0.4,1.2) +
        theme(plot.title = element_text(hjust = 0.5)) +
        ylab("Activity AUC") + xlab("Gene expression AUC") +
        ggtitle(regulon_name)

}


results_df <- reshape2::melt(results, measure.vars= c("gex", "epiregulon",  "dorothea", "humanbase", "get_collectri", "FigR", "scenicplus", "pando" ))
colnames(results_df) <-  c("TF", "methods", "auc")

results_df$source <- "database"
results_df$source[results_df$methods %in% c("epiregulon.pruned","FigR","scenicplus","pando")] <- "inferred"

library(dplyr)
results.summary <- results_df %>% group_by(methods) %>% summarize(median(auc, na.rm=TRUE))

# plot summary bar plots
library(ggbeeswarm)
results_df$methods <- factor(results_df$methods, levels= results.summary$methods[order(results.summary$`median(auc, na.rm = TRUE)`, decreasing = TRUE)])
barplots <- ggplot(results_df, aes(x=methods, y=auc, fill=source)) +
    geom_boxplot() +
    geom_quasirandom(width=0.3, size=0.5) +    #xlim(0,1)+
    #geom_text(aes(label=auc), vjust=-0.3, size=3.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = "none")


pdf("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/summary.plot.pdf", height=9)
marrangeGrob(grobs=summary_plots, ncol=2, nrow=3)
dev.off()

pdf("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/summary.plot2.pdf", height=3, width=3)
barplots
dev.off()



# plot an example gene
example_plot <- list()
example_plot[[1]] <- roc.plot$PRDM1$gex
example_plot[[2]] <- roc.plot$PRDM1$epiregulon
example_plot[[3]] <- roc.plot$PRDM1$dorothea
example_plot[[4]] <- violin.plot$PRDM1$gex
example_plot[[5]] <- violin.plot$PRDM1$epiregulon
example_plot[[6]] <- violin.plot$PRDM1$dorothea

pdf("/gstore/project/lineage/manuscript/epiregulon/PBMC_K562/OUTPUT/example.plot.pdf", height = 5, width=8)
marrangeGrob(grobs=example_plot, nrow=2, ncol=3)
dev.off()


#ggcorrplot::ggcorrplot(cor(as.matrix(t(score.combine$epiregulon)[1:100,])))
