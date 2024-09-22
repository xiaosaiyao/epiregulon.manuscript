rm(list = ls())
###### 09 construct GRN using correlation between expressions of TFs and GRN targets for benchmarking
library(tidyverse)
library(yardstick)
library(gridExtra)
library(ggplotify)
library(epiregulon)

# ######## source functions ############
# source("/gstore/project/lineage/GRN/analysis/calculateActivity.R")
# source("/gstore/project/lineage/analysis/function/functions.R")
# source("/gstore/project/lineage/analysis/function/function.3.R")

######## load regulon ###########
# regulon <- readRDS("/gstore/project/lineage/sam/heme_GRN/OUTPUT/Peak2GeneLinks/p2g_long.rds")
# regulon <- regulon %>% dplyr::select(TF, Gene, Correlation)
# colnames(regulon) <- c("tf","target","corr")

regulon.full <- readRDS("OUTPUT/heme.regulon.full.rds")

####### load filtered CRISPR data ####
sce <- readRDS("/gstore/project/lineage/sam/heme_GRN/OUTPUT/small.sce.rds")

######## trim regulon
tf_unique <- unique(c(sce$first_gRNA, sce$second_gRNA))
regulon <- regulon.full %>% dplyr::filter(tf %in% tf_unique)

####### calculate activity score using correlation between expressions ########
# ### log normalize data
# library(scater)
# seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
# table(colData(seRNA)$BioClassification)
# seRNA <- logNormCounts(seRNA)
# assayNames(sce)

#### Run addWeights function to compute correlation #######
#source("/gstore/project/lineage/GRN/analysis/addWeights.R")
#add weights to regulon based on expression matrix
regulon.w <- addWeights(regulon = regulon,
                        sce = sce,
                        cluster_factor = "guide_identity",
                        block_factor = NULL,
                        corr = TRUE,
                        MI = FALSE,
                        BPPARAM = BiocParallel::MulticoreParam()
)

####### calculate activity score ########

score.combine <- calculateActivity(sce = sce,
                                   regulon = regulon.w)

assay(sce, "logcounts") <- as(assay(sce, "logcounts"), Class = "dgCMatrix")

score.combine.aucell <- calculateActivity(sce = sce,
                                         regulon = regulon.w,
                                         method ="aucell",
                                         ncore = 4)

#################### plot heatmap
pheatmap::pheatmap(score.combine[,1:1000],show_colnames = F, scale = "row", breaks = seq(-3, 3, length.out=101))
pheatmap::pheatmap(score.combine.aucell[,1:1000],show_colnames = F, scale = "row", breaks = seq(-3, 3, length.out=101))
saveRDS(score.combine, "./OUTPUT/CRIPSR_activity_expression_CRISPR.rds")
saveRDS(score.combine.aucell, "./OUTPUT/CRIPSR_activity_expression_CRISPR.aucell.rds")

############## plot and generate ROC ##############
#select TF KO with Neg_control
NT_index <- which(sce$guide_identity == "NegCtrl0_NegCtrl0__NegCtrl0_NegCtrl0")

score.combine.names <- c("score.combine","score.combine.aucell")
pdf.names <- c("OUTPUT/roc.plot_expression_CRISPR.pdf","OUTPUT/roc.plot.aucell_expression_CRISPR.pdf")
matrix.names <- c("ROC.avg.matrix","ROC.AUCell.matrix")
matrix.filenames <- c("OUTPUT/roc.matrix.avg_expression_CRISPR.rds","OUTPUT/roc.matrix.aucell_expression_CRISPR.rds")

#?
#score.combine <- score.combine[-c(13),]; score.combine.aucell <- score.combine.aucell[-c(32),]

for (j in 1:2){

    pdf(pdf.names[j], width=18, height=9)
    score.combine.now <- get(score.combine.names[j])
    genes_to_plot <- rownames(score.combine.now)
    ROC.matrix <- data.frame(matrix(NA, nrow = length(rownames(score.combine.now)), ncol=6))
    rownames(ROC.matrix) <- rownames(score.combine.now)
    colnames(ROC.matrix) <- c("expression_Active","expression_NTC","activity_Active","activity_NTC","AUC_expression","AUC_activity")

    for (i in 1:length(genes_to_plot)){
        ##
        tf_name <- unname(genes_to_plot[i])

        #select TF being activated
        gRNA_name <- paste0(tf_name,"_NegCtrl0","__",tf_name,"_NegCtrl0")
        gRNA_name2 <- paste0("NegCtrl0_",tf_name,"__","NegCtrl0_",tf_name)

        #identify TF gene expression and activity respectively
        gene_expr <- logcounts(sce)[match(tf_name,rowData(sce)$symbol),]
        activity <- score.combine.now[tf_name,]
        df <- data.frame(gene_expr,activity)
        expr_vs_activity <- ggplot(df, aes(gene_expr, activity)) + geom_point() + theme_classic()
        gRNA_col <- rep(NA, ncol(sce))
        gRNA_col[NT_index] <- "NTC"
        gRNA_col[grep(paste0(gRNA_name, "|", gRNA_name2), sce$guide_identity)] <- "Active"

        ROC.matrix$expression_Active[i] <- mean(gene_expr[which(gRNA_col== "Active")],na.rm=T)
        ROC.matrix$expression_NTC[i] <- mean(gene_expr[which(gRNA_col== "NTC")],na.rm=T)

        ROC.matrix$activity_Active[i] <- mean(activity[which(gRNA_col== "Active")],na.rm=T)
        ROC.matrix$activity_NTC[i] <- mean(activity[which(gRNA_col== "NTC")],na.rm=T)


        #plot ROC of activity
        rocdata <- data.frame(truth = factor(as.character(gRNA_col),levels=c("Active","NTC")), estimate = activity)
        colnames(rocdata) <- c("truth","estimate")

        auc <- roc_auc(rocdata, truth, estimate, na_rm = T)
        print(paste(gRNA_name,"activity"))
        print(auc)

        roc <- roc_curve(rocdata,truth, estimate, na_rm = T)
        activity.roc.plot <- autoplot(roc) + ggtitle(paste("activity AUC=",round(auc$.estimate,4)))
        #print(activity.roc.plot)

        ROC.matrix$AUC_activity[i] <- auc$.estimate

        #plot ROC of gene expression
        rocdata <- data.frame(truth=factor(as.character(gRNA_col),levels=c("Active","NTC")), estimate=gene_expr)
        colnames(rocdata) <- c("truth","estimate")
        print(paste(gRNA_name,"expression"))
        auc <- roc_auc(rocdata,truth, estimate, na_rm = T)
        print(auc)

        roc <- roc_curve(rocdata,truth, estimate, na_rm = T)
        expr.roc.plot <- autoplot(roc) + ggtitle(paste("expression AUC=", round(auc$.estimate,4)))
        #print(expression.roc.plot)

        ROC.matrix$AUC_expression[i] <- auc$.estimate

        #plot umap plot of activity
        score.update <- data.frame(activity=activity, gene_expr, gRNA_col,UMAP1=reducedDim(sce,"UMAP")[,1], UMAP2=reducedDim(sce,"UMAP")[,2])
        plotscore <- ggplot(score.update, aes_string("UMAP1","UMAP2",color=activity)) + ggtitle(paste(tf_name,"activity")) +
            geom_point(size=0.3) + scale_color_gradient(low = "lightgrey", high = "red") +
            theme_classic() + theme(plot.title = element_text(hjust = 0.5)) #+
        #geom_text(data=centers.cluster, mapping = aes(x = mean_UMAP1, y = mean_UMAP2,label = 1:length(unique(perturb.gp$merged$cluster))),color = "black", size = 4)
        #limit=c(0,0.05),
        #print(plotscore)

        #plot umap plot of expression
        plotexpr <- ggplot(score.update, aes_string("UMAP1","UMAP2",color=gene_expr)) + ggtitle(paste(tf_name,"expression")) +
            geom_point(size=0.3) + scale_color_gradient(low = "lightgrey", high = "blue") +
            theme_classic() + theme(plot.title = element_text(hjust = 0.5)) #+
        #geom_text(data=centers.cluster, mapping = aes(x = mean_UMAP1, y = mean_UMAP2,label = 1:length(unique(perturb.gp$merged$cluster))),color = "black", size = 4)
        #limit=c(0,0.05),
        #print(plotexpr)

        #plot umap plot of KO
        plotKO <- ggplot(score.update, aes_string("UMAP1","UMAP2", color="gRNA_col")) + ggtitle(paste(tf_name,"guides")) +
            geom_point(size=0.3) + theme_classic() + theme(plot.title = element_text(hjust = 0.5))

        #print(plotKO)

        #plot violin plot of activity
        activity.violin.plot <- ggplot(score.update, aes_string(x = "gRNA_col", y = "activity", fill = "gRNA_col")) +
            geom_violin(scale="width") + theme_classic() + ggtitle(gRNA_name) +
            theme(plot.title = element_text(hjust = 0.5), legend.position="none",
                  axis.text.x = element_text(angle = 90, hjust=1)) + #+ #ylim(-0.2,0.2) +
            ylab("activity") + xlab ("cluster")
        #print(activity.violin.plot)

        #plot violin plot of expression
        expr.violin.plot <- ggplot(score.update, aes_string(x="gRNA_col", y="gene_expr", fill = "gRNA_col")) +
            geom_violin(scale="width") + theme_classic() + ggtitle(gRNA_name) +
            theme(plot.title = element_text(hjust = 0.5), legend.position="none",
                  axis.text.x = element_text(angle = 90, hjust=1)) + #+ #ylim(-0.2,0.2) +
            ylab("gene expression") + xlab ("cluster")
        #print(expr.violin.plot)

        grid.arrange(expr_vs_activity,  plotscore, plotexpr, plotKO, activity.roc.plot, expr.roc.plot, activity.violin.plot, expr.violin.plot, ncol=4)
    }
    saveRDS(ROC.matrix, matrix.filenames[j])
    assign(matrix.names[j],ROC.matrix)
    dev.off()
}

############## summary plot ########
print(mean(ROC.AUCell.matrix$AUC_expression))
print(mean(ROC.AUCell.matrix$AUC_activity))
print(mean(ROC.avg.matrix$AUC_expression))
print(mean(ROC.avg.matrix$AUC_activity))

ROC.avg.matrix$TF <- rownames(ROC.avg.matrix)
ROC.avg.matrix$logFC.activity <- log2(ROC.avg.matrix$activity_Active/ROC.avg.matrix$activity_NTC)
ROC.avg.matrix$logFC.expression <- log2(ROC.avg.matrix$expression_Active/ROC.avg.matrix$expression_NTC)

ROC.AUCell.matrix$TF <- rownames(ROC.AUCell.matrix)
ROC.AUCell.matrix$logFC.activity <- log2(ROC.AUCell.matrix$activity_Active/ROC.AUCell.matrix$activity_NTC)
ROC.AUCell.matrix$logFC.expression <- log2(ROC.AUCell.matrix$expression_Active/ROC.AUCell.matrix$expression_NTC)


pdf("OUTPUT/summary.plot_expression_CRISPR.pdf")
ggplot(ROC.avg.matrix, aes(AUC_expression, AUC_activity)) + geom_point() +geom_text(aes(label=TF,hjust=0, vjust=0)) +
    geom_abline(slope = 1, intercept = 0) + theme_classic() + ylim (0.4, 1) + xlim (0.4,1)
ggplot(ROC.AUCell.matrix, aes(AUC_expression, AUC_activity)) + geom_point() +geom_text(aes(label=TF,hjust=0, vjust=0)) +
    geom_abline(slope = 1, intercept = 0) + theme_classic() + ylim (0.4, 1) + xlim (0.4,1)

dev.off()

