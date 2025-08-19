source("/gstore/project/epigen/benchmark/src/utils.R")
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

#################
### Figure S2 ###
#################

source("/gstore/project/epigen/benchmark/src/utils.R")
# use a fixed logFC and p-value threshold
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/GeneExpressionMatrix.rds")
perturbed_tfs <- c("STAT6","ELK1", "GATA3", "JUN", "NFATC3", "NFKB1", "STAT3", "MAF")
all_genes <- unique(rownames(GeneExpressionMatrix))
TF_map <- read.csv("/gstore/project/epigen/PBMC/TF_to_dataset")
# select TFs
selected_tfs <- c()
for (tf in perturbed_tfs){
  print(tf)
  file_name <- TF_map[grep(tf, TF_map$TF),"Dataset_file"]
  file_path <- file.path("/gstore/project/epigen/PBMC/raw_data/KnockTF" ,file_name)
  knockout_data <- read.csv(file_path)
  KD_data_targets <- unique(knockout_data[(knockout_data$Corrected_P>0.05) & (abs(knockout_data$Log2FC)>0.5), "Target.Gene"])
  print(length(KD_data_targets))
  if(length(KD_data_targets)>1000) selected_tfs <- c(selected_tfs, tf)
}
perturbed_tfs <- selected_tfs

library(dplyr)

# check how the precision and recall depend on the epiregulon weights
df <- data.frame()
TF_map <- read.csv("/gstore/project/epigen/PBMC/TF_to_dataset")
package_output <- c(Epiregulon_chip_seq_merged = "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_merged_trimmed_no_clusters.rds",
                    Epiregulon_motifs = "/gstore/project/epigen/PBMC/OUTPUT/regulon.w_motifs.rds",
                    Epiregulon_chip_seq_tissue = "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_blood_trimmed_no_clusters.rds",
                    Epiregulon_motif_score_blood="/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_blood_trimmed_no_clusters_motifs.rds",
                    Epiregulon_motif_score_merged="/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_merged_trimmed_no_clusters_motifs.rds",
                    FigR = "/gstore/project/epigen/benchmark/PBMC/OUTPUT/FigR_GRN.rds",
                    Pando = "/gstore/project/epigen/PBMC/OUTPUT/Seurat_obj.rds",
                    scplus = '/gstore/project/epigen/PBMC/analysis/scenicplus/OUTPUT/scenic_GRN.csv',
                    GRaNIE = '/gstore/project/epigen/PBMC/OUTPUT/GRN_connections.all.rds'

)
for(package in names(package_output)){
  print(package)
  if(package == "scplus") grn <- read.csv('/gstore/project/epigen/PBMC/analysis/scenicplus/OUTPUT/scenic_GRN.csv')
  else grn <- readRDS(package_output[package])
  if(package=="Pando") grn <- coef(grn)
  for(tf in perturbed_tfs){
    new_row <- c()
    print(tf)
    if(!tf %in% switch(package, "Epiregulon_chip_seq_merged" = grn$tf,
                       "Epiregulon_motifs" = grn$tf,
                       "Pando" = grn$tf,
                       "FigR" = grn$Motif,
                       "scplus" = grn$TF,
                       "GRaNIE" = grn$TF.name)) next
    print("tf found")
    file_name <- TF_map[grep(tf, TF_map$TF),"Dataset_file"]
    file_path <- file.path("/gstore/project/epigen/PBMC/raw_data/KnockTF" ,file_name)
    knockout_data <- read.csv(file_path)
    KD_data_targets <- unique(knockout_data[(knockout_data$Corrected_P>0.05) & (abs(knockout_data$Log2FC)>0.5), "Target.Gene"])
    regulon_targets <- switch(package,
                              "Pando" = unique(unlist(grn[grn$tf==tf,"target"])),
                              "FigR" = unique(grn[grn$Motif==tf&grn$Score!=0,"DORC"]),
                              "scplus" = unique(unlist(grn[grn$TF==tf,"Gene"])),
                              "GRaNIE" = unique(unlist(grn[grn$TF.name==tf,"gene.name"])),
                              unique(grn[grn$tf==tf,"target"]))
    new_row = c(new_row, KD_data_targets = length(KD_data_targets))
    new_row = c(new_row, regulon_size = length(regulon_targets))
    new_row = c(new_row, precision = sum(regulon_targets %in% KD_data_targets)/length(regulon_targets))
    new_row = c(new_row, recall = sum(regulon_targets %in% KD_data_targets)/length(intersect(KD_data_targets, all_genes)))
    new_row = c(new_row, precision_random = length(intersect(all_genes, KD_data_targets))/length(all_genes))
    new_row = c(new_row, recall_random = length(regulon_targets)/length(all_genes))
    new_row = c(new_row, tf=tf, package = package)
    df <- rbind(df, as.list(new_row))
  }
}

df_plot_1 <- gather(df, "metric","value", "recall", "precision")
df_plot_1$value <- as.numeric(df_plot_1$value)
df_plot_2 <- gather(df, "metric","value", "recall_random", "precision_random")
df_plot_2 <- df_plot_2[,c("tf", "metric", "value", "package")]
df_plot_2$value <- as.numeric(df_plot_2$value)
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/PR_barplot_supplementary.pdf", width = 6, height = 14.5)
ggplot(data=df_plot_1, aes(fill = metric, y = value, x = tf))+
  geom_bar(stat = "identity", width = 0.5, position=position_dodge(width = 0.5))+
  geom_point(data=df_plot_2, aes(fill = metric, y = value, x = tf), shape=95, size=6, position=position_dodge(width = 0.5))+
  facet_wrap(~package, ncol=1)+
  theme_classic()
dev.off()


#################
### Figure 1g ###
#################

perturbed_tfs <- c("STAT6","ELK1", "GATA3", "JUN", "NFATC3", "NFKB1", "STAT3", "MAF")
df_copy <- df
df <- df[df$package %in% c("Epiregulon_chip_seq_merged", "Epiregulon_motifs", "FigR", "Pando", "scplus", "GRaNIE"),]
df_plot <- gather(df, "metric","value", "recall", "precision", "recall_random", "precision_random")
df_plot$random <- FALSE
df_plot$random[grep("random",df_plot$metric)] <- TRUE
df_plot$metric <- gsub("_random","",df_plot$metric)
df_plot$value <- as.numeric(df_plot$value)
df_plot <- df_plot[!df_plot$random,]
df_plot_copy <- df_plot
df_plot <- df_plot[df_plot$metric=="recall", ]
colors <- c(Epiregulon_chip_seq_merged = "red", Epiregulon_motifs ="red", FigR = "#888888", Pando = "blue", GRaNIE = "green1",
            scplus = "orange",
            "Gene expression" = "darkslategray2")

colors <- colors[unique(df_plot$package)]
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/PR_barplot_recall.pdf", width = 5, height = 3.5)
ggplot(df_plot, aes(x=package, y=value, fill=package)) +
  #geom_boxplot() +
  geom_bar(position = position_dodge(width = 0.90), fun = "median", stat = "summary")+
  #geom_point(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white") +
  #geom_jitter(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white")+
  theme_classic() +
  ylab("recall")+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggnewscale::new_scale_fill()+
  geom_quasirandom(dodge.width=0.75, size=1.2, color="white", shape = 21, fill="black")
dev.off()

df_plot %>% split(df_plot$package) %>% lapply(function(x) median(x$value))

df_plot <- df_plot_copy
df_plot <- df_plot[df_plot$metric=="precision", ]
library(ggbeeswarm)
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/PR_barplot_precision.pdf", width = 5, height = 3.5)
ggplot(df_plot, aes(x=package, y=value)) +
  #geom_boxplot() +
  geom_bar(aes(fill=package), position = position_dodge(width = 0.90), fun = "median", stat = "summary")+

  #geom_point(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white") +
  #geom_jitter(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white")+
  theme_classic() +
  ylab("precision")+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggnewscale::new_scale_fill()+
  geom_quasirandom(dodge.width=0.75, size=1.2, color="white", shape = 21, fill="black")
dev.off()

df_plot %>% split(df_plot$package) %>% lapply(function(x) median(x$value))


#################
### Figure 1d ###
#################


lineage_tfs <- c("TCF7", "GATA3", "BCL11B", "RUNX1", "FOXP3", "RUNX3", "IKZF1", "SPI1",
                 "CEBPA", "KLF4", "CEBPB", "EBF1", "PAX5", "POU2AF1", "TCF3", "EOMES",
                 "PRDM1", "TBX21", "IRF8", "TCF4")
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/GeneExpressionMatrix.rds")

GeneExpressionMatrix <- GeneExpressionMatrix[,GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "FCGR3A+ Mono", "Memory CD4+ T",
                                                                                    "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "CD14+ Mono")]


GeneExpressionMatrix <- GeneExpressionMatrix[match(lineage_tfs, rownames(GeneExpressionMatrix)),]
cell_type_colors <- c("B" = "#1F78B4", "CD14+ Mono" = "#B2DF8A", "Memory CD4+ T" = "#E31A1C", "DC" = "darkmagenta",
                      "Naive CD8+ T" = "#E31A1C", "FCGR3A+ Mono" = "#B2DF8A", "Monocytes" = "#B2DF8A",
                      "Naive CD4+ T" = "#E31A1C", "Memory CD8+ T" = "#E31A1C", "NK" = "wheat2")

tf_colors <- c("TCF7"= "#E31A1C", "GATA3"= "#E31A1C", "BCL11B"= "#E31A1C","RUNX1" = "#E31A1C",
               "FOXP3"= "#E31A1C", "RUNX3" = "#E31A1C", "IKZF1" ="#E31A1C","SPI1"= "#B2DF8A",
               "CEBPA"= "#B2DF8A", "KLF4"= "#B2DF8A","CEBPB"= "#B2DF8A","EBF1"= "#1F78B4",
               "PAX5"= "#1F78B4","POU2AF1"= "#1F78B4", "TCF3" = "#1F78B4","EOMES" = "wheat2",
               "PRDM1" = "wheat2", "TBX21" ="wheat2","IRF8"= "darkmagenta","TCF4"= "darkmagenta")


library(SummarizedExperiment)
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/master_regulators_expression_heatmap.pdf", width = 12, height = 6)
plotHeatmapActivity(activity=assay(GeneExpressionMatrix)[rownames(assay(GeneExpressionMatrix)) %in% lineage_tfs,],
                    sce=GeneExpressionMatrix,
                    tfs=lineage_tfs,
                    downsample=1000,
                    cluster_rows = FALSE,
                    cell_attributes="cell_type",
                    col_gap="cell_type",
                    name = "transcription factor expression",
                    column_title_rot = 45,
                    columns_col =list(cell_type = cell_type_colors),
                    row_col = list(transcription_factor = tf_colors))
dev.off()

#################
### Figure 1e ###
#################

regulon.w <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_merged_trimmed_no_clusters.rds")
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/GeneExpressionMatrix.rds")

GeneExpressionMatrix <- GeneExpressionMatrix[,GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "FCGR3A+ Mono", "Memory CD4+ T",
                                                                                    "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "CD14+ Mono")]

activity.matrix <- epiregulon::calculateActivity(expMatrix = GeneExpressionMatrix, exp_assay="normalizedCounts",
                                                 regulon = regulon.w)
activity.matrix <- activity.matrix[, GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "FCGR3A+ Mono", "Memory CD4+ T",
                                                                           "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "CD14+ Mono")]

activity.matrix <- activity.matrix[rownames(activity.matrix) %in% lineage_tfs, ]
activity.matrix <- activity.matrix[na.omit(match(lineage_tfs, rownames(activity.matrix))),]

set.seed(2091)
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/master_regulators_activity_heatmap.pdf", width = 12, height = 6)
plotHeatmapActivity(activity=as.matrix(activity.matrix),
                    sce=GeneExpressionMatrix,
                    tfs=lineage_tfs,
                    downsample=1000,
                    cell_attributes="cell_type",
                    col_gap="cell_type",
                    name = "transcription factor activity",
                    cluster_rows = FALSE,
                    column_title_rot = 45,
                    columns_col =list(cell_type = cell_type_colors),
                    row_col = list(transcription_factor = tf_colors)


)
dev.off()


#################
### Figure 1h ###
#################

library(tidyr)
library(ggplot2)
library(ggbeeswarm)
df <- get_runtime_results("/gstore/project/epigen/benchmark/resource_use.txt")
df_plot <- df
colnames(df_plot)[which(colnames(df_plot) == "JobName")] = "package"
colnames(df_plot)[which(colnames(df_plot) == "Elapsed")] = "run_time"
df_plot <- df_plot[grep("Epiregulon_main_PBMC|FigR_main_PBMC_250k|GRaNIE_main_PBMC|Pando_main_PBMC|scenic_plus_main_PBMC",df_plot$package),]
# exclude nh002 which seems to be corrpt
df_plot <- df_plot[(df_plot$node!="nh002"|df_plot$package!="FigR_main_PBMC_250k"),]
df_plot$package <- gsub("_main_PBMC.*","",df_plot$package)
df_plot <- df_plot[!(df_plot$package == "GRaNIE" & df_plot$requested_memory != "128G"),]
df_plot <- df_plot[!(df_plot$package == "Epiregulon" & df_plot$requested_memory != "64G"),]
df_plot <- df_plot[!(df_plot$package == "scenic_plus" & df_plot$requested_memory != "64G"),]
df_plot <- df_plot[!(df_plot$package == "Pando" & df_plot$requested_memory != "64G"),]
runtime <- strsplit(df_plot$run_time, "-")
library(lubridate)
for(i in 1:length(runtime)){
  if(length(runtime[[i]])==2){
    second_period = strsplit(runtime[[i]][2], ":")[[1]]
    second_period <- hours(as.numeric(second_period[1])) + minutes(as.numeric(second_period[2])) + seconds(as.numeric(second_period[3]))
    runtime[[i]] <- days(as.numeric(runtime[[i]][1])) + second_period
    runtime[[i]] <- as.numeric(period_to_seconds(runtime[[i]]))
  }
  else{
    period = strsplit(runtime[[i]][1], ":")[[1]]
    runtime[[i]] <- hours(as.numeric(period[1])) + minutes(as.numeric(period[2])) + seconds(as.numeric(period[3]))
    runtime[[i]] <- period_to_seconds(runtime[[i]])
  }
}

colors <- c(Epiregulon = "red", FigR = "#888888", Pando = "blue", GRaNIE = "green1",
            scenic_plus = "orange",
            "Gene expression" = "darkslategray2")

colors <- colors[unique(df_plot$package)]

df_plot$package <- factor(df_plot$package, levels = names(sort(tapply(df_plot$run_time, df_plot$package, median))))

df_plot$run_time <- unlist(runtime)
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/runtime.pdf", width = 5.5, height = 3.5)
ggplot(df_plot, aes(x=package, y=run_time, fill=package)) +
  #geom_boxplot() +
  geom_bar(position = position_dodge(width = 0.90), fun = "median", stat = "summary")+
  geom_quasirandom(dodge.width=0.75, size=0.8, color="black")+
  scale_fill_manual(values = colors)+
  scale_y_time(labels = duration_to_string, limits = c(0, max(df_plot$run_time)+10000), breaks = seq(0, max(df_plot$run_time+10000), by = 14400))+
  #geom_point(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white") +
  #geom_jitter(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Run time [hours]")
dev.off()

df_plot$MaxRSS <- gsub("G", "", df_plot$MaxRSS)
df_plot$MaxRSS <- as.numeric(df_plot$MaxRSS)
df_plot$package <- factor(df_plot$package, levels = names(sort(tapply(df_plot$MaxRSS, df_plot$package, median))))
dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/memory_use.pdf", width = 5.5, height = 3.5)
ggplot(df_plot, aes(x=package, y=MaxRSS, fill=package)) +
  #geom_boxplot() +
  geom_bar(position = position_dodge(width = 0.90), fun = "median", stat = "summary")+
  geom_quasirandom(dodge.width=0.75, size=0.8, color="black")+
  scale_fill_manual(values = colors)+
  #geom_point(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white") +
  #geom_jitter(position = position_dodge(width = 0.90), size=1.5, shape=21, stroke=1, fill="white")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Maximum memory use [GB]")
dev.off()



#################
### Figure 1f ###
#################

master_regulator_pairs <- list(c("B","EBF1"),
                               c("B","PAX5"),
                               c("B","POU2AF1"),
                               c("B","TCF3"),
                               c("NK","EOMES"),
                               c("NK","PRDM1"),
                               c("NK","TBX21"),
                               c("DC","IRF8"),
                               c("DC","TCF4"),
                               c("Monocytes","SPI1"),
                               c("Monocytes","CEBPA"),
                               c("Monocytes","CEBPB"),
                               c("Monocytes","KLF4"),
                               c("T cells","TCF7"),
                               c("T cells","GATA3"),
                               c("T cells","BCL11B"),
                               c("T cells","RUNX1"),
                               c("T cells","RUNX3"),
                               c("T cells","FOXP3"),
                               c("T cells","IKZF1"))

plot_df <- do.call(rbind, master_regulator_pairs)
colnames(plot_df) <- c("cell type", "transcription factor")
library(epiregulon)
regulon.w <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_merged_trimmed_no_clusters.rds")
master_regulators <- c("EBF1","PAX5","POU2AF1","TCF3","EOMES","PRDM1","TBX21","IRF8",
                       "TCF4","SPI1","CEBPA","CEBPB","KLF4","TCF7","GATA3","BCL11B",
                       "RUNX1","RUNX3","FOXP3","IKZF1")
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/GeneExpressionMatrix.rds")

GeneExpressionMatrix <- GeneExpressionMatrix[,GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "FCGR3A+ Mono", "Memory CD4+ T", "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "CD14+ Mono")]

GeneExpressionMatrix$cell_type <- gsub(".* T$","T cells",GeneExpressionMatrix$cell_type)
GeneExpressionMatrix$cell_type <- gsub(".*Mono.*","Monocytes",GeneExpressionMatrix$cell_type)
GeneExpressionMatrix <- GeneExpressionMatrix[match(master_regulators, rownames(GeneExpressionMatrix)),]
GeneExpressionMatrix <- GeneExpressionMatrix[!rownames(GeneExpressionMatrix) %in% c("POU2F2", "SPIB"),]
library(scran)
assay(GeneExpressionMatrix, "logcounts") <- log(assay(GeneExpressionMatrix, "normalizedCounts")+1)
marker.info <- scoreMarkers(GeneExpressionMatrix, groups = GeneExpressionMatrix$cell_type)

plot_df <- as.data.frame(plot_df)
plot_df$AUC <- NA
plot_df$cohenD <- NA
for(i in 1:nrow(plot_df)){
  plot_df$AUC[i] <- marker.info[[plot_df[["cell type"]][i]]][plot_df$`transcription factor`[i],"mean.AUC"]
  plot_df$cohenD[i] <- marker.info[[plot_df[["cell type"]][i]]][plot_df$`transcription factor`[i],"mean.logFC.cohen"]
}

plot_df$data <- "gene expression"
plot_df_activity <- plot_df

library(SingleCellExperiment)
# restore full GeneExpressionMatrix (with all target genes) to calculate activity
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/GeneExpressionMatrix.rds")

GeneExpressionMatrix <- GeneExpressionMatrix[,GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "FCGR3A+ Mono", "Memory CD4+ T",
                                                                                    "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "CD14+ Mono")]

GeneExpressionMatrix$cell_type <- gsub(".* T$","T cells",GeneExpressionMatrix$cell_type)
GeneExpressionMatrix$cell_type <- gsub(".*Mono.*","Monocytes",GeneExpressionMatrix$cell_type)
activity.matrix <- epiregulon::calculateActivity(expMatrix = GeneExpressionMatrix, exp_assay="normalizedCounts",
                                                 regulon = regulon.w)
activity.matrix <- activity.matrix[, GeneExpressionMatrix$cell_type %in% c("B", "NK", "DC", "Monocytes", "T cells")]

activity.matrix <- activity.matrix[na.omit(match(master_regulators, rownames(activity.matrix))),]
activity.sce <- SingleCellExperiment(assays = list(logcounts = activity.matrix), colData = colData(GeneExpressionMatrix))
library(scran)
marker.info.activity <- scoreMarkers(activity.sce, groups = activity.sce$cell_type)


plot_df_activity$AUC <- NA
plot_df_activity$cohenD <- NA
for(i in 1:nrow(plot_df)){
  plot_df_activity$AUC[i] <- marker.info.activity[[plot_df_activity[["cell type"]][i]]][plot_df_activity$`transcription factor`[i],"mean.AUC"]
  plot_df_activity$cohenD[i] <- marker.info.activity[[plot_df_activity[["cell type"]][i]]][plot_df_activity$`transcription factor`[i],"mean.logFC.cohen"]
}

plot_df_activity$data <- "TF activity"
plot_df <- rbind(plot_df, plot_df_activity)

library(ggplot2)
library(ggbeeswarm)

dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/AUC_differential_expression.pdf", width = 4, height = 6)
ggplot(plot_df, aes(y=AUC, x = data)) +
  geom_bar(position = "identity", fun = "median", stat="summary", fill="#86244f")+
  geom_quasirandom(size=1.4, aes(color = .data[["transcription factor"]]))+
  ylab("Mean AUC")+
  xlab("Data type")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =11, lineheight = 0.8, family = "sans"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(fill = "white", colour="black"))
dev.off()



#############################
### Response to reviewers ###
#############################

library(epiregulon.extra)
library(epiregulon)
library(dsassembly)
library(ArchR)
library(scran)
library(ggrepel)

# VERSION 1
topgenes <- 10
regulon.size.filter <- 35


regulon.w <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_corr_merged_no_clusters.rds")
GeneExpressionMatrix <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/GeneExpressionMatrix.rds")
activity_matrix <- epiregulon::calculateActivity(expMatrix = GeneExpressionMatrix, exp_assay="normalizedCounts",
                                                 regulon = regulon.w)





# differential activity
markers  <- findDifferentialActivity(activity_matrix = activity_matrix,
                                     clusters = GeneExpressionMatrix$cell_type ,
                                     pval.type = "some",
                                     direction = "up",
                                     log.p=TRUE,
)

markers <- lapply(markers, function(x) {names(x) <- gsub("^log\\.","",names(x)); x})

markers.sig <- getSigGenes(markers, direction = "up", topgenes =topgenes)

dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/actvity_all_tfs_v1.pdf", height = 9, width = 6)
plotBubble(activity_matrix = activity_matrix,
           tf = markers.sig$tf,
           clusters = GeneExpressionMatrix$cell_type,
           bubblesize = "summary.logFC",
           direction = "up")
dev.off()

data_matrix <- Seurat::Read10X_h5("/gstore/project/epigen/PBMC/raw_data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")[["Gene Expression"]]
colnames(data_matrix) <- paste0("PBMC_10k#", colnames(data_matrix))
data_matrix <- data_matrix[rownames(GeneExpressionMatrix),colnames(GeneExpressionMatrix)]
sce <- GeneExpressionMatrix
assay(sce) <- data_matrix
names(assays(sce)) <- "counts"
# prepare a downsampled gene expression matrix for computing logFC of target genes
sce <- logNormCounts(sce)

# order cell_types
sce$cell_type <- factor(sce$cell_type, levels = c("Memory CD4+ T",
                                                  "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T",
                                                  "Monocytes", "FCGR3A+ Mono", "CD14+ Mono","B", "NK", "DC"))

# add logFC onto regulons
regulon.w <- addLogFC(expMatrix = sce,
                      clusters = sce$cell_type,
                      regulon = regulon.w,
                      sig_type = "FDR")


# filter regulons based on logFC
regulon.w.filtered <- regulon.w[purrr::pmap_lgl(regulon.w,function(...) any(unlist(list(...)[grep("logFC",names(list(...)))]>0.5)*unlist(list(...)[grep("FDR",names(list(...)))]<0.05))), ]

# calculate altered regulon size
regulon.size <- table(regulon.w.filtered$tf)

# filter top TFs based on altered regulon size
markers.sig$regulon.size <-  as.vector(regulon.size[markers.sig$tf])
markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]

markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]


bubble <- plotBubble(activity_matrix = activity_matrix,
                     tf = markers.sig$tf,
                     clusters = GeneExpressionMatrix$cell_type ,
                     bubblesize = "summary.logFC", direction = "up")



# assemble a sce with colData info
.colData <- colData(GeneExpressionMatrix)
sce.dummy <- SingleCellExperiment(colData=.colData)
downsample_seq <- seq(from = 1, to = ncol(sce.dummy), by = floor(max(1,
                                                                     ncol(sce.dummy)/1000)))

activity_matrix_2 <- activity_matrix[match(unique(markers.sig$tf), rownames(activity_matrix)),]
sce.dummy <- sce.dummy[, downsample_seq]

cell_type_colors <- list(cell_type = c("B" = "#1F78B4", "CD14+ Mono" = "#B2DF8A", "Memory CD4+ T" = "#E31A1C", "DC" = "darkmagenta",
                      "Naive CD8+ T" = "#E31A1C", "FCGR3A+ Mono" = "#B2DF8A", "Monocytes" = "#B2DF8A",
                      "Naive CD4+ T" = "#E31A1C", "Memory CD8+ T" = "#E31A1C", "NK" = "wheat2"))
top_annotation <- data.frame(colData(sce.dummy)["cell_type"])
top_annotation$cell_type <- factor(top_annotation$cell_type, levels = c("Memory CD4+ T",
                                                                        "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "Monocytes", "FCGR3A+ Mono", "CD14+ Mono","B", "NK", "DC"))

row_annotation <- data.frame(transcription_factor=rownames(activity_matrix))

activity_matrix_2 <- activity_matrix_2[unique(markers.sig$tf), downsample_seq]
activity_matrix_2 <- Matrix::t(scale(Matrix::t(activity_matrix_2),
                                   scale = TRUE, center = TRUE))


col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
tfplot <- ComplexHeatmap::Heatmap(activity_matrix_2, col = col_fun,
                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation, col = cell_type_colors),
                                  #right_annotation = ComplexHeatmap::rowAnnotation(df = row_annotation),
                                  column_split = top_annotation["cell_type"], use_raster = TRUE,
                                  raster_quality = 10, cluster_rows = FALSE,
                                  cluster_columns = FALSE, border = TRUE,
                                  show_column_names = FALSE,
                                  column_title_rot = 45)




gex <- assay(sce, "logcounts")
gex <- as(gex,"dgCMatrix")
gex <- gex[unique(markers.sig$tf), downsample_seq]
row_annotation <- data.frame(transcription_factor=rownames(gex))
gex <- Matrix::t(scale(Matrix::t(gex), scale = TRUE, center = TRUE))

gexplot <- ComplexHeatmap::Heatmap(gex, col = col_fun,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation, col = cell_type_colors),
                                   # right_annotation = ComplexHeatmap::rowAnnotation(df = row_annotation),
                                   column_split = top_annotation["cell_type"], use_raster = TRUE,
                                   raster_quality = 10, cluster_rows = FALSE,
                                   cluster_columns = FALSE, border = TRUE,
                                   show_column_names = FALSE,
                                   column_title_rot = 45)

pdf(paste0("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/actvity_filtered_tfs_v1.pdf"), height = min(20, round(length(unique(markers.sig$tf))*0.2)+2))
print(bubble)
print(tfplot)
print(gexplot)
dev.off()


# VERSION 2

topgenes <- 15
regulon.size.filter <- 50

markers.sig <- getSigGenes(markers, direction = "up", topgenes =topgenes)

dev.new()
pdf("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/actvity_all_tfs_v2.pdf", height = 9, width = 6)
plotBubble(activity_matrix = activity_matrix,
           tf = markers.sig$tf,
           clusters = GeneExpressionMatrix$cell_type,
           bubblesize = "summary.logFC",
           direction = "up")
dev.off()

# filter top TFs based on altered regulon size
markers.sig$regulon.size <-  as.vector(regulon.size[markers.sig$tf])
markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]

markers.sig <- markers.sig[which(markers.sig$regulon.size > regulon.size.filter), ]


bubble <- plotBubble(activity_matrix = activity_matrix,
                     tf = markers.sig$tf,
                     clusters = GeneExpressionMatrix$cell_type ,
                     bubblesize = "summary.logFC", direction = "up")



# assemble a sce with colData info
.colData <- colData(GeneExpressionMatrix)
sce.dummy <- SingleCellExperiment(colData=.colData)
downsample_seq <- seq(from = 1, to = ncol(sce.dummy), by = floor(max(1,
                                                                     ncol(sce.dummy)/1000)))

activity_matrix_2 <- activity_matrix[match(unique(markers.sig$tf), rownames(activity_matrix)),]
sce.dummy <- sce.dummy[, downsample_seq]

top_annotation <- data.frame(colData(sce.dummy)["cell_type"])
top_annotation$cell_type <- factor(top_annotation$cell_type, levels = c("Memory CD4+ T",
                                                                        "Memory CD8+ T", "Naive CD4+ T", "Naive CD8+ T", "Monocytes", "FCGR3A+ Mono", "CD14+ Mono","B", "NK", "DC"))

row_annotation <- data.frame(transcription_factor=rownames(activity_matrix))

activity_matrix_2 <- activity_matrix_2[unique(markers.sig$tf), downsample_seq]
activity_matrix_2 <- Matrix::t(scale(Matrix::t(activity_matrix_2),
                                   scale = TRUE, center = TRUE))


col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
tfplot <- ComplexHeatmap::Heatmap(activity_matrix_2, col = col_fun,
                                  top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation, col = cell_type_colors),
                                  #right_annotation = ComplexHeatmap::rowAnnotation(df = row_annotation),
                                  column_split = top_annotation["cell_type"], use_raster = TRUE,
                                  raster_quality = 10, cluster_rows = FALSE,
                                  cluster_columns = FALSE, border = TRUE,
                                  show_column_names = FALSE,
                                  column_title_rot = 45)

gex <- assay(sce, "logcounts")
gex <- as(gex,"dgCMatrix")
gex <- gex[unique(markers.sig$tf), downsample_seq]
row_annotation <- data.frame(transcription_factor=rownames(gex))
gex <- Matrix::t(scale(Matrix::t(gex), scale = TRUE, center = TRUE))

gexplot <- ComplexHeatmap::Heatmap(gex, col = col_fun,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(df = top_annotation, col = cell_type_colors),
                                   #right_annotation = ComplexHeatmap::rowAnnotation(df = row_annotation),
                                   column_split = top_annotation["cell_type"], use_raster = TRUE,
                                   raster_quality = 10, cluster_rows = FALSE,
                                   cluster_columns = FALSE, border = TRUE,
                                   show_column_names = FALSE,
                                   column_title_rot = 45)

pdf(paste0("/gstore/project/epigen/benchmark/PBMC/OUTPUT/plots/actvity_filtered_tfs_v2.pdf"), height = min(20, round(length(unique(markers.sig$tf))*0.2)+2))
print(bubble)
print(tfplot)
print(gexplot)
dev.off()
