library(reshape2)
library(vegan)
library(ggplot2)
library(gridExtra)
library(depmap)
library(ggrepel)
library(pheatmap)


#get chronos score and sample metadata

chronos <- depmap::depmap_crispr()
sampledata <- depmap::depmap_metadata()

#get a list of TF factor
d.TF <- data.table::fread(
    file.path("/gstore/project/lineage/data/TF_data/DatabaseExtract_v_1.01.txt"), sep = "\t",header = TRUE)
colnames(d.TF) <- make.names(colnames(d.TF))
d.TF.ensg <- d.TF %>% filter(Is.TF. == "Yes")  %>% dplyr::select(Ensembl.ID) %>% pull
d.TF.symbol <- d.TF$HGNC.symbol[which(d.TF$Is.TF. == "Yes")]

#filter for TF
chronos.TF <- chronos[chronos$gene_name %in% d.TF.symbol,]

# merge sample data with chronos
chronos.TF <- merge(chronos.TF, sampledata, by = "depmap_id")

total_count <- chronos.TF %>% group_by (lineage) %>% dplyr::count(lineage, gene_name)
total_count_mat <- dcast(total_count, lineage ~ gene_name)

#filter out rare lineages
total_count_mat <- na.omit(total_count_mat[total_count_mat$ADNP >10, ])
rownames(total_count_mat) <- total_count_mat$lineage
total_count_mat <- total_count_mat[,-1]

# count total number of TF with ceres < -1
chronos.TF$chronos <- as.numeric(chronos.TF$dependency)

cutoff_count <- chronos.TF %>% dplyr::filter(dependency < -1.0) %>% group_by (lineage) %>%
    dplyr::count(lineage, gene_name)

cutoff_count_mat <- data.frame(matrix(0, nrow=nrow(total_count_mat), ncol=ncol(total_count_mat)))

colnames(cutoff_count_mat) <- colnames(total_count_mat)
rownames(cutoff_count_mat) <- rownames(total_count_mat)

cutoff_count <- na.omit(cutoff_count)
cutoff_count <- cutoff_count[cutoff_count$lineage %in% rownames(total_count_mat),]

for (i in 1:nrow(cutoff_count)){
    tissue <- as.character(cutoff_count[i,"lineage"])
    TF <- as.character(cutoff_count[i,"gene_name"])
    cutoff_count_mat[tissue,TF] <- cutoff_count[i,"n"]
}

fraction <- cutoff_count_mat/total_count_mat

#remove all 0 columns
colmaxs <- apply(fraction,2, max)
keep <- which(colmaxs > 0.05 )
fraction <- fraction[,keep]
diversity <- diversity(fraction, index = "shannon", MARGIN = 2)
diversity <- diversity[order(diversity)]
diversity_df <- data.frame(TF=names(diversity), diversity=diversity, rank=1:length(diversity))

#plot genes by shannon diversity score


#find min value of lineage factor
TF_lineage_summary <- chronos.TF %>% group_by(gene_name, lineage) %>%
    summarize(chronos_median = median(dependency,na.rm=T),
              chronos_min = min(dependency,na.rm=T),
              chronos_10percent=quantile(dependency, probs=0.1,na.rm=T) )

TF_summary <- TF_lineage_summary %>% group_by(gene_name) %>%
    summarize(chronos_median = min(chronos_median),
              chronos_min = min(chronos_min),
              chronos_10percent = min(chronos_10percent))

chronos.min <- TF_lineage_summary %>% group_by(gene_name) %>% slice(which.min(chronos_min))
chronos_10percent <- TF_lineage_summary %>% group_by(gene_name) %>% slice(which.min(chronos_10percent))

diversity_df <- cbind(diversity_df, TF_summary[match(diversity_df$TF, TF_summary$gene_name),])


pdf("OUTPUT/chronos.tissue.pdf", width=12)
#diversity_df_select = diversity_df %>% filter(TF %in% c("BCL6", "MYOD1","ESR1", "FOXA1","PAX8","SOX10","MYB","CENPA","CTCF","MYC"))
diversity_df_select = diversity_df %>% filter(chronos_10percent < -1.5)
ggplot(diversity_df, aes_string(x="rank", y="diversity",  color="chronos_10percent")) + ylim(0,5) +
    geom_point(size=1) +
    geom_text_repel(
        data = diversity_df_select,
        aes_string(label = "TF"),
        size = 4,
        nudge_x=5,
        nudge_y=1,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.7, "lines"),
        min.segment.length=unit(0, "lines"),
        max.overlaps = Inf
    ) +
    scale_colour_gradient(low="yellow", high="blue") +
    theme_classic() +
    ylab("shannon diversity") + xlab("least diverse                                                 most diverse") +
    theme(text = element_text(size=20),
          axis.text = element_text(size = 16),
          panel.spacing = unit(0.5, "lines"),
          plot.margin = margin(2, 2, 2, 2, "cm"))


diversity_df_select = diversity_df %>% filter(chronos_min < -1.5)
ggplot(diversity_df, aes_string(x="rank", y="diversity",  color="chronos_min")) + ylim(0,5) +
    geom_point(size=1) +
    geom_text_repel(
        data = diversity_df_select,
        aes_string(label = "TF"),
        size = 4,
        nudge_x=5,
        nudge_y=1,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.7, "lines"),
        min.segment.length=unit(0, "lines"),
        max.overlaps = Inf
    ) +
    scale_colour_gradient(low="yellow", high="blue") +
    theme_classic() +
    ylab("shannon diversity") + xlab("least diverse                                                 most diverse") +
    theme(text = element_text(size=20),
          axis.text = element_text(size = 16),
          panel.spacing = unit(0.5, "lines"),
          plot.margin = margin(2, 2, 2, 2, "cm"))


#plot heatmap
topTF <- fraction[,names(head(diversity,60))]

pheatmap(topTF,
         #annotation_col = annotation_col,
         #labels_col=labels,
         cluster_rows = F,
         cluster_cols = T,
         color = colorRampPalette(c("white", "blue"))(100),
         breaks = seq(0, 0.2, length.out = 101),
         main = "tissue-specific factor",
         fontsize_row = 12,
         treeheight_row=0,
         treeheight_col=0,
         legend = T

)

dev.off()

# plot examples

genes_to_plot=c("MYB","PAX8","FOXA1","MYCN","STAT3","FLI1")
chronosplot=list()
soft_tissue=c(rep("black",13), "red",rep("black",3))
breast=c(rep("black",4), "red",rep("black",12))
color_tissue=list(MYOD1=soft_tissue, MYOG=soft_tissue, PAX3=soft_tissue,
                  FOXA1=breast, ESR1=breast, TFAP2C=breast)


for (gene in genes_to_plot){

    chronosplot[[gene]]=ggplot(chronos.TF %>% filter(symbol==gene & lineage %in% rownames(total_count_mat)) , aes(x=lineage, y=chronos)) +
        geom_point(size=0.5) + coord_flip() + theme_classic() + ggtitle(gene) +
        geom_hline(yintercept=-1, linetype="dashed", color = "red", size=1) +
        theme(plot.title = element_text(hjust = 0.5)) #, axis.text.y = element_text( colour = color_tissue[[gene]]))
}


grid.arrange(grobs=chronosplot,nrow=2)
dev.off()
