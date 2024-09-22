
full.regulon <- readRDS("reprogramseq/preprocessing/OUTPUT/reprogram.seq.regulon.rds")
pruned.regulon <- readRDS("reprogramseq/preprocessing/OUTPUT/reprogram.seq.pruned.regulon.rds")


dim(full.regulon)
#743656
dim(pruned.regulon)
#124923

length(unique(full.regulon$tf))
#1410
length(unique(pruned.regulon$tf))
#1201


# add annotation
merged <- data.frame(readxl::read_xlsx("QC/data/epiregulon_chipseq_table.xlsx", sheet=3))

# add motif annotations
TF.motif <- merged$regulator[merged$X..peaks.with.motifs >0]
pruned.regulon$motif <- 0
pruned.regulon$motif[pruned.regulon$tf %in% TF.motif] <- 1

pruned.regulon.motif <- pruned.regulon[pruned.regulon$motif==1,]
pruned.regulon.nomotif <- pruned.regulon[pruned.regulon$motif==0,]
summary(as.numeric(table(pruned.regulon.motif$tf)))
summary(as.numeric(table(pruned.regulon.nomotif$tf)))
t.test(table(pruned.regulon.motif$tf), table(pruned.regulon.nomotif$tf))

# add activator/repressor annotations
activator <- merged$regulator[merged$annotation == "activator"]
repressor <- merged$regulator[merged$annotation == "repressor"]
pruned.regulon$mode <- NA
pruned.regulon$mode[pruned.regulon$tf %in% activator] <- "activator"
pruned.regulon$mode[pruned.regulon$tf %in% repressor] <- "repressor"

pruned.regulon.activator <- pruned.regulon[which(pruned.regulon$mode == "activator"),]
pruned.regulon.repressor <- pruned.regulon[which(pruned.regulon$mode == "repressor"),]
summary(as.numeric(table(pruned.regulon.activator$tf)))
summary(as.numeric(table(pruned.regulon.repressor$tf)))
t.test(table(pruned.regulon.activator$tf), table(pruned.regulon.repressor$tf))


#
source("/gstore/project/crc/summary/COAD/analysis/bulk.visualization.R")
plotTCGA("GATA6", output_name = "QC/OUTPUT/GATA6.TCGA.pdf")
