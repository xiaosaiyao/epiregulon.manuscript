GeneExpressionMatrix_GRaNIE <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/GeneExpressionMatrix.rds")
peakMatrix_GRaNIE <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/peakMatrix.rds")


# GRaNIE
library(GenomeInfoDb)
library(GRaNIE)
library(readr)
library(SummarizedExperiment)

rowData(GeneExpressionMatrix_GRaNIE)$name <- gsub("\\.[0-9]{1,2}$", "", rowData(GeneExpressionMatrix_GRaNIE)$name)
new_names <- ensembldb::mapIds(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys= as.vector(rowData(GeneExpressionMatrix_GRaNIE)$name), keytype = "SYMBOL", columns = c("GENEID"))
new_names <- new_names[!is.na(new_names)]
GeneExpressionMatrix_GRaNIE <- GeneExpressionMatrix_GRaNIE[rowData(GeneExpressionMatrix_GRaNIE)$name %in% names(new_names),]
rowData(GeneExpressionMatrix_GRaNIE)$name <- new_names[rowData(GeneExpressionMatrix_GRaNIE)$name]
rownames(GeneExpressionMatrix_GRaNIE) <- rowData(GeneExpressionMatrix_GRaNIE)$name
peakMatrix_peakIDs <- paste0(seqnames(peakMatrix_GRaNIE),":",start(peakMatrix_GRaNIE),"-",end(peakMatrix_GRaNIE))
peakMatrix_assay <- as.matrix(assay(peakMatrix_GRaNIE))
peakMatrix_df <- cbind(data.frame(peakID = peakMatrix_peakIDs), peakMatrix_assay)
GeneExpressionMatrix_assay <- as.matrix(assay(GeneExpressionMatrix_GRaNIE))
GeneExpressionMatrix_data <- tibble::tibble(cbind(data.frame(geneID = rowData(GeneExpressionMatrix_GRaNIE)$name), GeneExpressionMatrix_assay))
sampleMetadata <- cbind(data.frame(sample = colnames(GeneExpressionMatrix_GRaNIE)), data.frame(cell_type = GeneExpressionMatrix_GRaNIE$cell_type))

rm(list=ls()[sapply(ls(), function(x) !x %in% c("sampleMetadata", "GeneExpressionMatrix_data", "peakMatrix_df"))])
save.image(".RData")

