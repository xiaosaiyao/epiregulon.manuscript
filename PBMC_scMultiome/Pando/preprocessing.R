library(Signac)
library(Seurat)
library(zellkonverter)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(Pando)
library(foreach)

args <- commandArgs(trailingOnly = TRUE)
cores <- as.numeric(args[1])

path_to_fragments <- "/gstore/project/epigen/PBMC/raw_data/filtered_feature_bc_matrix.h5"

load("/gstore/project/epigen/benchmark/data/geneAnnoHg38.rda")
# Get motif data
data(motifs)
annotations <- geneAnnoHg38$genes
overlaps <- GenomicRanges::findOverlaps(geneAnnoHg38$genes, geneAnnoHg38$TSS)
annotations <- annotations[unique(overlaps@from),]
mcols(annotations)$gene_biotype <- "protein_coding" # other possible types are
mcols(annotations)$gene_name <- mcols(annotations)$symbol

counts <- Read10X_h5(path_to_fragments)
Seurat_obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay ="RNA")

Seurat_obj[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"),
                                             fragments = "/gstore/project/epigen/PBMC/raw_data/atac_fragments.tsv.gz",
                                             annotation = annotations)

geneExprMatrix <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/GeneExpressionMatrix.rds")
peakMatrix <- readRDS("/gstore/project/epigen/PBMC/OUTPUT/peakMatrix.rds")

colnames(geneExprMatrix) <- gsub("PBMC_10k#","",colnames(geneExprMatrix))
colnames(peakMatrix) <- gsub("PBMC_10k#","",colnames(peakMatrix))

Seurat_obj <- subset(Seurat_obj, cells = colnames(geneExprMatrix))
geneExprMatrix <- geneExprMatrix[,colnames(Seurat_obj)]

# add hash tag data
Seurat_obj <- AddMetaData(Seurat_obj, metadata = geneExprMatrix$cell_type,
                          col.name = "cell_type")

# add peak counts
peakCounts <- assay(peakMatrix)
peakRanges <- rowRanges(peakMatrix)
peakRanges <- paste(seqnames(peakRanges), ranges(peakRanges) ,sep = "-")
rownames(peakCounts) <- peakRanges
peakCounts <- as(peakCounts, "dgCMatrix")


# adjust cell names to those found in fragment files
#colnames(peakCounts) <- gsub("reprogram#", "", colnames(peakCounts))

Seurat_obj[["peaks"]] <- CreateChromatinAssay(counts = peakCounts,
                                              fragments = "/gstore/project/epigen/PBMC/raw_data/atac_fragments.tsv.gz",
                                              annotation = annotations)

missing_sequences <- setdiff(as.character(levels(Seurat_obj@assays$peaks@ranges@seqnames@values)),
                             names(BSgenome.Hsapiens.UCSC.hg38))

peaks_assay <- GetAssay(Seurat_obj, assay = "peaks")
Seurat_obj@assays$peaks <- subset(peaks_assay, features = rownames(peaks_assay)[!as.character(seqnames(Seurat_obj@assays$peaks@ranges)) %in% missing_sequences])


rm(list=ls()[sapply(ls(), function(x) !x %in% c("Seurat_obj"))])
save.image(".RData")


