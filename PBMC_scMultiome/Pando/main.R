library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Pando)
library(foreach)
library(Signac)

load(".RData")

data(motifs)

DefaultAssay(Seurat_obj) <- "RNA"
# data normalization
Seurat_obj <- SCTransform(Seurat_obj)

DefaultAssay(Seurat_obj) <- "peaks"

Seurat_obj <- RunTFIDF(Seurat_obj)

Seurat_obj <- FindTopFeatures(Seurat_obj, min.cutoff = "q0") # use all features
# complete LSI by performing single value decomposition
Seurat_obj <- RunSVD(Seurat_obj)

Seurat_obj <- Seurat::FindVariableFeatures(Seurat_obj, assay='SCT')

DefaultAssay(Seurat_obj) <- "SCT"

Seurat_obj <- as(Seurat_obj, "SeuratPlus")

missing_sequences <- setdiff(as.character(levels(Seurat_obj@assays$peaks@ranges@seqnames@values)),
                             names(BSgenome.Hsapiens.UCSC.hg38))
peaks_assay <- GetAssay(Seurat_obj, assay = "peaks")
Seurat_obj@assays$peaks <- subset(peaks_assay, features = rownames(peaks_assay)[!as.character(seqnames(Seurat_obj@assays$peaks@ranges)) %in% missing_sequences])

# Initiate GRN object and select candidate regions
Seurat_obj <- initiate_grn(Seurat_obj,
                           rna_assay = "SCT",
                           exclude_exons = FALSE)

# Scan candidate regions for TF binding motifs
Seurat_obj <- find_motifs(
    Seurat_obj,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer gene regulatory network (only variable features are taken)
Seurat_obj <- infer_grn(Seurat_obj,
                        parallel = TRUE)
