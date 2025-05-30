library(dsassembly)
library(scater)
library(ArchR)

number <-  commandArgs(trailingOnly=TRUE)
print(number)
part_no <- paste0('part_', number)
sce <- dsassembly::getExperiment('DS000016526', part_no)
assay(sce,1) <- as(assay(sce,1), "dgCMatrix")
rownames(sce) <- rowData(sce)$symbol

proj <- ArchR::loadArchRProject(paste0("/gstore/scratch/u/yaox19/TGI/No_866_ATAC/", part_no, "/output_TSS4_Frags1000"))
getAvailableMatrices(proj)

proj$meta_pieceID <- sapply(strsplit(as.character(proj@cellColData$Sample), "_"),"[",1)
ATAC_samples <- unique(proj$meta_pieceID)
RNA_samples <- unique(sce$meta_pieceID)

common_samples <- intersect(RNA_samples, ATAC_samples)

proj <- proj[proj$cellNames[which(proj$meta_pieceID %in% common_samples)],]
sce <-  sce[,which(colData(sce)$meta_pieceID %in% common_samples)]

generateGroupList <- function(meta_pieceID, sce, proj){
    cells <- list(ATAC=rownames(proj@cellColData)[which(proj$meta_pieceID == meta_pieceID)],
                  RNA=colnames(sce)[which(sce$meta_pieceID == meta_pieceID)])
    assign(meta_pieceID, cells)
}

groupList <- lapply(unique(proj$meta_pieceID), generateGroupList, sce, proj)
names(groupList) <- unique(proj$meta_pieceID)

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = sce,
    addToArrow = TRUE,
    force = TRUE,
    groupATAC = "meta_pieceID",
    groupRNA = "meta_pieceID",
    groupList = groupList,
    embeddingATAC = getEmbedding(proj),
    useImputation = FALSE
)

saveArchRProject(
   ArchRProj = proj,
   #outputDirectory = "/gstore/scratch/u/yaox19/TGI/ATAC_RNA/ccRCC",
   overwrite = TRUE,
   load = TRUE,
   dropCells = FALSE,
   logFile = createLogFile("saveArchRProject"),
   threads = getArchRThreads()
)

