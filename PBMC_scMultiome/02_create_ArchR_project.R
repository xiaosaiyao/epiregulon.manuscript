suppressPackageStartupMessages(library(ArchR))
addArchRGenome("hg38")
addArchRThreads(10)

#Get Input Fragment Files
inputFiles <- getInputFiles("/gstore/project/epigen/PBMC/raw_data/")
names(inputFiles) <- "PBMC_10k"

#Create Arrow Files (disabled here)
createArrowFiles(inputFiles, force = TRUE)

ArrowFiles <- "/gstore/project/epigen/PBMC/PBMC_10k.arrow"

#ArchRProject
proj <- ArchRProject(ArrowFiles)

saveArchRProject(ArchRProj = proj, outputDirectory = "saved_project", load = FALSE)

seRNA <- import10xFeatureMatrix(
  input = c("/gstore/project/epigen/PBMC/raw_data/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"),
  names = c("PBMC_10k")
)

# filter out genes which are expressed in less than 3 cells
seRNA <- seRNA[colSums(assay(seRNA))>2,]

proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

#Filter Cells
proj <- proj[proj$TSSEnrichment > 6 & proj$nFrags > 2500 & !is.na(proj$Gex_nUMI)]

#Doublet Filtration
proj <- addDoubletScores(proj)
proj <- filterDoublets(proj)

#LSI-ATAC
proj <- addIterativeLSI(
  ArchRProj = proj,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

#LSI-RNA
proj <- addIterativeLSI(
  ArchRProj = proj,
  clusterParams = list(
    resolution = 0.2,
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix",
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA"
)

#Combined Dims
proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

#UMAPs
proj <- addUMAP(proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)

proj <- addUMAP(proj, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)

proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

#Add Clusters
proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)
#proj <- addClusters(proj, reducedDims = "LSI_RNA", name = "Clusters_genes", resolution = 0.4, force = TRUE)

#Plot Embedding
p1 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_ATAC", size = 1.5, labelAsFactors=F, labelMeans=F)

p2 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_RNA", size = 1.5, labelAsFactors=F, labelMeans=F)

p3 <- plotEmbedding(proj, name = "Clusters", embedding = "UMAP_Combined", size = 1.5, labelAsFactors=F, labelMeans=F)

p1

p2

p3

plotPDF(p1, p2, p3, name = "UMAP-scATAC-scRNA-Combined", addDOC = FALSE)

saveArchRProject(ArchRProj = proj, outputDirectory = "saved_project", load = FALSE)
