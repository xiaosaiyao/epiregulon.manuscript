library(ArchR)

archr.proj <- ArchR::loadArchRProject("OUTPUT/ArchRProject/")

# filter doublets
archr.proj <- filterDoublets(
    ArchRProj = archr.proj,
    cutEnrich = 1,
    cutScore = -Inf,
    filterRatio = 1
)

# filter cells that do not contain rna
archr.proj <- archr.proj[!is.na(archr.proj$Gex_nUMI)]

# add reduced dims
archr.proj <- addIterativeLSI(
    ArchRProj = archr.proj,
    useMatrix = 'TileMatrix',
    name = 'IterativeLSI_TileMatrix',
    threads = 4,
    seed = 2,
    force = TRUE
)

archr.proj <- addIterativeLSI(
    ArchRProj = archr.proj,
    useMatrix = 'GeneExpressionMatrix',
    name = 'IterativeLSI_GeneExpressionMatrix',
    firstSelection = "variable",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    binarize = FALSE,
    threads = 4,
    seed = 2,
    force = TRUE
)

archr.proj <- addCombinedDims(
    archr.proj,
    reducedDims = c('IterativeLSI_TileMatrix', 'IterativeLSI_GeneExpressionMatrix'),
    name =  'IterativeLSI_Combined'
)


saveArchRProject(
    ArchRProj = archr.proj,
    overwrite = TRUE,
    logFile = createLogFile("saveArchRProject"),
    threads = getArchRThreads()
)




# add embeddings

archr.proj <- addUMAP(
    ArchRProj = archr.proj,
    reducedDims = 'IterativeLSI_Combined',
    name = 'UMAP_Combined',
    seed = 2,
    threads = 1
)

archr.proj <- addUMAP(
    ArchRProj = archr.proj,
    reducedDims = 'IterativeLSI_TileMatrix',
    name = 'UMAP_ATAC',
    seed = 2,
    threads = 1,
    force = TRUE
)


archr.proj <- addUMAP(
    ArchRProj = archr.proj,
    reducedDims = 'IterativeLSI_GeneExpressionMatrix',
    name = 'UMAP_RNA',
    seed = 2,
    threads = 1
)

saveArchRProject(
    ArchRProj = archr.proj,
    overwrite = TRUE,
    logFile = createLogFile("saveArchRProject"),
    threads = getArchRThreads()
)

# Plot UMAP embedding
p.umap.sample <- plotEmbedding(
    ArchRProj = archr.proj,
    embedding = 'UMAP_Combined',
    colorBy = 'cellColData',
    name = 'Sample',
    threads = 4
)
plotPDF(p.umap.sample, name = "GPSA-UMAP_Combined-samples.pdf", ArchRProj = archr.proj)

p.umap.doublet <- plotEmbedding(
    ArchRProj = archr.proj,
    embedding = 'UMAP_Combined',
    colorBy = 'cellColData',
    name = 'DoubletEnrichment',
    pal = c("grey", "#FB8861FF", "#B63679FF", "#51127CFF", "#000004FF"),
    threads = 4
)
plotPDF(p.umap.doublet, name = "GPSA-UMAP_Combined-DoubletEnrichment.pdf", ArchRProj = archr.proj)



# Peak calling

library(BSgenome.Hsapiens.Genentech.GRCh38)

archr.proj <- addGroupCoverages(
    ArchRProj = archr.proj,
    groupBy = 'hash_assignment2',
    threads = 4
)


archr.proj <- ArchR.helper::addReproduciblePeakSet(
    ArchRProj = archr.proj,
    groupBy = 'hash_assignment2',
    peakMethod = "MACSr",
    excludeChr = c('chrMT','chrY'),
    genomeSize = 2.7e9,
    threads = 4,
    force = TRUE
)

archr.proj <- addPeakMatrix(
    ArchRProj = archr.proj,
    binarize = FALSE,
    threads = 4,
    force = TRUE
)

saveArchRProject(
    ArchRProj = archr.proj,
    overwrite = TRUE,
    logFile = createLogFile("saveArchRProject"),
    threads = getArchRThreads()
)

# TF annotation
peaks.anno <- scMultiome::tfBinding()
archr.proj <- addPeakAnnotations(
    ArchRProj = archr.proj,
    regions = peaks.anno,
    name = 'TF_peaks',
    force = TRUE
)

archr.proj <- addDeviationsMatrix(
    ArchRProj = archr.proj,
    peakAnnotation = 'TF_peaks',
    matrixName = 'TFPeaksDeviationsMatrix',
    threads = 4,
    force = TRUE
)

saveArchRProject(
    ArchRProj = archr.proj,
    overwrite = TRUE,
    logFile = createLogFile("saveArchRProject"),
    threads = getArchRThreads()
)


# motif annotation
archr.proj  <- addMotifAnnotations(ArchRProj = archr.proj, motifSet = 'cisbp', name = 'Motif', species='Homo sapiens')

archr.proj <- addDeviationsMatrix(
    ArchRProj = archr.proj,
    peakAnnotation = 'Motif',
    threads = 1,
    force = TRUE
)

saveArchRProject(
    ArchRProj = archr.proj,
    overwrite = TRUE,
    logFile = createLogFile("saveArchRProject"),
    threads = getArchRThreads()
)

# add bigwigs

getGroupBW(
    ArchRProj = archr.proj,
    groupBy = 'hash_assignment2',
    normMethod = 'ReadsInTSS',
    threads = 1
)

