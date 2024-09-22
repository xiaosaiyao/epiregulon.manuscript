library(ArchR)
library(parallel)
outdir <- "/gstore/project/ar_ligands/NE/reprogram_seq/multiome_arrayed/OUTPUT/doubletremoved/"
setwd(outdir)
proj_clean <- loadArchRProject(path = outdir, showLogo = TRUE)

#Identify differential peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_clean,
    useMatrix = "PeakMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
# a summarizedExperiment of peak stats by clusters
markersPeaks

#Get statistically significant peaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = T)

# motif enrichment

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj_clean ,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.01 & Log2FC >= 1",
    logFile = "peakAnnoEnrichment"
)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, pal=paletteContinuous(set = "solarExtra"), n = 20, transpose = TRUE, logFile="plotEnrichHeatmap")

pdf("/gstore/project/lineage/manuscript/epiregulon/reprogramseq/OUTPUT/differential.peaks.motif.pdf", height=4)
print(heatmapEM)
dev.off()

