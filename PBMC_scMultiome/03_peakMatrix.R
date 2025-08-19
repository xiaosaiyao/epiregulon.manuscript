#system("pip install macs2")
proj <- loadArchRProject("/gstore/project/epigen/PBMC/saved_project")

pathToMacs2 <- "/mnt/site-library/r440-bioc319-20231108/MACSr/basilisk/env_macs/bin/macs3"
library(BSgenome.Hsapiens.UCSC.hg38)

proj<- addGroupCoverages(ArchRProj = proj, groupBy = "Clusters")

proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  force = TRUE,
  threads = 10
)

# add information about sequence motifs recognized by known transcriptions factors
proj <- addMotifAnnotations(ArchRProj = proj,
                               motifSet = "cisbp", name = "Motif")

proj <- addPeakMatrix(proj)

# add background peaks to be compared against during peak variation assessement
proj <- addBgdPeaks(proj)

# calculate per-cell devations of motif annotations
proj <- addDeviationsMatrix(proj, peakAnnotation = "Motif", force = TRUE)

# save project
saveArchRProject(ArchRProj = proj, outputDirectory = "saved_project", load = FALSE)

