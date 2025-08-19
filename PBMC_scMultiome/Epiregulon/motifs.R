library(epiregulon)
regulon.w <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_merged_trimmed_no_clusters.rds")
peakMatrix <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/peakMatrix.rds")

regulon.w <- addMotifScore(regulon.w,
                           species="human",
                           genome="hg38",
                           peaks = rowRanges(peakMatrix))

regulon.w$weight[regulon.w$motif==0,] <- 0

saveRDS(regulon.w, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_merged_trimmed_no_clusters_motifs.rds")

regulon.w <- readRDS("/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_blood_trimmed_no_clusters.rds")

regulon.w <- addMotifScore(regulon.w,
                           species="human",
                           genome="hg38",
                           peaks = rowRanges(peakMatrix))

regulon.w$weight[regulon.w$motif==0,] <- 0

saveRDS(regulon.w, "/gstore/project/epigen/benchmark/PBMC/OUTPUT/regulon.w_blood_trimmed_no_clusters_motifs.rds")
