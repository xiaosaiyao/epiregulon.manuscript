library(maw.utils)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

samplesheet <- read.delim("OUTPUT/samplesheet.csv")

extraCols_narrowPeak <- c(singnalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

samplesheet$Factor <- gsub("BRG1/","", samplesheet$Factor)
for (TF in unique(samplesheet$Factor)){
    peak_files <- samplesheet[samplesheet$Factor==TF, "Peaks"]
    peaks <- list()
    for(peak_file in peak_files){
        peaks <- c(peaks, rtracklayer::import.bed(peak_file, extraCols = extraCols_narrowPeak))
    }

    peaks <- lapply(peaks, GRanges)
    peaks <- lapply(peaks, function(x) x[x$qValue>5,])
    peaks <- do.call(c, peaks)
    peaks <- GenomicRanges::reduce(peaks)
    export.bed(peaks, paste0("OUTPUT/", TF, "_peaks.bed"))
}







