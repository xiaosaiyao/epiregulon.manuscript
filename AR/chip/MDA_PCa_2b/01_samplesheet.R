library(maw.utils)

# retrieve bam and peak info
FRSID <- c("FRS23816")

peaks <- list()
ctrl <- list()
bam <- list()
bigwigs <- list()

rep_num <- 2

for (id in FRSID) {
    info <- getFireDBResourceSetInfo(id)

    for (i in seq_len(rep_num)){
        peaks[[id]][[i]] <- sapply(info$files$uri,
                                   function(file){list.files(file.path(file, "peak", paste0("rep",i)),
                                                             pattern = "bfilt.narrowPeak.gz$",
                                                             full.names = TRUE,
                                                             recursive=FALSE)})
        ctrl[[id]][[i]] <- sapply(info$files$uri,
                                  function(file){list.files(file.path(file, "align", "ctl1"),
                                                            pattern = ".nodup.bam$",
                                                            recursive = TRUE,
                                                            full.names = TRUE)})
        bam[[id]][[i]] <- sapply(info$files$uri,
                                 function(file){list.files(file.path(file, "align", paste0("rep",i)),
                                                           pattern = ".nodup.bam$",
                                                           recursive = TRUE,
                                                           full.names = TRUE)})
        bigwigs[[id]][[i]]<- sapply(info$files$uri,
                                    function(file){list.files(file.path(file, "signal", paste0("rep",i)),
                                                              pattern = "nodup.fc.signal.bigwig$",
                                                              recursive = TRUE,
                                                              full.names = TRUE)})
    }
}

sample_info <- read.delim("data/sampleinfo.txt")

peaks <- sort(unlist(peaks))
peaks_id <- peaks |> basename() |>substr(start=1, stop=10)
matched_names <- sample_info$SampleID[sapply(as.list(peaks_id), grep, sample_info$fastq)]
names(peaks) <- matched_names

ctrl <- unlist(ctrl)
names(ctrl) <- matched_names

bam <- sort(unlist(bam))
names(bam) <- matched_names

bigwigs <- sort(unlist(bigwigs))
names(bigwigs) <- matched_names

bam_info <- data.frame(SampleID=names(peaks), Peaks=peaks, bamReads=bam, bamControl=ctrl, bigwigs=bigwigs)

# merge sample and bam info
info_merge <- merge(bam_info, sample_info, by= "SampleID")

info_merge$ScoreCol <- 9
info_merge$PeakFormat <- "regionPeak"
info_merge$PeakCaller <- "macs2"
#info_merge$Replicate <- as.numeric(as.factor(info_merge$Replicate))

write.csv(info_merge, "OUTPUT/sampleSheet.csv")

num_peaks <- c()
for (i in seq_along(peaks)) {
    regions <- rtracklayer::import(peaks[[i]])
    regions <- regions[regions$qValue >=5]
    num_peaks[i] <- length(regions)
}

names(num_peaks) <- matched_names
