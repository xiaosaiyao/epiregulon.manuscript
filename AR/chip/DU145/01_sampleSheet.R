library(maw.utils)
library(data.table)

# retrieve bam and peak info

FRSID <- c("FRS23855")

# download from Sunrise
sample_data <- read.csv("data/samplesheet.csv")

peaks <- c()
ctrl <- c()
bam <- c()
bigwigs <- c()

id <- FRSID[1]
info <- getFireDBResourceSetInfo(id)
info_files <- info$files

SampleID <- c()

for (i in 1:nrow(info_files)) {


    current_sample <- info_files[i,]

    id <- current_sample$sampleName

    peaks[id] <- sapply(current_sample$uri,
                          function(file){list.files(file.path(file, "peak/rep1"),
                                                    #pattern = ".xls$",
                                                    pattern = "*.bfilt.narrowPeak.gz",
                                                    full.names = TRUE)})


    sampleName <- gsub("(.*?)(SAM\\d*)(_R1.*)", "\\2", peaks[id])

    curr_ctrl <- sample_data[sample_data$experiment_sample_id==sampleName, "control_sample_id"]

    ctrl[id] <- sapply(current_sample$uri,
                         function(file){list.files(file.path(file, "align/ctl1/"),
                                                   pattern = paste0(".*",curr_ctrl, ".*srt.nodup.bam$"),
                                                   recursive = TRUE,
                                                   full.names = TRUE)})


    bam[id] <- list.files(file.path(current_sample$uri, "align/rep1/"), pattern = ".*srt.nodup.bam$",
                           full.names = TRUE)


    bigwigs[id] <- list.files(file.path(current_sample$uri, "signal/rep1"), pattern = ".*fc.signal.bigwig$",
                               full.names = TRUE)


    SampleID <- c(SampleID, sampleName)
}

bam_info <- data.frame(SampleID=SampleID, Peaks=unlist(peaks), bamReads=unlist(bam), bamControl=unlist(ctrl), bigwigs=unlist(bigwigs))
rownames(bam_info) <- SampleID



mapping_file <- "data/metadata.tsv"
mapping <- read.delim2(mapping_file)
samplesheet <- merge(mapping, bam_info, by.x="SampleID", by.y="SampleID")
samplesheet$ScoreCol <- 9
samplesheet$PeakFormat <- "macs2"
samplesheet$PeakCaller <- "macs2"
write.table(samplesheet, "OUTPUT/samplesheet.csv", sep="\t", quote=F, row.names=F)

num_peaks <- c()
for (i in seq_along(peaks)) {
    regions <- rtracklayer::import(peaks[[i]])
    regions <- regions[regions$qValue >=5]
    num_peaks[i] <- length(regions)
}

names(num_peaks) <- names(peaks)
