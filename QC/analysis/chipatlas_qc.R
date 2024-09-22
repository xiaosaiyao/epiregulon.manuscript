# chipatlas QC
chipatlas_metadata <- read.delim("/gne/data/genomics/external_sources/chipAtlas/chipAtlas/experimentList.tab", header=FALSE)
chipatlas_metadata <- chipatlas_metadata[chipatlas_metadata$V2 == "hg38",]
chipatlas_metadata <- chipatlas_metadata[chipatlas_metadata$V3 == "TFs and others",]
chipatlas_metadata_qc <- data.frame(do.call(rbind, strsplit(chipatlas_metadata$V8, split = ",")))
chipatlas_metadata_qc <- apply(chipatlas_metadata_qc, 2, as.numeric)
colnames(chipatlas_metadata_qc) <- c("Number of reads", "Percentage mapped", "Percentage duplicates", "Number of peaks q < 1E-05")


chipatlas_metadata_qc <- data.frame(TF = chipatlas_metadata$V4,
                                    data_source = "ChIP-Atlas",
                                    experiment_id = chipatlas_metadata$V1,
                                    sample = chipatlas_metadata$V6,
                                    tissue = chipatlas_metadata$V5,
                                    chipatlas_metadata_qc)

chipatlas_metadata_qc$unique.reads <- round(chipatlas_metadata_qc[,"Number.of.reads"]*
    chipatlas_metadata_qc[,"Percentage.mapped"]*
    (100-chipatlas_metadata_qc[,"Percentage.duplicates"])/10000)

# QC

chipatlas.pass <- chipatlas_metadata_qc[which(chipatlas_metadata_qc$unique.reads >= 10000000 &
                                                  chipatlas_metadata_qc$Percentage.mapped >= 70 &
                                                  chipatlas_metadata_qc$Number.of.peaks.q...1E.05 >= 100), ]

write.table(chipatlas.pass, "QC/output/chipatlas.pass.txt", sep="\t", quote = FALSE, row.names = FALSE)
chipatlas.fail <- chipatlas_metadata_qc[which(!chipatlas_metadata_qc$experiment_id %in% chipatlas.pass$experiment_id), ]
# stats after filtering
summary(chipatlas.pass$Number.of.peaks.q...1E.05)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 101    1370    4933   13240   18831  247614


length(unique(chipatlas.pass$TF)) #  1709
length(unique(chipatlas.pass$experiment_id)) #23745

# compare to before filtering
length(unique(chipatlas_metadata_qc$TF)) #1805
length(unique(chipatlas_metadata_qc$experiment_id)) #29178

summary(as.numeric(unname(table(chipatlas.pass$TF))))

############
encode_metadata <- read.delim("/gne/data/genomics/external_sources/chipAtlas/encode/human/raw/metadata.tsv", header=TRUE)
encode_metadata <- encode_metadata[encode_metadata$File.assembly == "GRCh38",]
encode_metadata <- encode_metadata[encode_metadata$Assay == "TF ChIP-seq",]


###########
merged <- data.frame(readxl::read_xlsx("QC/data/epiregulon_chipseq_table.xlsx", sheet=3))
summary(merged$coverage....of.genome.)
merged.motif <- merged[merged$peaks.with.motifs>0,]
summary(merged.motif$X..peaks.with.motifs)
