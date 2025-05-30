#load ArchR proj
library(ArchR)
seRNA_final <- readRDS("OUTPUT/seRNA_final.rds")
seRNA_final$hash_assignment2 <- paste0(seRNA_final$library, seRNA_final$hash_assignment)
proj <- ArchR::loadArchRProject("OUTPUT/ArchRProject/")
common <- intersect(proj$cellNames, colnames(seRNA_final))
proj <- proj[common,]

# add HTO information

for (row_data in colnames(colData(seRNA_final))){
    proj <- addCellColData(
        ArchRProj = proj,
        data = colData(seRNA_final)[common, row_data],
        cells = common,
        name = row_data,
        force = TRUE
    )
}


#add extra  cell information
sample_info <- read.csv("/gstore/project/ar_ligands/AR/scRNAseq/pipeline_nonpipeline/data/HTO_SAMID.csv")
sample_info$SAMID_HTO <- paste0(sample_info$SAMID, sample_info$HTO)

proj$TREATMENT <- sample_info$TREATMENT[match(proj$hash_assignment2, sample_info$SAMID_HTO)]
proj$Cell <- unlist(lapply(strsplit(proj$TREATMENT, split = "-"),"[",1))
proj$TEST_ARTICLE <- unlist(lapply(strsplit(proj$TREATMENT, split = "-"),"[",2))

proj$Cell[proj$Cell == "22RV1"] <- "22Rv1"

#filter out doublets and non-confident calls
proj <- proj[which(proj$Confident == TRUE & proj$Doublet == FALSE), ]

#filter out unwanted samples
proj <- proj[!proj$hash_assignment2 %in% c(paste0("SAM24425416HTO-",14:16),
                                           paste0("SAM24425417HTO-",1:10),
                                           paste0("SAM24428812HTO-",1:6)),]

proj <- proj[proj$TEST_ARTICLE %in% c("DMSO", "Enza","ARV110", "A9690"),]

saveArchRProject(proj, outputDirectory = "OUTPUT/ArchRProject/")





