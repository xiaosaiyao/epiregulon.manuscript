library(rtracklayer)
library(eulerr)
library(Vennerable)
library(ChIPpeakAnno)
library(ggplotify)

# plot venn diagrams to show overlap of AR targets
LNCaP_regulon_matched <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.LNCaP.chip.rds"))
VCaP_regulon_matched <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.VCaP.chip.rds"))
MDA_regulon_matched <- readRDS(paste0("OUTPUT/ArchRProject/Epiregulon/regulon.w.MDA.chip.rds"))

LNCaP_AR_targets <- unique(LNCaP_regulon_matched$target[LNCaP_regulon_matched$tf=="AR"])
VCaP_AR_targets <- unique(VCaP_regulon_matched$target[VCaP_regulon_matched$tf=="AR"])
MDA_AR_targets <- unique(MDA_regulon_matched$target[MDA_regulon_matched$tf=="AR"])
venn_diagram <- euler(combinations = list(VCaP=VCaP_AR_targets,
                                          MDA=MDA_AR_targets))
venn_diagram <- plot(venn_diagram, main = "AR target", quantities = list(type = c("percent", "counts")))

print(venn_diagram)


# plot venn diagrams to show overlap of AR binding sites
grl.all <- readRDS("/gne/data/genomics/external_sources/chipAtlas/chipAtlas/chip_peaks/grl.chipatlas.sample.hg38.qc.rds")
VCaP_AR <- grl.all[["VCaP"]]$AR
MDA_AR <- import("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allAR.DBA_DESEQ2.bed")

venn_cnt2venn <- function(venn_cnt){
    n <- which(colnames(venn_cnt)=="Counts") - 1
    SetNames=colnames(venn_cnt)[1:n]
    Weight=venn_cnt[,"Counts"]
    names(Weight) <- apply(venn_cnt[,1:n], 1, paste, collapse="")
    Venn(SetNames=SetNames, Weight=Weight)
}
binding_venn <- makeVennDiagram(list(VCaP_AR,MDA_AR),
                           NameOfPeaks=c("VCaP","MDA"))
binding_venn  <- venn_cnt2venn(binding_venn$vennCounts)

pdf("OUTPUT/ArchRProject/Epiregulon/VCaP.vs.MDA.pdf")
plot(venn_diagram)
plot(binding_venn)
dev.off()


