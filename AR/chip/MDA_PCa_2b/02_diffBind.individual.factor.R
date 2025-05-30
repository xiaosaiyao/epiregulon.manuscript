library(DiffBind)
library(rtracklayer)


sampleSheet <- read.csv("/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/sampleSheet.csv", row.names = 1)

for (factor in c("SMARCA4","AR","FOXA1")){
    sampleSheet.factor <- sampleSheet[sampleSheet$Factor==factor,]
    samples <- dba(sampleSheet=sampleSheet.factor, minOverlap=2, filter=5)
    samples <- dba.count(samples, bParallel=TRUE, minOverlap=0)
    samples


    # Normalize background regions using csaw gives best differential
    samples <- dba.normalize(samples,
                             library=DBA_LIBSIZE_BACKGROUND,
                             method=DBA_ALL_METHODS,
                             normalize=DBA_NORM_NATIVE,
                             background=TRUE)

    filename <- paste0("OUTPUT/diffbind/samples", factor, ".rds")
    saveRDS(samples, filename )
    samples <- readRDS(filename)


    contrast <- list()


    samples.block <- dba(samples)
    contrast <- dba.contrast(samples.block,
                             #design="~Treatment",
                             contrast=c("Treatment","A9690","DMSO"))



    contrast <- dba.analyze(contrast,
                            bBlacklist=F,
                            bGreylist=F,
                            method=DBA_ALL_METHODS)

    filename <- paste0("OUTPUT/diffbind/contrast.",factor, ".rds")
    saveRDS(contrast, filename)
    contrast <- readRDS(filename)

    for (diff_method in c("DBA_EDGER","DBA_DESEQ2")){
        # Plot MA plot
        pdf(paste0("OUTPUT/diffbind/MAplot.",
                   diff_method, ".", factor,".pdf"))

        if (diff_method == "DBA_EDGER"){
            dba.plotMA(contrast, th=0.1, yrange=c(-10,10),
                       fold=0, method=DBA_EDGER)

        }else if (diff_method == "DBA_DESEQ2"){
            dba.plotMA(contrast, th=0.1, yrange=c(-10,10),
                       fold=0, method=DBA_DESEQ2)

        }
        dev.off()

        # Generate report

        if (diff_method == "DBA_EDGER"){
            report <- dba.report(contrast, contrast=1, fold=0,
                                 th=0.1, method=DBA_EDGER)
        } else if (diff_method =="DBA_DESEQ2"){
            report <- dba.report(contrast, contrast=1, fold=0,
                                 th=0.1, method=DBA_DESEQ2)
        }

        if (!is.null(report)){
            report.gain <- report[report$Fold>0]
            report.lost <- report[report$Fold<0]
        }


        report.all <- dba.report(contrast, contrast=1, th=1)
        report.all$name <- paste0("peak_",1:length(report.all))

        if (exists("report.gain")){
            if (!is.null(report.gain)){
                message("exporting gained peaks")
                report.gain$name <-
                    paste0("peak_",1:length(report.gain))
                export.bed(report.gain,
                           paste0("OUTPUT/diffbind/diffbind.gain", factor, ".",
                                  diff_method,".bed")
                )
            }
        }

        if (exists("report.lost")) {
            if (!is.null(report.lost)){
                message("exporting lost peaks")
                report.lost$name <-
                    paste0("peak_",1:length(report.lost))
                export.bed(report.lost,
                           paste0("OUTPUT/diffbind/diffbind.lost",factor, ".",
                                  diff_method,".bed")
                )
            }
        }

        if (exists("report.all")) {
            message("exporting all peaks")
            head(report.all)
            export.bed(report.all,
                       paste0("OUTPUT/diffbind/diffbind.all", factor, ".",
                              diff_method,".bed")
            )
        }
    }
}

