#!/bin/bash
#! -n 4
#! -N 1
#! --mem=120G
#! --qos=medium

module load deeptools

BRG1_SMARCA4_DMSO=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/DMSO_BRG1/croo_output/signal/rep1/LIB6235914_SAM24452211_R1.srt.nodup_x_LIB6235913_SAM24452210_R1.srt.nodup.fc.signal.bigwig
FOSL1_DMSO=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/DMSO_FOSL1/croo_output/signal/rep1/LIB6235916_SAM24452213_R1.srt.nodup_x_LIB6235913_SAM24452210_R1.srt.nodup.fc.signal.bigwig
TEAD1_DMSO=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/DMSO_TEAD1/croo_output/signal/rep1/LIB6235917_SAM24452214_R1.srt.nodup_x_LIB6235913_SAM24452210_R1.srt.nodup.fc.signal.bigwig
BRG1_SMARCA4_A9690=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/A9690_BRG1/croo_output/signal/rep1/LIB6235919_SAM24452216_R1.srt.nodup_x_LIB6235918_SAM24452215_R1.srt.nodup.fc.signal.bigwig
FOSL1_A9690=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/A9690_FOSL1/croo_output/signal/rep1/LIB6235921_SAM24452218_R1.srt.nodup_x_LIB6235918_SAM24452215_R1.srt.nodup.fc.signal.bigwig
TEAD1_A9690=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/A9690_TEAD1/croo_output/signal/rep1/LIB6235922_SAM24452219_R1.srt.nodup_x_LIB6235918_SAM24452215_R1.srt.nodup.fc.signal.bigwig
GR_DMSO=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/DMSO_GR/croo_output/signal/rep1/LIB6235915_SAM24452212_R1.srt.nodup_x_LIB6235913_SAM24452210_R1.srt.nodup.fc.signal.bigwig
GR_A9690=/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/A9690_GR/croo_output/signal/rep1/LIB6235920_SAM24452217_R1.srt.nodup_x_LIB6235918_SAM24452215_R1.srt.nodup.fc.signal.bigwig

/gstore/data/genomics/congee_rest_runs/6745698daa83427048c4452e/A9690_TEAD/croo_output/signal/rep1/LIB6235922_SAM24452219_R1.srt.nodup_x_LIB6235918_SAM24452215_R1.srt.nodup.fc.signal.bigwig
computeMatrix reference-point --referencePoint center \
    -S  $BRG1_SMARCA4_DMSO $BRG1_SMARCA4_A9690  $FOSL1_DMSO $FOSL1_A9690 $TEAD1_DMSO $TEAD1_A9690 $GR_DMSO $GR_A9690\
    -R ./OUTPUT/SMARCA4_peaks.bed \
    -b 500 -a 500 \
    --skipZeros -p max/2 \
    -o ./OUTPUT/deeptools/merged_BGR1.gz

plotHeatmap -m ./OUTPUT/deeptools/merged_BGR1.gz \
    -out ./OUTPUT/deeptools/merged_BGR1.pdf \
    --samplesLabel  BRG1_SMARCA4_DMSO BRG1_SMARCA4_A9690 FOSL1_DMSO FOSL1_A9690 TEAD1_DMSO TEAD1_A9690 GR_DMSO GR_A9690 \
    --regionsLabel "merged BGR1/SMARCA4 DMSO and A9690 as reference points" --colorMap 'Greens' 'Greens' 'Purples' 'Purples' 'Oranges' 'Oranges' 'Blues' 'Blues' --missingDataColor 1



