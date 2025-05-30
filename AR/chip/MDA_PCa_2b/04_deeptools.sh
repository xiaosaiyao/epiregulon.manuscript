#!/bin/bash
#! -n 4
#! -N 1
#! --mem=120G
#! --qos=medium

module add deeptools
#source activate my_root
dir=/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/deeptools

cd $dir

# peaks
allBRG1=/gstore/project/ar_ligands/AR/chip/LAB11158_A9690/OUTPUT/diffbind/diffbind.allSMARCA4.DBA_DESEQ2.bed

# bigwigs
SMARCA4_DMSO_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/SMARCA4_DMSO/croo_output/signal/rep1/01_0FP6_0255Genen_DMSO-1a_BRG1_hs_i48_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
SMARCA4_DMSO_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/SMARCA4_DMSO/croo_output/signal/rep2/02_0FP7_0255Genen_DMSO-1b_BRG1_hs_i49_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
SMARCA4_A9690_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/SMARCA4_A9690/croo_output/signal/rep1/03_0FP8_0255Genen_A9690-2a_BRG1_hs_i50_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
SMARCA4_A9690_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/SMARCA4_A9690/croo_output/signal/rep2/04_0FP9_0255Genen_A9690-2b_BRG1_hs_i52_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
AR_DMSO_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_DMSO/croo_output/signal/rep1/05_0FP2_0255Genen_DMSO-3a_AR_hs_i37_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
AR_DMSO_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_DMSO/croo_output/signal/rep2/06_0FP3_0255Genen_DMSO-3b_AR_hs_i39_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
AR_A9690_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_A9690/croo_output/signal/rep1/07_0FP4_0255Genen_A9690-4a_AR_hs_i43_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
AR_A9690_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_A9690/croo_output/signal/rep2/08_0FP5_0255Genen_A9690-4b_AR_hs_i45_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
FOXA1_DMSO_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/FOXA1_DMSO/croo_output/signal/rep1/09_0FPA_0255Genen_DMSO-5a_FOXA1_hs_i56_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
FOXA1_DMSO_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/FOXA1_DMSO/croo_output/signal/rep2/10_0FPB_0255Genen_DMSO-5b_FOXA1_hs_i57_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
FOXA1_A9690_1=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/FOXA1_A9690/croo_output/signal/rep1/11_0FPC_0255Genen_A9690-6a_FOXA1_hs_i59_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
FOXA1_A9690_2=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/FOXA1_A9690/croo_output/signal/rep2/12_0FPD_0255Genen_A9690-6b_FOXA1_hs_i64_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig


computeMatrix reference-point --referencePoint center -S $SMARCA4_DMSO_1 $SMARCA4_A9690_1 $AR_DMSO_1 $AR_A9690_1 $FOXA1_DMSO_1 $FOXA1_A9690_1 -b 500 -a 500 -R $allBRG1 --skipZeros -o BRG1.all.gz -p max/2
plotHeatmap -m BRG1.all.gz -out BRG1.all.pdf --samplesLabel "SMARCA4_DMSO" "SMARCA4_A9690" "AR_DMSO" "AR_A9690" "FOXA1_DMSO" "FOXA1_A9690" "ATAC_DMSO" "ATAC_A9690" --regionsLabel "all binding sites" --zMax 10 10 15 15 20 20 0.1 0.1 --yMax 10 10 15 15 20 20 0.1 0.1 --colorMap 'Greens' 'Greens' 'Purples' 'Purples' 'Oranges' 'Oranges' 'Blues' 'Blues'

