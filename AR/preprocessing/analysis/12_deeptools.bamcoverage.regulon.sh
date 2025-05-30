#!/bin/bash
#SBATCH -n 4
#SBATCH -N 1
#SBATCH --mem=120G
#SBATCH --qos=medium

module add deeptools
dir=/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline-rerun/OUTPUT/ArchRProject/Deeptools
cd $dir
wigdir=/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline-rerun/data
beddir=/gstore/project/ar_ligands/AR/scRNAseq/nonpipeline-rerun/OUTPUT/ArchRProject/Epiregulon


declare -A AR_all
AR_all[MDA]=regulon.MDA.chip.AR.all.bed


MDA_AR_DMSO=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_DMSO/croo_output/signal/rep1/05_0FP2_0255Genen_DMSO-3a_AR_hs_i37_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig
MDA_AR_A9690=/gstore/data/genomics/congee_rest_runs/67422e60aa83427048c43a8a/AR_A9690/croo_output/signal/rep1/07_0FP4_0255Genen_A9690-4a_AR_hs_i43_R1.srt.nodup_x_00_0FPL_0255Genen_Pooled_Input_hs_i67_R1.srt.nodup.fc.signal.bigwig


for cell in MDA
do
computeMatrix reference-point --referencePoint center -S $MDA_AR_DMSO $MDA_AR_A9690 -b 500 -a 500 -R $beddir/${AR_all[$cell]}  --skipZeros -o $dir/AR_regulon_"$cell".all.gz -p max/2
plotHeatmap -m $dir/AR_regulon_"$cell".all.gz -out $dir/AR_regulon_"$cell".all.pdf  --colorMap  'Purples' 'Purples' --samplesLabel 'AR_DMSO' 'AR_A9690'   --regionsLabel 'all' --plotTitle "${cell}" --zMax 15 15 --yMax 15 15
done

