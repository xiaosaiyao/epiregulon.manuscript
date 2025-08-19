#!/bin/sh



# sbatch_output=$(sbatch peak_matrix.sh)
# job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')
#
# specification="1,64,/gstore/project/epigen/benchmark/PBMC/scenicplus,PBMC,scenicplus,peak_matrix,Rosalind,$(date +%Y-%m-%d_%H:%M)"
# if [ -f "/gstore/project/epigen/benchmark/resource_use.txt" ]; then
#     echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
# else
#     colnames="job_id,n_cpu,requested_memtory,path,dataset,tool,part,enviroment,timestamp"
#     echo -e "${colnames}\n${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
# fi
#
# sbatch_output=$(sbatch --dependency=afterok:$job_id adata.sh)
# job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')
#
# specification="1,64,/gstore/project/epigen/benchmark/PBMC/scenicplus,PBMC,scenicplus,adata,Rosalind,$(date +%Y-%m-%d_%H:%M)"
#
# echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
#
#
# ######## TOPICS
# out_dir="/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/"
# ml add Miniforge3
# conda activate scenicplus2
# mkdir -p "${out_dir}/qc"
# pycistopic tss get_tss \
#     --output /gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/qc/tss.bed \
#     --name "hsapiens_gene_ensembl" \
#     --to-chrom-source ucsc \
#     --ucsc hg38
#
# sbatch_output=$(sbatch --dependency=afterok:$job_id topics.sh)
# job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')
#
# specification="20,256,/gstore/project/epigen/benchmark/PBMC/scenicplus,PBMC,scenicplus,topics,Rosalind,$(date +%Y-%m-%d_%H:%M)"
#
# echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
#
#
# ######## MAIN ANALYSIS
ml Miniforge3
conda activate scenicplus2
mkdir -p scplus_pipeline
scenicplus init_snakemake --out_dir scplus_pipeline
cp config.yaml ./scplus_pipeline/Snakemake/config/config.yaml
#
# sbatch_output=$(sbatch --dependency=afterok:$job_id main.sh)

sbatch_output=$(sbatch main.sh)
job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')

specification="20,800,/gstore/project/epigen/benchmark/PBMC/scenicplus,PBMC,scenicplus,main,Rosalind,$(date +%Y-%m-%d_%H:%M)"

echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt

# python scplus_results_processing.py





