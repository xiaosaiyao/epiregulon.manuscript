#!/bin/bash
if [ -d "scplus_pipeline" ]; then
  rm -fr scplus_pipeline
  ml add Miniforge3
  conda activate scenicplus2
  mkdir -p scplus_pipeline
  scenicplus init_snakemake --out_dir scplus_pipeline
  cp config.yaml ./scplus_pipeline/Snakemake/config/config.yaml
fi

sbatch_output=$(sbatch main.sh)
job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')
file_path="/gstore/project/epigen/benchmark/PBMC/scenicplus/slurm-${job_id}.out"
export file_path=$file_path
export job_id=$job_id
sbatch --dependency=afterok:$job_id --export=ALL /gstore/project/epigen/benchmark/PBMC/save_stats.sh
specification="20,64,/gstore/project/epigen/benchmark/PBMC/scenicplus/,PBMC,scenicplus,main,Rosalind,$(date +%Y-%m-%d_%H:%M)"
echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt




