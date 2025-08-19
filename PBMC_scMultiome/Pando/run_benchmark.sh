#!/bin/sh
sbatch_output=$(sbatch preprocessing.sh)
job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')

specification="1,64,/gstore/project/epigen/benchmark/PBMC/Pando,PBMC,Pando,preprocessing,Rosalind,$(date +%Y-%m-%d_%H:%M)"
if [ -f "/gstore/project/epigen/benchmark/resource_use.txt" ]; then
    echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
else
    colnames="job_id,n_cpu,requested_memtory,job_id,path,dataset,tool,part,enviroment,timestamp"
    echo -e "${colnames}\n${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
fi

file_path="/gstore/project/epigen/benchmark/PBMC/Pando/slurm-${job_id}.out"
export file_path=$file_path
export job_id=$job_id
sbatch --dependency=afterok:$job_id --export=ALL /gstore/project/epigen/benchmark/PBMC/save_stats.sh

sbatch_output=$(sbatch --dependency=afterok:$job_id main.sh)
job_id=$(echo "$sbatch_output" | grep -o 'Submitted batch job [0-9]*' | awk '{print $NF}')
file_path="/gstore/project/epigen/benchmark/PBMC/Pando/slurm-${job_id}.out"
export file_path=$file_path
export job_id=$job_id


specification="20,64,/gstore/project/epigen/benchmark/PBMC/Pando,PBMC,Pando,main,Rosalind,$(date +%Y-%m-%d_%H:%M)"

echo "${job_id},${specification}" >> /gstore/project/epigen/benchmark/resource_use.txt
sbatch --dependency=afterok:$job_id --export=ALL /gstore/project/epigen/benchmark/PBMC/save_stats.sh
