#!/bin/bash
out_dir="/gstore/project/epigen/benchmark/PBMC/OUTPUT/"

ml CEDAR/2020.08
ml Python/scenic/2024_06_27

mkdir -p "${out_dir}/qc"
pycistopic tss get_tss \
    --output /gstore/project/epigen/benchmark/AR/OUTPUT/LNCaP/qc/tss.bed \
    --name "hsapiens_gene_ensembl" \
    --to-chrom-source ucsc \
    --ucsc hg38


# Start Time
start_time=$(date +%s.%N)
python cistopic.py $2
if [ $? -ne 0 ]; then
    exit 1
    echo "pycistopic error"
fi

# End Time
end_time=$(date +%s.%N)

timestamp=$(date +%Y-%m-%d_%H:%M)

# Calculate execution time
elapsed_time=$(echo "$end_time - $start_time" | bc)

output_file="/gstore/project/epigen/runtime.txt"
output="0,0,$elapsed_time,$1,$2,0,0,scenicplus,pycistopic_LNCaP,$timestamp,Rosalind"
echo "$output" >> "$output_file"
