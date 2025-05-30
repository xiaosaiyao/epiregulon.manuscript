#!/bin/bash

if [ -d "create_cisTarget_databases" ]; then
    echo "repository already cloned"
else
    git clone https://github.com/aertslab/create_cisTarget_databases
fi

# wget https://resources.aertslab.org/cistarget/programs/cbust

if [ ! -d "aertslab_motif_colleciton" ]; then
    mkdir -p aertslab_motif_colleciton
    wget -O aertslab_motif_colleciton/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip
    cd aertslab_motif_colleciton; unzip -q v10nr_clust_public.zip
    cd ..
fi


ml add BEDTools/2.30.0-GCCcore-6.3.0

# genome fasta downloadedd from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
# chromosome sizes downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
REGION_BED="/gstore/project/epigen/benchmark/AR/OUTPUT/peakRegions_VCaP.bed"
GENOME_FASTA="/gstore/project/epigen/PBMC/analysis/scenicplus/hg38.fa"
CHROMSIZES="/gstore/project/epigen/PBMC/analysis/scenicplus/hg38.chrom.sizes"
DATABASE_PREFIX="AR_VCaP"
SCRIPT_DIR=$(pwd)
SCRIPT_DIR="${SCRIPT_DIR}/create_cisTarget_databases"

# Start Time
start_time=$(date +%s.%N)

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.AR_VCaP_1kb_bg_padding.fa \
        1000 \
        yes

if [ $? -ne 0 ]; then
    echo "Error create_fasta_with_padded_bg_from_bed"
    exit 1
fi
#!/bin/bash
ml add Miniforge3
conda activate scenicplus

ls aertslab_motif_colleciton/v10nr_clust_public/singletons > motifs.txt

OUT_DIR=""${PWD}""
CBDIR="${OUT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="${OUT_DIR}/hg38.AR_VCaP_1kb_bg_padding.fa"
MOTIF_LIST="${OUT_DIR}/motifs.txt"
DATABASE_PREFIX="AR_VCaP"
SCRIPT_DIR="/gstore/project/epigen/PBMC/analysis/scenicplus/create_cisTarget_databases"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 500 \
    -t 26 \
    -c /gstore/project/epigen/PBMC/analysis/scenicplus/cbust

if [ $? -ne 0 ]; then
    echo "Error cbust"
    exit 1
fi
# End Time
end_time=$(date +%s.%N)

timestamp=$(date +%Y-%m-%d_%H:%M)

# Calculate execution time
elapsed_time=$(echo "$end_time - $start_time" | bc)

output_file="/gstore/project/epigen/runtime.txt"
output="0,0,$elapsed_time,0,0,$1,$2,scenicplus,cistarget_VCaP,$timestamp,Rosalind"
echo "$output" >> "$output_file"
