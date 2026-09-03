#!/bin/bash
export USER_ID="n026t"
source /home/${USER_ID}/.bashrc

# Explicitly initialize and activate conda in this script
eval "$(conda shell.bash hook)"
conda activate dynavelo

export READ_LOCATION=/omics/groups/OE0132/internal/metzj/projects/Multiome/CellRanger_Output/Exp4/ctr/outs
export OUTPUT_LOCATION=/omics/groups/OE0132/tandem/nschmidt/velocyto_output

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_LOCATION}

# Create timestamped log filename
LOG_FILE="${OUTPUT_LOCATION}/velocyto_$(date +%Y%m%d_%H%M%S).log"

# Log the job start
{
  echo "=========================================="
  echo "Job started at $(date)"
  echo "Command: velocyto run"
  echo "=========================================="
  
  velocyto run \
    -b ${READ_LOCATION}/filtered_feature_bc_matrix/barcodes.tsv.gz \
    -o ${OUTPUT_LOCATION} \
    /omics/groups/OE0132/internal/nschmidt/velocyto_data/cellsorted_gex_possorted_bam.bam \ 
    /home/${USER_ID}/mm10_genes.gtf
  
  EXIT_STATUS=$?
  
  echo "=========================================="
  echo "Job completed at $(date)"
  echo "Exit status: $EXIT_STATUS"
  echo "=========================================="
  exit $EXIT_STATUS
} | tee ${LOG_FILE}