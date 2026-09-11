#!/bin/bash
export USER_ID="n026t"
source /home/${USER_ID}/.bashrc

# Explicitly initialize and activate conda in this script
eval "$(conda shell.bash hook)"
conda activate archvelo

export OUTPUT_LOCATION=/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_LOCATION}

# Create timestamped log filename
LOG_FILE="${OUTPUT_LOCATION}/archvelo_$(date +%Y%m%d_%H%M%S).log"

# Log the job start
{
  echo "=========================================="
  echo "Job started at $(date)"
  echo "Command: archvelo run"
  echo "=========================================="
  
  python /omics/groups/OE0132/tandem/nschmidt/DynaVelo_testing/velocity_experiments/archvelo_server.py
  
  EXIT_STATUS=$?
  
  echo "=========================================="
  echo "Job completed at $(date)"
  echo "Exit status: $EXIT_STATUS"
  echo "=========================================="
  exit $EXIT_STATUS
} | tee ${LOG_FILE}