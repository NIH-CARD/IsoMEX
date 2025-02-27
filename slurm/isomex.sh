#!/bin/bash
#SBATCH --job-name=IsoMEX
#SBATCH -o run_logs/slurm.%J.%x.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --array=1-$(wc -l < samples.txt)  # Dynamically set array size

#submit via
# filtered_feature_bc_matrix

# Load the sample list into an array
SAMPLES=($(cat samples.txt))

# Determine the sample for this specific array task
SAMPLE_NAME=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}

# Define paths
GENE_MAP="gene_map.txt"
TRANSCRIPT_MAP="transcript_map.txt"
OUTDIR="isomex/${SAMPLE_NAME}"

# Print the sample name
echo "Processing sample: $SAMPLE_NAME"

# Ensure output directory exists
mkdir -p ${OUTDIR}
cd ${OUTDIR}

# Execute the IsoMEX script for the selected sample.
python3 /path/to/isomex.py "${SAMPLE_NAME}" \
    --gene_map ${GENE_MAP} \
    --transcript_map ${TRANSCRIPT_MAP} \
    --output_dir filtered_feature_bc_matrix

echo "Processing of sample $SAMPLE_NAME complete."
