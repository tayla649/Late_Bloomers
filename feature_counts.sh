#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --output=featurecounts.out
#SBATCH --error=featurecounts.err
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4

set -euo pipefail

# Move to counts directory
cd /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/counts/

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/health_sciences/bms/biochemistry/brownfield_lab/conda_envs/rna_seq

# Find only valid paired BAMs
paired_bams=()
for bam in ../mapping/*_*_*.bam; do
    filename=$(basename "$bam")
    # Exclude junk files (.err.bam, .out.bam, .sh.bam, .html.bam)
    if [[ ! "$filename" =~ \.(err|out|sh|html)\.bam$ ]]; then
        paired_bams+=("$bam")
    fi
done

# Check if we found any BAMs
if [ ${#paired_bams[@]} -eq 0 ]; then
    echo "No valid paired BAM files found! Exiting."
    exit 1
fi

echo "Counting features for ${#paired_bams[@]} BAM files..."

# Run featureCounts
featureCounts -T ${SLURM_CPUS_PER_TASK} \
  -p \
  -t exon \
  -g gene_id \
  -a ../genome/2.6.1_genome.gtf \
  -o ryegrass_exon_counts.txt \
  "${paired_bams[@]}"

echo "Feature counting completed!"
