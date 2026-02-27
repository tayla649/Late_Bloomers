#!/bin/bash
#SBATCH --job-name=readalign
#SBATCH --output=readalign.out
#SBATCH --error=readalign.err
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=4

#current seff 3621883

set -euo pipefail

cd /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/

module load hisat2

mkdir -p mapping

for r1 in raw/*_1.fq.gz
do
base=$(basename "$r1" _1.fq.gz)
r2="raw/${base}_2.fq.gz"

if [[ ! -f "$r2" ]]; then
echo "Missing pair for $base - skipping"
continue
fi

echo "Aligning $base"

hisat2 -p ${SLURM_CPUS_PER_TASK} \
-x genome/2.6.1_genome \
-1 "$r1" \
-2 "$r2" \
--summary-file mapping/${base}_summary.txt \
| samtools view -bS - > mapping/${base}.bam

echo "[$(date)] Finished $base"
done