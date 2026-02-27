#!/bin/bash
#SBATCH --job-name=mapping_multiqc
#SBATCH --output=mapping_multiqc.out
#SBATCH --error=mapping_multiqc.err
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

module load samtools
module load multiqc

BAM_DIR=/projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/mapping
LOG_DIR=/projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/scripts

mkdir -p /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/qc/samtools /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/qc/hisat2 multiqc

for bam in ${BAM_DIR}/*.bam; do
base=$(basename "$bam" .bam)
samtools flagstat "$bam" > qc/samtools/${base}.flagstat.txt
done

for log in ${LOG_DIR}/*.err; do
awk '
/reads; of these:/ {p=1}
p {print}
/overall alignment rate/ {p=0}
' "$log" > qc/hisat2/$(basename "$log").hisat2.txt
done

multiqc qc -o multiqc -f
