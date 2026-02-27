#!/bin/bash
#SBATCH --job-name=genomeanno
#SBATCH --output=genomeanno.out
#SBATCH --error=genomeanno.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256GB

set -ueo pipefail

cd /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/genome/

module load hisat2

hisat2_extract_splice_sites.py 2.6.1_genome.gtf > 2.6.1_genome.ss
hisat2_extract_exons.py 2.6.1_genome.gtf > 2.6.1_genome.exon

hisat2-build -p 8 \
--ss 2.6.1_genome.ss \
--exon 2.6.1_genome.exon \
2.6.1_genome.fasta \
2.6.1_genome