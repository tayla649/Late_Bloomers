#!/bin/bash
#SBATCH --job-name=gffread
#SBATCH --output=gffread.out
#SBATCH --error=gffread.err
#SBATCH --time=00:01:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB

cd /projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/genome/

gffread 2.6.1_genome.gff3 -T -o 2.6.1_genome.gtf