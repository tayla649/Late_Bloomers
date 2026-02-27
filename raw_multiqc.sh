#!/bin/bash
#SBATCH --job-name=lachlan_multiqc
#SBATCH --output=multiqc.out
#SBATCH --error=multiqc.err
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

set -euo pipefail
QC_DIR="/projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/qc"
OUT_DIR="/projects/health_sciences/bms/biochemistry/brownfield_lab/lachlan_rna_seq/multiqc"

mkdir -p "$OUT_DIR"
MULTIQC="$HOME/.local/bin/multiqc"

echo "Running MultiQC"
echo "QC_DIR = $QC_DIR"
echo "OUT_DIR = $OUT_DIR"
echo "MULTIQC = $MULTIQC"

"$MULTIQC" --version
"$MULTIQC" "QC_DIR" -o "$OUT_DIR" --filename multiqc_report.html

