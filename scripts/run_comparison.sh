#!/bin/bash
# SLURM job running scripts/mt_sample_comparison.py on the pipeline results.
# Copy to the analysis folder, edit the settings below and submit with
#   sbatch run_comparison.sh
#SBATCH -A <account>
#SBATCH -p <partition>
#SBATCH -J mt_comparison
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G
#SBATCH -t 12:00:00
#SBATCH -o mt_comparison_%j.log

set -euo pipefail

# ---- settings ----
MTOOLBOX_DIR=/path/to/MToolBox_snakemake
ANALYSIS_DIR=/path/to/analysis/dir
CONDA=/path/to/conda/bin/conda
REF_MT=rCRS                  # ref_genome_mt in data/reference_genomes.tab
REF_N=GRCh38                 # ref_genome_n
MT_FASTA=data/genomes/rCRS.fasta
NUMT_BED=numts.bed           # leave empty to skip the NUMT allele list
UNMASKED_GENOME=             # UNMASKED nuclear genome fasta, needed with NUMT_BED
SAMPLES=                     # tsv with columns sample, individual, type; empty:
                             # split sample names at the last "_" (P01_tumour)
EXCLUDE=""                   # samples to leave out, e.g. "P01_cfDNA P02_normal"
OUT=comparison/mt_comparison
# ------------------

eval "$("$CONDA" shell.bash hook)"
conda activate mtoolbox
cd "$ANALYSIS_DIR"
COMPARE="python $MTOOLBOX_DIR/scripts/mt_sample_comparison.py"

NUMT_ARGS=()
if [[ -n "$NUMT_BED" ]]; then
    if [[ ! -s numt_alleles.tsv ]]; then
        if [[ -z "$UNMASKED_GENOME" ]]; then
            echo "NUMT_BED is set but UNMASKED_GENOME is not" >&2
            exit 1
        fi
        echo "== Listing NUMT alleles ($(date))"
        $COMPARE numt-alleles --bed "$NUMT_BED" --genome "$UNMASKED_GENOME" \
            --mt-fasta "$MT_FASTA" --gmap-db "$REF_MT" \
            --threads "${SLURM_CPUS_PER_TASK:-4}" --out numt_alleles.tsv
    else
        echo "== Using the existing numt_alleles.tsv"
    fi
    NUMT_ARGS=(--numt-alleles numt_alleles.tsv)
fi

SAMPLE_ARGS=()
if [[ -n "$SAMPLES" ]]; then
    SAMPLE_ARGS=(--samples "$SAMPLES")
fi

EXCLUDE_ARGS=()
if [[ -n "$EXCLUDE" ]]; then
    EXCLUDE_ARGS=(--exclude $EXCLUDE)
fi

echo "== Comparing samples ($(date))"
mkdir -p "$(dirname "$OUT")"
$COMPARE compare --vcf "results/vcf/${REF_MT}_${REF_N}.annotated.vcf" \
    "${NUMT_ARGS[@]+"${NUMT_ARGS[@]}"}" "${SAMPLE_ARGS[@]+"${SAMPLE_ARGS[@]}"}" \
    "${EXCLUDE_ARGS[@]+"${EXCLUDE_ARGS[@]}"}" \
    --out "$OUT"

echo "== Done ($(date)). Sample concordance:"
column -t "${OUT}_concordance.tsv" 2>/dev/null || cat "${OUT}_concordance.tsv"
