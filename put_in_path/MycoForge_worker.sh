#!/bin/bash -l
set -euo pipefail

# # # Args # # #
INPUT_FILE="$1"
LINE_NUM="$2"

if [[ -z "${INPUT_FILE:-}" || -z "${LINE_NUM:-}" ]]; then
  echo "Usage: MycoAssembler_worker.sh <metadata.tsv> <line_number>"
  exit 1
fi

# # # Get headers/make files # # #
INPUT_CLEAN=$(sed 's/\r$//' "$INPUT_FILE")

LR_DIR=$(sed -n 's/^#Location of Long Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
SR_DIR=$(sed -n 's/^#Location of Short Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
RNA_DIR=$(sed -n 's/^#Location of RNAseq Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)

LOGDIR="$OUTPUT_DIR/01_logs"
WORKDIR="$OUTPUT_DIR/02_processed_small_reads"
ASM_DIR="$OUTPUT_DIR/03_assembled_genomes"
ANNO_DIR="$OUTPUT_DIR/04_funannotate_genomes"
BUSCO_DIR="$OUTPUT_DIR/05_busco_output"

mkdir -vp "$LOGDIR" "$WORKDIR" "$ASM_DIR" "$ANNO_DIR" "$BUSCO_DIR"

# # # Set threads/memory # # #
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  MAX_MEMORY=$((SLURM_MEM_PER_NODE / 1024))
else
  MAX_MEMORY=128
fi

# # # Extract metadata line # # #
LINE=$(grep -vE '^\s*#|^\s*$' <<< "$INPUT_CLEAN" | sed -n "${LINE_NUM}p")

if [[ -z "$LINE" ]]; then
  echo "Error: No metadata found for line $LINE_NUM"
  exit 1
fi

IFS=$'\t' read -r \
  LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER GENOME_SIZE BUSCO_TYPE BUSCO_LIN \
  <<< "$LINE"

BASE="${DES//[^a-zA-Z0-9._-]/_}"

echo
echo "============================================"
echo " Processing $BASE (line $LINE_NUM)"
echo "============================================"
echo

LOGFILE="$LOGDIR/worker.${BASE}.line${LINE_NUM}.log"
exec > >(tee -a "$LOGFILE") 2>&1

# # # Load functions # # #
source MycoForge_helper_functions.sh

# # # Run pipeline # # #
echo "Step 1: Preprocessing short reads"
module load AAFTF
preprocess_reads "$BASE" "$R1" "$R2" "$RNA1" "$RNA2"

echo
echo "Step 2: Assembly"
module unload java || true
module load canu Flye
assemble_genome "$BASE" "$ASSEMBLER" "$LR" "$R1" "$R2" "$GENOME_SIZE"

echo
echo "Step 3: Polishing & assessment"
module load AAFTF samtools bcftools medaka
polish_and_assess "$BASE" "$ASSEMBLER" "$LR" "$R1" "$R2"

echo
echo "Step 4: Annotation"
module load funannotate
export PASACONF="$PASAHOME/pasa_conf/alignAssembly.sqlite.config" # this will be removed before update
export PATH="${CODINGQUARRY_PATH}:$PATH"
annotate_genome "$BASE" "$ASSEMBLER" "$RNA1" "$RNA2" "$SPECIES"

echo
echo "Step 5: BUSCO"
module load busco
run_busco "$BASE" "$ASSEMBLER" "$SPECIES" "$BUSCO_LIN" "$BUSCO_TYPE"

echo
echo "============================================"
echo " Finished $BASE"
echo "============================================"
