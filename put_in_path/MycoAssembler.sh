#!/bin/bash -l
set -euo pipefail

# # # PARSE ARGS # # #
INPUT_FILE=""
SLURM_OPTS_FILE=""
CODINGQUARRY_PATH=""

show_help() {
  echo "Usage: $0 --input metadata.tsv [--slurm-opts slurm.opts] [--codingquarry /path/to/CodingQuarry]"
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT_FILE="$2"; shift 2 ;;
    --codingquarry) CODINGQUARRY_PATH="$2"; shift 2 ;;
    --slurm-opts) SLURM_OPTS_FILE="$2"; shift 2 ;;
    -h|--help) show_help ;;
    *) echo "Unknown option $1"; show_help ;;
  esac
done

# Check required args
if [[ -z "$INPUT_FILE" ]]; then
  echo "Error: --input <metadata.tsv> is required."
  show_help
fi

if [[ -z "$CODINGQUARRY_PATH" ]]; then
  echo "Error: --codingquarry is required"
  exit 1
fi

if [[ -n "$SLURM_OPTS_FILE" && ! -f "$SLURM_OPTS_FILE" ]]; then
  echo "Error: SLURM opts file not found: $SLURM_OPTS_FILE"
  exit 1
fi

if [[ -n "$CODINGQUARRY_PATH" ]]; then
  if [[ ! -d "$CODINGQUARRY_PATH" ]]; then
    echo "Error: CodingQuarry path not found: $CODINGQUARRY_PATH"
    exit 1
  fi
  export CODINGQUARRY_PATH
fi

# (CRLF -> LF)
INPUT_CLEAN=$(sed 's/\r$//' "$INPUT_FILE")

# # # EXTRACT HEADERS FOR PATHS AND SET UP WORKING DIRECTORY # # #
LR_DIR=$(sed -n 's/^#Location of Long Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
SR_DIR=$(sed -n 's/^#Location of Short Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
RNA_DIR=$(sed -n 's/^#Location of RNAseq Reads=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)
OUTPUT_DIR=$(sed -n 's/^#Desired Output Location=//p' <<< "$INPUT_CLEAN" | head -n1 | xargs)

for header in \
  "#Location of Long Reads=" \
  "#Location of Short Reads=" \
  "#Location of RNAseq Reads=" \
  "#Desired Output Location="; do
    count=$(grep -c "^${header}" <<< "$INPUT_CLEAN")
    [[ "$count" -ne 1 ]] && {
      echo "Error: Expected exactly one '${header}' header, found $count"
      exit 1
    }
done

# # # CHECK DES COLLISIONS # # #
check_des_collisions() {
  declare -A MAP
  while IFS=$'\t' read -r LR R1 R2 RNA1 RNA2 DES SPECIES ASSEMBLER GENOME_SIZE BUSCO_TYPE BUSCO_LIN; do
    [[ "$LR" =~ ^# || -z "$LR" ]] && continue
    BASE="${DES//[^a-zA-Z0-9._-]/_}"
    [[ -n "${MAP[$BASE]:-}" ]] && {
      echo "Duplicate designations: $DES vs ${MAP[$BASE]}"
      exit 1
    }
    MAP[$BASE]="$DES"
  done < <(printf "%s\n" "$INPUT_CLEAN")
}
check_des_collisions

# # # FINAL VALIDATIONS # # #
for f in MycoAssembler_helper_functions.sh MycoAssembler_eggnog.py MycoAssembler_worker.sh; do
  printf "%-30s : " "$f"
  if command -v "$f" >/dev/null 2>&1 && [[ -x "$(command -v "$f")" ]]; then
    echo "OK ($(command -v "$f"))"
  else
    echo "MISSING or NOT EXECUTABLE"
  fi
done

for v in LR_DIR SR_DIR RNA_DIR OUTPUT_DIR; do
  if [[ -z "${!v}" ]]; then
    echo "Error: Required header $v missing from metadata file"
    exit 1
  fi
done

[[ "$LR_DIR" != "NA" && ! -d "$LR_DIR" ]] && { echo "Missing LR_DIR: $LR_DIR"; exit 1; }
[[ "$SR_DIR" != "NA" && ! -d "$SR_DIR" ]] && { echo "Missing SR_DIR: $SR_DIR"; exit 1; }
[[ "$RNA_DIR" != "NA" && ! -d "$RNA_DIR" ]] && { echo "Missing RNA_DIR: $RNA_DIR"; exit 1; }

source MycoAssembler_helper_functions.sh
validate_metadata

# # # BUILD SBATCH OPTIONS # # #
build_sbatch_args() {
  local file="$1"
  local args=()

  while IFS='=' read -r key value; do
    [[ -z "$key" || "$key" =~ ^# ]] && continue
    args+=( "--$key=$value" )
  done < "$file"

  echo "${args[@]}"
}

SBATCH_ARGS=()

if [[ -n "$SLURM_OPTS_FILE" ]]; then
  SBATCH_ARGS+=( $(build_sbatch_args "$SLURM_OPTS_FILE") )
else
  SBATCH_ARGS+=(
    --cpus-per-task=32
    --mem=128G
    --time=14-00:00:00
  )
fi

# # # PROCESS # # #
DATA_LINES=$(grep -vE '^\s*#|^\s*$' <<< "$INPUT_CLEAN")
LINE_NUM=0

while IFS= read -r line; do
  LINE_NUM=$((LINE_NUM+1))

  sbatch \
    --job-name=fass_${LINE_NUM} \
    "${SBATCH_ARGS[@]}" \
    MycoAssembler_worker.sh "$INPUT_FILE" "$LINE_NUM"

done <<< "$DATA_LINES"
