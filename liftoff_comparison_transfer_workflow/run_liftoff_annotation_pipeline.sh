#!/usr/bin/env bash
# Run the four-stage Hyalorbilia Liftoff annotation-transfer workflow.
#
# Usage:
#   bash run_liftoff_annotation_pipeline.sh genomes.tsv output_directory [threads]
#
# Manifest columns (tab-delimited):
#   genome_id  fasta  gff3  target  donor
#
# The 20 final genomes should have target=1 and donor=1. Alternative donor-only
# assemblies such as HsImV25B and HsImV27B should have target=0 and donor=1.
# Relative FASTA/GFF3 paths are resolved relative to the manifest directory.
#
# This script deliberately does not use `set -e` and does not return a nonzero
# status by default, making it safe to run from an interactive allocation. It
# stops launching downstream stages after a failure and prints the failed stage.
# Set PIPELINE_STRICT=1 to return a nonzero status for batch jobs or CI.

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
MANIFEST=${1:-}
OUTPUT_ROOT=${2:-}
THREADS=${3:-8}
PYTHON=${PYTHON:-python3}
LIFTOFF=${LIFTOFF:-liftoff}
GFFREAD=${GFFREAD:-gffread}

if [ -z "$MANIFEST" ] || [ -z "$OUTPUT_ROOT" ]; then
    printf 'Usage: bash %s genomes.tsv output_directory [threads]\n' "$0" >&2
    if [ "${PIPELINE_STRICT:-0}" = "1" ]; then
        exit 2
    fi
    return 0 2>/dev/null || exit 0
fi

mkdir -p "$OUTPUT_ROOT"

run_stage() {
    STAGE_NAME=$1
    shift

    printf '\n===== %s =====\n' "$STAGE_NAME"
    "$@"
    STAGE_STATUS=$?

    if [ "$STAGE_STATUS" -ne 0 ]; then
        printf 'PIPELINE_STAGE_FAILED\t%s\tstatus=%s\n' \
            "$STAGE_NAME" "$STAGE_STATUS" >&2
        FAILED_STAGE=$STAGE_NAME
        return "$STAGE_STATUS"
    fi

    printf 'PIPELINE_STAGE_PASSED\t%s\n' "$STAGE_NAME"
    return 0
}

FAILED_STAGE=""

run_stage "01_run_liftoff" \
    "$PYTHON" "$SCRIPT_DIR/01_run_liftoff.py" \
    --manifest "$MANIFEST" \
    --output-dir "$OUTPUT_ROOT/01_liftoff" \
    --threads "$THREADS" \
    --liftoff "$LIFTOFF"
STATUS=$?

if [ "$STATUS" -eq 0 ]; then
    run_stage "02_filter_liftoff_models" \
        "$PYTHON" "$SCRIPT_DIR/02_filter_liftoff_models.py" \
        --input-dir "$OUTPUT_ROOT/01_liftoff/raw_gff" \
        --output-dir "$OUTPUT_ROOT/02_threshold_filtered" \
        --coverage-min 0.70 \
        --identity-min 0.95
    STATUS=$?
fi

if [ "$STATUS" -eq 0 ]; then
    run_stage "03_validate_complete_cds" \
        "$PYTHON" "$SCRIPT_DIR/03_validate_complete_cds.py" \
        --manifest "$MANIFEST" \
        --input-dir "$OUTPUT_ROOT/02_threshold_filtered" \
        --output-dir "$OUTPUT_ROOT/03_cds_validated" \
        --gffread "$GFFREAD"
    STATUS=$?
fi

if [ "$STATUS" -eq 0 ]; then
    run_stage "04_merge_liftoff_models" \
        "$PYTHON" "$SCRIPT_DIR/04_merge_liftoff_models.py" \
        --manifest "$MANIFEST" \
        --input-dir "$OUTPUT_ROOT/03_cds_validated" \
        --output-dir "$OUTPUT_ROOT/04_merged" \
        --novel-overlap-fraction 0.50
    STATUS=$?
fi

if [ "$STATUS" -eq 0 ]; then
    printf '\nPIPELINE_RESULT\tPASS\n'
    printf 'Merged annotations: %s\n' "$OUTPUT_ROOT/04_merged"
else
    printf '\nPIPELINE_RESULT\tFAIL\tstage=%s\tstatus=%s\n' \
        "$FAILED_STAGE" "$STATUS" >&2
fi

if [ "${PIPELINE_STRICT:-0}" = "1" ]; then
    exit "$STATUS"
fi

# Default interactive behavior: report failure but do not terminate the parent
# shell or an interactive Slurm allocation with a nonzero status.
return 0 2>/dev/null || exit 0
