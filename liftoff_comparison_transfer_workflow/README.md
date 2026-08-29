# Hyalorbilia Liftoff annotation-transfer pipeline

This workflow reproduces the corrected gene-model transfer procedure used for
the Hyalorbilia genome project. Twenty final assemblies are targets and donors;
the alternative HsImV25B and HsImV27B assemblies are donors only.

## Selection procedure

1. Run each donor annotation against every final target with Liftoff 1.6.3.
2. Retain mappings with `coverage >= 0.70` and `sequence_ID >= 0.95`.
3. Run `gffread -J -F --keep-genes` and remove gene records left without a
   retained transcript.
4. Require `valid_ORFs=1` and reject lifted genes overlapping a native gene on
   the same strand.
5. Rank mappings by `coverage * sequence_ID`. Retain the best mapping of each
   source gene, then greedily resolve lifted-model conflicts in descending
   score order. A candidate is rejected when at least 50% of its span overlaps
   an already accepted lifted model on the same strand.

Opposite-strand overlaps are intentionally retained. Later submission-specific
conflict resolution is not part of this pipeline.

## Requirements

- Python 3.9 or later
- Liftoff 1.6.3
- gffread 0.12.7

The Python scripts use only the standard library. A Conda environment is
provided in `environment.yml`.

## Manifest

Create a tab-delimited manifest with these columns:

```text
genome_id  fasta  gff3  target  donor
```

- `target=1`: produce a merged annotation for this assembly.
- `donor=1`: use this assembly as an annotation donor.
- Relative paths are resolved relative to the manifest file.
- Donor-only assemblies use `target=0` and `donor=1`.

See `genomes.example.tsv` for the expected genome rows.

## Run

```bash
conda env create -f environment.yml
conda activate hyalorbilia-liftoff

bash run_liftoff_annotation_pipeline.sh \
    genomes.tsv \
    liftoff_results \
    8
```

The wrapper does not use `set -e`. It reports a failed stage and, by default,
returns status zero to avoid closing an interactive Slurm allocation. For batch
jobs and continuous integration, request conventional nonzero failure status:

```bash
PIPELINE_STRICT=1 bash run_liftoff_annotation_pipeline.sh \
    genomes.tsv liftoff_results 8
```

Existing nonempty raw Liftoff results are skipped. Add `--overwrite` when
running `01_run_liftoff.py` directly if they must be regenerated.

## Outputs

```text
liftoff_results/
├── 01_liftoff/
│   ├── raw_gff/
│   ├── unmapped/
│   └── logs/
├── 02_threshold_filtered/
│   └── 02_filter_audit.tsv
├── 03_cds_validated/
│   └── 03_cds_validation_audit.tsv
└── 04_merged/
    ├── GENOME.merged.gff3
    └── 04_merge_audit.tsv
```

Every filtering and merge decision is recorded in an audit TSV. The merged
GFF3 files retain donor product names and other attributes present in the
source annotations and add `liftoff_*` provenance attributes.

## Reproducibility note

Tie-breaking is deterministic: score, sequence identity, coverage, donor ID,
and gene ID are considered in that order. This removes the filesystem-order
dependence present in the historical exploratory scripts.
