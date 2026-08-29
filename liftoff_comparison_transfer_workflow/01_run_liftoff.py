#!/usr/bin/env python3
"""Run all requested donor-to-target Liftoff mappings.

The manifest is tab-delimited and must contain:
    genome_id, fasta, gff3, target, donor

Rows with target=1 are final assemblies to annotate. Rows with donor=1 are
annotation donors. A donor-only alternative assembly therefore uses target=0
and donor=1. Self comparisons are skipped.
"""

import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path


TRUE = {"1", "true", "yes", "y"}


def enabled(value):
    return value.strip().lower() in TRUE


def read_manifest(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "fasta", "gff3", "target", "donor"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: " + ", ".join(sorted(missing)))
        rows = list(reader)

    seen = set()
    base = path.parent
    for row in rows:
        genome_id = row["genome_id"].strip()
        if not genome_id or genome_id in seen:
            raise ValueError(f"Missing or duplicate genome_id: {genome_id!r}")
        seen.add(genome_id)
        row["genome_id"] = genome_id
        for field in ("fasta", "gff3"):
            item = Path(row[field])
            row[field] = item if item.is_absolute() else (base / item).resolve()
            if not row[field].is_file() or row[field].stat().st_size == 0:
                raise ValueError(f"Missing or empty {field} for {genome_id}: {row[field]}")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--liftoff", default="liftoff")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    rows = read_manifest(args.manifest.resolve())
    targets = sorted((r for r in rows if enabled(r["target"])), key=lambda r: r["genome_id"])
    donors = sorted((r for r in rows if enabled(r["donor"])), key=lambda r: r["genome_id"])
    if not targets or not donors:
        raise ValueError("Manifest must define at least one target and one donor")

    gff_dir = args.output_dir / "raw_gff"
    unmapped_dir = args.output_dir / "unmapped"
    log_dir = args.output_dir / "logs"
    for directory in (gff_dir, unmapped_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    failures = []
    completed = skipped = 0
    for target in targets:
        for donor in donors:
            if target["genome_id"] == donor["genome_id"]:
                continue
            stem = f"{target['genome_id']}_vs_{donor['genome_id']}"
            output = gff_dir / f"{stem}.gff"
            unmapped = unmapped_dir / f"{stem}.unmapped.txt"
            log = log_dir / f"{stem}.log"
            if output.is_file() and output.stat().st_size > 0 and not args.overwrite:
                print(f"SKIP\t{stem}\texisting output")
                skipped += 1
                continue

            command = [
                args.liftoff,
                str(target["fasta"]),
                str(donor["fasta"]),
                "-g", str(donor["gff3"]),
                "-o", str(output),
                "-u", str(unmapped),
                "-p", str(args.threads),
            ]
            print("RUN\t" + stem + "\t" + shlex.join(command))
            with log.open("w") as handle:
                result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
            if result.returncode or not output.is_file() or output.stat().st_size == 0:
                failures.append((stem, result.returncode, str(log)))
                print(f"FAIL\t{stem}\tstatus={result.returncode}\tlog={log}", file=sys.stderr)
            else:
                completed += 1

    print(f"completed={completed}\tskipped={skipped}\tfailed={len(failures)}")
    if failures:
        failure_file = args.output_dir / "01_liftoff_failures.tsv"
        with failure_file.open("w") as out:
            out.write("comparison\tstatus\tlog\n")
            for row in failures:
                out.write("\t".join(map(str, row)) + "\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
