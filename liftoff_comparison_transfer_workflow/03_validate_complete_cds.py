#!/usr/bin/env python3
"""Validate Liftoff CDS models with gffread and remove orphan genes.

For each filtered target-vs-donor GFF3, this stage runs:
    gffread -J -F --keep-genes -g TARGET_FASTA

The historical --keep-genes option leaves gene records whose transcripts were
discarded by -J. This script removes those orphan records and writes an audit.
"""

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


TRUE = {"1", "true", "yes", "y"}
NAME_RE = re.compile(r"^(.+)_vs_(.+)\.gff$")


def parse_attrs(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def target_fastas(manifest):
    result = {}
    base = manifest.parent
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "fasta", "target"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            if row["target"].strip().lower() not in TRUE:
                continue
            fasta = Path(row["fasta"])
            fasta = fasta if fasta.is_absolute() else (base / fasta).resolve()
            if not fasta.is_file() or fasta.stat().st_size == 0:
                raise ValueError(f"Missing target FASTA for {row['genome_id']}: {fasta}")
            result[row["genome_id"]] = fasta
    return result


def orphan_genes(path):
    genes = set()
    parented = set()
    with path.open() as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(columns) != 9:
                continue
            feature_attrs = parse_attrs(columns[8])
            if columns[2] == "gene" and feature_attrs.get("ID"):
                genes.add(feature_attrs["ID"])
            elif columns[2] in {"mRNA", "transcript"}:
                parented.update(x for x in feature_attrs.get("Parent", "").split(",") if x)
    return genes - parented, genes


def remove_orphans(source, destination, orphans):
    with source.open() as incoming, destination.open("w") as outgoing:
        for line in incoming:
            columns = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(columns) == 9 and columns[2] == "gene":
                if parse_attrs(columns[8]).get("ID") in orphans:
                    continue
            outgoing.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gffread", default="gffread")
    args = parser.parse_args()
    fastas = target_fastas(args.manifest.resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temporary = args.output_dir / "gffread_intermediate"
    temporary.mkdir(exist_ok=True)
    inputs = sorted(args.input_dir.glob("*.gff"))
    if not inputs:
        print(f"No *.gff files found in {args.input_dir}", file=sys.stderr)
        return 1

    failures = 0
    audit_path = args.output_dir / "03_cds_validation_audit.tsv"
    with audit_path.open("w") as audit:
        audit.write("file\ttarget\tgenes_after_gffread\torphan_genes_removed\tresult\n")
        for source in inputs:
            match = NAME_RE.match(source.name)
            if not match or match.group(1) not in fastas:
                print(f"SKIP\t{source.name}\tcannot identify target FASTA", file=sys.stderr)
                audit.write(f"{source.name}\t\t0\t0\tTARGET_LOOKUP_FAILED\n")
                failures += 1
                continue
            target = match.group(1)
            intermediate = temporary / f"{source.name}.filtered"
            destination = args.output_dir / f"{source.name}.filtered"
            command = [
                args.gffread, "-J", "-F", "--keep-genes",
                "-g", str(fastas[target]), "-o", str(intermediate), str(source),
            ]
            result = subprocess.run(command)
            if result.returncode or not intermediate.is_file():
                print(f"FAIL\t{source.name}\tgffread_status={result.returncode}", file=sys.stderr)
                audit.write(f"{source.name}\t{target}\t0\t0\tGFFREAD_FAILED\n")
                failures += 1
                continue
            orphans, genes = orphan_genes(intermediate)
            remove_orphans(intermediate, destination, orphans)
            audit.write(f"{source.name}\t{target}\t{len(genes)}\t{len(orphans)}\tPASS\n")
            print(f"{source.name}\tgenes={len(genes)}\torphans_removed={len(orphans)}\tPASS")
    print(f"files={len(inputs)}\tfailures={failures}\taudit={audit_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
