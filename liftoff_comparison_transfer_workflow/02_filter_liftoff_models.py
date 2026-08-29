#!/usr/bin/env python3
"""Filter Liftoff GFF3 files by nucleotide coverage and identity.

Defaults reproduce the corrected historical filter: coverage >= 0.70 and
sequence_ID >= 0.95. Complete feature hierarchies are retained using ID and
Parent relationships. Missing or malformed metrics cause rejection.
"""

import argparse
import csv
import sys
from pathlib import Path


def attrs(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def number(attributes, key):
    try:
        return float(attributes[key])
    except (KeyError, TypeError, ValueError):
        return None


def filter_file(source, destination, coverage_min, identity_min, audit):
    lines = source.read_text().splitlines(keepends=True)
    keep_genes = set()
    feature_parent = {}
    metrics = {}

    for line in lines:
        columns = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(columns) != 9:
            continue
        feature_attrs = attrs(columns[8])
        feature_id = feature_attrs.get("ID")
        parents = [x for x in feature_attrs.get("Parent", "").split(",") if x]
        if feature_id:
            feature_parent[feature_id] = parents
        if columns[2] == "gene" and feature_id:
            coverage = number(feature_attrs, "coverage")
            identity = number(feature_attrs, "sequence_ID")
            keep = (
                coverage is not None and identity is not None
                and coverage >= coverage_min and identity >= identity_min
            )
            metrics[feature_id] = (coverage, identity, keep)
            if keep:
                keep_genes.add(feature_id)

    def belongs_to_retained_gene(feature_attrs):
        pending = [x for x in feature_attrs.get("Parent", "").split(",") if x]
        visited = set()
        while pending:
            parent = pending.pop()
            if parent in keep_genes:
                return True
            if parent not in visited:
                visited.add(parent)
                pending.extend(feature_parent.get(parent, []))
        return False

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as out:
        for line in lines:
            columns = line.rstrip("\n").split("\t")
            if line.startswith("#") or not line.strip():
                out.write(line)
            elif len(columns) == 9:
                feature_attrs = attrs(columns[8])
                if columns[2] == "gene":
                    if feature_attrs.get("ID") in keep_genes:
                        out.write(line)
                elif belongs_to_retained_gene(feature_attrs):
                    out.write(line)

    for gene_id, (coverage, identity, keep) in sorted(metrics.items()):
        reasons = []
        if coverage is None:
            reasons.append("missing_coverage")
        elif coverage < coverage_min:
            reasons.append("coverage_below_threshold")
        if identity is None:
            reasons.append("missing_sequence_ID")
        elif identity < identity_min:
            reasons.append("identity_below_threshold")
        audit.writerow({
            "file": source.name,
            "gene": gene_id,
            "coverage": "" if coverage is None else coverage,
            "sequence_ID": "" if identity is None else identity,
            "decision": "KEEP" if keep else "REMOVE",
            "reason": "thresholds_passed" if keep else ",".join(reasons),
        })
    return len(metrics), len(keep_genes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--coverage-min", type=float, default=0.70)
    parser.add_argument("--identity-min", type=float, default=0.95)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    inputs = sorted(args.input_dir.glob("*.gff"))
    if not inputs:
        print(f"No *.gff files found in {args.input_dir}", file=sys.stderr)
        return 1
    audit_path = args.output_dir / "02_filter_audit.tsv"
    with audit_path.open("w", newline="") as handle:
        fields = ["file", "gene", "coverage", "sequence_ID", "decision", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        total = kept = 0
        for source in inputs:
            destination = args.output_dir / source.name
            n_total, n_kept = filter_file(
                source, destination, args.coverage_min, args.identity_min, writer
            )
            total += n_total
            kept += n_kept
            print(f"{source.name}\ttotal={n_total}\tkept={n_kept}\tremoved={n_total-n_kept}")
    print(f"total={total}\tkept={kept}\tremoved={total-kept}\taudit={audit_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
