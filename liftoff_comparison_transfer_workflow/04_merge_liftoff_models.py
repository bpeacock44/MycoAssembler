#!/usr/bin/env python3
"""Merge validated Liftoff models into native target annotations.

This implements the corrected historical selection rules:
  1. require Liftoff valid_ORFs=1;
  2. reject any gene-span overlap with a native gene on the same strand;
  3. retain the highest coverage*sequence_ID mapping per source gene;
  4. process candidates by descending score and reject a candidate when at
     least 50% of its span overlaps an accepted lifted model on the same strand.

Opposite-strand overlaps are intentionally allowed. Ties are resolved
deterministically by identity, coverage, donor ID, and gene ID. An audit TSV
records every decision.
"""

import argparse
import bisect
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


TRUE = {"1", "true", "yes", "y"}
FILE_RE = re.compile(r"^(.+)_vs_(.+)\.gff\.filtered$")


def parse_attrs(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def manifest_targets(path):
    result = {}
    base = path.parent
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "gff3", "target"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Manifest missing columns: " + ", ".join(sorted(missing)))
        for row in reader:
            if row["target"].strip().lower() not in TRUE:
                continue
            gff3 = Path(row["gff3"])
            gff3 = gff3 if gff3.is_absolute() else (base / gff3).resolve()
            if not gff3.is_file() or gff3.stat().st_size == 0:
                raise ValueError(f"Missing native GFF3 for {row['genome_id']}: {gff3}")
            result[row["genome_id"]] = gff3
    return result


class StaticIntervalIndex:
    """Efficient overlap queries against an immutable set of closed intervals."""

    def __init__(self, intervals):
        ordered = sorted(intervals)
        self.starts = [x[0] for x in ordered]
        self.prefix_max_end = []
        maximum = -1
        for _, end in ordered:
            maximum = max(maximum, end)
            self.prefix_max_end.append(maximum)

    def overlaps(self, start, end):
        index = bisect.bisect_right(self.starts, end) - 1
        return index >= 0 and self.prefix_max_end[index] >= start


def read_native(path):
    lines = path.read_text().splitlines(keepends=True)
    if any(line.startswith("##FASTA") for line in lines):
        raise ValueError(f"Embedded FASTA is not supported in native GFF3: {path}")
    intervals = defaultdict(list)
    genes = 0
    for line in lines:
        columns = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(columns) != 9 or columns[2] != "gene":
            continue
        intervals[(columns[0], columns[6])].append((int(columns[3]), int(columns[4])))
        genes += 1
    return lines, {key: StaticIntervalIndex(value) for key, value in intervals.items()}, genes


def descendants_by_gene(lines):
    """Return complete gene blocks using ID/Parent ancestry, not line order."""
    records = []
    parents = {}
    gene_ids = set()
    gene_record = {}
    for order, line in enumerate(lines):
        columns = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(columns) != 9:
            continue
        feature_attrs = parse_attrs(columns[8])
        feature_id = feature_attrs.get("ID")
        parent_ids = [x for x in feature_attrs.get("Parent", "").split(",") if x]
        if feature_id:
            parents[feature_id] = parent_ids
        record = (order, line, columns, feature_attrs, parent_ids)
        records.append(record)
        if columns[2] == "gene" and feature_id:
            gene_ids.add(feature_id)
            gene_record[feature_id] = record

    ancestry_cache = {}

    def gene_ancestors(node):
        if node in ancestry_cache:
            return ancestry_cache[node]
        found = set()
        pending = [node]
        visited = set()
        while pending:
            item = pending.pop()
            if item in visited:
                continue
            visited.add(item)
            if item in gene_ids:
                found.add(item)
            else:
                pending.extend(parents.get(item, []))
        ancestry_cache[node] = found
        return found

    blocks = {gene_id: [gene_record[gene_id]] for gene_id in gene_ids}
    for record in records:
        if record[2][2] == "gene":
            continue
        ancestors = set()
        for parent in record[4]:
            ancestors.update(gene_ancestors(parent))
        for gene_id in ancestors:
            blocks[gene_id].append(record)
    return {gene_id: sorted(block, key=lambda x: x[0]) for gene_id, block in blocks.items()}


def candidate_file(path, target, donor):
    lines = path.read_text().splitlines(keepends=True)
    candidates = []
    for gene_id, records in descendants_by_gene(lines).items():
        gene_record = next((r for r in records if r[2][2] == "gene"), None)
        if gene_record is None:
            continue
        columns = gene_record[2]
        attributes = gene_record[3]
        try:
            coverage = float(attributes["coverage"])
            identity = float(attributes["sequence_ID"])
            valid_orf = int(attributes.get("valid_ORFs", "0"))
        except (KeyError, TypeError, ValueError):
            coverage = identity = 0.0
            valid_orf = 0
        candidates.append({
            "gene": gene_id,
            "target": target,
            "donor": donor,
            "seqid": columns[0],
            "start": int(columns[3]),
            "end": int(columns[4]),
            "strand": columns[6],
            "coverage": coverage,
            "identity": identity,
            "valid_orf": valid_orf,
            "score": coverage * identity,
            "records": records,
            "file": path.name,
        })
    return candidates


def native_overlap(candidate, indexes):
    index = indexes.get((candidate["seqid"], candidate["strand"]))
    return bool(index and index.overlaps(candidate["start"], candidate["end"]))


def overlap_fraction(candidate, accepted):
    overlap = min(candidate["end"], accepted["end"]) - max(candidate["start"], accepted["start"]) + 1
    if overlap <= 0:
        return 0.0
    return overlap / (candidate["end"] - candidate["start"] + 1)


def annotated_lines(candidate):
    additions = {
        "liftoff_source_genome": candidate["donor"],
        "liftoff_target_genome": candidate["target"],
        "liftoff_coverage": f"{candidate['coverage']:.3f}",
        "liftoff_sequence_ID": f"{candidate['identity']:.3f}",
        "liftoff_valid_ORF": str(candidate["valid_orf"]),
        "liftoff_score": f"{candidate['score']:.3f}",
        "liftoff_status": "novel",
    }
    output = []
    for _, line, columns, attributes, _ in candidate["records"]:
        attributes.update(additions)
        columns = list(columns)
        columns[8] = ";".join(f"{key}={value}" for key, value in attributes.items()) + ";"
        output.append("\t".join(columns) + "\n")
    return output


def rank(candidate):
    return (-candidate["score"], -candidate["identity"], -candidate["coverage"],
            candidate["donor"], candidate["gene"], candidate["seqid"], candidate["start"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--novel-overlap-fraction", type=float, default=0.50)
    args = parser.parse_args()
    targets = manifest_targets(args.manifest.resolve())
    args.output_dir.mkdir(parents=True, exist_ok=True)

    files_by_target = defaultdict(list)
    for path in sorted(args.input_dir.glob("*.gff.filtered")):
        match = FILE_RE.match(path.name)
        if match:
            files_by_target[match.group(1)].append((path, match.group(2)))

    audit_path = args.output_dir / "04_merge_audit.tsv"
    fields = ["target", "donor", "gene", "seqid", "start", "end", "strand",
              "coverage", "sequence_ID", "score", "decision", "reason"]
    with audit_path.open("w", newline="") as audit_handle:
        audit = csv.DictWriter(audit_handle, fieldnames=fields, delimiter="\t")
        audit.writeheader()
        for target, native_path in sorted(targets.items()):
            native_lines, native_indexes, native_count = read_native(native_path)
            candidates = []
            for path, donor in files_by_target.get(target, []):
                candidates.extend(candidate_file(path, target, donor))

            decisions = {}
            eligible = []
            for candidate in candidates:
                key = id(candidate)
                if candidate["valid_orf"] != 1:
                    decisions[key] = ("REMOVE", "invalid_orf")
                elif native_overlap(candidate, native_indexes):
                    decisions[key] = ("REMOVE", "same_strand_native_overlap")
                else:
                    eligible.append(candidate)

            grouped = defaultdict(list)
            for candidate in eligible:
                grouped[candidate["gene"]].append(candidate)
            best_per_gene = []
            for versions in grouped.values():
                ordered = sorted(versions, key=rank)
                best_per_gene.append(ordered[0])
                for candidate in ordered[1:]:
                    decisions[id(candidate)] = ("REMOVE", "lower_score_same_source_gene")

            accepted = []
            accepted_by_location = defaultdict(list)
            for candidate in sorted(best_per_gene, key=rank):
                conflict = any(
                    overlap_fraction(candidate, previous) >= args.novel_overlap_fraction
                    for previous in accepted_by_location[(candidate["seqid"], candidate["strand"])]
                )
                if conflict:
                    decisions[id(candidate)] = ("REMOVE", "overlap_with_higher_ranked_lifted_model")
                else:
                    accepted.append(candidate)
                    accepted_by_location[(candidate["seqid"], candidate["strand"])].append(candidate)
                    decisions[id(candidate)] = ("KEEP", "accepted")

            output = args.output_dir / f"{target}.merged.gff3"
            with output.open("w") as out:
                out.writelines(native_lines)
                if native_lines and not native_lines[-1].endswith("\n"):
                    out.write("\n")
                out.write("# Novel genes added from Liftoff\n")
                for candidate in sorted(accepted, key=lambda x: (x["seqid"], x["start"], x["end"], x["gene"])):
                    out.writelines(annotated_lines(candidate))

            for candidate in sorted(candidates, key=lambda x: (x["donor"], x["gene"], x["seqid"], x["start"])):
                decision, reason = decisions[id(candidate)]
                audit.writerow({
                    "target": target, "donor": candidate["donor"], "gene": candidate["gene"],
                    "seqid": candidate["seqid"], "start": candidate["start"], "end": candidate["end"],
                    "strand": candidate["strand"], "coverage": candidate["coverage"],
                    "sequence_ID": candidate["identity"], "score": candidate["score"],
                    "decision": decision, "reason": reason,
                })
            print(f"{target}\tnative={native_count}\tcandidates={len(candidates)}\taccepted={len(accepted)}\toutput={output}")
    print(f"audit={audit_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
