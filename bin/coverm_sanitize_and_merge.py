#!/usr/bin/env python3

"""
Fix CoverM depth table headers and merge into a single matrix.

Usage:
    coverm_sanitize_and_merge.py -i <dir> -o <merged.tsv>
    coverm_sanitize_and_merge.py -o <merged.tsv> <file1> <file2> ...

Inputs:
    - One or more CoverM depth tables (*.depth.tsv or *.tsv).
    - Optional input directory with glob pattern(s) to discover tables.

Outputs:
    - A merged contig-by-sample TSV.
    - Input tables with standardized header ("Contig", "TPM").

Requirements:
    - Python 3.

The script standardizes headers to "Contig" and "TPM" before merging.
"""


import argparse
import csv
import sys
from pathlib import Path


def parse_patterns(pattern_text):
    return [p.strip() for p in pattern_text.split(",") if p.strip()]


def collect_input_files(input_dir, patterns):
    input_dir = Path(input_dir)
    files = []
    for pattern in patterns:
        files.extend(sorted(input_dir.glob(pattern)))
    unique = []
    seen = set()
    for path in files:
        key = str(path.resolve())
        if key in seen:
            continue
        seen.add(key)
        unique.append(path)
    return unique


def derive_sample_name(path):
    name = Path(path).name
    for suffix in (".depth.tsv", ".tsv", ".txt"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_coverm"):
        name = name[: -len("_coverm")]
    return name


def standardize_header(path, contig_label="Contig", metric_label="TPM"):
    with open(path, "r") as handle:
        lines = handle.readlines()
    
    if not lines:
        print(f"WARNING: File is empty: {path}", file=sys.stderr)
        return False
    
    header_fields = lines[0].rstrip("\n").split("\t")
    if len(header_fields) < 2:
        print(f"WARNING: Header has fewer than two columns: {path}", file=sys.stderr)
        return False
    
    header_fields[0] = contig_label
    header_fields[1] = metric_label
    new_header = "\t".join(header_fields) + "\n"
    
    if lines[0] != new_header:
        lines[0] = new_header
        with open(path, "w") as handle:
            handle.writelines(lines)
    return True


def read_depth_tsv(path, metric=None):
    with open(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header in {path}")
        fieldnames = reader.fieldnames
        field_map = {name.lower(): name for name in fieldnames}
        
        contig_col = field_map.get("contig", fieldnames[0])
        
        metric_col = None
        if metric:
            metric_col = field_map.get(metric.lower())
        else:
            metric_col = field_map.get("tpm") or field_map.get("rpkm")
        
        if not metric_col:
            metric_col = next((c for c in fieldnames if c != contig_col), None)
        
        if not metric_col or metric_col not in fieldnames:
            raise ValueError(f"Metric column '{metric_col}' not found in {path}")
        
        data = {}
        order = []
        for row in reader:
            contig = row.get(contig_col)
            if not contig:
                continue
            if contig not in data:
                order.append(contig)
            data[contig] = row.get(metric_col, "")
        return data, order, metric_col


def main():
    parser = argparse.ArgumentParser(description="Fix CoverM headers and merge contig depth tables")
    parser.add_argument("input_files", nargs="*", help="CoverM depth tables (.tsv)")
    parser.add_argument("-i", "--input-dir", default=None, help="Directory with CoverM tables")
    parser.add_argument("--pattern", default="*_coverm*.tsv", help="Glob pattern(s) for --input-dir (comma-separated)")
    parser.add_argument("-o", "--output", required=True, help="Output merged TSV")
    parser.add_argument("--metric", default=None, help="Metric column to extract (default: tpm if present)")
    parser.add_argument("--sample-names", default=None, help="Comma-separated sample names in input order")
    parser.add_argument("--fill", default="0", help="Fill value for missing contigs")
    args = parser.parse_args()
    
    files = []
    if args.input_dir:
        patterns = parse_patterns(args.pattern)
        files.extend(collect_input_files(args.input_dir, patterns))
        
    files.extend([Path(p) for p in args.input_files])
    
    if not files:
        print("ERROR: No input files found", file=sys.stderr)
        sys.exit(2)
    
    for path in files:
        if not path.exists():
            print(f"ERROR: Input file not found: {path}", file=sys.stderr)
            sys.exit(2)
    
    for path in files:
        standardize_header(path)
    
    if args.sample_names:
        samples = [s.strip() for s in args.sample_names.split(",") if s.strip()]
        if len(samples) != len(files):
            print("ERROR: --sample-names count must match number of input files", file=sys.stderr)
            sys.exit(2)
    else:
        samples = [derive_sample_name(p) for p in files]
    
    merged = {}
    contig_order = []
    metric = args.metric
    
    for path, sample in zip(files, samples):
        data, order, metric_col = read_depth_tsv(path, metric)
        if metric is None:
            metric = metric_col
        merged[sample] = data
        for contig in order:
            if contig not in contig_order:
                contig_order.append(contig)
    
    for sample in samples:
        for contig in merged[sample].keys():
            if contig not in contig_order:
                contig_order.append(contig)
    
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["Contig"] + samples)
        for contig in contig_order:
            row = [contig]
            for sample in samples:
                row.append(merged[sample].get(contig, args.fill))
            writer.writerow(row)


if __name__ == "__main__":
    main()
