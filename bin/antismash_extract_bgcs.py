#!/usr/bin/env python3

"""
Extract BGC FASTA sequences from contigs using antiSMASH regions TSVs.

Usage:
    For each sample, slice contig sequences to region coordinates and emit a
    BGC_FASTAs/<sample>_regions.fa file that can be used as a CoverM reference.

Inputs:
    --sample-list   Text file with one sample ID per line.
    --megahit-dir   Base directory with <sample>/final.contigs.fa.
    --antismash-dir Base directory with <sample>/<sample>_regions.tsv.

Outputs:
    <antismash-dir>/<sample>/BGC_FASTAs/<sample>_regions.fa

Requirements:
    Python 3 and Biopython (SeqIO).
"""


import argparse
import csv
from pathlib import Path
from Bio import SeqIO

def load_regions(regions_tsv):
    regions = []
    with open(regions_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                regions.append({
                    "sample": row["sample"],
                    "record_id": row["record_id"],
                    "region_number": row["region_number"],
                    "start": int(row["region_start"]),
                    "end": int(row["region_end"]),
                    "product": row.get("product", "NA"),
                })
            except Exception as e:
                print(f"[WARN] Skipping row in {regions_tsv}: {e}")
    return regions

def extract_regions_fasta(sample, megahit_dir, antismash_dir, verbose=False):
    contigs_fa = Path(megahit_dir) / sample / "final.contigs.fa"
    regions_tsv = Path(antismash_dir) / sample / f"{sample}_regions.tsv"
    out_dir = Path(antismash_dir) / sample / "BGC_FASTAs"

    # Ensure parent and BGC_FASTAs folder exist
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(exist_ok=True)

    if not contigs_fa.exists():
        print(f"[WARN] Missing contigs fasta: {contigs_fa}")
        return
    if not regions_tsv.exists():
        print(f"[WARN] Missing regions.tsv: {regions_tsv}")
        return

    contigs = SeqIO.to_dict(SeqIO.parse(str(contigs_fa), "fasta"))
    regions = load_regions(regions_tsv)

    out_fa = out_dir / f"{sample}_regions.fa"
    with open(out_fa, "w") as out:
        for r in regions:
            rid = r["record_id"]
            if rid not in contigs:
                if verbose:
                    print(f"[WARN] Record {rid} not found in {contigs_fa}")
                continue
            seq = contigs[rid].seq[r["start"]-1:r["end"]]
            header = f">{sample}|{rid}|region{r['region_number']}|{r['start']}-{r['end']}|{r['product']}"
            out.write(header + "\n")
            out.write(str(seq) + "\n")

    print(f"[INFO] {sample}: wrote {len(regions)} regions -> {out_fa}")

def main():
    parser = argparse.ArgumentParser(description="Extract BGC FASTAs for Healthy samples")
    parser.add_argument("--sample-list", "-s", default=None)
    parser.add_argument("--megahit-dir", "-m", default=None)
    parser.add_argument("--antismash-dir", "-a", default=None)
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    with open(args.sample_list) as f:
        samples = [ln.strip() for ln in f if ln.strip()]

    for sample in samples:
        extract_regions_fasta(sample, args.megahit_dir, args.antismash_dir, verbose=args.verbose)

if __name__ == "__main__":
    main()
