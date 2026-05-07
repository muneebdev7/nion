#!/usr/bin/env python3

"""
Process antiSMASH .gbk results into per-sample TSV summaries.

Usage:
    Parse antiSMASH GBK files per sample and generate region/gene/summary tables.

Inputs:
    --sample-list   Text file with one sample ID per line.
    --output-base   Base directory containing antismash_output/<sample>/*.gbk.

Outputs (written into antismash_output/<sample>/):
    <sample>_regions.tsv  Region-level metadata.
    <sample>_genes.tsv    Gene-level metadata.
    <sample>_summary.tsv  Region summary.
    <sample>_merged.tsv   Gene table merged with region/summary fields.

Requirements:
    Python 3 and Biopython (SeqIO).
"""


from __future__ import print_function
import argparse
import csv
import sys
from pathlib import Path
from Bio import SeqIO


def safe_first(x, default="NA"):
    if x is None:
        return default
    if isinstance(x, (list, tuple)):
        return x[0] if x else default
    return x


def extract_domains_from_qual(q):
    if not isinstance(q, dict):
        return "NA"
    doms = []
    for k in ("db_xref", "Dbxref", "dbxref"):
        if k in q:
            vals = q.get(k) or []
            for v in vals:
                s = str(v)
                if "PFAM" in s.upper() or "Pfam" in s or "pfam" in s:
                    doms.append(s)
    if "note" in q:
        notes = q.get("note") or []
        if isinstance(notes, (list, tuple)):
            for n in notes:
                s = str(n)
                if "pfam" in s.lower() or "domain" in s.lower() or "hmm" in s.lower():
                    doms.append(s)
    return ";".join(doms) if doms else "NA"


def parse_one_sample(sample_dir, sample, verbose=False):
    sd = Path(sample_dir)
    region_rows, gene_rows, summary_rows = [], [], []

    if not sd.exists() or not sd.is_dir():
        if verbose:
            print("WARN: sample dir missing:", sd, file=sys.stderr)
        return region_rows, gene_rows, summary_rows

    gbk_files = sorted(sd.glob("*.region*.gbk"))
    if not gbk_files:
        gbk_files = sorted(sd.glob("*.gbk"))

    if not gbk_files:
        if verbose:
            print("WARN: no .gbk files found in", sd, file=sys.stderr)
        return region_rows, gene_rows, summary_rows

    for gbk in gbk_files:
        try:
            for rec in SeqIO.parse(str(gbk), "genbank"):
                rec_id = rec.id or rec.name or gbk.stem
                rec_desc = (rec.description or "").strip()
                found_region_for_record = False

                # Regions
                for feat in rec.features:
                    if feat.type.lower() == "region":
                        found_region_for_record = True
                        q = feat.qualifiers
                        region_number = safe_first(q.get("region_number"), "NA")
                        try:
                            region_start = int(feat.location.start) + 1
                            region_end = int(feat.location.end)
                        except Exception:
                            region_start, region_end = 1, len(rec.seq)
                        contig_edge = safe_first(q.get("contig_edge"), "NA")
                        product = safe_first(q.get("product"), "NA")
                        notes_list = []
                        for kk in ("note", "notes", "similarity", "similarities"):
                            if kk in q:
                                v = q.get(kk) or []
                                if isinstance(v, (list, tuple)):
                                    notes_list.extend([str(x) for x in v])
                                else:
                                    notes_list.append(str(v))
                        kcb = safe_first(q.get("knownclusterblast"), "NA")
                        kcb_acc = safe_first(q.get("knownclusterblast_acc") or q.get("knownclusterblast_accession"), "NA")
                        kcb_sim = safe_first(q.get("knownclusterblast_sim") or q.get("knownclusterblast_similarity") or q.get("similarity"), "NA")
                        mibig = safe_first(q.get("mibig_accession") or q.get("mibig_id") or q.get("mibig"), "NA")

                        region_rows.append({
                            "sample": sample,
                            "record_id": rec_id,
                            "gbk_file": gbk.name,
                            "region_number": region_number,
                            "region_start": region_start,
                            "region_end": region_end,
                            "contig_edge": contig_edge,
                            "product": product,
                            "notes": ";".join(notes_list) if notes_list else "NA",
                            "kcb_hit": kcb,
                            "kcb_acc": kcb_acc,
                            "kcb_sim": kcb_sim,
                            "mibig_hit": mibig,
                        })
                        summary_rows.append({
                            "file": gbk.name,
                            "record_id": rec_id,
                            "region": region_number,
                            "start": region_start,
                            "end": region_end,
                            "contig_edge": contig_edge,
                            "product": product,
                            "KCB_hit": kcb,
                            "KCB_acc": kcb_acc,
                            "KCB_%": kcb_sim,
                            "record_desc": rec_desc or "NA",
                        })

                if not found_region_for_record:
                    region_rows.append({
                        "sample": sample,
                        "record_id": rec_id,
                        "gbk_file": gbk.name,
                        "region_number": "NA",
                        "region_start": 1,
                        "region_end": len(rec.seq),
                        "contig_edge": "NA",
                        "product": rec.description or "NA",
                        "notes": "NA",
                        "kcb_hit": "NA",
                        "kcb_acc": "NA",
                        "kcb_sim": "NA",
                        "mibig_hit": "NA",
                    })
                    summary_rows.append({
                        "file": gbk.name,
                        "record_id": rec_id,
                        "region": "NA",
                        "start": 1,
                        "end": len(rec.seq),
                        "contig_edge": "NA",
                        "product": rec.description or "NA",
                        "KCB_hit": "NA",
                        "KCB_acc": "NA",
                        "KCB_%": "NA",
                        "record_desc": rec_desc or "NA",
                    })

                # Genes
                for feat in rec.features:
                    if feat.type.lower() == "cds":
                        q = feat.qualifiers
                        try:
                            gstart = int(feat.location.start) + 1
                            gend = int(feat.location.end)
                        except Exception:
                            gstart, gend = "NA", "NA"
                        gstrand = feat.location.strand
                        locus = safe_first(q.get("locus_tag") or q.get("gene") or q.get("locus"), "NA")
                        protein_id = safe_first(q.get("protein_id"), "NA")
                        gene_id = locus if locus != "NA" else (protein_id or "NA")
                        annotation = safe_first(q.get("product") or q.get("note") or q.get("function"), "NA")
                        domains = extract_domains_from_qual(q)
                        gene_rows.append({
                            "sample": sample,
                            "record_id": rec_id,
                            "gbk_file": gbk.name,
                            "gene_id": gene_id,
                            "gene_start": gstart,
                            "gene_end": gend,
                            "strand": gstrand,
                            "annotation": annotation,
                            "domains": domains
                        })
        except Exception as e:
            if verbose:
                print("WARN: failed to parse", gbk, ":", e, file=sys.stderr)

    return region_rows, gene_rows, summary_rows


def write_tsv_dicts(path, rows, header):
    p = Path(path)
    with p.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow([r.get(c, "NA") for c in header])


def merge_and_write(sample_dir, sample, region_rows, gene_rows, summary_rows, verbose=False):
    region_map = {(r.get("sample"), r.get("record_id"), r.get("gbk_file")): r for r in region_rows}
    summary_map = {(sample, s.get("record_id")): s for s in summary_rows}

    merged_rows = []
    out_header = [
        "sample", "record_id", "gbk_file",
        "gene_id", "gene_start", "gene_end", "strand", "annotation", "domains",
        "region_start", "region_end", "product", "notes",
        "KCB_hit", "KCB_acc", "KCB_%", "record_desc"
    ]

    for g in gene_rows:
        key = (g.get("sample"), g.get("record_id"), g.get("gbk_file"))
        r = region_map.get(key, {})
        s = summary_map.get((g.get("sample"), g.get("record_id")), {})
        merged = {
            "sample": g.get("sample"),
            "record_id": g.get("record_id"),
            "gbk_file": g.get("gbk_file"),
            "gene_id": g.get("gene_id"),
            "gene_start": g.get("gene_start"),
            "gene_end": g.get("gene_end"),
            "strand": g.get("strand"),
            "annotation": g.get("annotation"),
            "domains": g.get("domains"),
            "region_start": r.get("region_start", "NA"),
            "region_end": r.get("region_end", "NA"),
            "product": r.get("product", "NA"),
            "notes": r.get("notes", "NA"),
            "KCB_hit": s.get("KCB_hit", r.get("kcb_hit", "NA")),
            "KCB_acc": s.get("KCB_acc", r.get("kcb_acc", "NA")),
            "KCB_%": s.get("KCB_%", r.get("kcb_sim", "NA")),
            "record_desc": s.get("record_desc", "NA")
        }
        merged_rows.append(merged)

    merged_path = Path(sample_dir) / f"{sample}_merged.tsv"
    write_tsv_dicts(merged_path, merged_rows, out_header)
    if verbose:
        print("WROTE merged rows:", len(merged_rows), "->", merged_path)
    return merged_path, len(merged_rows)


# Main driver
def main():
    parser = argparse.ArgumentParser(description="Process all samples: extract regions/genes/summary from .gbk and merge per sample")
    parser.add_argument("--sample-list", "-s", default=None, help="Path to sample list (one sample per line)")
    parser.add_argument("--output-dir", "-o", default=None, help=" Antismash output dir containing sample subdirs")
    parser.add_argument("--logs-dir", "-l", default=None, help="Optional logs dir (defaults to <base_dir>/logs)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    sample_list = Path(args.sample_list)
    output_directory = Path(args.output_dir)
    logs_dir = Path(args.logs_dir) if args.logs_dir else output_directory / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    if not sample_list.exists():
        print("ERROR: sample list not found:", sample_list, file=sys.stderr)
        sys.exit(2)

    with sample_list.open() as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]

    for sample in samples:
        sample_dir = output_directory / sample
        log_file = logs_dir / f"{sample}.log"
        try:
            print("\n[INFO] Processing sample:", sample)
            if not sample_dir.exists():
                print("WARN: sample dir missing:", sample_dir)
                with log_file.open("w") as L:
                    L.write("MISSING_SAMPLE_DIR\n")
                continue

            region_rows, gene_rows, summary_rows = parse_one_sample(sample_dir, sample, verbose=args.verbose)

            write_tsv_dicts(sample_dir / f"{sample}_regions.tsv", region_rows, [
                "sample", "record_id", "gbk_file", "region_number", "region_start", "region_end",
                "contig_edge", "product", "notes", "kcb_hit", "kcb_acc", "kcb_sim", "mibig_hit"
            ])
            write_tsv_dicts(sample_dir / f"{sample}_genes.tsv", gene_rows, [
                "sample", "record_id", "gbk_file", "gene_id", "gene_start", "gene_end", "strand", "annotation", "domains"
            ])
            write_tsv_dicts(sample_dir / f"{sample}_summary.tsv", summary_rows, [
                "file", "record_id", "region", "start", "end", "contig_edge", "product", "KCB_hit", "KCB_acc", "KCB_%", "record_desc"
            ])

            merge_and_write(sample_dir, sample, region_rows, gene_rows, summary_rows, verbose=args.verbose)

            with log_file.open("w") as L:
                L.write("sample\tregions\tgenes\tmerged_rows\n")
                L.write(f"{sample}\t{len(region_rows)}\t{len(gene_rows)}\t{len(gene_rows)}\n")

            print("[DONE] sample:", sample, "regions:", len(region_rows), "genes:", len(gene_rows), "merged_rows:", len(gene_rows))
        except Exception as e:
            print("ERROR processing sample", sample, ":", e, file=sys.stderr)
            with log_file.open("w") as L:
                L.write("ERROR: " + str(e) + "\n")
            continue


if __name__ == "__main__":
    main()
