#!/usr/bin/env python3

"""
compare_metrics.py

Compares mapped read counts from samtools flagstat and idxstats outputs.

Usage:
    python scripts/compare_metrics.py \
        --flagstat results/example_flagstat.txt \
        --idxstats results/example_idxstats.txt \
        --out results/example_comparison.tsv

Arguments:
    --flagstat   Path to TSV file from compute_flagstat.sh
    --idxstats   Path to TSV file from compute_idxstats.sh
    --out        Output file path for the comparison result

Author: Samson Olofinsae
"""


import csv
import argparse

def read_flagstat(path):
    data = {}
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["Sample"]
            data[sample] = {
                "total_reads": int(row["Total_Reads"]),
                "mapped_reads_flagstat": int(row["Mapped_Reads"]),
                "mapped_percent": float(row["Mapped_Percent"]),
                "properly_paired": int(row["Properly_Paired"]),
            }
    return data

def read_idxstats(path):
    data = {}
    with open(path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample = row["Sample"]
            data[sample] = int(row["Mapped_Reads_Idxstats"])
    return data

def compare_metrics(flagstat_data, idxstats_data):
    rows = []
    for sample in flagstat_data:
        if sample not in idxstats_data:
            continue
        fs = flagstat_data[sample]
        idx = idxstats_data[sample]
        abs_diff = fs["mapped_reads_flagstat"] - idx
        pct_diff = (abs(abs_diff) / fs["mapped_reads_flagstat"]) * 100 if fs["mapped_reads_flagstat"] else 0.0

        rows.append({
            "Sample": sample,
            "Flagstat_Mapped": fs["mapped_reads_flagstat"],
            "Idxstats_Mapped": idx,
            "Absolute_Diff": abs_diff,
            "Percent_Diff": f"{pct_diff:.2f}",
        })
    return rows

def write_output(rows, path):
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Sample", "Flagstat_Mapped", "Idxstats_Mapped", "Absolute_Diff", "Percent_Diff"
        ], delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--idxstats", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    flagstat_data = read_flagstat(args.flagstat)
    idxstats_data = read_idxstats(args.idxstats)
    results = compare_metrics(flagstat_data, idxstats_data)
    write_output(results, args.out)

if __name__ == "__main__":
    main()
