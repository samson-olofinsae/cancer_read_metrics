"""
test_compare_metrics.py

Unit test for compare_metrics.py, which compares mapped read counts between
samtools flagstat and idxstats output files.

Purpose:
    This test ensures that the compare_metrics() function correctly:
    - Matches samples between flagstat and idxstats TSV files
    - Computes the absolute and percentage differences in mapped read counts
    - Handles simple, valid inputs without error

Usage:
    Run the test with pytest:
        pytest tests/test_compare_metrics.py

Author:
    Samson Olofinsae
"""

import os
import sys
import tempfile
import csv

# Fix: Add the project root to the Python module path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from scripts.compare_metrics import read_flagstat, read_idxstats, compare_metrics

def create_temp_tsv(headers, rows):
    tmp = tempfile.NamedTemporaryFile(mode='w+', delete=False, newline='')
    writer = csv.DictWriter(tmp, fieldnames=headers, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    tmp.close()
    return tmp.name

def test_compare_metrics_match_exact():
    # Simulate flagstat and idxstats TSVs with matching mapped read counts
    flagstat_headers = ["Sample", "Total_Reads", "Mapped_Reads", "Mapped_Percent", "Properly_Paired"]
    idxstats_headers = ["Sample", "Mapped_Reads_Idxstats"]

    flagstat_rows = [{
        "Sample": "testsample",
        "Total_Reads": "100",
        "Mapped_Reads": "90",
        "Mapped_Percent": "90.00",
        "Properly_Paired": "45"
    }]

    idxstats_rows = [{
        "Sample": "testsample",
        "Mapped_Reads_Idxstats": "90"
    }]

    flagstat_file = create_temp_tsv(flagstat_headers, flagstat_rows)
    idxstats_file = create_temp_tsv(idxstats_headers, idxstats_rows)

    try:
        fs_data = read_flagstat(flagstat_file)
        idx_data = read_idxstats(idxstats_file)
        result = compare_metrics(fs_data, idx_data)

        assert len(result) == 1
        row = result[0]
        assert row["Sample"] == "testsample"
        assert row["Flagstat_Mapped"] == 90
        assert row["Idxstats_Mapped"] == 90
        assert row["Absolute_Diff"] == 0
        assert float(row["Percent_Diff"]) == 0.00
    finally:
        os.remove(flagstat_file)
        os.remove(idxstats_file)
