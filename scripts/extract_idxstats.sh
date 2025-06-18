#!/bin/bash

# extract_idxstats.sh
# --------------------
# Extracts mapped read count from a BAM file using samtools idxstats.
# Outputs a tab-separated summary line for appending to a master TSV file.
#
# Usage:
#   ./extract_idxstats.sh input.bam output.tsv
#
# Requirements:
#   - samtools must be in PATH
#   - BAM index (.bai) must exist

set -euo pipefail

# Input arguments
BAM_FILE="$1"
OUTPUT_TSV="$2"

# Get sample name
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

# Run idxstats and sum mapped reads (exclude '*' which lists unmapped)
MAPPED_READS=$(samtools idxstats "$BAM_FILE" | awk '$1 != "*" {sum += $3} END {print sum}')

# Output header if file doesn’t exist
if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo -e "Sample\tMapped_Reads_Idxstats" > "$OUTPUT_TSV"
fi

# Append result
echo -e "${SAMPLE_NAME}\t${MAPPED_READS}
