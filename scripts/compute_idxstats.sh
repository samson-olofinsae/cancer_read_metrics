#!/bin/bash

# compute_idxstats.sh
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

#!/bin/bash
set -euo pipefail

BAM_FILE="$1"
OUTPUT_TSV="$2"
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file '$BAM_FILE' not found." >&2
    exit 1
fi

MAPPED_READS=$(samtools idxstats "$BAM_FILE" | awk '$1 != "*" {sum += $3} END {print sum}')

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo -e "Sample\tMapped_Reads_Idxstats" > "$OUTPUT_TSV"
fi

echo -e "${SAMPLE_NAME}\t${MAPPED_READS}" >> "$OUTPUT_TSV"
