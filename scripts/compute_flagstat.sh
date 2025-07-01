#!/bin/bash
set -euo pipefail

BAM_FILE="$1"
OUTPUT_TSV="$2"
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

if [[ ! -f "$BAM_FILE" ]]; then
    echo "Error: BAM file '$BAM_FILE' not found." >&2
    exit 1
fi

FLAGSTAT_OUTPUT=$(samtools flagstat "$BAM_FILE")
TOTAL_READS=$(echo "$FLAGSTAT_OUTPUT" | grep "in total" | cut -d ' ' -f 1)
MAPPED_READS=$(echo "$FLAGSTAT_OUTPUT" | grep -E '^[0-9]+ \+ [0-9]+ mapped \(' | cut -d ' ' -f 1)
PROPERLY_PAIRED=$(echo "$FLAGSTAT_OUTPUT" | grep "properly paired (" | cut -d ' ' -f 1)
MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" 'BEGIN {printf "%.2f", (m/t)*100}')

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Percent\tProperly_Paired" > "$OUTPUT_TSV"
fi

echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAPPED_PCT}\t${PROPERLY_PAIRED}" >> "$OUTPUT_TSV"
