#!/bin/bash

# compute_flagstat.sh
# --------------------
# Extracts mapped read metrics from a BAM file using samtools flagstat.
# Outputs a tab-separated summary line for appending to a master TSV file.
#
# Usage:
#   ./extract_flagstat.sh input.bam output.tsv
#
# Requirements:
#   - samtools must be in PATH
#   - BAM index (.bai) should exist alongside BAM

set -euo pipefail

# Input arguments
BAM_FILE="$1"
OUTPUT_TSV="$2"

# Get sample name from BAM filename
SAMPLE_NAME=$(basename "$BAM_FILE" .bam)

# Run flagstat
FLAGSTAT_OUTPUT=$(samtools flagstat "$BAM_FILE")

# Extract total reads and mapped reads
TOTAL_READS=$(echo "$FLAGSTAT_OUTPUT" | grep "in total" | cut -d ' ' -f 1)
MAPPED_READS=$(echo "$FLAGSTAT_OUTPUT" | grep "mapped (" | cut -d ' ' -f 1)
PROPERLY_PAIRED=$(echo "$FLAGSTAT_OUTPUT" | grep "properly paired (" | cut -d ' ' -f 1)

# Calculate % mapped (as float)
MAPPED_PCT=$(awk -v m="$MAPPED_READS" -v t="$TOTAL_READS" 'BEGIN {printf "%.2f", (m/t)*100}')

# Output header if file does not exist
if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo -e "Sample\tTotal_Reads\tMapped_Reads\tMapped_Percent\tProperly_Paired" > "$OUTPUT_TSV"
fi

# Append results
echo -e "${SAMPLE_NAME}\t${TOTAL_READS}\t${MAPPED_READS}\t${MAPPED_PCT}\t${PROPERLY_PAIRED}" >> "$OUTPUT_TSV"
