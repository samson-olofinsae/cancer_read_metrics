#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/compute_flagstat.sh input.bam output.txt

BAM_FILE="$1"
OUTPUT_TXT="$2"

if [[ ! -f "${BAM_FILE}" ]]; then
    echo "Error: BAM file '${BAM_FILE}' not found." >&2
    exit 1
fi

mkdir -p "$(dirname "${OUTPUT_TXT}")"

samtools flagstat "${BAM_FILE}" > "${OUTPUT_TXT}"

echo "Flagstat written to ${OUTPUT_TXT}"