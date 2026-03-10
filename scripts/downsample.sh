#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash scripts/downsample.sh input.bam fraction output.bam
#
# Example:
#   bash scripts/downsample.sh sample.bam 0.1 sample_10pct.bam

INPUT_BAM="$1"
FRACTION="$2"
OUTPUT_BAM="$3"

if [[ ! -f "${INPUT_BAM}" ]]; then
    echo "Error: input BAM '${INPUT_BAM}' not found." >&2
    exit 1
fi

# Ensure output directory exists
mkdir -p "$(dirname "${OUTPUT_BAM}")"

# Downsample using samtools
samtools view -@ 2 -s "${FRACTION}" -b "${INPUT_BAM}" -o "${OUTPUT_BAM}"

# Index output BAM
samtools index "${OUTPUT_BAM}"

echo "Downsampled ${INPUT_BAM} -> ${OUTPUT_BAM} at fraction ${FRACTION}"