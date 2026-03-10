#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# Cancer Read Metrics Pipeline
# Full analysis runner
# -----------------------------

RAW_DIR="data/raw_bams"
DOWNSAMPLED_DIR="data/downsampled"
RESULTS_DIR="results"
FLAGSTAT_DIR="${RESULTS_DIR}/flagstat"
IDXSTATS_DIR="${RESULTS_DIR}/idxstats"

METRICS_TABLE="${RESULTS_DIR}/metrics_table.tsv"
SUMMARY_TABLE="${RESULTS_DIR}/discrepancy_summary.tsv"

DOWNSAMPLE_SCRIPT="scripts/downsample.sh"
FLAGSTAT_SCRIPT="scripts/compute_flagstat.sh"
IDXSTATS_SCRIPT="scripts/compute_idxstats.sh"

# Downsampling fractions
FRACTIONS=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

# -----------------------------
# Checks
# -----------------------------
echo "Checking required files and directories..."

[[ -d "${RAW_DIR}" ]] || { echo "ERROR: ${RAW_DIR} not found"; exit 1; }
[[ -f "${DOWNSAMPLE_SCRIPT}" ]] || { echo "ERROR: ${DOWNSAMPLE_SCRIPT} not found"; exit 1; }
[[ -f "${FLAGSTAT_SCRIPT}" ]] || { echo "ERROR: ${FLAGSTAT_SCRIPT} not found"; exit 1; }
[[ -f "${IDXSTATS_SCRIPT}" ]] || { echo "ERROR: ${IDXSTATS_SCRIPT} not found"; exit 1; }

mkdir -p "${DOWNSAMPLED_DIR}" "${RESULTS_DIR}" "${FLAGSTAT_DIR}" "${IDXSTATS_DIR}"

# -----------------------------
# Initialise output tables
# -----------------------------
echo -e "sample\tdepth_fraction\tdepth_percent\tdownsampled_bam\tflagstat_total\tflagstat_mapped\tidxstats_total\tidxstats_mapped\tmapped_diff\tmapped_diff_pct" > "${METRICS_TABLE}"

# -----------------------------
# Main loop
# -----------------------------
echo "Starting analysis..."

shopt -s nullglob
BAM_FILES=("${RAW_DIR}"/*.bam)

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${RAW_DIR}"
    exit 1
fi

for BAM in "${BAM_FILES[@]}"; do
    SAMPLE=$(basename "${BAM}" .bam)
    echo "Processing sample: ${SAMPLE}"

    # Ensure BAM index exists
    if [[ ! -f "${BAM}.bai" ]]; then
        echo "Index not found for ${BAM}. Creating index..."
        samtools index "${BAM}"
    fi

    for FRACTION in "${FRACTIONS[@]}"; do
        DEPTH_PERCENT=$(python3 - <<PY
fraction = float("${FRACTION}")
print(int(fraction * 100))
PY
)

        DOWNSAMPLED_BAM="${DOWNSAMPLED_DIR}/${SAMPLE}_${DEPTH_PERCENT}pct.bam"
        FLAGSTAT_OUT="${FLAGSTAT_DIR}/${SAMPLE}_${DEPTH_PERCENT}pct.flagstat.txt"
        IDXSTATS_OUT="${IDXSTATS_DIR}/${SAMPLE}_${DEPTH_PERCENT}pct.idxstats.txt"

        echo "  Downsampling ${SAMPLE} to ${DEPTH_PERCENT}%..."

        bash "${DOWNSAMPLE_SCRIPT}" "${BAM}" "${FRACTION}" "${DOWNSAMPLED_BAM}"

        if [[ ! -f "${DOWNSAMPLED_BAM}.bai" ]]; then
            samtools index "${DOWNSAMPLED_BAM}"
        fi

        echo "  Computing flagstat..."
        bash "${FLAGSTAT_SCRIPT}" "${DOWNSAMPLED_BAM}" "${FLAGSTAT_OUT}"

        echo "  Computing idxstats..."
        bash "${IDXSTATS_SCRIPT}" "${DOWNSAMPLED_BAM}" "${IDXSTATS_OUT}"

        # -----------------------------
        # Parse flagstat
        # -----------------------------
        FLAGSTAT_TOTAL=$(awk '/in total/ {print $1 + $3; exit}' "${FLAGSTAT_OUT}")
        FLAGSTAT_MAPPED=$(awk '/ mapped \(/ && !seen {print $1 + $3; seen=1}' "${FLAGSTAT_OUT}")

        # -----------------------------
        # Parse idxstats
        # total = mapped + unmapped across contigs
        # -----------------------------
        read -r IDXSTATS_TOTAL IDXSTATS_MAPPED < <(
            awk 'BEGIN{mapped=0; unmapped=0}
                 {mapped += $3; unmapped += $4}
                 END{print mapped + unmapped, mapped}' "${IDXSTATS_OUT}"
        )

        # -----------------------------
        # Compute differences
        # -----------------------------
        read -r MAPPED_DIFF MAPPED_DIFF_PCT < <(
            python3 - <<PY
flag_mapped = int("${FLAGSTAT_MAPPED}")
idx_mapped = int("${IDXSTATS_MAPPED}")
diff = flag_mapped - idx_mapped
pct = (diff / flag_mapped * 100) if flag_mapped != 0 else 0
print(f"{diff} {pct:.6f}")
PY
        )

        # -----------------------------
        # Append row
        # -----------------------------
        echo -e "${SAMPLE}\t${FRACTION}\t${DEPTH_PERCENT}\t${DOWNSAMPLED_BAM}\t${FLAGSTAT_TOTAL}\t${FLAGSTAT_MAPPED}\t${IDXSTATS_TOTAL}\t${IDXSTATS_MAPPED}\t${MAPPED_DIFF}\t${MAPPED_DIFF_PCT}" >> "${METRICS_TABLE}"
    done
done

# -----------------------------
# Build summary table
# -----------------------------
echo "Generating summary table..."

python3 - <<PY
import pandas as pd

metrics = pd.read_csv("${METRICS_TABLE}", sep="\t")

summary = (
    metrics.groupby("depth_percent")
    .agg(
        n_samples=("sample", "count"),
        mean_flagstat_mapped=("flagstat_mapped", "mean"),
        mean_idxstats_mapped=("idxstats_mapped", "mean"),
        mean_mapped_diff=("mapped_diff", "mean"),
        mean_mapped_diff_pct=("mapped_diff_pct", "mean"),
        sd_mapped_diff_pct=("mapped_diff_pct", "std"),
        min_mapped_diff_pct=("mapped_diff_pct", "min"),
        max_mapped_diff_pct=("mapped_diff_pct", "max"),
    )
    .reset_index()
)

summary.to_csv("${SUMMARY_TABLE}", sep="\t", index=False)
print("Summary written to ${SUMMARY_TABLE}")
PY

echo "Analysis complete."
echo "Metrics table: ${METRICS_TABLE}"
echo "Summary table: ${SUMMARY_TABLE}"