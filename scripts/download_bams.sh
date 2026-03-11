#!/usr/bin/env bash
set -euo pipefail

MANIFEST="docs/bam_manifest.tsv"
OUTDIR="data/raw_bams"

mkdir -p "${OUTDIR}"

if [[ ! -f "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found: ${MANIFEST}" >&2
    exit 1
fi

echo "Reading manifest: ${MANIFEST}"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r sample_id study_accession cancer_type assay_type source file_url md5 size_gb notes
do
    [[ -z "${sample_id}" ]] && continue
    [[ -z "${file_url}" ]] && continue

    outfile="${OUTDIR}/${sample_id}.bam"

    echo "----------------------------------------"
    echo "Sample: ${sample_id}"
    echo "Study: ${study_accession}"
    echo "Cancer: ${cancer_type}"
    echo "Assay: ${assay_type}"
    echo "Source: ${source}"
    echo "URL: ${file_url}"
    echo "Output: ${outfile}"

    if [[ -f "${outfile}" ]]; then
        echo "File already exists, skipping download: ${outfile}"
    else
        wget -c -O "${outfile}" "${file_url}"
    fi

    if [[ -n "${md5}" && "${md5}" != "NA" ]]; then
        echo "${md5}  ${outfile}" | md5sum -c -
    else
        echo "No MD5 provided; skipping checksum."
    fi
done

echo "Download step complete."