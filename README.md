# Cancer Read Metrics: A Comparative Evaluation of `flagstat` vs `idxstats`

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

This project investigates discrepancies between `samtools flagstat` and `samtools idxstats` in evaluating mapped read counts in cancer BAM files. Our goal is to assess how these differences affect downstream analyses such as downsampling accuracy and variant allele frequency (VAF) detection in clinical cancer genomics.

---

## Motivation

Accurate estimation of mapped reads is essential for:

- Downsampling tumor/normal BAMs to matched coverage  
- Calibrating Variant Allele Frequency (VAF) thresholds  
- Ensuring consistency in cancer genomics pipelines  

This project explores:

- The divergence between `flagstat` and `idxstats` read counts  
- Whether discrepancies affect clinical metrics like VAF detection  
- Proposing a standardised approach or validation script  

---

## Folder Structure

```
cancer_read_metrics/
├── CHANGELOG.md
├── LICENSE
├── README.md
├── data/
│   └── sample_metadata.tsv
├── figures/
│   ├── vaf_shift_plot.png
│   └── flagstat_vs_idxstats.pdf
├── notebooks/
│   └── plot_vaf_metrics.ipynb
├── scripts/
│   ├── compute_flagstat.sh
│   ├── compute_idxstats.sh
│   └── compare_metrics.py
├── results/
│   └── metrics_comparison.tsv
└── tests/
    └── example.sam
```

---

## Tools and Dependencies

- `samtools` (v1.22 or above)
- `python` (3.10+)
- `pandas`, `matplotlib`, `seaborn` (for analysis and plotting)
- `jupyterlab` (for notebook-based exploration)

---

## How to Use

### 1. Compute Read Metrics

Run both tools on each BAM file:
```bash
bash scripts/compute_flagstat.sh input.bam > results/input_flagstat.txt
bash scripts/compute_idxstats.sh input.bam > results/input_idxstats.txt
```

### 2. Compare Outputs

```bash
python scripts/compare_metrics.py --flagstat results/input_flagstat.txt --idxstats results/input_idxstats.txt --out results/input_comparison.tsv
```

### 3. Visualise Differences

Use the provided Jupyter notebook to explore discrepancies:
```bash
jupyter notebook notebooks/plot_vaf_metrics.ipynb
```

---

## Research Questions

- Which tool produces higher/lower mapped read counts?
- Are differences consistent across tumour vs normal BAMs?
- Do discrepancies meaningfully shift VAF calculations in low-frequency variant calling?

---

## Preliminary Results

Early data suggest `flagstat` often reports **slightly higher mapped reads** than `idxstats`, particularly for files with soft-clipped regions or secondary alignments. This may lead to subtle downstream effects in:

- Coverage-matching during tumour/normal pairing  
- VAF thresholds in somatic variant detection  
- LoD (Limit of Detection) validation strategies  

---

## Future Directions

- Evaluate ≥100 public cancer BAMs (e.g., TCGA, ICGC)
- Cross-reference with downsampling ratios
- Publish as a short paper and release validation tools on GitHub

---

## Testing and Validation

A minimal SAM file (`tests/example.sam`) is included for development and testing. It allows quick validation of `samtools` parsing tools (`flagstat`, `idxstats`) without requiring large BAM files.

This file contains dummy reads aligned to chr1 and chr2:

```
@HDVN:1.0SO:unsorted
@SQSN:chr1LN:248956422
@SQSN:chr2LN:242193529
read0010chr11006010M*00ACTGACTGAAIIIIIIIIII
read0020chr21506010M*00TGACTGACTGIIIIIIIIII
read0034*00**00NNNNNNNNNNIIIIIIIIII
```

---

## Citation and Reuse

Feel free to reuse, cite, or adapt this project. Cite as:

> Olofinsae, S. (2025). *Cancer Read Metrics: Discrepancies in `flagstat` vs `idxstats` and Their Impact on VAF Detection.* GitHub. https://github.com/samson-olofinsae/cancer_read_metrics

---

## License

This project is licensed under the [MIT License](./LICENSE) © 2025 Samson Olofinsae.  
You are free to use, modify, and distribute this software with proper attribution.

---

## Training Context

This project supports my Clinical Bioinformatics training by evaluating discrepancies between mapped read metrics (`flagstat` vs `idxstats`) in cancer BAM files.

Its aims include:

- Validating downsampling workflows for tumour/normal pair matching
- Assessing impact on Variant Allele Frequency (VAF) thresholds
- Laying groundwork for a reproducible, citable validation tool

The work aligns with STP Equivalence Domains 1 (Clinical Care), 2 (Scientific Practice), and 4 (Research & Innovation).

---

## Ethical Statement
This project is intended strictly for educational and training purposes. It does not represent, replicate, or leak any internal pipeline, data, or intellectual property of the NHS, WMRGL, or affiliated institutions. All workflows and configurations are independently developed based on public documentation. No patient data was used. This project is released openly to promote transparency, validation, and reproducibility in clinical bioinformatics training.

## Environment & Reproducibility

All analyses in this project were performed using a reproducible Conda environment defined in [`environment.yml`](./environment.yml). Key software tools were version-pinned to ensure consistency across runs and reproducibility of results.

**Core dependencies include:**
- `samtools v1.22` – for all BAM-level metric extraction (`flagstat`, `idxstats`)
- `python v3.10` – for scripting and data processing
- `pandas`, `matplotlib`, `seaborn` – for downstream analysis and visualization
- `jupyterlab` – for interactive data exploration and reproducibility

This environment reflects stable versions commonly used in clinical and cancer genomics workflows.

To recreate the environment:
```bash
conda env create -f environment.yml
conda activate cancer-metrics-env
```



<!-- Test commit for contribution verification -->
