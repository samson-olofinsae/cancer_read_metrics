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
│   └── analysis_summary.ipynb
├── scripts/
│   ├── compute_flagstat.sh
│   ├── compute_idxstats.sh
│   └── compare_metrics.py
└── results/
    └── metrics_comparison.tsv
```

---

## Tools and Dependencies

- `samtools` (v1.10 or above)
- `python` (3.8+)
- `pandas`, `matplotlib`, `seaborn` (for analysis and plotting)
- `jupyter` (for notebook-based exploration)

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
jupyter notebook notebooks/analysis_summary.ipynb
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

## Citation and Reuse

Feel free to reuse, cite, or adapt this project. Cite as:

> Olofinsae, S. (2025). *Cancer Read Metrics: Discrepancies in `flagstat` vs `idxstats` and Their Impact on VAF Detection.* GitHub. https://github.com/samson-olofinsae/cancer_read_metrics

---

## License

This project is licensed under the [MIT License](./LICENSE).


---

## 📄 License

This project is licensed under the [MIT License](./LICENSE) © 2025 Samson Olofinsae.  
You are free to use, modify, and distribute this software with proper attribution.

