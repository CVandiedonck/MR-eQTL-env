# Practical: eQTL and Mendelian Randomization
## Master ST4Health — Human Genetics — 2025-2026

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YOUR_USERNAME/YOUR_REPO/main)

---

## Overview

This repository contains materials for a 2-hour practical on **Mendelian Randomization (MR)**,
preceded by a 15-20 min lecture introduction.

**Biological question:** Does LDL-cholesterol causally increase the risk of coronary artery disease (CAD)?

**Structure:**
- Part I  — eQTL as an instrumental variable (~20 min)
- Part II  — Exploring MR summary statistics data (~15 min)
- Part III — IVW by hand: Wald ratios, weighted mean, WLS regression (~25 min)
- Part IV  — Sensitivity analyses (~30 min)

---

## Two versions of the practical

| File | When | Dependencies |
|---|---|---|
| `practical_eQTL_MR_baseR.ipynb` | **March 2025** (current env) | base R + ggplot2 only |
| `practical_eQTL_MR_MRpackage.ipynb` | **2026 onwards** (new env) | `MendelianRandomization` package |

Both notebooks are identical for Parts I–III.  
Part IV differs: `baseR` implements MR-Egger and Weighted Median from scratch;  
`MRpackage` uses `mr_allmethods()`, `mr_plot()`, `mr_forest()`.

---

## Data

### Provided (simulated, ready to use)

| File | Description |
|---|---|
| `data/eqtl_HMGCR.csv` | 250 individuals, SNP rs12916, HMGCR expression (GTEx-like, Liver) |
| `data/MR_IVs_LDL_CAD.csv` | **PLACEHOLDER** — 77 simulated IVs (see below to replace with real data) |

> **Note on data:** The eQTL data are simulated from published GTEx v8 values  
> (rs12916, beta = -0.35, GTEx Consortium 2020).  
> The MR IVs are based on real rsIDs and published effect sizes  
> (Willer et al. *Nat Genet* 2013; Nikpay et al. *Nat Genet* 2015)  
> but were generated computationally, not downloaded from the original GWAS.

### Reproducibility: generate real data from scratch

To replace the placeholder CSV with real summary statistics:

```bash
# Step 1: download GWAS summary statistics from EBI
bash scripts/01_download_data.sh

# Step 2: filter, clump (PLINK 1.9 required) and harmonise
# Edit PLINK_BIN and REF_PANEL paths in the script header first
Rscript scripts/02_prepare_MR_IVs.R
```

**Requirements for Step 2:**
- PLINK 1.9 (`plink --version`)
- 1000 Genomes EUR reference panel, build GRCh37, PLINK format  
  Download: `wget https://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz`

---

## Environment

The `binder/` folder configures the computational environment via
[repo2docker](https://repo2docker.readthedocs.io/) (used by Plasma JupyterHub).

| File | Role |
|---|---|
| `binder/environment.yml` | conda/mamba: Python 3.10, R 4.3, rpy2, PLINK 1.9, genal-python |
| `binder/apt.txt` | system libraries (libcurl, cairo, etc.) |
| `binder/postBuild` | R packages from CRAN, Jupyter kernel registration |

**To rebuild the environment on Plasma:**
```bash
chmod +x binder/postBuild   # required by repo2docker
git add binder/postBuild
git commit -m "make postBuild executable"
git push
# Then trigger rebuild on Plasma admin interface
```

---

## Repository structure

```
.
├── binder/
│   ├── environment.yml
│   ├── apt.txt
│   └── postBuild
├── data/
│   ├── eqtl_HMGCR.csv
│   └── MR_IVs_LDL_CAD.csv        ← replace with real data (see scripts/)
├── scripts/
│   ├── 01_download_data.sh
│   └── 02_prepare_MR_IVs.R
├── practical_eQTL_MR_baseR.ipynb       ← March 2025, no external R package
├── practical_eQTL_MR_MRpackage.ipynb   ← 2026+, uses MendelianRandomization
└── README.md
```

---

## References

- Willer CJ et al. *Discovery and refinement of loci associated with lipid levels.*  
  **Nat Genet** 2013;45:1274–1283. (LDL GWAS — GCST90468080)
- Nikpay M et al. *A comprehensive 1000 Genomes-based genome-wide association meta-analysis of coronary artery disease.*  
  **Nat Genet** 2015;47:1121–1130. (CAD GWAS — GCST005194)
- Bowden J et al. *Mendelian randomization with invalid instruments.*  
  **Int J Epidemiol** 2015;44:512–525. (MR-Egger)
- Bowden J et al. *Consistent estimation in Mendelian randomization with some invalid instruments.*  
  **Genet Epidemiol** 2016;40:304–314. (Weighted Median)
- GTEx Consortium. *The GTEx Consortium atlas of genetic regulatory effects across human tissues.*  
  **Science** 2020;369:1318–1330. (eQTL data)
