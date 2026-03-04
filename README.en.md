# Mendelian Randomization and eQTL Practical

🇬🇧 English version | [🇫🇷 Version française](README.md)

Practical course in statistical genetics for the European Master's in Genetics (Université Paris Cité).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/CVandiedonck/MR-eQTL-env/main?urlpath=lab)

---

## 📚 Course Content

**Duration**: 2h  
**Level**: Master's (M1)  
**Format**: Jupyter notebook with R (via rpy2)

### Part I: eQTL as Instrumental Variable
- cis-eQTL concept
- Instrument strength (F-statistic)
- Example: SNP rs12916 and HMGCR expression (statin target gene)

### Part II: Exploring MR Data
- 280 real SNPs (LDL cholesterol and coronary artery disease GWAS)
- Gene annotation (HMGCR, PCSK9, LDLR, APOE, APOB, CETP, LPA...)
- Instrument visualization

### Part III: Inverse-Variance Weighted (IVW)
- Manual calculation of Wald ratio per SNP
- Forest plot
- Weighted average (IVW)
- Verification via WLS regression

### Part IV: Sensitivity Analyses
- MR-Egger (manually coded)
- Weighted Median
- Outlier detection
- Method comparison

**Pedagogical approach**: All MR methods are manually coded in base R for deep understanding (no MendelianRandomization package).

---

## 🚀 Usage

### mybinder (recommended)

Click the badge above or use this link:
```
https://mybinder.org/v2/gh/CVandiedonck/MR-eQTL-env/main?urlpath=lab
```

The environment builds automatically (~5 min first time).

### Local installation

```bash
# Clone the repo
git clone https://github.com/CVandiedonck/MR-eQTL-env.git
cd MR-eQTL-env

# Create conda environment
conda env create -f binder/environment.yml
conda activate mr-env

# Install rpy2
pip install rpy2==3.4.5

# Launch Jupyter
jupyter lab
```

Open `practical_eQTL_MR_baseR.ipynb`.

---

## 📊 Data

### Repository Structure

```
MR-eQTL-env/
├── binder/                          # Environment configuration
│   ├── environment.yml              # Conda/pip packages
│   ├── postBuild                    # rpy2 installation
│   └── apt.txt                      # System dependencies
├── data/                            # Course data
│   ├── eqtl_HMGCR.csv              # Simulated eQTL data
│   └── MR_IVs_LDL_CAD.csv          # 280 annotated MR SNPs
├── data_raw/                        # Raw GWAS (not versioned)
│   ├── LDL.tsv.gz                  # 372M - LDL GWAS
│   ├── CAD.tsv.gz                  # 261M - CAD GWAS
│   ├── LDL_md5sum.txt              # Checksum
│   └── [1KG_file.zip]              # 1000 Genomes panel (optional)
├── scripts/                         # Generation scripts
│   ├── 00_generate_eQTL_data.py    # eQTL data generation (simulated)
│   ├── 01_download_data.sh         # GWAS download
│   ├── 02_prepare_MR_IVs.R         # R/PLINK pipeline (doc)
│   ├── 03_generate_MR_data.py      # IV generation (genal)
│   ├── 04_annotate_genes.R         # GENCODE annotation
│   └── tmp_GENAL/                  # Temp files (not versioned)
├── practical_eQTL_MR_baseR.ipynb   # Main notebook (base R)
├── practical_eQTL_MR_tidyverse.ipynb # Tidyverse version (archive)
├── practical_eQTL_MR_MRpackage.ipynb # MR package version (archive)
├── README.md                        # French version
├── README.en.md                     # This file (English)
├── README_DATA_GENERATION.md        # Data reproduction pipeline
├── INSTALL.md                       # Local installation
├── ETAT.md                          # Development status
├── CREDITS.md                       # Acknowledgments
└── LICENSE                          # CC BY 4.0

Files not versioned (.gitignore):
├── data_raw/                        # Raw GWAS (372M + 261M)
├── gencode.v19.annotation.gtf.gz    # Annotations (40M, auto-download)
├── venv/                            # Local Python environment
├── tmp_GENAL/                       # Genal temp files (root)
└── scripts/tmp_GENAL/               # Genal temp files (scripts)
```

### eQTL Data (simulated)
- **File**: `data/eqtl_HMGCR.csv`
- **Script**: `scripts/00_generate_eQTL_data.py`
- **Source**: Based on GTEx v8 (Liver, rs12916, HMGCR)
- **Content**: 250 individuals with:
  - rs12916 genotypes (MAF=0.42, HW equilibrium)
  - HMGCR expression (effect β=-0.28/T allele)
  - Covariates: age, sex, PC1, PC2
  - Pedagogical variables: BMI, LDL_cholesterol, CAD_status

### MR Data (real)
- **File**: `data/MR_IVs_LDL_CAD.csv`
- **GWAS sources**:
  - LDL cholesterol: GCST90468080 (7.9M SNPs, build37)
  - Coronary artery disease: GCST005194 (13.3M SNPs, build37)
- **Pipeline**: Clumping via `genal-python` (p<5e-8, r²<0.01, 10Mb)
- **Result**: 280 independent SNPs annotated (GENCODE v19)

### Key Genes Identified
HMGCR, PCSK9, LDLR, APOE, APOB, CETP, LPA, LIPC, LIPG

---

## 🛠 Data Reproduction

See `README_DATA_GENERATION.md` for the complete MR data generation pipeline.

**Available scripts**:
- `scripts/00_generate_eQTL_data.py`: Simulated eQTL data generation
- `scripts/01_download_data.sh`: GWAS download
- `scripts/03_generate_MR_data.py`: Clumping and harmonization
- `scripts/04_annotate_genes.R`: Gene annotation (GenomicRanges + GENCODE)

**Note**: The script `03_generate_MR_data.py` uses relative paths (`../data_raw/`, `../data/`) and must be run from the `scripts/` directory.

---

## 📝 License

This project is licensed under [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE).

You are free to:
- Share and adapt this material
- Provided you credit the original author

---

## 👥 Credits

**Author**: Claire Vandiedonck (Université Paris Cité)  
**GitHub**: [@CVandiedonck](https://github.com/CVandiedonck)  
**Email**: claire.vandiedonck@u-paris.fr  
**Year**: 2026

See [CREDITS.md](CREDITS.md) for detailed acknowledgments.

---

## 🔗 Useful Links

- **Python MR Course (Marie Verbanck)**: https://github.com/CVandiedonck/MR_course_Python_Verbanck
- **genal documentation**: https://genal.readthedocs.io/
- **GENCODE**: https://www.gencodegenes.org/
- **GTEx Portal**: https://gtexportal.org/
- **GWAS Catalog**: https://www.ebi.ac.uk/gwas/

---

## 📧 Contact

**Claire Vandiedonck**  
claire.vandiedonck@u-paris.fr  
Université Paris Cité

---

**Enjoy the practical!** 🧬
