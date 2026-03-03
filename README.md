# TP Randomisation Mendélienne et eQTL

[🇬🇧 English version](README.en.md) | 🇫🇷 Version française

Travaux pratiques de génétique statistique pour le Magistère Européen de Génétique (Université Paris Cité).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/CVandiedonck/MR-eQTL-env/main?urlpath=lab)

---

## 📚 Contenu du TP

**Durée** : 2h  
**Niveau** : M1  
**Format** : Jupyter notebook avec R (via rpy2)

### Partie I : eQTL comme variable instrumentale
- Concept de cis-eQTL
- Force de l'instrument (F-statistic)
- Exemple : SNP rs12916 et expression de HMGCR (gène cible des statines)

### Partie II : Exploration données MR
- 280 SNPs réels (GWAS LDL cholesterol et maladie coronarienne)
- Annotation des gènes (HMGCR, PCSK9, LDLR, APOE, APOB, CETP, LPA...)
- Visualisation des instruments

### Partie III : Inverse-Variance Weighted (IVW)
- Calcul manuel du ratio de Wald par SNP
- Forest plot
- Moyenne pondérée (IVW)
- Vérification via régression WLS

### Partie IV : Analyses de sensibilité
- MR-Egger (codé manuellement)
- Weighted Median
- Détection des outliers
- Comparaison des méthodes

**Approche pédagogique** : Toutes les méthodes MR sont codées à la main en R base pour une compréhension approfondie (pas de package MendelianRandomization).

---

## 🚀 Utilisation

### mybinder (recommandé)

Cliquez sur le badge ci-dessus ou utilisez ce lien :
```
https://mybinder.org/v2/gh/CVandiedonck/MR-eQTL-env/main?urlpath=lab
```

L'environnement se construit automatiquement (~5 min la première fois).

### Installation locale

```bash
# Clone le repo
git clone https://github.com/CVandiedonck/MR-eQTL-env.git
cd MR-eQTL-env

# Crée l'environnement conda
conda env create -f binder/environment.yml
conda activate mr-env

# Installe rpy2
pip install rpy2==3.4.5

# Lance Jupyter
jupyter lab
```

Ouvrez `practical_eQTL_MR_baseR.ipynb`.

---

## 📊 Données

### Structure du dépôt

```
MR-eQTL-env/
├── binder/                          # Configuration environnement
│   ├── environment.yml              # Packages conda/pip
│   ├── postBuild                    # Installation rpy2
│   └── apt.txt                      # Dépendances système
├── data/                            # Données du TP
│   ├── eqtl_HMGCR.csv              # Données eQTL simulées
│   └── MR_IVs_LDL_CAD.csv          # 280 SNPs MR annotés
├── data_raw/                        # GWAS bruts (non versionnés)
│   ├── LDL.tsv.gz                  # 372M - GWAS LDL
│   ├── CAD.tsv.gz                  # 261M - GWAS CAD
│   ├── LDL_md5sum.txt              # Checksum
│   └── [fichier_1KG.zip]           # Panel 1000 Genomes (optionnel)
├── scripts/                         # Scripts de génération
│   ├── 01_download_data.sh         # Téléchargement GWAS
│   ├── 02_prepare_MR_IVs.R         # Pipeline R/PLINK (doc)
│   ├── 03_generate_MR_data.py      # Génération IVs (genal)
│   ├── 04_annotate_genes.R         # Annotation GENCODE
│   └── tmp_GENAL/                  # Fichiers temp (non versionnés)
├── practical_eQTL_MR_baseR.ipynb   # Notebook principal (R base)
├── practical_eQTL_MR_tidyverse.ipynb # Version tidyverse (archive)
├── practical_eQTL_MR_MRpackage.ipynb # Version package MR (archive)
├── README.md                        # Ce fichier (français)
├── README.en.md                     # Version anglaise
├── README_DATA_GENERATION.md        # Pipeline reproduction données
├── INSTALL.md                       # Installation locale
├── ETAT.md                          # État développement
├── CREDITS.md                       # Remerciements
└── LICENSE                          # CC BY 4.0

Fichiers non versionnés (.gitignore) :
├── data_raw/                        # GWAS bruts (372M + 261M)
├── gencode.v19.annotation.gtf.gz    # Annotations (40M, auto-download)
├── venv/                            # Environnement Python local
├── tmp_GENAL/                       # Fichiers temporaires genal (racine)
└── scripts/tmp_GENAL/               # Fichiers temporaires genal (scripts)
```

### Données eQTL (simulées)
- **Fichier** : `data/eqtl_HMGCR.csv`
- **Source** : Basé sur GTEx v8 (Liver, rs12916, HMGCR)
- 250 individus simulés avec génotypes, expression, covariables

### Données MR (réelles)
- **Fichier** : `data/MR_IVs_LDL_CAD.csv`
- **Sources GWAS** :
  - LDL cholesterol : GCST90468080 (7.9M SNPs, build37)
  - Maladie coronarienne : GCST005194 (13.3M SNPs, build37)
- **Pipeline** : Clumping via `genal-python` (p<5e-8, r²<0.01, 10Mb)
- **Résultat** : 280 SNPs indépendants annotés (GENCODE v19)

### Gènes clés identifiés
HMGCR, PCSK9, LDLR, APOE, APOB, CETP, LPA, LIPC, LIPG

---

## 🛠 Reproduction des données

Voir `README_DATA_GENERATION.md` pour le pipeline complet de génération des données MR.

**Scripts disponibles** :
- `scripts/01_download_data.sh` : Téléchargement GWAS
- `scripts/03_generate_MR_data.py` : Clumping et harmonisation
- `scripts/04_annotate_genes.R` : Annotation des gènes (GenomicRanges + GENCODE)

**Note** : Le script `03_generate_MR_data.py` utilise des chemins relatifs (`../data_raw/`, `../data/`) et doit être exécuté depuis le dossier `scripts/`.

---

## 📝 Licence

Ce projet est sous licence [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE).

Vous êtes libre de :
- Partager et adapter ce matériel
- À condition de créditer l'auteur original

---

## 👥 Crédits

**Auteur** : Claire Vandiedonck (Université Paris Cité)  
**GitHub** : [@CVandiedonck](https://github.com/CVandiedonck)  
**Email** : claire.vandiedonck@u-paris.fr  
**Année** : 2026

Voir [CREDITS.md](CREDITS.md) pour les remerciements détaillés.

---

## 🔗 Liens utiles

- **TP Python MR (Marie Verbanck)** : https://github.com/CVandiedonck/MR_course_Python_Verbanck
- **Documentation genal** : https://genal.readthedocs.io/
- **GENCODE** : https://www.gencodegenes.org/
- **GTEx Portal** : https://gtexportal.org/
- **GWAS Catalog** : https://www.ebi.ac.uk/gwas/

---

## 📧 Contact

**Claire Vandiedonck**  
claire.vandiedonck@u-paris.fr  
Université Paris Cité

---

**Bon TP !** 🧬
