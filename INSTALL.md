# Installation et Configuration

Guide d'installation des dépendances pour reproduire les données et exécuter le TP.

---

## Pour utiliser le TP (utilisateurs)

**Aucune installation nécessaire** si vous utilisez mybinder ou l'environnement Adenine configuré.

**URL mybinder** : https://mybinder.org/v2/gh/CVandiedonck/MR-eQTL-env/main?urlpath=lab

---

## Pour reproduire les données (développeurs)

### 1. Environnement Python

#### Créer un environnement virtuel

```bash
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# ou
venv\Scripts\activate  # Windows
```

#### Installer les dépendances Python

```bash
pip install numpy pandas genal-python
```

**Packages requis** :
- `numpy` : Calculs numériques (génération données eQTL)
- `pandas` : Manipulation de données (tous les scripts)
- `genal-python` : Clumping et harmonisation MR

### 2. Environnement R

#### Packages R de base

```r
install.packages(c("ggplot2", "dplyr", "data.table"))
```

#### Packages Bioconductor (pour annotation gènes)

```r
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

### 3. Outils système

**Linux/Mac** :
```bash
# wget pour télécharger les GWAS
sudo apt-get install wget  # Debian/Ubuntu
# ou
brew install wget  # macOS
```

**Windows** :
- Installer wget via Chocolatey ou utiliser PowerShell `Invoke-WebRequest`

---

## Génération des données

### Données eQTL (simulées)

```bash
# Depuis la racine du repo
python scripts/00_generate_eQTL_data.py
```

**Sortie** : `data/eqtl_HMGCR.csv`

**Dépendances** : numpy, pandas

### Données MR (réelles)

#### Étape 1 : Télécharger les GWAS (optionnel)

```bash
bash scripts/01_download_data.sh
```

**Taille** : ~633 MB  
**Dépendances** : wget

#### Étape 2 : Générer les IVs MR

```bash
python scripts/03_generate_MR_data.py
```

**Sortie** : `data/MR_IVs_LDL_CAD.csv` (280 SNPs, sans annotations)

**Dépendances** : genal-python, pandas

#### Étape 3 : Annoter les gènes

```bash
Rscript scripts/04_annotate_genes.R
```

**Sortie** : `data/MR_IVs_LDL_CAD.csv` (280 SNPs annotés avec gènes)

**Dépendances** : GenomicRanges, rtracklayer

---

## Configuration mybinder

**Pour tester l'environnement mybinder en local** :

```bash
# Installer repo2docker
pip install jupyter-repo2docker

# Builder l'image
jupyter-repo2docker --no-run .

# Ou directement exécuter
jupyter-repo2docker .
```

---

## Dépendances complètes

### Python (données + notebook)

**Requis pour génération données** :
- numpy >= 1.24
- pandas >= 2.0
- genal-python >= 1.4

**Requis pour notebook** :
- rpy2 == 3.4.5 (installé via postBuild sur mybinder)

**Optionnel** :
- jupyter-repo2docker (pour tester mybinder localement)

### R (analyses + annotation)

**Requis pour notebook** :
- R >= 4.3
- ggplot2
- (optionnel : dplyr, data.table si version tidyverse)

**Requis pour annotation gènes** :
- GenomicRanges (Bioconductor)
- rtracklayer (Bioconductor)

### Système

- wget (téléchargement GWAS)
- git (clonage repo)

---

## Installation complète (exemple Ubuntu/Debian)

```bash
# 1. Système
sudo apt-get update
sudo apt-get install -y python3 python3-venv python3-pip r-base wget git

# 2. Clone repo
git clone https://github.com/CVandiedonck/MR-eQTL-env.git
cd MR-eQTL-env

# 3. Python
python3 -m venv venv
source venv/bin/activate
pip install numpy pandas genal-python

# 4. R packages
R -e "install.packages(c('ggplot2', 'BiocManager'), repos='https://cloud.r-project.org')"
R -e "BiocManager::install(c('GenomicRanges', 'rtracklayer'))"

# 5. Générer données
python scripts/00_generate_eQTL_data.py
python scripts/03_generate_MR_data.py
Rscript scripts/04_annotate_genes.R

# 6. Vérification
ls -lh data/
# Devrait afficher eqtl_HMGCR.csv et MR_IVs_LDL_CAD.csv
```

---

## Versions testées

**Environnement de développement** :
- Python 3.10
- R 4.3.3
- numpy 1.26
- pandas 2.2
- genal-python 1.4.3
- rpy2 3.4.5

**Environnement mybinder** :
- Python 3.10
- R 4.3
- rpy2 3.4.5 (via postBuild)
- Packages R via conda

---

## Troubleshooting

### `ModuleNotFoundError: No module named 'numpy'`

**Solution** : Activer le venv et installer les dépendances
```bash
source venv/bin/activate
pip install numpy pandas
```

### `rpy2.rinterface_lib.openrlib.ValueError: r_home is None`

**Solution** : Sur Plasma, définir R_HOME manuellement (voir PLASMA_TECHNICAL.md)

### Packages R manquants

**Solution** :
```r
install.packages("nom_du_package", repos="https://cloud.r-project.org")
```

### GENCODE download timeout

**Solution** : Le script réessaie automatiquement. Si échec répété, télécharger manuellement :
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
```

---

## Fichiers exclus du repo

Les fichiers suivants sont générés localement et exclus du repo (`.gitignore`) :

- `data_raw/` : GWAS bruts (633 MB)
- `gencode.v19.annotation.gtf.gz` : Annotations (40 MB)
- `venv/` : Environnement Python
- `tmp_GENAL/` : Fichiers temporaires

---

**Dernière mise à jour** : 5 mars 2026
