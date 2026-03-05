# Génération des données - Guide complet

Ce document explique comment reproduire les données du TP à partir des sources GWAS publiques.

---

## Vue d'ensemble

Le TP utilise deux types de données :

1. **Données eQTL** (simulées) : `data/eqtl_HMGCR.csv`
2. **Données MR** (réelles) : `data/MR_IVs_LDL_CAD.csv`

---

## Partie 1 : Données eQTL (simulées)

### Script de génération

**Fichier** : `scripts/00_generate_eQTL_data.py`

Ce script génère 250 individus avec :
- Génotypes rs12916 (MAF=0.42, équilibre Hardy-Weinberg)
- Expression HMGCR corrélée avec le génotype (β = -0.28/allèle T)
- Covariables : age, sex, PC1, PC2
- Variables pédagogiques : BMI, LDL_cholesterol, CAD_status

### Prérequis

```bash
# Créer environnement virtuel
python3 -m venv venv
source venv/bin/activate

# Installer dépendances
pip install numpy pandas
```

### Exécution

```bash
# Depuis la racine du repo
python scripts/00_generate_eQTL_data.py
```

**Sortie** : `data/eqtl_HMGCR.csv` (~10 KB)

### Paramètres de simulation

**Basé sur GTEx** (rs12916 cis-eQTL de HMGCR) :

- **Génotype rs12916** :
  - Allèle T (majeur) : fréquence = 0.60 (protecteur : diminue HMGCR et LDL)
  - Allèle C (mineur) : MAF = 0.40 (risque : augmente HMGCR et LDL) [source: dbSNP]
  - Équilibre Hardy-Weinberg
  - **Codage : nombre d'allèles C (mineur)** : 0=TT, 1=CT, 2=CC
  - Standard en génétique : coder l'allèle mineur
  
- **Expression HMGCR** :
  - Baseline : 8.5 (log2-normalized)
  - Effet par allèle C : β = +0.28 (augmente l'expression)
  - Effet par allèle T : β = -0.28 (diminue l'expression, effet protecteur)
  - Ajusté sur : age, sex, PC1, PC2
  - Résiduel SD : 0.45

- **Covariables** :
  - `age` : 20-70 ans (moyenne 45, SD 12)
  - `sex` : 0=F, 1=M (50/50)
  - `PC1` : Premier axe ACP (population structure)
  - `PC2` : Deuxième axe ACP (variance réduite, corrélation 0.15 avec PC1)

- **Variables pédagogiques** (à ne PAS ajuster dans le modèle eQTL) :
  - `BMI` : 18-40 kg/m² (indépendant du génotype)
  - `LDL_cholesterol` : 1.5-7.0 mmol/L (médiateur : C → HMGCR↑ → LDL↑)
  - `CAD_status` : 0/1 (outcome : C → LDL↑ → CAD↑)

### Vérifications

Le script affiche automatiquement :
- Équilibre Hardy-Weinberg (chi²)
- Effet génotype observé vs attendu
- Corrélations PC1-PC2
- Chaîne causale génotype → LDL → CAD
- Statistiques descriptives complètes

---

## Partie 2 : Données MR (réelles)

### Étape 1 : Téléchargement des GWAS

**Sources** :
- **LDL Cholesterol** : GCST90468080 (7.9M SNPs, build37)
- **Coronary Artery Disease** : GCST005194 (13.3M SNPs, build37)

**Script** : `scripts/01_download_data.sh` (documentation)

```bash
# LDL
wget -P data_raw/ https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468080/GCST90468080.tsv.gz

# CAD
wget -P data_raw/ https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005194/CAD_META.gz
mv data_raw/CAD_META.gz data_raw/CAD.tsv.gz
```

**Taille totale** : ~633 MB

### Étape 2 : Génération des IVs MR

**Script** : `scripts/03_generate_MR_data.py`

**Prérequis** :

```bash
# Installer genal-python
pip install genal-python pandas
```

**Exécution** :

```bash
python scripts/03_generate_MR_data.py
```

**Pipeline** :
1. Charge les GWAS (LDL et CAD)
2. Preprocessing avec genal (Fill_delete)
3. Clumping (p < 5e-8, r² < 0.01, window 10Mb)
4. Panel référence : EUR_37 (1000 Genomes, via API genal)
5. Harmonisation exposure/outcome
6. Export CSV

**Sortie** : `data/MR_IVs_LDL_CAD.csv` (280 SNPs, sans annotations)

### Étape 3 : Annotation des gènes

**Script** : `scripts/04_annotate_genes.R`

**Prérequis** :

```r
install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "rtracklayer"))
```

**Exécution** :

```bash
Rscript scripts/04_annotate_genes.R
```

**Pipeline** :
1. Télécharge GENCODE v19 GTF (1ère fois, 40 MB)
2. Charge 57,820 gènes
3. Étend fenêtres gènes ±500kb
4. FindOverlaps SNPs ↔ gènes
5. Priorise protein-coding genes
6. Choisit gène le plus proche si multiples

**Sortie** : `data/MR_IVs_LDL_CAD.csv` (280 SNPs annotés avec colonne GENE)

### Gènes clés identifiés

| Gène | SNPs | Fonction |
|------|------|----------|
| HMGCR | 1 | Cible des statines |
| PCSK9 | 2 | Régulation récepteur LDL |
| LDLR | 2 | Récepteur LDL |
| APOE | 1 | Apolipoprotéine E |
| APOB | 6 | Apolipoprotéine B |
| CETP | 2 | Transfert cholestérol |
| LPA | 4 | Lipoprotéine(a) |
| LIPC | 2 | Lipase hépatique |
| LIPG | 2 | Lipase endothéliale |

---

## Workflow complet (tout régénérer)

```bash
# 1. Configuration environnement
python3 -m venv venv
source venv/bin/activate
pip install numpy pandas genal-python

# 2. Données eQTL
python scripts/00_generate_eQTL_data.py

# 3. Télécharger GWAS (optionnel si déjà fait)
# bash scripts/01_download_data.sh

# 4. Générer IVs MR
python scripts/03_generate_MR_data.py

# 5. Annoter gènes (nécessite R + Bioconductor)
Rscript scripts/04_annotate_genes.R

# 6. Vérification
ls -lh data/
# Devrait afficher :
# eqtl_HMGCR.csv (~10 KB)
# MR_IVs_LDL_CAD.csv (~47 KB)
```

---

## Fichiers exclus du repo

**`.gitignore`** exclut :
- `data_raw/` : GWAS bruts (633 MB)
- `gencode.v19.annotation.gtf.gz` : Annotations (40 MB)
- `venv/` : Environnement Python local
- `tmp_GENAL/` : Fichiers temporaires genal

Ces fichiers sont auto-téléchargés ou générés par les scripts.

---

## Dépendances complètes

**Python** :
- numpy
- pandas
- genal-python

**R / Bioconductor** :
- GenomicRanges
- rtracklayer

**Système** :
- wget (téléchargement GWAS)

---

## Notes techniques

### Build génome

**Tous les fichiers utilisent GRCh37** :
- GWAS LDL : build37
- GWAS CAD : build37
- Panel 1000G : EUR_37
- GENCODE : v19 (GRCh37)

### Panel 1000 Genomes

Le panel EUR_37 est téléchargé automatiquement par l'API genal lors du clumping. 
Pas besoin de le télécharger manuellement.

### GENCODE

Le GTF est téléchargé automatiquement par `04_annotate_genes.R` la première fois.
Ensuite il est mis en cache localement.

---

## Troubleshooting

**Problème** : `ModuleNotFoundError: No module named 'numpy'`  
**Solution** : Activer le venv et installer les dépendances

**Problème** : `biomaRt` timeout  
**Solution** : Utiliser `04_annotate_genes.R` (GenomicRanges + GENCODE local)

**Problème** : Panel 1000G download failed  
**Solution** : genal gère automatiquement via API, rien à faire

---

**Dernière mise à jour** : 5 mars 2026
