# Génération des données MR réelles

## Données sources

- **LDL** : GCST90468080 (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468080/)
- **CAD** : GCST005194 (https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005194/)

## Procédure de génération

### Prérequis
```bash
python3 -m venv venv
source venv/bin/activate
pip install genal-python pandas
```

### Téléchargement GWAS
```bash
mkdir -p data_raw
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468080/GCST90468080.tsv.gz \
     -O data_raw/LDL.tsv.gz
wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005194/CAD_META.gz \
     -O data_raw/CAD.tsv.gz
```

### Génération des IVs
```bash
python3 generate_MR_data.py
```

**Sortie** : `data/MR_IVs_LDL_CAD.csv` (~280 SNPs)

### Résultat

- 280 instruments génétiques indépendants
- F-stat médiane : ~80
- Build : GRCh37
- Panel de référence : 1000G EUR (via genal)

## Fichier généré le

28 février 2026 - 280 SNPs
