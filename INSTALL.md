# Installation et dépendances

## Environnement validé

- **Python** : 3.10+
- **R** : 4.3+
- **OS** : Ubuntu 22.04+ / macOS / Windows (WSL2)

## Dépendances Python
```bash
python3 -m venv venv
source venv/bin/activate  # Linux/Mac
# venv\Scripts\activate   # Windows

pip install genal-python pandas
```

### Installation automatique de PLINK2

genal installe automatiquement PLINK 2.0 au premier usage :
```python
import genal
genal.install_plink()  # Installation dans ~/.genal/plink2
```

Ceci est fait automatiquement par le script `scripts/03_generate_MR_data.py`.

## Dépendances R (pour notebooks)
```r
install.packages(c("ggplot2", "dplyr", "data.table"))
```

## Validation environnement
```bash
# Python
python3 -c "import genal; import pandas; print('Python OK')"

# R
R -e "library(ggplot2); library(dplyr); cat('R OK\n')"
```

## Fichiers téléchargés automatiquement

Au premier usage, genal télécharge :
- **Panel 1000G EUR** : `~/.genal/Reference_files/EUR_37/` (~970 Mo)
- **PLINK 2.0** : `~/.genal/plink2` (~10 Mo)

Ces fichiers sont réutilisés entre sessions (pas besoin de re-télécharger).

## Générer les données MR réelles

Voir `README_DATA_GENERATION.md` pour la procédure complète.
