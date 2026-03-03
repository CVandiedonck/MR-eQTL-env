# État actuel - 3 mars 2026

## Commit actuel
Dernière version stable (postBuild sans /etc/environment)

## Statut plateformes
- ✅ mybinder : fonctionne
- ❌ Plasma/Uracile : nécessite workaround R_HOME

## Workaround Plasma
Ajouter en 1ère cellule du notebook :
```python
import os
os.environ['R_HOME'] = '/srv/conda/envs/notebook/lib/R'
```

## Notes
- postBuild essayait d'écrire dans /etc/environment sans permission
- Revenu à version simple (sans tentative fix R_HOME)
