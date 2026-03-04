# Note pédagogique à ajouter au notebook après Question 1c

## Markdown cell à insérer

```markdown
---

### ⚠️ Note importante : Choix des covariables

Le tableau `eqtl` contient plusieurs colonnes supplémentaires : `BMI`, `LDL_cholesterol`, et `CAD_status`.

**Pourquoi ne pas les inclure dans le modèle de régression eQTL ?**

#### 1. `LDL_cholesterol` - Variable MÉDIATRICE ❌

**Chaîne causale** :
```
rs12916 (T allele) → HMGCR expression ↓ → LDL production ↓ → LDL sanguin ↓
```

- Le génotype **affecte** le taux de LDL via l'expression de HMGCR
- **Ajuster sur LDL casserait la chaîne causale** qu'on cherche à étudier !
- On perdrait l'information sur l'effet du génotype qui transite par HMGCR

**Analogie** : C'est comme vouloir mesurer l'effet d'un médicament sur la guérison, 
mais en ajustant sur le taux du médicament dans le sang. On efface l'effet qu'on veut observer !

---

#### 2. `CAD_status` - Variable OUTCOME (collider) ❌

**Chaîne causale** :
```
rs12916 → HMGCR expression → LDL ↓ → Risque CAD ↓
```

- Le statut CAD est une **conséquence** du génotype (via LDL)
- **Ajuster sur CAD créerait un biais de sélection** (collider bias)
- On sélectionnerait un sous-groupe non représentatif

**Exemple** : Parmi les personnes avec CAD, celles avec le génotype protecteur TT 
ont probablement d'autres facteurs de risque élevés. Ajuster sur CAD compare des "pommes et des oranges".

---

#### 3. `BMI` - Acceptable mais non nécessaire ✓

- Le BMI est **indépendant du génotype** (randomisation mendélienne)
- Pas de biais si on l'inclut, mais **pas nécessaire** dans un contexte génétique
- Les variants génétiques sont distribués aléatoirement à la conception

---

### 📝 Règle générale en analyse eQTL

**Variables à ajuster** : 
- Facteurs techniques (batch, plateforme)
- Structure de population (PC1, PC2, ...)
- Covariables démographiques non-médiées (âge, sexe)

**Variables à NE PAS ajuster** :
- Médiateurs (entre génotype et expression)
- Outcomes/conséquences (colliders)
- Variables post-exposition

---

### 🎓 Exercice bonus (optionnel)

**Question** : Que se passe-t-il si on ajuste sur `LDL_cholesterol` dans le modèle eQTL ?

```r
%%R
# Modèle INCORRECT (à ne pas faire!)
model_wrong <- lm(HMGCR_expression ~ genotype_rs12916 + age + sex + PC1 + PC2 + LDL_cholesterol, 
                  data = eqtl)
summary(model_wrong)$coefficients["genotype_rs12916", ]
```

**Comparez** l'effet du génotype avec et sans ajustement sur LDL. Que remarquez-vous ?

<details>
<summary>Cliquez pour voir la réponse</summary>

L'effet du génotype **diminue** ou **disparaît** ! 

**Pourquoi ?** On "explique" l'effet du génotype par le LDL, alors que le LDL 
est justement la conséquence de HMGCR. On efface l'effet causal qu'on cherche à détecter.

C'est un exemple classique de **sur-ajustement** (over-adjustment bias).

</details>

---
```

## Placement dans le notebook

Insérer cette note **après la Question 1c** (analyse eQTL avec régression), 
avant de passer à la Question 1d (F-statistic).

Cela permet de :
1. ✅ Consolider la compréhension du modèle eQTL
2. ✅ Introduire les concepts de médiateur/collider
3. ✅ Préparer la logique MR pour les parties II-IV

---

## Alternative (version plus courte)

Si tu veux une version plus concise :

```markdown
### ⚠️ Choix des covariables

**Variables dans le tableau** : `BMI`, `LDL_cholesterol`, `CAD_status`

**Pourquoi ne pas les ajuster ?**
- **`LDL_cholesterol`** : Médiateur (rs12916 → HMGCR → LDL). Ajuster casserait la chaîne causale.
- **`CAD_status`** : Outcome (collider). Ajuster créerait un biais de sélection.
- **`BMI`** : Acceptable mais non nécessaire (indépendant du génotype).

**Règle** : En analyse eQTL, ajuster uniquement sur les confondants (âge, sexe, PCs), 
pas sur les médiateurs ou outcomes.
```

---

Voilà ! Utilise la version qui te convient le mieux selon le niveau de détail souhaité.
