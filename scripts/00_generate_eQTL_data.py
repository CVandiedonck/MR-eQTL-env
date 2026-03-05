#!/usr/bin/env python3
"""
Generate simulated eQTL data for HMGCR gene
Based on GTEx (rs12916 cis-eQTL of HMGCR)

Includes additional covariates for pedagogical purposes:
- BMI (confounding variable)
- LDL_cholesterol (mediator variable - NOT to adjust for in eQTL model)
- CAD_status (outcome variable - NOT to adjust for in eQTL model)

Allele frequencies (from dbSNP): T (major, freq=0.60), C (minor, MAF=0.40)
Genotypes coded as number of MINOR allele C (0=TT, 1=CT, 2=CC)
Effect: C allele increases HMGCR expression (β=+0.28 per C allele)
       T allele (protective) decreases HMGCR and LDL

Output: data/eqtl_HMGCR.csv
"""

import numpy as np
import pandas as pd
from pathlib import Path

# Set seed for reproducibility
np.random.seed(42)

# Parameters
N = 250  # Number of individuals
FREQ_C = 0.40  # Frequency of C allele (minor allele, MAF from dbSNP)
FREQ_T = 0.60  # Frequency of T allele (major allele)

print("=" * 60)
print("Generating eQTL data for HMGCR")
print("=" * 60)

# 1. Generate genotypes under Hardy-Weinberg equilibrium
print("\n1. Generating genotypes (rs12916, T=major, C=minor, MAF=0.40 from dbSNP)...")

# HWE frequencies (genotypes coded as number of T alleles: 0, 1, 2)
p_TT = FREQ_T ** 2  # 2 copies of T (major)
p_CT = 2 * FREQ_T * FREQ_C  # 1 copy of T
p_CC = FREQ_C ** 2  # 0 copies of T (2 copies of C, minor)

# Generate genotypes (0, 1, 2 copies of MINOR allele C)
# Standard practice: code minor allele
genotypes = np.random.choice(
    [0, 1, 2],  # 0=TT, 1=CT, 2=CC (number of C alleles)
    size=N, 
    p=[p_TT, p_CT, p_CC]  # IMPORTANT: order changed for minor allele coding
)

# Verify HWE
obs_counts = np.bincount(genotypes, minlength=3)
print(f"   Genotype counts: TT (0C)={obs_counts[0]}, CT (1C)={obs_counts[1]}, CC (2C)={obs_counts[2]}")
print(f"   Observed MAF(C): {(obs_counts[1] + 2*obs_counts[2]) / (2*N):.3f} (expected: {FREQ_C:.3f})")
print(f"   Observed freq(T): {(obs_counts[1] + 2*obs_counts[0]) / (2*N):.3f} (expected: {FREQ_T:.3f})")

# 2. Generate basic covariates
print("\n2. Generating basic covariates...")

# Age (20-70 years, normal distribution)
age = np.random.normal(45, 12, N).clip(20, 70)

# Sex (0=F, 1=M, ~50/50)
sex = np.random.binomial(1, 0.5, N)

# PC1 (first principal component of genotype matrix)
# Simulates population structure, mean=0, sd=1
PC1 = np.random.normal(0, 1, N)

# PC2 (second principal component)
# Lower variance than PC1, slight correlation with PC1
PC2_base = np.random.normal(0, 0.6, N)  # Reduced variance
PC2 = PC2_base + 0.15 * PC1  # Slight correlation with PC1
PC2 = (PC2 - PC2.mean()) / PC2.std()  # Standardize

print(f"   Age: mean={age.mean():.1f}, sd={age.std():.1f}")
print(f"   Sex: {sex.sum()} males, {N-sex.sum()} females")
print(f"   PC1: mean={PC1.mean():.3f}, sd={PC1.std():.3f}")
print(f"   PC2: mean={PC2.mean():.3f}, sd={PC2.std():.3f}")
print(f"   Correlation PC1-PC2: {np.corrcoef(PC1, PC2)[0,1]:.3f}")

# 3. Generate HMGCR expression
print("\n3. Generating HMGCR expression...")

# Parameters based on GTEx (rs12916 cis-eQTL of HMGCR)
# T allele (major) decreases HMGCR expression
# We code number of C alleles (minor), so C increases HMGCR

baseline_expression = 8.5  # Mean log2-normalized expression
beta_genotype = 0.28  # Effect size per C allele (positive: C increases expression)
beta_age = 0.008  # Small age effect
beta_sex = 0.15  # Males slightly higher
beta_PC1 = 0.05  # Small population structure effect
beta_PC2 = 0.02  # Very small PC2 effect
residual_sd = 0.45  # Residual variation

# Generate expression
HMGCR_expression = (
    baseline_expression +
    beta_genotype * genotypes +
    beta_age * (age - age.mean()) +
    beta_sex * sex +
    beta_PC1 * PC1 +
    beta_PC2 * PC2 +
    np.random.normal(0, residual_sd, N)
)

print(f"   Mean expression: {HMGCR_expression.mean():.3f}")
print(f"   SD expression: {HMGCR_expression.std():.3f}")

# Verify genotype effect
expr_by_geno = [HMGCR_expression[genotypes == g].mean() for g in [0, 1, 2]]
print(f"   Expression by genotype:")
print(f"      TT (0C): {expr_by_geno[0]:.3f}")
print(f"      CT (1C): {expr_by_geno[1]:.3f}")
print(f"      CC (2C): {expr_by_geno[2]:.3f}")
print(f"   Effect per C allele: {(expr_by_geno[2] - expr_by_geno[0]) / 2:.3f} (should be ~+0.28)")

# 4. Generate additional covariates (for pedagogical purposes)
print("\n4. Generating additional covariates (for pedagogical purposes)...")
print("   ⚠️  These should NOT be included in the eQTL regression model!")

# BMI (confounding variable, independent of genotype)
BMI = np.random.normal(25, 4, N).clip(18, 40)

# LDL cholesterol (MEDIATOR: genotype → HMGCR expression → LDL)
# C allele increases HMGCR → increases LDL production
# Mean LDL ~3.5 mmol/L, effect ~+0.15 mmol/L per C allele
LDL_baseline = 3.5
LDL_beta_genotype = 0.15  # Per C allele (positive: C raises LDL)
LDL_beta_age = 0.01  # Age effect
LDL_beta_sex = 0.2  # Males slightly higher
LDL_beta_BMI = 0.05  # BMI effect

LDL_cholesterol = (
    LDL_baseline +
    LDL_beta_genotype * genotypes +
    LDL_beta_age * (age - age.mean()) +
    LDL_beta_sex * sex +
    LDL_beta_BMI * (BMI - BMI.mean()) +
    np.random.normal(0, 0.6, N)
).clip(1.5, 7.0)

# CAD status (OUTCOME: genotype → LDL → CAD)
# Lower LDL = lower CAD risk
# Logistic model: P(CAD) depends on LDL, age, sex, BMI
logit_CAD = (
    -3.5 +  # Baseline low prevalence
    0.4 * (LDL_cholesterol - 3.5) +  # LDL increases risk
    0.03 * (age - 45) +  # Age increases risk
    0.3 * sex +  # Males higher risk
    0.05 * (BMI - 25)  # BMI increases risk
)
prob_CAD = 1 / (1 + np.exp(-logit_CAD))
CAD_status = np.random.binomial(1, prob_CAD)

print(f"   BMI: mean={BMI.mean():.1f}, sd={BMI.std():.1f}")
print(f"   LDL cholesterol: mean={LDL_cholesterol.mean():.2f} mmol/L, sd={LDL_cholesterol.std():.2f}")
print(f"   CAD prevalence: {CAD_status.sum()}/{N} ({100*CAD_status.mean():.1f}%)")

# Verify causal relationships
print(f"\n   Causal chain verification:")
print(f"   - LDL by genotype: TT (0C)={LDL_cholesterol[genotypes==0].mean():.2f}, "
      f"CT (1C)={LDL_cholesterol[genotypes==1].mean():.2f}, "
      f"CC (2C)={LDL_cholesterol[genotypes==2].mean():.2f}")
print(f"     → C allele increases LDL (via increased HMGCR)")
print(f"   - CAD by genotype: TT (0C)={100*CAD_status[genotypes==0].mean():.1f}%, "
      f"CT (1C)={100*CAD_status[genotypes==1].mean():.1f}%, "
      f"CC (2C)={100*CAD_status[genotypes==2].mean():.1f}%")
print(f"     → C allele increases CAD risk (via LDL)")

# 5. Create DataFrame
print("\n5. Creating DataFrame...")

eqtl_data = pd.DataFrame({
    'individual_id': range(1, N + 1),
    'genotype_rs12916': genotypes,  # Number of minor allele C (0=TT, 1=CT, 2=CC)
    'HMGCR_expression': HMGCR_expression,
    'age': age,
    'sex': sex,
    'PC1': PC1,
    'PC2': PC2,
    'BMI': BMI,
    'LDL_cholesterol': LDL_cholesterol,
    'CAD_status': CAD_status
})

print(eqtl_data.head())
print(f"\nShape: {eqtl_data.shape}")

# 6. Save to file
print("\n6. Saving to data/eqtl_HMGCR.csv...")

output_dir = Path("data")
output_dir.mkdir(exist_ok=True)

output_file = output_dir / "eqtl_HMGCR.csv"
eqtl_data.to_csv(output_file, index=False)

print(f"✓ Saved to {output_file}")
print(f"  File size: {output_file.stat().st_size / 1024:.1f} KB")

# 7. Summary statistics
print("\n" + "=" * 60)
print("Summary Statistics")
print("=" * 60)

print("\nGenotype distribution:")
print(eqtl_data['genotype_rs12916'].value_counts().sort_index())

print("\nExpression statistics:")
print(eqtl_data.groupby('genotype_rs12916')['HMGCR_expression'].agg(['count', 'mean', 'std']))

print("\nLDL cholesterol by genotype:")
print(eqtl_data.groupby('genotype_rs12916')['LDL_cholesterol'].agg(['mean', 'std']))

print("\nCAD prevalence by genotype:")
cad_by_geno = eqtl_data.groupby('genotype_rs12916')['CAD_status'].agg(['sum', 'count', 'mean'])
cad_by_geno['percentage'] = 100 * cad_by_geno['mean']
print(cad_by_geno[['sum', 'count', 'percentage']])

print("\nCorrelation matrix (key variables):")
corr_vars = ['genotype_rs12916', 'HMGCR_expression', 'LDL_cholesterol', 'CAD_status', 
             'age', 'sex', 'BMI', 'PC1', 'PC2']
corr = eqtl_data[corr_vars].corr()
print(corr.round(3))

print("\n" + "=" * 60)
print("Pedagogical note:")
print("=" * 60)
print("""
The dataset includes variables that should NOT be adjusted for in the eQTL model:

1. LDL_cholesterol: MEDIATOR variable
   - Causal chain: rs12916 → HMGCR expression → LDL production
   - Adjusting for LDL would break the causal chain we want to study!

2. CAD_status: OUTCOME variable (collider)
   - Causal chain: rs12916 → LDL → CAD
   - Adjusting for CAD creates selection bias

3. BMI: Acceptable but unnecessary
   - Independent of genotype (Mendelian randomization assumption)
   - Not a confounder in this genetic context

Valid covariates for eQTL model: age, sex, PC1, PC2
""")

print("\n" + "=" * 60)
print("✓ eQTL data generation complete!")
print("=" * 60)
