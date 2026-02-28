#!/usr/bin/env python3
"""
Generate real MR IVs (LDL → CAD) from GWAS summary statistics
Using genal package (same as Marie's notebook)
"""

import pandas as pd
import genal

print("=" * 60)
print("Génération des IVs MR réels (LDL → CAD)")
print("=" * 60)

# 1. Load LDL GWAS
print("\n[1/5] Chargement LDL GWAS...")
ldl = pd.read_csv("../data_raw/LDL.tsv.gz", sep="\t", compression='gzip')
print(f"  → {len(ldl):,} SNPs chargés")

E_Geno = genal.Geno(
    ldl, 
    CHR="chromosome", 
    POS="base_pair_location", 
    EA="effect_allele", 
    NEA="other_allele", 
    BETA="beta", 
    SE="standard_error", 
    P="p_value", 
    EAF="effect_allele_frequency",
    keep_columns=False
)

print("[2/5] Preprocessing LDL...")
E_Geno.preprocess_data(preprocessing='Fill_delete', reference_panel="EUR_37")

# 2. Load CAD GWAS
print("[3/5] Chargement CAD GWAS...")
cad = pd.read_csv("../data_raw/CAD.tsv.gz", sep="\t", compression='gzip')
print(f"  → {len(cad):,} SNPs chargés")

O_Geno = genal.Geno(
    cad, 
    CHR="CHR", 
    POS="BP", 
    EA="Allele1", 
    NEA="Allele2",
    BETA="Effect", 
    SE="StdErr", 
    P="P-value", 
    EAF="Freq1", 
    keep_columns=False
)

print("[4/5] Preprocessing CAD...")
O_Geno.preprocess_data(preprocessing="Fill_delete", reference_panel="EUR_37")

# 3. Clump + harmonise (comme Marie)
print("[5/5] Clumping + harmonisation...")
print("  (ceci peut prendre 2-3 min, genal télécharge le panel 1000G EUR)")
E_clumped = E_Geno.clump(p1=5e-8, r2=0.01, kb=10000, reference_panel="EUR_37")
E_clumped.query_outcome(O_Geno, proxy=False)

# 4. Export CSV
E_df = pd.DataFrame(E_clumped.MR_data[0])
O_df = pd.DataFrame(E_clumped.MR_data[1])
merged = pd.merge(E_df, O_df, on="SNP", suffixes=("_LDL", "_CAD"))

# Calculate F-statistic
merged['F_stat'] = (merged['BETA_LDL'] / merged['SE_LDL']) ** 2

merged.to_csv("../data/MR_IVs_LDL_CAD.csv", index=False)

print("\n" + "=" * 60)
print(f"✓ SUCCÈS ! Généré {len(merged)} IVs")
print(f"  Fichier : data/MR_IVs_LDL_CAD.csv")
print(f"  F-stat médiane : {merged['F_stat'].median():.1f}")
print(f"  F-stat min/max : {merged['F_stat'].min():.1f} / {merged['F_stat'].max():.1f}")
print("=" * 60)
