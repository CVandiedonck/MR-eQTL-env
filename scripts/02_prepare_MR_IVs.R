#!/usr/bin/env Rscript
# Prepare MR instrumental variables: LDL -> CAD
# Inputs : data_raw/LDL.tsv.gz, data_raw/CAD.tsv.gz
# Output : data/MR_IVs_LDL_CAD.csv (replaces placeholder)

library(data.table)
library(dplyr)

# ── PARAMETERS — EDIT THESE ─────────────────────────────────────────────────
PLINK_BIN <- "plink"       # or full path e.g. "/usr/local/bin/plink"
REF_PANEL <- "ref/EUR"     # path to EUR.bed/bim/fam WITHOUT extension
P_THRESH  <- 5e-8
R2_THRESH <- 0.01
KB_WINDOW <- 10000
OUTFILE   <- "data/MR_IVs_LDL_CAD.csv"
# ────────────────────────────────────────────────────────────────────────────

cat("=== Loading LDL ===\n")
ldl_raw <- fread("data_raw/LDL.tsv.gz")
ldl <- ldl_raw %>%
    rename(SNP=rs_id, CHR=chromosome, POS=base_pair_location,
           EA=effect_allele, NEA=other_allele,
           BETA=beta, SE=standard_error, P=p_value, EAF=effect_allele_frequency) %>%
    filter(!is.na(SNP), !is.na(BETA), !is.na(SE), !is.na(P),
           nchar(EA)==1, nchar(NEA)==1) %>%
    mutate(EA=toupper(EA), NEA=toupper(NEA))
cat("After QC:", nrow(ldl), "SNPs\n")

ldl_sig <- ldl %>% filter(P < P_THRESH)
cat("Significant (p <", P_THRESH, "):", nrow(ldl_sig), "SNPs\n")
fwrite(ldl_sig %>% select(SNP, P), "tmp_plink_input.txt", sep="\t")

cat("\n=== PLINK clumping ===\n")
plink_cmd <- paste(PLINK_BIN,
    "--bfile", REF_PANEL,
    "--clump tmp_plink_input.txt",
    "--clump-p1", P_THRESH,
    "--clump-p2 0.01",
    "--clump-r2", R2_THRESH,
    "--clump-kb", KB_WINDOW,
    "--clump-field P --clump-snp-field SNP",
    "--out tmp_clump")
system(plink_cmd)

if (!file.exists("tmp_clump.clumped")) stop("PLINK output not found. Check PLINK_BIN and REF_PANEL.")
clumped <- fread("tmp_clump.clumped")
cat("Independent IVs after clumping:", nrow(clumped), "\n")
ldl_iv <- ldl %>% filter(SNP %in% clumped$SNP)

cat("\n=== Loading CAD ===\n")
cad_raw <- fread("data_raw/CAD.tsv.gz")
cad <- cad_raw %>%
    rename(SNP=oldID, CHR=CHR, POS=BP, EA=Allele1, NEA=Allele2,
           BETA=Effect, SE=StdErr, P=`P-value`, EAF=Freq1) %>%
    mutate(EA=toupper(EA), NEA=toupper(NEA)) %>%
    filter(!is.na(SNP), !is.na(BETA), !is.na(SE))

cat("\n=== Harmonisation ===\n")
merged <- inner_join(
    ldl_iv %>% select(SNP, CHR, POS, EA, NEA, BETA, SE, P, EAF),
    cad    %>% select(SNP, EA, NEA, BETA, SE, P),
    by="SNP", suffix=c("_LDL","_CAD"))
cat("SNPs after merge:", nrow(merged), "\n")

ok   <- merged$EA_LDL == merged$EA_CAD
flip <- merged$EA_LDL == merged$NEA_CAD & merged$NEA_LDL == merged$EA_CAD
ambiguous <- !ok & !flip
cat("  Concordant:", sum(ok), "| Flipped:", sum(flip), "| Removed:", sum(ambiguous), "\n")

merged$BETA_CAD[flip] <- -merged$BETA_CAD[flip]
merged <- merged %>%
    filter(!ambiguous) %>%
    mutate(BETA_CAD = BETA_CAD * sign(BETA_LDL),
           BETA_LDL = abs(BETA_LDL))

merged$F_stat <- (merged$BETA_LDL / merged$SE_LDL)^2

final <- merged %>%
    select(SNP, CHR, POS, EA=EA_LDL, NEA=NEA_LDL,
           BETA_LDL, SE_LDL, P_LDL=P_LDL, EAF,
           BETA_CAD, SE_CAD, P_CAD=P_CAD, F_stat)

fwrite(final, OUTFILE)
cat("\n=== Summary ===\n")
cat("Final IVs:", nrow(final), "\n")
cat("F-stat: min =", round(min(final$F_stat),1),
    " median =", round(median(final$F_stat),1),
    " max =", round(max(final$F_stat),1), "\n")
cat("Written to:", OUTFILE, "\n")

file.remove(c("tmp_plink_input.txt","tmp_clump.clumped","tmp_clump.log"))
cat("Done.\n")
