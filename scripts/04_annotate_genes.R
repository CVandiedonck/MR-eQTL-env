#!/usr/bin/env Rscript
# Annotation avec GenomicRanges + GENCODE (GRCh37) - version simplifiée

cat("============================================================\n")
cat("Annotation des SNPs avec GenomicRanges + GENCODE\n")
cat("============================================================\n\n")

# Install packages if needed
if (!requireNamespace("GenomicRanges", quietly=TRUE)) {
  cat("[!] Installation de GenomicRanges...\n")
  if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)
  }
  BiocManager::install("GenomicRanges", update=FALSE, ask=FALSE)
}

if (!requireNamespace("rtracklayer", quietly=TRUE)) {
  cat("[!] Installation de rtracklayer...\n")
  BiocManager::install("rtracklayer", update=FALSE, ask=FALSE)
}

if (!requireNamespace("data.table", quietly=TRUE)) {
  install.packages("data.table", repos="https://cloud.r-project.org", quiet=TRUE)
}

library(GenomicRanges)
library(rtracklayer)
library(data.table)

# Load data
cat("[1/4] Chargement des données MR...\n")
df <- fread("data/MR_IVs_LDL_CAD.csv")
cat(sprintf("  → %d SNPs chargés\n\n", nrow(df)))

# Download GENCODE GRCh37 GTF if not present
gtf_file <- "gencode.v19.annotation.gtf.gz"
if (!file.exists(gtf_file)) {
  cat("[2/4] Téléchargement GENCODE v19 (GRCh37)...\n")
  cat("  (fichier ~40 Mo)\n")
  download.file(
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    destfile = gtf_file,
    method = "auto",
    quiet = FALSE
  )
  cat("  → Téléchargement terminé\n\n")
} else {
  cat("[2/4] GENCODE v19 déjà présent\n\n")
}

# Load GENCODE annotations
cat("[3/4] Chargement des annotations GENCODE...\n")
gtf <- import(gtf_file)

# Keep only genes and expand boundaries by 500kb
genes <- gtf[gtf$type == "gene"]

# Expand gene ranges by 500kb upstream/downstream
genes <- resize(genes, width(genes) + 1000000, fix="center")

cat(sprintf("  → %d gènes chargés (fenêtre ±500kb)\n", length(genes)))

# Separate protein-coding from others for prioritization
protein_coding <- genes$gene_type == "protein_coding"
cat(sprintf("  → %d protein-coding, %d autres types\n\n", 
            sum(protein_coding), sum(!protein_coding)))

# Check chromosome naming
chr_style <- as.character(seqnames(genes)[1])
use_chr_prefix <- grepl("^chr", chr_style)

# Create GRanges for SNPs (match GTF naming)
cat("[4/4] Annotation des SNPs par overlap...\n")
if (use_chr_prefix) {
  snp_seqnames <- paste0("chr", df$CHR_LDL)
} else {
  snp_seqnames <- as.character(df$CHR_LDL)
}

snp_ranges <- GRanges(
  seqnames = snp_seqnames,
  ranges = IRanges(start = df$POS_LDL, end = df$POS_LDL),
  strand = "*"
)

# Find overlaps (SNP within expanded gene boundaries)
overlaps <- findOverlaps(snp_ranges, genes)

# Annotate - for multiple overlaps, choose closest gene
df$GENE <- "intergenic"

# Group overlaps by SNP
snp_gene_pairs <- data.frame(
  snp_idx = queryHits(overlaps),
  gene_idx = subjectHits(overlaps)
)

for (snp_idx in unique(snp_gene_pairs$snp_idx)) {
  # Get all genes overlapping this SNP
  gene_indices <- snp_gene_pairs$gene_idx[snp_gene_pairs$snp_idx == snp_idx]
  
  if (length(gene_indices) == 1) {
    # Single gene - easy
    df$GENE[snp_idx] <- genes$gene_name[gene_indices]
  } else {
    # Multiple genes - prioritize protein-coding
    is_protein_coding <- genes$gene_type[gene_indices] == "protein_coding"
    
    if (sum(is_protein_coding) > 0) {
      # Choose closest protein-coding gene
      gene_indices_pc <- gene_indices[is_protein_coding]
      snp_pos <- df$POS_LDL[snp_idx]
      gene_centers <- (start(genes[gene_indices_pc]) + end(genes[gene_indices_pc])) / 2
      distances <- abs(gene_centers - snp_pos)
      closest_idx <- gene_indices_pc[which.min(distances)]
      df$GENE[snp_idx] <- genes$gene_name[closest_idx]
    } else {
      # No protein-coding - choose closest overall
      snp_pos <- df$POS_LDL[snp_idx]
      gene_centers <- (start(genes[gene_indices]) + end(genes[gene_indices])) / 2
      distances <- abs(gene_centers - snp_pos)
      closest_idx <- gene_indices[which.min(distances)]
      df$GENE[snp_idx] <- genes$gene_name[closest_idx]
    }
  }
}

n_annotated <- sum(df$GENE != "intergenic")
n_intergenic <- sum(df$GENE == "intergenic")

cat(sprintf("  → %d SNPs dans des gènes\n", n_annotated))
cat(sprintf("  → %d SNPs intergéniques\n\n", n_intergenic))

# Reorder columns: GENE after SNP
cols <- names(df)
if ("GENE" %in% cols) {
  cols <- cols[cols != "GENE"]
}
snp_idx <- which(cols == "SNP")
new_order <- c(cols[1:snp_idx], "GENE", cols[(snp_idx+1):length(cols)])
setcolorder(df, new_order)

# Save
fwrite(df, "data/MR_IVs_LDL_CAD.csv")

cat("============================================================\n")
cat("✓ Annotation terminée avec succès\n")
cat(sprintf("  SNPs dans des gènes: %d / %d\n", n_annotated, nrow(df)))
cat(sprintf("  Intergéniques: %d\n\n", n_intergenic))

cat("Gènes clés du métabolisme lipidique :\n")
for (g in c("HMGCR","PCSK9","LDLR","APOE","APOB","CETP","LPA","SORT1","LIPC","LIPG")) {
  n <- sum(df$GENE == g, na.rm=TRUE)
  if (n > 0) cat(sprintf("  - %s: %d SNP(s)\n", g, n))
}

if (n_annotated > 0) {
  cat("\nTop 10 gènes (hors intergenic) :\n")
  top_genes <- head(sort(table(df$GENE[df$GENE != "intergenic"]), decreasing=TRUE), 10)
  for (i in seq_along(top_genes)) {
    cat(sprintf("  %2d. %s: %d SNP(s)\n", i, names(top_genes)[i], top_genes[i]))
  }
}
cat("============================================================\n")
