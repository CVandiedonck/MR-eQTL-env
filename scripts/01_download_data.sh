#!/bin/bash
# Download GWAS summary statistics from EBI GWAS Catalog
mkdir -p data_raw

LDL_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90468001-GCST90469000/GCST90468080/GCST90468080.tsv.gz"
CAD_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005194/CAD_META.gz"

if [ ! -f data_raw/LDL.tsv.gz ]; then
    echo "Downloading LDL summary stats..."
    wget "$LDL_URL" -O data_raw/LDL.tsv.gz
    echo "LDL: $(zcat data_raw/LDL.tsv.gz | wc -l) lines"
else
    echo "LDL already present, skipping."
fi

if [ ! -f data_raw/CAD.tsv.gz ]; then
    echo "Downloading CAD summary stats..."
    wget "$CAD_URL" -O data_raw/CAD.tsv.gz
    echo "CAD: $(zcat data_raw/CAD.tsv.gz | wc -l) lines"
else
    echo "CAD already present, skipping."
fi

echo "Done."
