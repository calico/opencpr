## code to download and prepare `example_gtex_data` dataset
# This is Median gene-level TPM by tissue from GTEx Analysis V8 data
# (dbGaP Accession phs000424.v8.p2)
# Median expression was calculated (by GTEX) from the file
# GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz.
# see also documentation for "example_gtex_data" in "/R/data.R"

suppressPackageStartupMessages(library(tidyverse))

GTEX_URL <- "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
GTEX_REL_PATH <- "data-raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"

# download data from GTEX
# this is from GTEx Analysis V8 RNA-Seq data (dbGaP Accession phs000424.v8.p2)
if(!file.exists(GTEX_REL_PATH)){
  download.file(GTEX_URL, destfile = GTEX_REL_PATH)
}

gtex_expression_data <- load_and_clean_gtex_data(GTEX_REL_PATH)

# trim to just a few tissues to reduce dataset size
example_gtex_data <- gtex_expression_data %>%
  dplyr::select(
    ensembl_gene_id,
    ensembl_geneTranscript_id,
    Description,
    Liver,
    Spleen,
    Lung,
    Testis,
    Ovary,
    Whole.Blood,
    Kidney...Cortex,
    Kidney...Medulla
    )

# further trim to genes present in the example string network
example_gtex_data <- example_gtex_data %>%
  filter(ensembl_gene_id %in% rcpr::example_string_graph$genes$gene)

# update file in /data
usethis::use_data(example_gtex_data, overwrite = TRUE)
