## code to download and prepare `protein_subcellular_localizations` dataset
# This is Subcellular location data from The Protein Atlas, available at:
# https://www.proteinatlas.org/about/download
#
# see also documentation for "example_protein_subcellular_localizations" in "/R/data.R"

library(tidyverse)

# download subcellular localization data from the protein atlas

HPA_BASE_URL <- "https://www.proteinatlas.org/download/tsv"
FILE_NAME <- "subcellular_location.tsv.zip"
OVERWRITE <- FALSE

hpa_download_url <- glue::glue("{HPA_BASE_URL}/{FILE_NAME}")
hpa_save_path <- file.path("data-raw", FILE_NAME)


if(!file.exists(hpa_save_path)){
  download.file(hpa_download_url, destfile = hpa_save_path)
}

example_protein_subcellular_localizations <- load_and_clean_hpa_data(hpa_save_path)

# update file in /data
usethis::use_data(example_protein_subcellular_localizations, overwrite = TRUE)
