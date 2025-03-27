#' Bioconductor Org Package Prefix
#'
#' Find the prefix of the bioconductor package and its functions for a species
#'
#' @param species Species name
#'
#' @return character prefix
#'
#' @examples
#' bioconductor_org_package_prefix("Homo sapiens")
#' @export
bioconductor_org_package_prefix <- function(species) {
  stopifnot(class(species) == "character", length(species) == 1)

  if (!stringr::str_detect(species, "^[A-Z][a-z]+ [A-Za-z]+$")) {
    stop(species, " is not a valid species name")
  }

  split_species_name <- stringr::str_split(species, pattern = " ")[[1]]
  species_abbrev <- paste0(
    stringr::str_split(split_species_name[1], pattern = "")[[1]][1],
    stringr::str_split(split_species_name[2], pattern = "")[[1]][1]
  )

  governing_org <- dplyr::case_when(
    species_abbrev == "Sc" ~ "sgd",
    .default = "eg" # ensembl genomes
  )


  prefix <- glue::glue("org.{species_abbrev}.{governing_org}")

  return(prefix)
}

#' Bioconductor Org Library
#'
#' Load the organism bioconductor package for the specified species
#'
#' @inheritParams bioconductor_org_package_prefix
#'
#' @return 0 invisibly
#'
#' @examples
#'
#' if (interactive()) {
#'   bioconductor_org_library("Homo sapiens")
#' }
bioconductor_org_library <- function(species) {
  stopifnot(class(species) == "character", length(species) == 1)

  prefix <- bioconductor_org_package_prefix(species)

  package_name <- paste0(prefix, ".db")

  if (!package_name %in% utils::installed.packages()) {
    stop(glue::glue("{package_name} is not installed, install with:

      if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

      BiocManager::install('{package_name}')

      (If this doesn't work then package may follow a non-standard naming convention.)"))
  }

  suppressPackageStartupMessages(do.call(library, list(package = package_name)))

  return(invisible(0))
}

#' Bioconductor Org Function
#'
#' Call a function from one of the bioconductor org packages
#'
#' @param object_type Type of function to call
#' @inheritParams bioconductor_org_package_prefix
#'
#' @return a tibble for function calls which can be coerced to a tabular result or the raw results otherwise
#'
#' @examples
#' \dontrun{
#' bioconductor_org_function("ENZYME", "Homo sapiens")
#' }
#' @export
bioconductor_org_function <- function(object_type, species) {
  stopifnot(class(object_type) == "character", length(object_type) == 1)
  stopifnot(class(species) == "character", length(species) == 1)

  bioconductor_org_library(species)

  prefix <- bioconductor_org_package_prefix(species)
  package_name <- paste0(prefix, ".db")

  pkg_functions <- tibble::tibble(full_object = ls(glue::glue("package:{package_name}"))) %>%
    dplyr::mutate(object_type = stringr::str_replace(full_object, prefix, ""))

  full_object <- pkg_functions["full_object"][pkg_functions["object_type"] == object_type]

  if (length(full_object) != 1) {
    print(pkg_functions %>% glue::glue_data("Use {object_type} for {full_object}"))
    stop(glue::glue("{object_type} is not a valid function type, all valid types have been printed"))
  }

  org_object <- get(full_object)

  result <- try(toTable(org_object), silent = TRUE)
  if ("try-error" %in% class(result)) {
    print("Object could not be converted to a table, returning raw result")
    return(org_object)
  } else {
    return(result %>%
      tibble::as_tibble())
  }
}
