#' Interactive Initialization Wrapper
#'
#' Determine whether object created by interactive_initialization() are
#'   are present in the global environment, and if not call
#'   interactive_initialization()
#'
#' @param testing_env Logical indicating whether to load assets for testing (TRUE)
#'   or, prime-time pathway representations (FALSE)
#'
#' @returns 0, invisibly
#'
#' @export
interactive_initialization_wrapper <- function(testing_env = FALSE) {
  required_globals <- c(
    "cpr",
    "consensus_model",
    "consensus_graph",
    "precomputed_distances",
    "species_names",
    "species_identifiers"
  )

  if (!all(required_globals %in% ls(envir = globalenv()))) {
    interactive_initialization(testing_env)
  }

  return(invisible(0))
}

#' Interactive Initialization
#'
#' Add objects which support MechNet to the global environment. This can be
#'   used to run the full app or test its interactive and non-interactive
#'   components.
#'
#' @param testing_env Work off a simple environment
#'
#' Note that this approach breaks lots of conventions (assignment to global
#' environment; sourcing packages within functions, ...).
interactive_initialization <- function(testing_env = FALSE) {

  checkmate::assertLogical(testing_env, len = 1)

  MY_LOCAL_VENV <- "<<PATH_TO_VENV>>"

  MY_CPR_ASSETS <- tibble::tribble(
    ~ type, ~ object, ~local_path,
    "model", "consensus_nocofactors.pkl", "<<PATH_TO_ASSETS>>/output_sbml_dfs.pkl",
    "identifiers", "species_identifiers.tsv", "<<PATH_TO_ASSETS>>/species_identifiers.tsv",
    "directed regulatory", "cpr_graph_weighted.pkl", "<<PATH_TO_ASSETS>>/cpr_graph_weighted.pkl",
    "directed regulatory distances", "precomputed_distances.pkl", "<<PATH_TO_ASSETS>>/precomputed_distances.pkl"
  )

  MY_LOCAL_CPR_PATH <- "<<PATH_TO_CPR>>/calicolabs-open-cpr-py"

  checkmate::assertDirectoryExists(MY_LOCAL_VENV)
  checkmate::assertDirectoryExists(MY_LOCAL_CPR_PATH)
  reticulate::use_virtualenv(MY_LOCAL_VENV, required = TRUE)

  # bind python packages - note the global assignment of "cpr"
  assign("cpr", reticulate::import("cpr"), envir = globalenv())

  # load R packages
  suppressPackageStartupMessages(library(dplyr))
  library(ggplot2)

  # trigger local debugging messages (you'll see these sprinkled throughout the code)
  debugr::debugr_switchOn()

  if (testing_env) {
    cpr_list <- load_assets_for_testing(MY_LOCAL_CPR_PATH)
  } else {
    cpr_list <- load_assets_from_files(MY_CPR_ASSETS)
  }

  species_names <- cpr_list$consensus_model$species %>%
    dplyr::select(-s_Identifiers, -s_Source)

  identifiers <- cpr_list$species_identifiers %>%
    tidyr::nest(ontology_ids = -ontology) %>%
    dplyr::mutate(n = purrr::map_int(ontology_ids, nrow)) %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(ontology = factor(ontology, levels = ontology))

  # assign to global environment
  assign("consensus_model", cpr_list$consensus_model, envir = globalenv())
  assign("consensus_graph", cpr_list$consensus_graph, envir = globalenv())
  assign("precomputed_distances", cpr_list$precomputed_distances, envir = globalenv())
  assign("species_identifiers", cpr_list$species_identifiers, envir = globalenv())
  assign("identifiers", cpr_list$species_identifiers, envir = globalenv())
  assign("species_names", species_names, envir = globalenv())

  return(invisible(0))
}

load_assets_for_testing <- function(cpr_path) {

  os <- reticulate::import("os")

  test_data <- file.path(cpr_path, "src", "tests", "test_data")
  checkmate::assertDirectoryExists(test_data)

  pw_index <- cpr$indices$PWIndex(os$path$join(test_data, "pw_index_metabolism.tsv"))
  sbml_dfs_dict <- cpr$consensus$construct_sbml_dfs_dict(pw_index)
  sbml_dfs <- cpr$consensus$construct_consensus_model(sbml_dfs_dict, pw_index)
  cpr_graph <- cpr$network$net_create$process_cpr_graph(sbml_dfs, graph_type="surrogate")
  precomputed_distances <- cpr$network$precompute$precompute_distances(cpr_graph)

  species_identifiers <- sbml_dfs$get_identifiers("species")
  species_identifiers <- cpr$sbml_dfs_core$filter_to_characteristic_species_ids(species_identifiers)

  output <- list(
    consensus_model = sbml_dfs,
    consensus_graph = cpr_graph,
    precomputed_distances = precomputed_distances,
    species_identifiers = species_identifiers
  )

  return(output)
}

load_assets_from_files <- function(cpr_assets) {

  cpr_asset_paths <- pluck_required_cpr_assets(
    cpr_assets,
    c("model", "identifiers", "directed regulatory", "directed regulatory distances")
  )

  debugr::dwatch(msg = "Unpickling network model [cpr<mechnet.R>::interactive_initialization]")

  # download consensus model with pickle
  # if there are weird import errors then the cpr python package probably needs to be reinstalled
  consensus_model <- cpr$utils$load_pickle(cpr_asset_paths[["model"]])
  consensus_model$validate()

  # read a python igraph network
  debugr::dwatch(msg = "Unpickling network graph [cpr<mechnet.R>::interactive_initialization]")
  consensus_graph <- cpr$utils$load_pickle(cpr_asset_paths[["directed regulatory"]])

  debugr::dwatch(msg = "Unpickling pre-computed distances [cpr<mechnet.R>::interactive_initialization]")
  precomputed_distances <- cpr$utils$load_pickle(cpr_asset_paths[["directed regulatory distances"]])

  debugr::dwatch(msg = "Organizing identifiers [cpr<mechnet.R>::interactive_initialization]")
  species_identifiers <- readr::read_tsv(cpr_asset_paths[["identifiers"]], show_col_types = FALSE)

  output <- list(
    consensus_model = consensus_model,
    consensus_graph = consensus_graph,
    precomputed_distances = precomputed_distances,
    species_identifiers = species_identifiers
  )

  return(output)
}

validate_cpr_assets <- function(cpr_assets) {

  checkmate::assertDataFrame(cpr_assets)

  # check for a misformatted DF
  required_cpr_asset_vars <- c("type", "object", "local_path")
  missing_vars <- setdiff(required_cpr_asset_vars, colnames(cpr_assets))

  if (length(missing_vars) > 0) {
    cli::cli_abort(
      "The following required column{?s} {?is/are} missing from the CPR assets data frame: {.field {missing_vars}}"
    )
  }

  # check for missing assets
  missing_assets <- cpr_assets %>%
    dplyr::filter(!purrr::map_lgl(local_path, file.exists))

  if (nrow(missing_assets) > 0) {
    cli::cli_abort(
      "The following {nrow(missing_assets)} asset{?s} {?is/are} missing from the local file system: {.field {missing_assets$object}}.
        They were expected at {.file {missing_assets$local_path}}"
    )
  }

  # check for duplicated assets
  nonunique_assets <- unique(cpr_assets$type[duplicated(cpr_assets$type)])
  if (length(nonunique_assets) > 0) {
    cli::cli_abort(
      "The following asset type{?s} {?is/are} duplicated in the CPR assets data frame: {.field {nonunique_assets}}"
    )
  }

  return(invisible(0))

}

pluck_required_cpr_assets <- function(cpr_assets, required_assets) {

  checkmate::assertDataFrame(cpr_assets)
  validate_cpr_assets(cpr_assets)

  checkmate::assertCharacter(required_assets, min.len = 1)
  missing_asset_types <- setdiff(required_assets, cpr_assets$type)
  if (length(missing_asset_types) > 0) {
    cli::cli_abort(
      "The following required asset type{?s} {?is/are} missing from the CPR assets data frame: {.field {missing_asset_types}}"
    )
  }

  named_asset_paths <- cpr_assets %>%
    dplyr::filter(type %in% required_assets) %>%
    {stats::setNames(.$local_path, .$type)}

  return(named_asset_paths)
}
