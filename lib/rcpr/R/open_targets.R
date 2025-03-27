#' Shiny Open Targets Test
#'
#' A stand-alone app to test querying open targets using ensembl IDs
shiny_opentargets_test <- function() {
  shiny::shinyApp(
    ui = shiny::fluidPage(
      shiny::textInput(
        "ensembl",
        "Provide an ensembl gene ID",
        value = "ENSG00000119718"
      )
    ),
    server = function(input, output, session) {
      ot_results <- shiny::reactive({
        req(input$ensembl)

        summarize_open_targets_targets(
          input$ensembl
        )
      })

      observe({
        print(ot_results())
      })
    }
  )
}

#' Summarize Diseases
#'
#' Using a molecular neighborhood, target and indication summaries for all
#'   genes from Open Targets.
#'
#' @inheritParams create_neighborhood_table
#' @param neighborhood_summary_table tibble produced by
#'   \code{create_neighborhood_summary_table}
#' @inheritParams selectEntityServer
#' @inheritParams query_open_targets_targets
#'
#' @returns A list containing targets and indications
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'
#'   species_id <- rownames(consensus_model$species)[
#'     consensus_model$species$s_name == "KRAS:GTP"
#'   ]
#'   neighborhood_summary_table <- create_neighborhood_summary_table(
#'     create_neighborhood_table(
#'       species_id,
#'       consensus_model,
#'       consensus_graph,
#'       cpr = cpr,
#'       max_steps = 7L
#'     )
#'   )
#'   summarize_diseases(
#'     consensus_model,
#'     neighborhood_summary_table,
#'     species_identifiers
#'     )
#' }
#' @export
summarize_diseases <- function(
    consensus_model,
    neighborhood_summary_table,
    species_identifiers
    ) {
  checkmate::assertClass(consensus_model, "cpr.sbml.SBML_dfs")
  checkmate::assertDataFrame(neighborhood_summary_table)
  checkmate::assertDataFrame(species_identifiers)

  # extract a tibble of ensembl IDs for neighbors
  neighborhood_ensembl_ids <- format_neighborhood_ensembl_ids(
    neighborhood_summary_table,
    species_identifiers
  )

  if (nrow(neighborhood_ensembl_ids) == 0) {
    return(list(
      targets = tibble(),
      indications = tibble()
    ))
  }

  open_targets_results <- summarize_open_targets_targets(
    neighborhood_ensembl_ids$identifier
  )

  if ("targets" %in% names(open_targets_results)) {
    open_targets_results$targets <- open_targets_results$targets %>%
      dplyr::left_join(
        neighborhood_ensembl_ids %>%
          dplyr::select(ensembl_id = identifier, path_weight:sc_name),
        by = "ensembl_id"
      )
  }

  return(open_targets_results)
}

format_neighborhood_ensembl_ids <- function(
    neighborhood_summary_table,
    species_identifiers
    ) {
  checkmate::assertDataFrame(neighborhood_summary_table)

  if (!("ensembl_gene" %in% species_identifiers$ontology)) {
    stop(
      "no entries for \"ensembl_gene\" among identifiers.
     ensembl gene IDs are required for open targets's API"
    )
  }

  ensembl_ids <- species_identifiers$ontology_ids[species_identifiers$ontology == "ensembl_gene"][[1]]

  included_sids <- neighborhood_summary_table %>%
    dplyr::filter(node_type == "species")

  # filter ensembl IDs to entries in neighborhoods and
  # add path length / path weights
  neighborhood_ensembl_ids <- ensembl_ids %>%
    dplyr::inner_join(
      included_sids %>%
        dplyr::select(s_id, path_weight, path_length, sc_name),
      by = "s_id"
    ) %>%
    dplyr::arrange(path_weight, path_length) %>%
    dplyr::group_by(identifier) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  if (nrow(neighborhood_ensembl_ids) == 0) {
    message("No neighbors with ensembl IDs")
  } else {
    message(glue::glue("Query open targets with {nrow(neighborhood_ensembl_ids)} ensembl IDs from {nrow(included_sids)} species"))
  }

  return(neighborhood_ensembl_ids)
}

#' Summarize Open Targets - Targets
#'
#' Using Open Target's GraphQL API to read and format indications and
#'   tractability information for a set of genes.
#'
#' @inheritParams query_open_targets_targets
#'
#' @returns a list containing:
#' \describe{
#'   \item{targets}{tibble summarizing each target and their known chemical
#'     modulators}
#'   \item{indications}{tibble summarizing disease/phenotypic associations
#'     scores related to each target. NULL if no scores are present}
#' }
#'
#' @examples
#' target_ensembl_gene_ids <- c(
#'   "ENSG00000999999",
#'   "ENSG00000119718",
#'   "ENSG00000175354",
#'   "ENSG00000091831"
#'   )
#' summarize_open_targets_targets(target_ensembl_gene_ids)
#' @export
summarize_open_targets_targets <- function(target_ensembl_gene_ids) {

  # query open target's graphQL API
  target_results <- query_open_targets_targets(
    target_ensembl_gene_ids
  )

  # check for ids with no open targets data
  if (length(query_open_targets_targets) == 0) {
    warning(
      "No results were found in Open Targets for the provided ensembl IDs:",
      paste(target_ensembl_gene_ids, collapse = ", ")
    )

    return(list(
      targets = tibble::tibble(),
      indications = tibble::tibble()
    ))
  }

  found_target_ensembl_gene_ids <- purrr::map_chr(
    target_results,
    function(x) {
      x$id
    }
  )

  if (length(found_target_ensembl_gene_ids) == 0) {
    warning(glue::glue(
      "No results were returned from Open Targets"
    ))

    return(list(
      targets = tibble::tibble(),
      indications = tibble::tibble()
    ))
  }

  missing_target_ensembl_gene_ids <- setdiff(
    target_ensembl_gene_ids,
    found_target_ensembl_gene_ids
  )

  if (length(missing_target_ensembl_gene_ids) != 0) {
    warning(glue::glue(
      "No results for {length(missing_target_ensembl_gene_ids)} ensembl IDs
      were found in Open Targets:
      {paste(missing_target_ensembl_gene_ids, collapse = ', ')}"
    ))
  }

  # format entries as one-row per target
  targets_spread <- tibble::tibble(
    ensembl_id = found_target_ensembl_gene_ids
  ) %>%
    dplyr::mutate(target_results = target_results) %>%
    tidyr::unnest_wider(target_results)

  # summarize how readily druggable each target seems to be

  if (!("tractability" %in% colnames(targets_spread))) {
    tractability <- tibble(
      approvedSymbol = NA_character_,
      ab_tractability = NA_character_,
      sm_tractability = NA_character_,
      n_sms = 0
    )
  } else {
    tractability <- targets_spread %>%
      dplyr::select(approvedSymbol, tractability) %>%
      dplyr::filter(!purrr::map_lgl(tractability, is.null))

    if (nrow(tractability) == 0) {
      tractability <- tibble(
        approvedSymbol = NA_character_,
        ab_tractability = NA_character_,
        sm_tractability = NA_character_,
        n_sms = 0
      )
    } else {
      tractability <- tractability %>%
        dplyr::mutate(tractability = purrr::map(
          tractability,
          format_open_targets_tractability
        )) %>%
        tidyr::unnest(tractability)
    }
  }

  # find known drugs and chemical probes for each target
  drugs <- format_open_targets_drug_list(targets_spread)
  probes <- format_open_targets_probe_list(targets_spread)
  conservation <- format_open_targets_genetic_conservation(targets_spread)

  # generate a final target-level summary
  target_summary <- targets_spread %>%
    dplyr::select(ensembl_id, approvedSymbol, approvedName) %>%
    dplyr::left_join(tractability, by = "approvedSymbol")

  if (!is.null(drugs)) {
    target_summary <- target_summary %>%
      dplyr::left_join(drugs, by = "approvedSymbol")
  }
  if (!is.null(probes)) {
    target_summary <- target_summary %>%
      dplyr::left_join(probes, by = "approvedSymbol")
  }
  if (!is.null(conservation)) {
    target_summary <- target_summary %>%
      dplyr::left_join(conservation, by = "approvedSymbol")
  }

  target_summary <- target_summary %>%
    dplyr::mutate_if(is.character, stringr::str_wrap, width = 40, exdent = 2)

  # find all top indications for each target (max 25 per target)
  indications <- targets_spread %>%
    dplyr::mutate(diseases = purrr::map(
      associatedDiseases,
      format_open_targets_diseases
    )) %>%
    dplyr::select(approvedSymbol, diseases) %>%
    tidyr::unnest(diseases)

  if (!is.null(indications)) {
    data_types_present <- setdiff(
      colnames(indications),
      c(
        "approvedSymbol",
        "disease_name",
        "score",
        "therapeutic_area",
        "description"
        )
    ) %>%
      sort()

    indications <- indications %>%
      dplyr::arrange(dplyr::desc(score)) %>%
      dplyr::mutate(description = stringr::str_wrap(description, exdent = 2)) %>%
      dplyr::select(
        approvedSymbol,
        disease_name,
        score,
        !!!rlang::syms(data_types_present),
        therapeutic_area,
        description
      ) %>%
      dplyr::arrange(dplyr::desc(score))
  } else {
    indications <- tibble::tibble()
  }

  output <- list(
    targets = target_summary,
    indications = indications
  )

  return(output)
}

#' Query Open Targets - Targets
#'
#' Using Open Target's GraphQL API return gene-metadata and disease/phenotype
#'   associations for a set of genes.
#'
#' @param target_ensembl_gene_ids a character vector of ensembl gene names
#' @param post_n_genes separate queries into posts of this number of genes
#'
#' @returns a list of target attributes
#'
#' @examples
#' target_ensembl_gene_ids <- c("ENSG00000119718", "ENSG00000175354")
#' query_open_targets_targets(target_ensembl_gene_ids)
#' @export
query_open_targets_targets <- function(
  target_ensembl_gene_ids,
  post_n_genes = 10
  ) {

  checkmate::assertCharacter(
    target_ensembl_gene_ids,
    pattern = "^ENSG[0-9]{11}$",
    min.len = 1,
    max.len = 1000
  )

  checkmate::assertNumber(post_n_genes)

  # split large queries into multiple posts
  max_posts <- ceiling(length(target_ensembl_gene_ids) / post_n_genes)

  post_sets <- tibble::tibble(
    gene = target_ensembl_gene_ids,
    set = ceiling(seq(
      0.0001,
      max_posts,
      length.out = length(target_ensembl_gene_ids)
      ))
  ) %>%
    tidyr::nest(genes = -set) %>%
    {
      purrr::map(.$genes, function(x) {
        x$gene
      })
    } %>%
    # post a query for each set of genes
    purrr::map(
      post_open_targets_query
    )

  post_sets <- do.call(c, post_sets)

  return(post_sets)
}

post_open_targets_query <- function(target_ensembl_gene_ids) {
  query_string <- "
    query targets(
    $ensemblIds: [String!]!,
    ){
    targets(ensemblIds: $ensemblIds){
        id
        approvedSymbol,
        approvedName,
        associatedDiseases {
          count,
          rows {
            score,
            datatypeScores {
              id,
              score
            },
            disease {
              id,
              name,
              description,
              therapeuticAreas {
                name
              }
            }
          }
        },
    		tractability {
        	label,
        	modality,
        	value
      	},
        geneticConstraint{
          constraintType
          exp
          obs
          score
          oe
          oeLower
          oeUpper
        },
        chemicalProbes {
          id,
          control,
          drugId,
          mechanismOfAction,
          isHighQuality,
          urls {
            url
          }
        },
        knownDrugs {
          count,
          rows {
            drugId,
            drug {
              name,
              tradeNames,
              drugType,
              yearOfFirstApproval
              }
            }
          }
        }
      }
  "

  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # Set variables object of arguments to be passed to endpoint
  variables <- list(
    "ensemblIds" = target_ensembl_gene_ids
  )
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)

  # Perform POST request
  r <- httr::POST(url = base_url, body = post_body, encode = "json")
  validate_open_targets_post(r)
  target_results <- httr::content(r)$data$targets

  return(target_results)
}

format_open_targets_diseases <- function(diseases) {
  checkmate::assertList(diseases, len = 2, names = "unique")
  stopifnot(names(diseases) == c("count", "rows"))

  if (diseases$count == 0) {
    return(NULL)
  }

  indication_scores <- purrr::map(diseases$rows, function(x) {
    spread_datatype_scores <- purrr::map(
      x$datatypeScores,
      function(y) {
        tibble::tibble(id = y$id, score = y$score)
      }
    ) %>%
      dplyr::bind_rows() %>%
      tidyr::spread(id, score)

    # x$disease$therapeuticAreas

    tibble::tibble(
      score = x$score,
      disease_name = x$disease$name,
      disease_id = x$disease$id,
      therapeutic_area = ifelse(
        is.null(x$disease$therapeuticAreas),
        "",
        purrr::map_chr(x$disease$therapeuticAreas, function(x) {
          x$name
        }) %>%
          paste(collapse = ", + ")
      ),
      description = ifelse(
        is.null(x$disease$description),
        "",
        x$disease$description
      )
    ) %>%
      dplyr::bind_cols(spread_datatype_scores)
  }) %>%
    dplyr::bind_rows()

  # rearrange and round scores
  data_types_present <- setdiff(
    colnames(indication_scores),
    c("disease_id", "disease_name", "score", "description")
  ) %>%
    sort()

  reformated_scores <- indication_scores %>%
    dplyr::mutate_if(is.numeric, round, 2) %>%
    dplyr::select(disease_name, score, !!!rlang::syms(data_types_present), description) %>%
    dplyr::arrange(dplyr::desc(score))

  return(reformated_scores)
}

format_open_targets_probe_list <- function(targets_spread) {
  checkmate::assertDataFrame(targets_spread)

  if (!("chemicalProbes" %in% colnames(targets_spread)) || !inherits(targets_spread$chemicalProbes, "list")) {
    return(NULL)
  }

  valid_probes <- targets_spread %>%
    dplyr::select(approvedSymbol, chemicalProbes) %>%
    dplyr::filter(!purrr::map_lgl(chemicalProbes, is.null)) %>%
    tidyr::unnest(chemicalProbes) %>%
    dplyr::filter(purrr::map_int(chemicalProbes, length) != 0)

  if (nrow(valid_probes) == 0) {
    # some targets have a probes entry which doesn't include a named probe
    return(NULL)
  }

  out <- valid_probes %>%
    tidyr::unnest_wider(chemicalProbes) %>%
    dplyr::group_by(approvedSymbol) %>%
    dplyr::summarize(
      `high quality probes` = paste(id[isHighQuality], collapse = ", "),
      `N probes` = length(id)
    )

  return(out)
}

format_open_targets_tractability <- function(tractability) {
  checkmate::assertList(tractability)

  # combine tractability summaries
  tractability_df <- tibble::tibble(dat = tractability) %>%
    tidyr::unnest_wider(dat) %>%
    dplyr::filter(value == TRUE) %>%
    dplyr::group_by(modality) %>%
    dplyr::summarize(label = paste(unique(sort(label)), collapse = ", "))

  # add missing categories
  missing_tractability_df <- tibble::tibble(
    modality = c("AB", "PR", "SM"),
    label = NA_character_
  ) %>%
    dplyr::anti_join(tractability_df, by = "modality")

  out <- dplyr::bind_rows(tractability_df, missing_tractability_df) %>%
    dplyr::mutate(
      modality = dplyr::case_when(
        modality == "AB" ~ "antibody tractability",
        modality == "PR" ~ "protac tractability",
        modality == "SM" ~ "small molecule tractability"
      ),
      modality = factor(modality, levels = c(
        "small molecule tractability",
        "antibody tractability",
        "protac tractability"
      ))
    ) %>%
    tidyr::spread(modality, label)

  return(out)
}


format_open_targets_drug_list <- function(targets_spread) {
  checkmate::assertDataFrame(targets_spread)

  # if there are no records then knownDrugs will just be an NA

  if (!("knownDrugs" %in% colnames(targets_spread)) || !inherits(targets_spread$knownDrugs, "list")) {
    return(NULL)
  }

  valid_drugs <- targets_spread %>%
    dplyr::select(approvedSymbol, knownDrugs) %>%
    dplyr::filter(!purrr::map_lgl(knownDrugs, is.null)) %>%
    tidyr::unnest_wider(knownDrugs) %>%
    dplyr::filter(!is.na(rows))

  if (nrow(valid_drugs) == 0) {
    return(NULL)
  }

  valid_drugs <- valid_drugs %>%
    tidyr::unnest(rows) %>%
    tidyr::unnest_wider(rows) %>%
    tidyr::unnest_wider(drug)

  if (!("yearOfFirstApproval" %in% colnames(valid_drugs))) {
    valid_drugs <- valid_drugs %>%
      dplyr::mutate(yearOfFirstApproval = NA)
  }

  valid_drugs <- valid_drugs %>%
    dplyr::distinct(approvedSymbol, name, yearOfFirstApproval) %>%
    dplyr::mutate(
      name = stringr::str_to_title(name),
      label = dplyr::case_when(
        is.na(yearOfFirstApproval) ~ name,
        TRUE ~ as.character(glue::glue("{name} ({yearOfFirstApproval})"))
      )
    ) %>%
    dplyr::arrange(dplyr::desc(yearOfFirstApproval)) %>%
    dplyr::group_by(approvedSymbol) %>%
    dplyr::summarize(drugs = stringr::str_c(label, collapse = "; "))

  return(valid_drugs)
}


format_open_targets_genetic_conservation <- function(targets_spread) {
  checkmate::assertDataFrame(targets_spread)

  if (!("geneticConstraint" %in% colnames(targets_spread))) {
    return(NULL)
  }

  out <- targets_spread %>%
    dplyr::select(approvedSymbol, geneticConstraint) %>%
    dplyr::filter(!purrr::map_lgl(geneticConstraint, is.null)) %>%
    dplyr::mutate(dat = purrr::map(geneticConstraint, function(x) {
      tibble::tibble(x) %>%
        tidyr::unnest_wider(x)
    })) %>%
    dplyr::select(-geneticConstraint) %>%
    tidyr::unnest(dat) %>%
    dplyr::mutate(
      var = dplyr::case_when(
        constraintType == "syn" ~ "o/e synonymous",
        constraintType == "mis" ~ "o/e non-synonymous",
        constraintType == "lof" ~ "o/e loss-of-function"
      ),
      var = factor(
        var,
        levels = c(
          "o/e synonymous",
          "o/e non-synonymous",
          "o/e loss-of-function"
        )
      )
    ) %>%
    dplyr::select(approvedSymbol, var, oe) %>%
    tidyr::spread(var, oe)

  return(out)
}


#' Summarize Indications
#'
#' Using a molecular neighborhood, find all Open Target scores for a specified disease.
#'
#' @inheritParams query_open_targets_indications
#' @inheritParams summarize_diseases
#'
#' @returns A tibble of species names and disease scores
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'
#'   species_id <- rownames(consensus_model$species)[
#'     consensus_model$species$s_name == "pyruvate kinase tetramer"
#'   ]
#'   neighborhood_summary_table <- create_neighborhood_summary_table(
#'     create_neighborhood_table(
#'       species_id,
#'       consensus_model,
#'       consensus_graph,
#'       cpr = cpr,
#'       max_steps = 7L
#'     )
#'   )
#'   disease_id <- "EFO_0000400" # diabetes
#'   summarize_indication(disease_id, consensus_model, neighborhood_summary_table, species_identifiers)
#' }
#' @export
summarize_indication <- function(
  disease_id,
  consensus_model,
  neighborhood_summary_table,
  species_identifiers
  ) {

  checkmate::assertString(disease_id)
  checkmate::assertClass(consensus_model, "cpr.sbml_dfs_core.SBML_dfs")
  checkmate::assertDataFrame(neighborhood_summary_table)
  checkmate::assertDataFrame(species_identifiers)

  # extract a tibble of ensembl IDs for neighbors
  neighborhood_ensembl_ids <- format_neighborhood_ensembl_ids(
    neighborhood_summary_table = neighborhood_summary_table,
    species_identifiers = species_identifiers
  )

  if (nrow(neighborhood_ensembl_ids) == 0) {
    return(
      tibble::tibble(
        s_id = NA_character_,
        id = NA_character_,
        score = NA_real_,
        approvedSymbol = NA_character_,
        approvedName = NA_character_,
        url = NA_character_
      ) %>%
        dplyr::slice(-1)
    )
  }

  # call the graphQL API repeatedly since it often says
  # "try again in 30 seconds"
  max_retries <- 2
  query_tries <- 0
  valid_query <- FALSE
  while (query_tries < max_retries && !valid_query) {

    open_targets_results <- query_open_targets_indications(
      neighborhood_ensembl_ids$identifier,
      disease_id
    )

    if (!is.null(open_targets_results)) {
      valid_query <- TRUE
    }
    query_tries <- query_tries + 1
    Sys.sleep(30)
  }

  if (!valid_query) {
    error_message <- glue::glue("The Open Targets query failed {max_retries} times, try again later")
  }

  if (valid_query && open_targets_results$target_scores$count > 0) {
    # at least one neighbor has a score for disease
    indications_df <- tibble::tibble(dat = open_targets_results$target_scores$rows) %>%
      tidyr::unnest_wider(dat) %>%
      tidyr::unnest_wider(target) %>%
      dplyr::mutate(url = glue::glue("https://platform.opentargets.org/evidence/{id}/{disease_id}"))
  } else {
    out <- tibble::tibble(
      s_id = NA_character_,
      id = NA_character_,
      score = NA_real_,
      approvedSymbol = NA_character_,
      approvedName = NA_character_,
      url = NA_character_
    ) %>%
      dplyr::slice(-1)

    if (!valid_query) {
      # smuggle an error message out of the query
      attr(out, 'ot_error') <- error_message
    }

    return(out)
  }

  indication_scores <- neighborhood_ensembl_ids %>%
    dplyr::select(s_id, id = identifier) %>%
    dplyr::left_join(indications_df, by = "id") %>%
    dplyr::mutate(score = ifelse(is.na(score), 0, score)) %>%
    dplyr::group_by(s_id) %>%
    dplyr::arrange(dplyr::desc(score)) %>%
    # retain only the top score when multiple ensembl IDs are associated with the same node
    dplyr::slice(1) %>%
    dplyr::ungroup()

  return(indication_scores)
}

#' Query Open Targets - Indications
#'
#' Using Open Target's GraphQL API return indication score for a
#'   disease/phenotype of interest across a set of genes.
#'
#' @param target_ensembl_gene_ids a character vector of ensembl gene names
#' @param disease_id A disease to search from the EFO ontology
#'
#' @returns a list of target attributes
#'
#' @examples
#' target_ensembl_gene_ids <- c("ENSG00000090104", "ENSG00000163599")
#' disease_id <- "EFO_0001060" # celiac's disease
#' query_open_targets_indications(target_ensembl_gene_ids, disease_id)
#' @export
query_open_targets_indications <- function(target_ensembl_gene_ids,
                                           disease_id) {
  checkmate::assertCharacter(
    target_ensembl_gene_ids,
    pattern = "^ENSG[0-9]{11}$",
    min.len = 1,
    max.len = 1000
  )
  checkmate::assertString(disease_id)

  query_string <- "
    query disease (
      $ensemblIds: [String!]!,
      $efoId: String!
    )
    {
    disease(efoId: $efoId){
      id,
    	name,
    	target_scores: associatedTargets (Bs: $ensemblIds) {
    	  count,
        rows {
          score,
          target {
            id,
            approvedSymbol,
            approvedName
            }
          }
        }
      }
    }
  "

  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

  # Set variables object of arguments to be passed to endpoint
  variables <- list(
    "efoId" = disease_id,
    "ensemblIds" = target_ensembl_gene_ids
  )
  # Construct POST request body object with query string and variables
  post_body <- list(query = query_string, variables = variables)

  # Perform POST request
  r <- httr::POST(url = base_url, body = post_body, encode = "json")
  validate_open_targets_post(r)

  disease_results <- httr::content(r)$data$disease

  return(disease_results)
}

validate_open_targets_post <- function(r) {
  errors <- httr::content(r)$errors[[1]]
  if (!is.null(errors)) {
    stop(errors$message)
  }

  if (r$status_code %in% c(400, 502)) {
    html_error <- httr::content(r) %>%
      rvest::html_element("body") %>%
      rvest::html_element("p") %>%
      rvest::html_text() %>%
      stringr::str_trim()

    if (r$status_code == 400) {
      stop(
        "The GraphQL query failed with the following error message:\n",
        html_error
      )
    } else if (r$status_code == 502) {
      warning(
        "The GraphQL query failed with the following error message:\n",
        html_error
      )
    } else {
      stop("undefined status code")
    }
  }

  return(invisible(0))
}
