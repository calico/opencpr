#' Compartmentalize STRING
#'
#' Apply netcontextr to compartmentalize the STRING interaction network
#' based on sub-cellular localization data from the human proteome atlas
#'
#' @inheritParams createStringGraph
#' @inheritParams load_and_clean_hpa_data
#'
#' @returns
#' a named list:
#' \describe{
#'   \item{compartmentalized_genes}{a tibble containing ensembl_gene and the comparments where the gene is expressed and interacts}
#'   \item{compartmentalized_edges}{an edgelist of interactions between pairs of genes and the compartments where these interactions could occur}
#'   \item{compartments}{a tibble of compartment meta-data}
#'   }
#'
#' @examples
#' if (interactive()) {
#'
#'   stringNetwork = "/group/cpr/external_pathways/string/9606.protein.links.full.v11.5.txt"
#'   aliases = "/group/cpr/external_pathways/string/9606.protein.aliases.v11.5.txt.gz"
#'   hpa_localization_data = "/group/cpr/external_pathways/protein_atlas/subcellular_location.tsv"
#'
#'   compartmentalize_string(
#'     stringNetwork,
#'     aliases,
#'     hpa_localization_data
#'     )
#' }
#' @export
compartmentalize_string <- function (
  stringNetwork,
  aliases,
  hpa_localization_data
) {

  # load STRING and add aliases
  full_string_graph <- createStringGraph(
    stringNetwork = stringNetwork,
    aliases = aliases
  )

  # load protein localization
  hpa_localization <- load_and_clean_hpa_data(hpa_localization_data)

  # localize string genes
  localized_full_string_graph <- annotate_genes(full_string_graph, hpa_localization, "compartment")

  # localize edges
  compartmentalized_string_graph <- trim_network_by_gene_attribute(
    localized_full_string_graph, "compartment",
    addAttributeToInteractions = TRUE
    )

  # return a named list of genes, interactions and compartments which
  # will map to a dict of pd.DataFrames

  # genes only metadata is compartment which can be derived from interactions
  # so there is no need to work with genes table directly

  compartmentalized_interactions <- compartmentalized_string_graph$interactions %>%
    dplyr::select(.data$protein1, .data$protein2, compartment, combined_score) %>%
    tidyr::unnest(compartment)

  stopifnot(compartmentalized_interactions %>% dplyr::distinct(protein1, protein2, compartment) %>% nrow() == nrow(compartmentalized_interactions))

  # format compartments
  compartments <- compartmentalized_interactions %>%
    dplyr::distinct(compartment) %>%
    dplyr::mutate(
      go_id = stringr::str_extract(compartment, "GO:[0-9]+"),
      c_name = stringr::str_extract(compartment, "^[a-zA-Z ]+"),
      c_name = stringr::str_trim(c_name)
    ) %>%
    # drop unusable compartment (rods & rings which has no GO ID)
    dplyr::filter(!is.na(go_id))

  compartmentalized_edges <- compartmentalized_interactions %>%
    dplyr::semi_join(compartments, by = "compartment")

  compartmenalized_genes <- dplyr::bind_rows(
    compartmentalized_edges %>%
      dplyr::distinct(ensembl_gene = protein1, compartment),
    compartmentalized_edges %>%
      dplyr::distinct(ensembl_gene = protein2, compartment)
  ) %>%
    dplyr::distinct(ensembl_gene, compartment)

  out <- list(
    compartmenalized_genes = compartmenalized_genes,
    compartmentalized_edges = compartmentalized_edges,
    compartments = compartments
  )

  return(out)
}

#' Trim a netcontextr graph by edge weight
#'
#' `trim_network_by_edge_weight` trims a netcontextr graph object to weight >=
#' "field_value" as specified in the "field_name" column of the "interactions"
#' element in the netcontextr graph. For string interaction data, weights
#' range from 0-1 (see `createStringGraph`).
#'
#' Confidence limits are as follows:
#' low confidence - 0.15 (or better) | or 150 or better in protein links files
#' medium confidence - 0.4 | or 400 or better in protein links files
#' high confidence - 0.7 | or 700 or better in protein links files
#' highest confidence - 0.9 | or 900 or better in protein links files
#'
#' @param netcontextr_graph a netcontextr graph object with "genes" and
#'  "interactions" elements
#' @param field_name the specific column in the "interactions" element of the
#' netcontextr graph to use for trimming interactions. This column should
#' contain some sort of numerical weight.
#' @param field_value the threshold for filtering (i.e. weight >= "field_value")
#'
#' @returns modified netcontextr graph object, with "interactions" filtered to
#'  those with weight > "field_value" and genes filtered to those contained
#'  within the remaining interactions
#'
#' @examples
#' # trim the edges of a netcontextr graph based on weight
#' high_confidence_string_graph <- trim_network_by_edge_weight(rcpr::example_string_graph,
#'   "combined_score", 0.7)
#' @export
trim_network_by_edge_weight <- function(netcontextr_graph, field_name, field_value){
  
  # Confirm valid input. Eventually should check that values in field_name match
  # the value type of field_value.
  checkmate::assertChoice(field_name, colnames(netcontextr_graph$interactions))
  checkmate::assertNumeric(field_value)
  
  # trim edges based on weight
  print("trimming interactions")
  netcontextr_graph$interactions <- dplyr::filter(netcontextr_graph$interactions, get({{field_name}}) >= {{field_value}})
  
  # trim network genes accordingly
  print("trimming gene list")
  genesLeft <- unique(c(netcontextr_graph$interactions$protein1, netcontextr_graph$interactions$protein2))
  netcontextr_graph$genes <- dplyr::filter(netcontextr_graph$genes, .data$gene %in% genesLeft)
  
  return(netcontextr_graph)
}


#' Trim a netcontextr graph based on gene metadata
#'
#' This function trims a netcontextr graph object to interactions where both
#' genes/nodes share the same value(s) in "field_name". Optionally, the user can
#' also specify a "field_value" to trim to interactions where both nodes have
#' that specific value in "field_name" only.
#'
#' @param netcontextr_graph a netcontextr graph object with "genes" and
#'  "interactions" elements
#' @param field_name the specific column in the "genes" element of the
#' netcontextr graph to use for trimming interactions. This column should be a
#' list of tibbles (i.e. as derived from `annotate_genes`)
#' @param field_value (optional) a specific value within "field_name" to filter
#'  interactions on. "field_value" can be a vector of length one or greater.
#'  If it is not specified (default) interactions are trimmed to those where
#'  both nodes have overlapping values in "field_name". If "field_value" is
#'  specified, then interactions are trimmed to those where both nodes have that
#'  specific "field_value" within "field_name". The latter is useful for
#'  trimming to a specific tissue (e.g. filter to interactions where both nodes
#'  are expressed in a given tissue or tissues).
#' @param addAttributeToInteractions (optional) a logical (default = FALSE)
#'  specifying whether to add field_name information as an interaction
#'  attribute (memory intensive) or, remove the field_name column so it does
#'  not become an interaction attribute. The default is to remove this column
#'  because all of this information is already stored in the genes part of the
#'  netcontextr graph object, and therefore, this information is redundant
#'  (and results in the netcontextr graph object becoming quite large).
#'
#' @return modified netcontextr graph object, with "interactions" filtered to
#'  those where both genes/nodes share the same value(s) in "field_name". If
#'  "field_value" is specified, the "interactions" will be filtered to only
#'  those where both nodes have "field_value(s)" in "field_column". "Genes" are
#'  filtered to those contained within the remaining interactions. If
#'  addAttributeToInteractions is TRUE, field_name information will be added as
#'  an interaction attribute.
#'
#' @examples
#' # annotate the genes in a netcontextr graph with subcellular localization(s)
#' # add this information to a new field called "compartment"
#' annotated_string_graph <- annotate_genes(
#'   rcpr::example_string_graph,
#'   rcpr::example_protein_subcellular_localizations,
#'   "compartment"
#'   )
#'   
#' # now trim this network to interactions where nodes are found in overlapping compartments
#' compartmentalized_string_graph <- trim_network_by_gene_attribute(
#'   annotated_string_graph,
#'   "compartment"
#'   )
#'
#' # Alternative example, where "field_value" is specified
#' # annotate the genes with expression data from GTEX
#' # first determine gene expression by tissue using zFPKM
#' gtex_expression_data_zFPKM <- gene_expression_by_tissue(rcpr::example_gtex_data)
#' # second annotate genes with tissue expression(s)
#' # add this information to a new field called "tissue"
#' annotated_string_graph <- annotate_genes(rcpr::example_string_graph,
#'   gtex_expression_data_zFPKM, "tissue")
#' # now trim network to interactions where both nodes are expressed in the liver
#' contextualized_string_graph <- trim_network_by_gene_attribute(
#'   annotated_string_graph,
#'   "tissue",
#'   "Liver"
#'   )
#' @export
trim_network_by_gene_attribute <- function(
  netcontextr_graph,
  field_name,
  field_value = NA,
  addAttributeToInteractions = FALSE
  ){
  
  # Confirm valid input. May eventually want to confirm that field_value exists,
  # though this is somewhat handled by assertInteger as part of filtering the
  # geneTable
  checkmate::assertChoice(field_name, colnames(netcontextr_graph$genes))
  
  # do joins to get interactions with overlapping values in field
  # select step removes other extraneous metadata so they don't get added
  # to interactions during join step
  cli::cli_alert_info("Trimming to interactions with nodes with overlapping {.field {field_name}}")
  geneTable <- netcontextr_graph$genes %>%
    tidyr::unnest(cols = c({{field_name}})) %>%
    dplyr::select(.data$gene, {{field_name}})
  
  # if field value is specified, filter on this value
  if(!is.na(field_value)){
    geneTable <- dplyr::filter(geneTable, get({{field_name}}) %in% {{field_value}})
    # Ensure that filtering doesn't reduce the number of genes to zero
    checkmate::assertInteger(nrow(geneTable), lower = 1)
  }
  
  # now do joins
  firstJoin <- netcontextr_graph$interactions %>%
    dplyr::inner_join(geneTable, by = c("protein1" = "gene"))
  secondJoin <- firstJoin %>%
    dplyr::semi_join(
      geneTable,
      by = c("protein2" = "gene", field_name),
      na_matches = "never"
      )
  
  # Add back to netcontextr graph object. Based on user input, either add
  #field_name information as an interaction attribute (memory intensive) or,
  # remove the field_name column so it does not become an interaction attribute.
  # NB the default is to remove this column because all of this information is
  # already stored in the genes part of the netcontextr graph object and so the
  # object starts becoming quite large.
  if(addAttributeToInteractions){
    netcontextr_graph$interactions <- secondJoin %>%
      tidyr::nest("{field_name}" := {{field_name}})
  } else {
    netcontextr_graph$interactions <- secondJoin %>%
      tidyr::nest("{field_name}" := {{field_name}}) %>% dplyr::select(-{{field_name}})
  }
  
  # trim network genes accordingly
  genesLeft <- unique(c(secondJoin$protein1, secondJoin$protein2))
  netcontextr_graph$genes <- netcontextr_graph$genes %>%
    dplyr::filter(.data$gene %in% genesLeft)
  
  return(netcontextr_graph)
}


#' Trim Network by Nodes
#' 
#' Trim a netcontextr graph object by specifying the nodes to keep
#'
#' @param graph netcontextr graph object with "genes" and "interactions"
#' @param nodesToKeep character vector specifying the IDs of nodes to keep
#' @param requireBothNodes logical specifying whether to trim to interactions
#'  where both nodesToKeep are present.
#'
#' @returns netcontextr graph object, with nodes and edges trimmed to the
#'  specified nodes
#'
#' @examples
#' # trim example string graph to first 1000 proteins with subcellular localization data
#' # set requireBothNodes = TRUE (default) in order to restrict to interactions
#' # where both nodes are present within the first 1000 proteins
#' trimmed_graph <- trim_network_by_nodes(
#'   rcpr::example_string_graph,
#'   rcpr::example_protein_subcellular_localizations$ensembl_gene_id[1:1000],
#'   requireBothNodes = TRUE
#'   )
#'
#' # trim example string graph to first 1000 proteins with subcellular localization data
#' # set requireBothNodes = FALSE in order to restrict to interactions where one
#' # or more nodes are present within the first 1000 proteins (this can be
#' # useful when elaborating the network of a key set of genes e.g. GWAS hits).
#' trimmed_graph <- trim_network_by_nodes(
#'   rcpr::example_string_graph,
#'   rcpr::example_protein_subcellular_localizations$ensembl_gene_id[1:1000],
#'   requireBothNodes = FALSE
#'   )
#' @export
trim_network_by_nodes <- function(graph, nodesToKeep, requireBothNodes = TRUE){
  
  if(requireBothNodes){
    # Filter genes, then filter edges based on these genes
    graph$genes <- dplyr::filter(graph$genes, graph$genes[["gene"]] %in% nodesToKeep)
    # Filter out all edges that don't contain remaining genes as both nodes
    # First create key for filtering
    filteringKey <- graph$interactions[["protein1"]] %in% graph$genes[["gene"]] & graph$interactions[["protein2"]] %in% graph$genes[["gene"]]
    # Now filter
    graph$interactions <- dplyr::filter(graph$interactions, filteringKey)
  } else {
    # Filter edges first, then use union of nodes to filter genes
    # Filter out all edges that don't contain remaining genes in one or both nodes
    # First create key for filtering edges
    filteringKey <- graph$interactions[["protein1"]] %in% nodesToKeep | graph$interactions[["protein2"]] %in% nodesToKeep
    # Now filter
    graph$interactions <- dplyr::filter(graph$interactions, filteringKey)
    # Finally, filter genes in netcontextr graph
    genesKept <- c(graph$interactions[["protein1"]], graph$interactions[["protein2"]]) %>% unique
    graph$genes <- dplyr::filter(graph$genes, graph$genes[["gene"]] %in% genesKept)
  }
  
  return(graph)
}

#' Trim Reactions By Gene Attribute
#'
#' @param netcontextr_reactions a netcontextr reactions object with "genes" and
#'  "reactions" elements
#' @inheritParams trim_network_by_gene_attribute
#'
#' @examples
#' netcontextr_reactions <- create_mock_reactions(100)
#' annotated_reactions <- annotate_genes(
#'   netcontextr_reactions,
#'   rcpr::example_protein_subcellular_localizations,
#'   "compartment"
#'   )
#'  
#' trim_reactions_by_gene_attribute(annotated_reactions, "compartment")
#' @export
trim_reactions_by_gene_attribute <- function(
  netcontextr_reactions,
  field_name,
  field_value = NA
){
  
  # validate netcontextr_reactions
  validate_netcontextr_reactions(netcontextr_reactions)
  
  checkmate::assertString(field_name)
  valid_contexts <- setdiff(colnames(netcontextr_reactions$genes), "gene")
  if (length(valid_contexts) == 0) {
    stop ("No contexts exist for genes please add a context first using annotate_genes()")
  }
  checkmate::assertChoice(field_name, valid_contexts)
  
  genes_w_contexts <- netcontextr_reactions$genes %>%
    tidyr::unnest(!!rlang::sym(field_name))
  
  if (!is.na(field_value)) {
    genes_w_contexts <- genes_w_contexts %>%
      dplyr::filter(!!rlang::sym(field_name) %in% field_value)
  }
  
  nested_reaction_species <- netcontextr_reactions$reactions %>%
    tidyr::nest(reaction_species = -reaction_id)
  
  nested_contextualized_reaction_species <- netcontextr_reactions$reactions %>%
    dplyr::left_join(genes_w_contexts, by = "gene", relationship = "many-to-many") %>%
    # create a nest of all reaction species existing in each compartment
    tidyr::nest(contextualized_reaction_species = -c(!!rlang::sym(field_name), reaction_id))
  
  contextualized_reactions <- nested_reaction_species %>%
    dplyr::left_join(
      nested_contextualized_reaction_species,
      by = "reaction_id"
    ) %>%
    # determine whether each reaction contains the reaction species
    # required to make it occur
    dplyr::mutate(valid_contextualized_rspecies = purrr::map2(
      reaction_species,
      contextualized_reaction_species,
      evaluate_reaction_context_sufficiency
    )) %>%
    # filter reactions which would not occur in the context
    dplyr::filter(!purrr::map_lgl(valid_contextualized_rspecies, is.null))
  
  # create output netcontextr reaction
  contextualized_netcontextr_reactions <- netcontextr_reactions
  netcontextr_reactions$reactions <- contextualized_reactions %>%
    dplyr::select(reaction_id, !!rlang::sym(field_name), valid_contextualized_rspecies) %>%
    tidyr::unnest(valid_contextualized_rspecies)
  
  netcontextr_reactions$genes <- netcontextr_reactions$genes %>%
    dplyr::filter(gene %in% netcontextr_reactions$reactions$gene)
  
  return(netcontextr_reactions)
}

evaluate_reaction_context_sufficiency <- function (
  reaction_species,
  contextualized_reaction_species
) {
  
  checkmate::assertDataFrame(reaction_species)
  checkmate::assertDataFrame(contextualized_reaction_species)
  checkmate::assertSubset(c("gene", "stoi", "role"), colnames(reaction_species))
  checkmate::assertSubset(c("gene", "stoi", "role"), colnames(contextualized_reaction_species))
  
  # compare reactions to controlled vocabulary of roles
  known_roles <- c("substrate", "product", "catalyst", "interactor", "activator", "inhibitor")
  reaction_species <- reaction_species %>%
    dplyr::mutate(role_fct = factor(role, levels = known_roles))
  
  if (any(is.na(reaction_species$role_fct))) {
    unknown_reaction_species <- reaction_species %>%
      dplyr::filter(is.na(role_fct)) %>%
      dplyr::mutate(glue_str = glue::glue("{gene} ({role})"))
    
    warning(glue::glue(
      "{nrow(unknown_reaction_species)} reaction species have unrecognized roles:
    {paste(unknown_reaction_species$glue_str, collapse = ', ')}
    known roles are {paste(known_roles, collapse = ', ')}
    "))
  }
  
  # determine whether the genes present in a given context are sufficient
  # to allow the reaction to occur in a context
  
  # some reactions may produce or consume complexes which if large may not
  # make sense to be annotated with with genic identifiers
  # because only the subset of species which are elemental genes will be
  # evaluated using this function we don't validate that all reactions contain
  # substrates & products thought this must be true
  
  if (any(reaction_species$role == "interactor")) {
    if (any(reaction_species$role != "interactor") || nrow(reaction_species) != 2) {
      warning("A reaction was defined which included at least one interactor,
  and was not defined using a pair of interactor. Interactors are
  generally used for compatibility with an 'edgelist' view of a pair
  of interacting molecules so an alternative role may be more appropriate")
    }
  }
  
  required_rspecies <- reaction_species %>%
    dplyr::mutate(
      requirement = dplyr::case_when(
        role_fct == "substrate" ~ "required", # otherwise the reaction can't occur
        role_fct == "product" ~ "required", # otherwise the reaction does not occur
        role_fct == "catalyst" ~ "1+ required if exists", # otherwise the reaction would not occur
        role_fct == "interactor" ~ "required", # for pairwise interactions all species are required
        TRUE ~ "optional" # other modifiers
      )) %>%
    dplyr::left_join(
      contextualized_reaction_species %>%
        dplyr::mutate(exists_in_context = TRUE),
      by = c("gene", "role")
    ) %>%
    dplyr::mutate(exists_in_context = ifelse(
      is.na(exists_in_context),
      FALSE,
      exists_in_context
    ))
  
  # determine whether the reaction occurs in the context
  missing_required_species <- required_rspecies %>%
    dplyr::filter(requirement == "required", !exists_in_context)
  
  missing_semi_required_species <- required_rspecies %>%
    dplyr::filter(requirement == "1+ required if exists") %>%
    dplyr::filter(all(exists_in_context == FALSE))
  
  if (nrow(missing_required_species) == 0 && nrow(missing_semi_required_species) == 0) {
    return(contextualized_reaction_species)
  } else {
    return(NULL)
  }
}

#' Create mock reactions
#'
#' Create a \code{netcontextr_reactions} list for use in
#' \code{trim_reactions_by_gene_attribute} by generating FAKE reactions from
#' the genes in STRING.
#'
#' @param n_reactions Number of fake reactions to create
#'
#' @returns a tibble of reactions and their reaction species
#'
#' @examples
#' mock_reactions <- create_mock_reactions(10)
#'
#' @export
create_mock_reactions <- function (n_reactions = 1000) {
  
  checkmate::assertIntegerish(n_reactions, len = 1)
  
  example_genes <- rcpr::example_string_graph$genes$gene
  
  reactions = tibble::tibble(reaction_id = seq(n_reactions)) %>%
    dplyr::mutate(reaction_members = purrr::map(reaction_id, ~ create_mock_reaction(example_genes))) %>%
    tidyr::unnest(reaction_members)
  
  mock_reactions <- list(
    genes = rcpr::example_string_graph$genes,
    reactions = reactions
  )
  
  return(mock_reactions)
}

create_mock_reaction <- function (example_genes) {
  
  checkmate::assertCharacter(example_genes, unique = TRUE)
  
  substrates <- sample(example_genes, pmax(stats::rpois(1, 2), 1)) # >= 1
  products <- sample(example_genes, pmax(stats::rpois(1, 2), 1)) # >= 1
  catalysts <- sample(example_genes, stats::rpois(1, 1)) # >= 0
  activators <- sample(example_genes, stats::rpois(1, 0.2)) # >= 0
  inhibitors <- sample(example_genes, stats::rpois(1, 0.2)) # >= 0
  
  mock_reaction_species <- tibble::tribble(
    ~ role, ~ gene, ~ stoi,
    "substrate", substrates, -1,
    "product", products, 1,
    "catalyst", catalysts, 0,
    "activator", activators, 0,
    "inhibitor", inhibitors, 0
  ) %>%
    tidyr::unnest(gene)
  
  return(mock_reaction_species)
}

validate_netcontextr_reactions <- function (netcontextr_reactions) {
  
  checkmate::assertList(netcontextr_reactions)
  required_names <- c("genes", "reactions")
  missing_required_names <- setdiff(required_names, names(netcontextr_reactions))
  if (length(missing_required_names) > 0) {
    stop (glue::glue("netcontextr_reactions must be a list containing named tables: {paste(required_names, collapse = ', ')}"))
  }
  
  checkmate::assertDataFrame(netcontextr_reactions$genes)
  required_gene_names <- c("gene")
  missing_required_names <- setdiff(required_gene_names, colnames(netcontextr_reactions$genes))
  if (length(missing_required_names) > 0) {
    stop (glue::glue("netcontextr_reactions's genes table is missing the required variables: {paste(missing_required_names, collapse = ', ')}"))
  }
  
  checkmate::assertDataFrame(netcontextr_reactions$reactions)
  required_gene_names <- c("gene", "reaction_id", "role")
  missing_required_names <- setdiff(required_gene_names, colnames(netcontextr_reactions$reactions))
  if (length(missing_required_names) > 0) {
    stop (glue::glue("netcontextr_reactions's reactions table is missing the required variables: {paste(missing_required_names, collapse = ', ')}"))
  }
  
  # check whether reaction species are represented among the genes
  
  unrecognized_genes <- setdiff(
    netcontextr_reactions$reactions$gene,
    netcontextr_reactions$genes$gene
  )
  
  if (length(unrecognized_genes) > 0) {
    warning(glue::glue(
      "{length(unrecognized_genes)} genes were present in the reactions but
    were missing from the genes tables. These reactions will be missing from
    all contexts where the missing gene are a required species (e.g., a substrate):
    Examples: {paste(unrecognized_genes[1:pmin(length(unrecognized_genes), 5)], collapse = ', ')}
    "))
  }
  
  return(invisible(0))
}

#' Discretize Expression
#'
#' Identify the relevant fields in \code{tall_expression_data} for applying \code{\link{tall_data_detect_unexpressed}} and then apply this function. This provides a flexible means of handling multiple datasets including those where expression context is jointly defined by multiple variables.
#'
#' @inheritParams tall_data_detect_unexpressed
#' @param valid_abundance_vars valid variables which can define the expression level; One and only one match is expected.
#' @param valid_feature_vars valid variables which can define a unique feature (e.g., gene); One and only one match is expected.
#' @param valid_context_vars valid variables which can define a unique expression context. 1+ matches are expected. If multiple matches are present then will be combined into a single context before applying \code{\link{tall_data_detect_unexpressed}}.
#'
#' @returns tall_expression_data with two variables added:
#' \describe{
#'   \item{zFPKM_score}{
#'   z-score of context \code{abundance_var} relative to other contexts.
#'   },
#'   \item{is_expressed}{
#'   Thresholded zFPKM scores indicating whether a gene is expressed (1) or
#'   not expressed (0) in this context.
#'   }
#' }
#'
#' @export
discretize_expression <- function (
  tall_expression_data,
  valid_abundance_vars,
  valid_feature_vars,
  valid_context_vars,
  zFPKM_cutoff = -3
) {
  
  checkmate::assertDataFrame(tall_expression_data)
  checkmate::assertCharacter(valid_abundance_vars)
  checkmate::assertCharacter(valid_feature_vars)
  checkmate::assertCharacter(valid_context_vars)
  checkmate::assertNumber(zFPKM_cutoff)
  
  # is a single abundance value present?
  abundance_var_matches <- intersect(colnames(tall_expression_data), valid_abundance_vars)
  n_abundance_var_matches <- length(abundance_var_matches)
  if (n_abundance_var_matches != 1) {
    cli::cli_abort(
      "{n_abundance_var_matches} valid abundance variables were found in {.emph tall_expression_data}. Exactly one match is required.\n
      {.strong tall_expression_data} contains the following vars: {.field {colnames(tall_expression_data)}}\n
      {.strong valid_abundance_vars} are: {.field {valid_abundance_vars}}"
    )
  }
  
  # were valid feature-defining variables found?
  feature_var_matches <- intersect(colnames(tall_expression_data), valid_feature_vars)
  n_feature_var_matches <- length(feature_var_matches)
  if (n_feature_var_matches == 0) {
    cli::cli_abort(
      "{n_feature_var_matches} valid feature variables were found in {.emph tall_expression_data}. One or more matches are required.\n
      {.strong tall_expression_data} contains the following vars: {.field {colnames(tall_expression_data)}}\n
      {.strong valid_feature_vars} are: {.field {valid_feature_vars}}"
    )
  }
  
  if (n_feature_var_matches > 1) {
    # provide an alert if multiple feature-defining variables are found
    cli::cli_alert_info(
      "{n_feature_var_matches} valid feature-defining variables were found in {.emph tall_expression_data}. {.field {feature_var_matches}} will be combined to create a combined feature definition."
    )
    
    tall_expression_data <- tall_expression_data %>%
      tidyr::unite(".combined_feature", !!!rlang::syms(feature_var_matches), remove = FALSE)
    
    feature_var_matches <- ".combined_feature"
  }
  
  # were valid contexts found?
  context_var_matches <- intersect(colnames(tall_expression_data), valid_context_vars)
  n_context_var_matches <- length(context_var_matches)
  if (n_context_var_matches == 0) {
    cli::cli_abort(
      "{n_context_var_matches} valid context variables were found in {.emph tall_expression_data}. One or more matches are required.\n
      {.strong tall_expression_data} contains the following vars: {.field {colnames(tall_expression_data)}}\n
      {.strong valid_context_vars} are: {.field {valid_context_vars}}"
    )
  }
  
  if (n_context_var_matches > 1) {
    # provide an alert if multiple contexts are found
    cli::cli_alert_info(
      "{n_context_var_matches} valid context variables were found in {.emph tall_expression_data}. {.field {context_var_matches}} will be combined to create a combined context."
    )
    
    tall_expression_data <- tall_expression_data %>%
      tidyr::unite(".combined_context", !!!rlang::syms(context_var_matches), remove = FALSE)
    
    context_var_matches <- ".combined_context"
  }
  
  # add zFPKM_score and is_expressed
  tall_thresholded_expression <- tall_data_detect_unexpressed(
    tall_expression_data,
    feature_var = feature_var_matches,
    context_var = context_var_matches,
    abundance_var = abundance_var_matches,
    zFPKM_cutoff = zFPKM_cutoff
  )
  
  # cleanup combined variables
  if (feature_var_matches == ".combined_feature") {
    tall_thresholded_expression <- tall_thresholded_expression %>%
      dplyr::select(-.combined_feature)
  }
  
  if (context_var_matches == ".combined_context") {
    tall_thresholded_expression <- tall_thresholded_expression %>%
      dplyr::select(-.combined_context)
  }
  
  return (tall_thresholded_expression)
}

#' Filter Unexpressed
#'
#' Filter a tabular table of feature (rows) by contexts (columns) using
#'   zFPKM::zFPKM to generate context-specificity scores.
#'
#' @param expression_data A tabular dataset with rows defining distinct features
#'   and columns as samples or contexts.
#' @param feature_defining_var A variable in \code{expression_data} defining
#'   a distinct feature
#' @param invalid_vars Variables in \code{expression_data} which should be dropped
#'   before applying zFPKM because they are either non-numeric or invalid contexts.
#' @param zFPKM_cutoff zFPKM values about this cutoff will be considered as
#'   meaningful expression.
#'
#' @returns A list containing:
#' \describe{
#'   \item{expression_data_df}{\code{expression_data} formatted as a data.frame.}
#'   \item{zFPKM_scores}{\code{zFPKM_scores} \code{expression_data_df} with zFPKM applied.}
#'   \item{is_expressed}{Boolean data.frame of \code{zFPKM_scores} thresholded by \code{zFPKM_cutoff}.}
#' }
detect_unexpressed <- function (
  expression_data,
  feature_defining_var,
  invalid_vars = NULL,
  zFPKM_cutoff = -3
) {
  
  # validate inputs
  checkmate::assertDataFrame(expression_data)
  checkmate::assertChoice(feature_defining_var, colnames(expression_data))
  # if invalid_vars are provided they must be present in expression_data
  purrr::walk(invalid_vars, checkmate::assertChoice, colnames(expression_data))
  checkmate::assertNumber(zFPKM_cutoff)
  
  expression_data_df <- expression_data %>%
    dplyr::select(!dplyr::any_of(invalid_vars)) %>%
    # convert to data.frame with feature_defining_var creating rownames
    tibble::column_to_rownames(var = feature_defining_var)
  
  # calculate zFPKM scores for each gene in each tissue
  zFPKM_scores <- zFPKM::zFPKM(expression_data_df)
  
  # trim to gene expressed in tissue of interest only
  # defined by zFPKM >= -3 and median expression > 0
  is_expressed <- zFPKM_scores >= zFPKM_cutoff & expression_data_df > 0
  
  stopifnot(rownames(expression_data_df) == rownames(zFPKM_scores))
  stopifnot(colnames(expression_data_df) == colnames(zFPKM_scores))
  stopifnot(rownames(expression_data_df) == rownames(is_expressed))
  stopifnot(colnames(expression_data_df) == colnames(is_expressed))
  
  out <- list(
    expression_data_df = expression_data_df,
    zFPKM_scores = zFPKM_scores,
    is_expressed = is_expressed
  )
  
  return(out)
}

#' Tall Data Filter Unexpressed
#'
#' Format a "tall" table of features, contexts, and expression intensity and
#'   identify contexts where each feature is expressed based on zFPKM.
#'
#' @param tall_expression_data Table with one row per feature x context pair.
#' @param feature_var Variable describing a unique feature (e.g., a gene ID).
#' @param context_var Variable describing a unique context (e.g., a tissue or sample).
#' @param abundance_var Expression level which should be thresholded.
#' @inheritParams detect_unexpressed
#'
#' @returns tall_expression_data with two variables added:
#' \describe{
#'   \item{zFPKM_score}{
#'   Z-score of context \code{abundance_var} relative to other contexts.
#'   }
#'   \item{is_expressed}{
#'   Thresholded zFPKM scores indicating whether a gene is expressed (1) or
#'   not expressed (0) in this context.
#'   }
#' }
#'
#' @examples
#' if (interactive()) {
#'   gtex_file <- "data-raw/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
#'   gtex_expression_data <- load_and_clean_gtex_data(gtex_file)
#'
#'   gtex_tall_expression_data <- gtex_expression_data %>%
#'     dplyr::slice(1:200) %>%
#'     tidyr::gather(tissue, abundance, -c(ensembl_gene_id:Description))
#'
#'   tall_data_detect_unexpressed(
#'     gtex_tall_expression_data, "ensembl_gene_id", "tissue", "abundance"
#'   )
#' }
#' @export
tall_data_detect_unexpressed <- function (
  tall_expression_data,
  feature_var,
  context_var,
  abundance_var,
  zFPKM_cutoff = -3
) {
  
  checkmate::assertDataFrame(tall_expression_data)
  checkmate::assertChoice(feature_var, colnames(tall_expression_data))
  checkmate::assertChoice(context_var, colnames(tall_expression_data))
  checkmate::assertFALSE(feature_var == context_var)
  checkmate::assertChoice(abundance_var, colnames(tall_expression_data))
  checkmate::assertTRUE(any(c("character", "factor") %in% class(tall_expression_data[[feature_var]])))
  checkmate::assertTRUE(any(c("character", "factor") %in% class(tall_expression_data[[context_var]])))
  checkmate::assertTRUE(any(c("numeric", "integer") %in% class(tall_expression_data[[abundance_var]])))
  
  redundant_measurements <- tall_expression_data %>%
    dplyr::group_by(!!rlang::sym(feature_var), !!rlang::sym(context_var)) %>%
    dplyr::filter(dplyr::n() > 1)
  
  if (nrow(redundant_measurements) != 0) {
    
    example_redundant_measurements <- redundant_measurements %>%
      dplyr::arrange(!!rlang::sym(feature_var), !!rlang::sym(context_var)) %>%
      dplyr::ungroup() %>%
      dplyr::slice(1:pmin(10, nrow(redundant_measurements)))
    
    example_redundant_measurements %>%
      knitr::kable() %>%
      kableExtra::kable_material()
    
    cli::cli_abort(
      "{.val {nrow(redundant_measurements)}} features included multiple pairs with
      the same combination of {.val {feature_var}} and {.val {context_var}}."
    )
  }
  
  wide_expression_data <- tall_expression_data %>%
    dplyr::select(dplyr::any_of(c(feature_var, context_var, abundance_var))) %>%
    tidyr::spread(context_var, abundance_var)
  
  # remove features which are undefined in some contexts (NAs)
  filtered_wide_expression_data <- tall_data_detect_unexpressed_filter_missing_values(
    wide_expression_data
  )
  
  thresholded_expression_data <- detect_unexpressed(
    filtered_wide_expression_data,
    feature_var,
    zFPKM_cutoff = zFPKM_cutoff
  )
  
  thresholded_expression_data <- thresholded_expression_data[
    names(thresholded_expression_data) %in% c("zFPKM_scores", "is_expressed")
  ]
  
  tall_thresholded_expression <- tibble::tibble(
    summary_type = names(thresholded_expression_data),
    dat = thresholded_expression_data
  ) %>%
    dplyr::mutate(
      dat = purrr::map(dat, as.data.frame),
      tall_dat = purrr::map(dat,
                            ~ .x %>%
                              tibble::rownames_to_column(var = feature_var) %>%
                              tibble::as_tibble() %>%
                              tidyr::gather(
                                !!rlang::sym(context_var),
                                "summary_var",
                                -!!rlang::sym(feature_var)
                              )
      )) %>%
    dplyr::select(-dat) %>%
    tidyr::unnest(tall_dat) %>%
    tidyr::spread(summary_type, summary_var)
  
  # add back features with missing values
  
  # drop entries from original data that were removed due to missing values
  valid_tall_expression_data <- tall_expression_data %>%
    dplyr::semi_join(filtered_wide_expression_data, by = feature_var)
  
  if (nrow(valid_tall_expression_data) != nrow(tall_thresholded_expression)) {
    cli::cli_abort(
      "Thresholded expression has {.val {nrow(tall_thresholded_expression)}} rows,
      but input had {.val {nrow(valid_tall_expression_data)}} rows")
  }
  
  out <- tall_expression_data %>%
    dplyr::left_join(tall_thresholded_expression, by = c(feature_var, context_var))
  
  return(out)
}

tall_data_detect_unexpressed_filter_missing_values <- function (wide_expression_data) {
  
  # check for missing values
  missing_n_vec <- rowSums(is.na(wide_expression_data))
  
  if (any(missing_n_vec != 0)) {
    missing_value_counts <- missing_n_vec %>%
      table() %>%
      tibble::as_tibble() %>%
      dplyr::rename(n_missing = ".") %>%
      dplyr::arrange(n_missing) %>%
      dplyr::mutate(report_val = glue::glue("{n_missing} [{n}]"))
    
    cli::cli_alert_warning(
      "Among {.strong {ncol(wide_expression_data)}} distinct expression contexts of {.strong {nrow(wide_expression_data)}} features
      {sum(missing_n_vec != 0)} features are missing in 1+ contexts and will be removed from consideration.
      Counts of missing values by feature are as follows {.field {missing_value_counts$report_val}}"
    )
  }
  
  filtered_wide_expression_data <- wide_expression_data[missing_n_vec == 0,]
  
  return(filtered_wide_expression_data)
}
