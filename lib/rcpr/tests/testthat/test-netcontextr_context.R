test_that("Test potential reaction contexts in evaluate_reaction_context_sufficiency()", {
  
  reaction_w_regulators <- tibble::tribble(
    ~ gene, ~ stoi, ~ role,
    "A", -1, "substrate",
    "B", -1, "substrate",
    "C", 1, "product",
    "D", 1, "product",
    "E", 0, "catalyst",
    "F", 0, "catalyst",
    "G", 0, "activator"
  )
  
  # missing substrates -> inoperable
  testthat::expect_null(
    rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(-1)
    ))
  
  # missing one enzyme -> operable
  testthat::expect_equal(
    nrow(rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(-5))
    ), 6)
  
  # missing one product -> inoperable
  testthat::expect_null(
    rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(-3)
    ))
  
  # missing all enzymes -> inoperable
  testthat::expect_null(
    rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(1:4)
    ))
  
  # missing regulators -> operable
  testthat::expect_equal(
    nrow(rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(-7))
    ), 6)
  
  weird_reaction_w_interactors <- tibble::tribble(
    ~ gene, ~ stoi, ~ role,
    "A", -1, "interactor",
    "B", 1, "interactor",
    "C", 1, "interactor"
  )
  
  # throw a warning for a weird interaction format
  testthat::expect_warning(
    rcpr:::evaluate_reaction_context_sufficiency(
      weird_reaction_w_interactors,
      weird_reaction_w_interactors
    ), regexp = "interactor")
  
  reaction_w_interactors <- tibble::tribble(
    ~ gene, ~ stoi, ~ role,
    "A", -1, "interactor",
    "B", 1, "interactor",
  )
  
  testthat::expect_null(
    rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_regulators,
      reaction_w_regulators %>% dplyr::slice(1:4)
    ))
  
  # missing one interactor
  testthat::expect_null(
    rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_interactors,
      reaction_w_interactors %>% dplyr::slice(-1)
    ))
  
  # identify function looks good
  testthat::expect_equal(
    nrow(rcpr:::evaluate_reaction_context_sufficiency(
      reaction_w_interactors,
      reaction_w_interactors
    )), 2)

})

test_that("trim_reactions_by_gene_attribute() and trim_network_by_gene_attribute() produce the same results for interactions", {
  
  annotated_string_graph <- annotate_genes(
    example_string_graph,
    example_protein_subcellular_localizations,
    "compartment"
  )
  
  # filter to a subset of interaction for a quick eval
  
  set.seed(1234)
  annotated_string_graph$interactions <- annotated_string_graph$interactions %>%
    dplyr::sample_n(100)
  
  compartmentalized_string_graph <- trim_network_by_gene_attribute(
    annotated_string_graph,
    "compartment",
    addAttributeToInteractions = TRUE
  )
  
  # convert the string graph to a reaction format
  annotated_string_graph_reactions <- annotated_string_graph
  annotated_string_graph_reactions$interactions <- NULL
  annotated_string_graph_reactions$reactions <- annotated_string_graph$interactions %>%
    dplyr::select(protein1, protein2) %>%
    dplyr::mutate(reaction_id = 1:dplyr::n()) %>%
    tidyr::gather(gene_partner, gene, -reaction_id) %>%
    dplyr::mutate(stoi = ifelse(gene_partner == "protein2", 1, -1)) %>%
    dplyr::select(-gene_partner) %>%
    dplyr::mutate(role = "interactor")
  
  compartmentalized_string_graph_reactions <- trim_reactions_by_gene_attribute(
    annotated_string_graph_reactions,
    field_name = "compartment"
  )
  
  spread_reactions <- compartmentalized_string_graph_reactions$reactions %>%
    dplyr::distinct(reaction_id, gene, stoi) %>%
    tidyr::spread(stoi, gene)
  
  testthat::expect_equal(
    nrow(spread_reactions),
    nrow(compartmentalized_string_graph$interactions)
  )
})