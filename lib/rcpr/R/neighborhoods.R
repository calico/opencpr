#' Summarize Neighborhoods
#'
#' @param neighborhood_table Neighborhood table created by
#'   \code{create_neighborhood_table}
#' @inheritParams create_neighborhood_table
#' @inheritParams plot_one_neighborhood
#' @param ... extra arguments passed to \code{\link{plot_one_neighborhood}}
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper(TRUE)
#'
#'   species_id <- rownames(consensus_model$species)[
#'     consensus_model$species$s_name == "cAMP"
#'     ]
#'   species_id <- sample(rownames(consensus_model$species), 1)
#'
#'   neighborhood_table <- create_neighborhood_table(
#'     species_id = species_id,
#'     consensus_model = consensus_model,
#'     consensus_graph = consensus_graph,
#'     cpr = cpr,
#'     max_steps = 3L,
#'     max_neighbors = 40L,
#'   )
#'
#'  # score_overlay <- summarize_indication(
#'  #   disease_id = "EFO_0000400",
#'  #   consensus_model,
#'  #   create_neighborhood_summary_table(neighborhood_table),
#'  #   species_identifiers
#'  # )
#'
#'  score_overlay <- neighborhood_table %>%
#'    dplyr::select(vertices) %>%
#'    tidyr::unnest(vertices) %>%
#'    dplyr::filter(node_type == "species") %>%
#'    dplyr::distinct(s_id) %>%
#'    dplyr::sample_frac(0.5) %>%
#'    dplyr::mutate(score = rnorm(n()))
#'
#'   neighborhood_summaries <- plot_neighborhoods(
#'     neighborhood_table,
#'     consensus_graph,
#'     score_overlay,
#'     score_label = "diabetes mellitus",
#'     score_palette = "log2 fold-change"
#'   )
#' }
#' @export
plot_neighborhoods <- function(
  neighborhood_table,
  consensus_graph,
  score_overlay = NULL,
  score_label = NULL,
  ...
  ) {

  checkmate::assertDataFrame(neighborhood_table)
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertDataFrame(score_overlay, null.ok = TRUE)
  checkmate::assertString(score_label, null.ok = TRUE)

  debugr::dwatch(msg = "Starting plot_one_neighborhood [cpr<neighborhoods.R>::plot_neighborhoods]")

  neighborhood_table <- neighborhood_table %>%
    # plot each neighborhood
    dplyr::mutate(neighborhood_grob = purrr::pmap(
      list(vertices, edges, edge_sources, sc_id, sc_name),
      plot_one_neighborhood,
      consensus_graph = consensus_graph,
      score_overlay = score_overlay,
      score_label = score_label,
      ...
    ))

  debugr::dwatch(msg = "Aggregating neighborhood plots [cpr<neighborhoods.R>::plot_neighborhoods]")

  if (nrow(neighborhood_table) > 4) {
    message(glue::glue(
      "{neighborhood_table$s_id[1]} was present in {nrow(neighborhood_table)} compartments, the 4 largest compartment-specific neighborhoods will be shown, excluding
      {paste(neighborhood_table$sc_name[5:nrow(neighborhood_table)], collapse = ', ')}"
    ))
    n_plots <- 4
  } else {
    n_plots <- nrow(neighborhood_table)
  }

  grid_params <- list(ncol = 1)
  arranged_neighborhood_plots <- try(
    do.call(
      gridExtra::grid.arrange,
      append(neighborhood_table$neighborhood_grob[1:n_plots], grid_params)
    ),
    silent = TRUE
  )

  if ("try-error" %in% class(arranged_neighborhood_plots)) {
    arranged_neighborhood_plots <- ggplot(
      data.frame(x = 0, y = 0), aes(x = x, y = y)
    ) +
      geom_text(
        label = glue::glue("Plotting failed :("),
        size = 10
      ) +
      theme(text = element_blank(), line = element_blank())
  }

  return(arranged_neighborhood_plots)
}

#' Create Neighborhood Table
#'
#' @param species_id species identifier for focal node
#' @param consensus_model network reconstruction
#' @param consensus_graph igraph model of network
#' @param precomputed_distances A table of precomputed distances between pairs of
#'   molecular species. Or, NULL to calculate this on-the-fly (default).
#' @param network_type what type of neighborhood should be formed (ignored
#'   if \code{consensus_graph} is undirected).
#'   \describe{
#'     \item{downstream}{descendants of the focal node}
#'     \item{upstream}{ancestors of the focal node}
#'     \item{hourglass}{descendants and ancestors of focal node}
#'   }
#' @param max_steps number of steps away from focal node allowed
#' @param max_neighbors prune to this number of upstream regulators and
#'   downstream targets
#' @param cpr binding to cpr python library
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'   species_id <- "S00000061"
#'
#'   create_neighborhood_table(
#'     species_id,
#'     consensus_model,
#'     consensus_graph,
#'     network_type = "hourglass",
#'     cpr = cpr
#'   )
#'
#'  # calculating distances on-the-fly
#'   create_neighborhood_table(
#'     species_id,
#'     consensus_model,
#'     consensus_graph,
#'     precomputed_distances = precomputed_distances,
#'     network_type = "hourglass",
#'     cpr = cpr
#'   )
#' }
#' @export
create_neighborhood_table <- function(
  species_id,
  consensus_model,
  consensus_graph,
  precomputed_distances = NULL,
  network_type = "downstream",
  max_steps = 3L,
  max_neighbors = 10L,
  cpr
  ) {

  checkmate::assertCharacter(species_id, len = 1)
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertChoice(network_type, c("downstream", "upstream", "hourglass"))
  checkmate::assertInteger(max_steps, len = 1, lower = 1)
  checkmate::assertInteger(max_neighbors, len = 1, lower = 1)
  checkmate::assertClass(cpr, "python.builtin.module")

  debugr::dwatch(msg = "Starting create_neighborhood_table [cpr<neighborhoods.R>::create_neighborhood_table]")

  compartmentalized_species <- cpr$network$net_utils$compartmentalize_species(consensus_model, species_id)
  compartmentalized_species_ids <- compartmentalized_species[["sc_id"]]

  compartmentalized_species_attrs <- consensus_model$compartmentalized_species %>%
    dplyr::filter(rownames(.) %in% compartmentalized_species_ids) %>%
    dplyr::mutate(
      sc_id = rownames(.),
      sc_id = factor(sc_id, levels = compartmentalized_species_ids)
    ) %>%
    dplyr::as_tibble() %>%
    dplyr::arrange(sc_id)

  if (
    nrow(compartmentalized_species_attrs) !=
      length(compartmentalized_species_ids)
  ) {
    cli::cli_abort(
      "compartmentalized species attributes were not aligned to IDs,
      this is unexpected behavior"
    )
  }

  debugr::dwatch(
    msg = "Species -> compartmentalized species finished [cpr<neighborhoods.R>::create_neighborhood_table]"
  )

  cli::cli_alert_info(
    "Neighborhoods will be created for {.field {length(compartmentalized_species_ids)}} compartment{?s}: {compartmentalized_species_ids}"
  )

  neighborhoods <- cpr$network$neighborhoods$find_and_prune_neighborhoods(
    consensus_model,
    consensus_graph,
    compartmentalized_species = compartmentalized_species_ids,
    precomputed_distances = precomputed_distances,
    order = max_steps,
    network_type = network_type,
    top_n = max_neighbors
  )

  debugr::dwatch(msg = "Extracting neighborhood attributes [cpr<neighborhoods.R>::create_neighborhood_table]")

  neighborhood_table <- compartmentalized_species_attrs %>%
    dplyr::mutate(
      neighborhoods = neighborhoods,
      vertices = purrr::map(
        neighborhoods,
        extract_neighborhood_df,
        "vertices",
        is_null_valid = FALSE
      ),
      edges = purrr::map(
        neighborhoods,
        extract_neighborhood_df,
        "edges",
        is_null_valid = FALSE
      ),
      edge_sources = purrr::map(
        neighborhoods,
        extract_neighborhood_df,
        "edge_sources"
      )
    ) %>%
    dplyr::select(-neighborhoods) %>%
    # sort by neighborhood size
    dplyr::arrange(dplyr::desc(purrr::map_int(vertices, nrow)))

  debugr::dwatch(msg = "Completed create_neighborhood_table [cpr<neighborhoods.R>::create_neighborhood_table]")

  return(neighborhood_table)
}

create_neighborhood_summary_table <- function(neighborhood_table) {
  neighborhood_summary_table <- neighborhood_table %>%
    dplyr::select(sc_name, vertices) %>%
    tidyr::unnest(vertices) %>%
    dplyr::rename(
      neighbor_sc_id = name,
      neighbor_sc_name = node_name
    )

  return(neighborhood_summary_table)
}


extract_neighborhood_df <- function(entry_list, entry_field, is_null_valid = TRUE) {
  chosen_table <- entry_list[[entry_field]]

  if (is.null(chosen_table)) {
    if (is_null_valid) {
      return(NULL)
    } else {
      stop("entry_list was NULL and is_null_valid was TRUE")
    }
  }

  reticulate::py_to_r(chosen_table) %>%
    dplyr::as_tibble()
}

#' Plot One Neighborhood
#'
#' @param vertices table of species and reactions, produced by \link{create_neighborhood_table}
#' @param edges table of connections between species and reactions, produced by \link{create_neighborhood_table}
#' @param edge_sources table describing the model(s) each reaction comes from, produced by \link{create_neighborhood_table}
#' @param sc_id compartmentalized species identifier of focal node
#' @param sc_name name of focal node
#' @inheritParams create_neighborhood_table
#' @param score_overlay optional, scores of neighbors for an indication from \link{summarize_indication}.
#' @param score_label optional, name of disease being overlaid.
#' @param score_palette optional, color palette for scores.
#' @param max_labeled_species maximum number of species to label (to avoid overplotting)
#' @param network_layout method to used for creating a network layout (e.g., fr, kk, drl)
#' @param edge_width width of edges on graph
#'
#' @returns a ggplot2 grob
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper(TRUE)
#'
#'   #species_id <- rownames(consensus_model$species)[
#'   #  consensus_model$species$s_name == "pyruvate kinase tetramer"
#'   #  ]
#'   species_id <- sample(rownames(consensus_model$species), 1)
#'
#'   neighborhood_table <- create_neighborhood_table(
#'     species_id,
#'     consensus_model,
#'     consensus_graph,
#'     network_type = "hourglass",
#'     cpr = cpr,
#'     max_neighbors = 30L,
#'     max_steps = 15L
#'   )
#'   entry <- 1
#'   vertices <- neighborhood_table$vertices[[entry]]
#'   edges <- neighborhood_table$edges[[entry]]
#'   edge_sources <- neighborhood_table$edge_sources[[entry]]
#'   sc_id <- neighborhood_table$sc_id[entry]
#'   sc_name <- neighborhood_table$sc_name[entry]
#'
#'   score_overlay <- vertices %>%
#'     dplyr::filter(node_type == "species") %>%
#'     dplyr::distinct(s_id) %>%
#'     dplyr::sample_frac(0.5) %>%
#'     dplyr::mutate(score = rnorm(n()))
#'
#'   # score_overlay <- summarize_indication(
#'   #   disease_id = "EFO_0000400",
#'   #   consensus_model,
#'   #   create_neighborhood_summary_table(neighborhood_table),
#'   #   species_identifiers
#'   #   )
#'
#'   plot_one_neighborhood(
#'     vertices,
#'     edges,
#'     edge_sources,
#'     sc_id,
#'     sc_name,
#'     consensus_graph,
#'     score_overlay = NULL
#'   )
#'
#'   plot_one_neighborhood(
#'     vertices,
#'     edges,
#'     edge_sources,
#'     sc_id,
#'     sc_name,
#'     consensus_graph,
#'     score_overlay = score_overlay,
#'     score_palette = "log2 fold-change"
#'   )
#' }
#' @export
plot_one_neighborhood <- function(
  vertices,
  edges,
  edge_sources,
  sc_id,
  sc_name,
  consensus_graph,
  score_overlay = NULL,
  score_label = NULL,
  score_palette = NULL,
  max_labeled_species = 30L,
  network_layout = "fr",
  edge_width = 0.1
  ) {

  checkmate::assertDataFrame(vertices)
  checkmate::assertDataFrame(edges)
  stopifnot(class(sc_id) %in% c("factor", "character"), length(sc_id) == 1)
  stopifnot(class(sc_name) %in% c("factor", "character"), length(sc_name) == 1)
  checkmate::assertDataFrame(score_overlay, null.ok = TRUE)
  checkmate::assertString(score_label, null.ok = TRUE)
  checkmate::assertString(score_palette, null.ok = TRUE)
  checkmate::assertInteger(max_labeled_species, len = 1, min = 1)
  checkmate::assertString(network_layout)
  checkmate::assertNumeric(edge_width, len = 1)

  debugr::dwatch(msg = "Starting plot_one_neighborhood [cpr<utils.R>::plot_one_neighborhood]")

  if (nrow(edges) == 0) {
    invalid_plot <- ggplot(
      data.frame(x = 0, y = 0), aes(x = x, y = y)
    ) +
      geom_text(label = glue::glue("{sc_name} has zero neighbors"), size = 10) +
      theme(text = element_blank(), line = element_blank())
    return(invalid_plot)
  }

  # add labels
  vertices <- plot_one_neighborhood_label_vertices(
    vertices %>% dplyr::arrange(path_weight),
    max_labeled_species
  )

  # add scores if they are present
  if (!is.null(score_overlay)) {
    checkmate::assertDataFrame(score_overlay)
    stopifnot(c("s_id", "score") %in% colnames(score_overlay))

    vertices <- vertices %>%
      dplyr::left_join(
        score_overlay %>%
          dplyr::select(s_id, score),
        by = "s_id"
      )
  } else {
    # ignore label even if one was passed
    score_label <- NULL
  }

  # add pathway sources to help organize layout
  if (!("NULL" %in% class(edge_sources))) {
    edges <- edges %>%
      dplyr::bind_rows(
        edge_sources %>%
          dplyr::select(from = "r_id", to = "pathway_id")
      )

    vertices <- vertices %>%
      dplyr::bind_rows(edge_sources %>% dplyr::distinct(name = pathway_id))
  }

  steps_colors <- vertices %>%
    dplyr::distinct(path_length) %>%
    dplyr::filter(!is.na(path_length)) %>%
    dplyr::arrange(path_length) %>%
    dplyr::mutate(path_color = grDevices::colorRampPalette(c("white", "dodgerblue"))(n()))

  neighborhood_network <- igraph::graph_from_data_frame(
    edges,
    directed = consensus_graph$is_directed(),
    vertices = vertices
  )

  plot_one_neighborhood_v2(
    neighborhood_network = neighborhood_network,
    edge_sources = edge_sources,
    score_label = score_label,
    score_palette = score_palette,
    network_layout = network_layout,
    edge_width = edge_width
  )
}

plot_one_neighborhood_v2 <- function(
  neighborhood_network,
  edge_sources,
  score_label = NULL,
  score_palette = NULL,
  network_layout = "fr",
  edge_width = 0.1
  ) {

  gg_network_layout <- layout_with_reaction_sources(neighborhood_network, network_layout)
  neighborhood_grob <- ggraph::ggraph(graph = gg_network_layout)

  # find obscured labels due to overplotting
  obscured_labels <- find_obscured_labels(gg_network_layout)

  graph_height <- diff(range(gg_network_layout$y))
  vertices_df <- tibble::as_tibble(gg_network_layout) %>%
    dplyr::mutate(
      ggtext_label = ifelse(ggtext_label %in% obscured_labels, NA, ggtext_label),
      n_lines = stringr::str_count(ggtext_label, "<br>"),
      label_y = y + ((n_lines + 3) * graph_height / 120)
    )

  focal_species <- vertices_df %>%
    dplyr::filter(node_orientation == "focal")
  stopifnot(nrow(focal_species) == 1)
  neighbors <- vertices_df %>%
    dplyr::filter(node_orientation != "focal")

  # title setup
  n_hops <- max(vertices_df$path_length)
  focal_species_name <- stringr::str_replace(focal_species["species_name"][[1]], "<br>", " ")
  focal_compartment <- stringr::str_trim(focal_species[["compartment"]])
  compartment_label <- ifelse(focal_compartment == "", "", glue::glue(" ({focal_compartment})"))
  plot_title <- glue::glue(
    "**{n_hops}** hops neighborhood of<br>
     **{focal_species_name}**{compartment_label}"
  )

  if (!is.null(score_label)) {
    plot_title <- plot_title + glue::glue("<br>overlaying **{score_label}**")
  }

  if (!("NULL" %in% class(edge_sources))) {
    # add reaction sources if available
    pathway_coords <- layout_pathway_sources(gg_network_layout, edge_sources)

    neighborhood_grob <- neighborhood_grob +
      geom_rect(
        data = pathway_coords,
        aes(
          xmin = x_min - 0.05,
          xmax = x_max + 0.05,
          ymin = y_min - 0.05,
          ymax = y_max + 0.05,
          fill = label_wrap
        ),
        alpha = 0.2
      ) +
      scale_fill_viridis_d("Pathway")
  }

  if ("score" %in% colnames(vertices_df)) {
    color_by <- rlang::sym("score")

    # add the score palette
    score_palette_obj <- select_score_overlay_palette(score_palette)
    neighborhood_grob <- neighborhood_grob + score_palette_obj
  } else {
    color_by <- rlang::quo(factor(node_orientation))

    neighborhood_grob <- neighborhood_grob +
      scale_color_manual(values = c(
        "upstream" = "gray25",
        "downstream" = "gray25",
        "focal" = "white"
      )) +
      guides(color = "none")
  }

  neighborhood_grob <- neighborhood_grob +
    # add edges - irreversible
    ggraph::geom_edge_link(
      data = ggraph_get_edges_by_reversibility(FALSE),
      color = "gray25",
      arrow = grid::arrow(type = "closed", length = unit(0.15, "inches")),
      edge_width = edge_width
    ) +
    # add edges - reversible
    ggraph::geom_edge_link(
      data = ggraph_get_edges_by_reversibility(TRUE),
      color = "gray25",
      edge_width = edge_width
    ) +
    # add focal node border
    ggraph::geom_node_point(
      data = focal_species,
      aes(x = x, y = y),
      size = 6.5,
      color = "black"
    ) +
    # add focal nodes
    ggraph::geom_node_point(aes(
      shape = factor(node_type),
      color = !!color_by,
      size = 6
    )) +
    # add number of steps
    geom_text(
      data = focal_species,
      aes(x = x, y = y, label = path_length),
      color = "gray10"
    ) +
    geom_text(
      data = neighbors %>% dplyr::filter(path_length <= 9),
      aes(
        x = x,
        y = y,
        label = path_length
      ),
      color = "gray70"
    ) +
    # add node labels
    ggtext::geom_richtext(
      data = vertices_df %>%
        dplyr::filter(!is.na(ggtext_label)),
      aes(x = x, y = label_y, label = ggtext_label),
      size = 3
    ) +
    scale_shape_manual("Type", values = c("reaction" = 15, "species" = 19)) +
    scale_size_identity() +
    guides(
      size = "none",
      shape = guide_legend(override.aes = list(size = 10), byrow = TRUE)
    ) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white"),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      legend.position = "bottom",
      plot.title.position = "plot",
      plot.title = ggtext::element_markdown(size = 11, lineheight = 1.2)
    ) +
    labs(title = plot_title)

  return(neighborhood_grob)
}

ggraph_get_edges_by_reversibility <- function (
  is_reversible,
  format = "short",
  collapse = "none",
  ...) {

  # modification of ggraph::get_edges to include filtering to subsets of edges

  checkmate::assertLogical(is_reversible, len = 1)

  if (!collapse %in% c("none", "all", "direction")) {
    cli::cli_abort("{.arg collapse} must be either {.val none}, {.val all} or {.val direction}")
  }
  dots <- rlang::enquos(...)

  function(layout) {

    edges <- ggraph::collect_edges(layout)

    # filter by reversibility
    if (!("r_isreversible" %in% colnames(edges))) {
      cli::cli_abort("Edges must have a column named 'r_isreversible'")
    }

    edges <- edges %>%
      dplyr::filter(r_isreversible == is_reversible)

    edges <- switch(
      collapse,
      none = edges,
      all = ggraph::collapse_all_edges(edges),
      direction = ggraph::collapse_dir_edges(edges)
    )
    edges <- switch(
      format,
      short = ggraph:::format_short_edges(edges, layout),
      long = ggraph:::format_long_edges(edges, layout),
      cli::cli_abort("Unknown {.arg format}. Use either {.val short} or {.val long}")
    )

    extra_data <- lapply(dots, function(x) {
      val <- rlang::eval_tidy(x, edges)
      rep(val, length.out = nrow(edges))
    })
    if (length(extra_data) > 0) {
      edges <- cbind(edges, ggraph:::data_frame0(!!!extra_data))
    }
    attr(edges, "type_ggraph") <- "edge_ggraph"
    edges
  }
}

plot_one_neighborhood_label_vertices <- function(
  vertices,
  max_labeled_species
  ) {

  nodes_to_label <- vertices %>%
    dplyr::filter(node_type == "species") %>%
    dplyr::arrange(path_weight) %>%
    dplyr::slice(1:max_labeled_species) %>%
    {
      .$name
    }

  vertices <- vertices %>%
    dplyr::mutate(
      selected_node_name = ifelse(name %in% nodes_to_label, node_name, NA_character_)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      # extract name and compartment from species name
      compartment = stringr::str_match(selected_node_name, "\\[([A-Za-z0-9 -]+)\\]$")[2],
      species_name = stringr::str_split(selected_node_name, " \\[([A-Za-z0-9 -]+)\\]$")[[1]][1]
    ) %>%
    dplyr::ungroup()

  focal_compartment <- vertices$compartment[vertices$node_orientation == "focal"]

  vertices <- vertices %>%
    dplyr::mutate(
      # ignore compartment if its the same as the focal species compartment
      compartment = ifelse(compartment == focal_compartment & node_orientation != "focal", NA, compartment),
      compartment = stringr::str_c("\n", stringr::str_trunc(compartment, 20)),
      # replace with string so glue works
      compartment = ifelse(is.na(compartment), "", compartment),
      species_name = stringr::str_wrap(stringr::str_trunc(species_name, 40), 25),
      label = ifelse(!is.na(species_name), glue::glue("{species_name}{compartment}"), NA),
      ggtext_label = ifelse(
        !is.na(species_name),
        glue::glue("**{species_name}**{compartment}") %>%
          stringr::str_replace_all(" ?\\n ?", "<br>"),
        NA
      )
    )

  return(vertices)
}

find_obscured_labels <- function(gg_network_layout) {
  plot_char_width <- 100
  plot_char_height <- 50
  graph_width <- diff(range(gg_network_layout$x))
  graph_height <- diff(range(gg_network_layout$y))

  # format labels ordered by plotting priority
  labels <- gg_network_layout %>%
    as.data.frame() %>%
    dplyr::select(x, y, ggtext_label) %>%
    dplyr::filter(!is.na(ggtext_label)) %>%
    dplyr::mutate(n_lines = stringr::str_count(ggtext_label, "<br>") + 1) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      max_width = max(purrr::map_int(
        stringr::str_split(stringr::str_replace_all(ggtext_label, "\\*\\*", ""), "<br>")[[1]],
        stringr::str_length
      ))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      row_id = 1:dplyr::n(),
      # adding characters to width and height to represent box
      w = graph_width * ((max_width + 2) / plot_char_width),
      h = graph_height * ((n_lines + 1) / plot_char_height),
      x_min = x - w / 2,
      x_max = x + w / 2,
      y_min = y - h / 2,
      y_max = y + h / 2
    )

  # define all x all combinations
  overlapping_labels <- tidyr::expand_grid(
    row_1 = labels$row_id,
    row_2 = labels$row_id
  ) %>%
    dplyr::left_join(
      labels %>%
        dplyr::select(
          row_1 = row_id,
          x_min_1 = x_min,
          x_max_1 = x_max,
          y_min_1 = y_min,
          y_max_1 = y_max
        ),
      by = "row_1"
    ) %>%
    dplyr::left_join(
      labels %>%
        dplyr::select(
          row_2 = row_id,
          x_min_2 = x_min,
          x_max_2 = x_max,
          y_min_2 = y_min,
          y_max_2 = y_max
        ),
      by = "row_2"
    ) %>%
    # filter to lower-diagonal comparison
    dplyr::filter(row_1 < row_2) %>%
    # find overlaps
    dplyr::mutate(
      # calculate overlaps in x, and y separately
      # two ranges overlap if the end point of one end of a range is between
      # the min and max of the other
      x_overlap = dplyr::case_when(
        x_min_1 < x_max_2 & x_min_1 > x_min_2 ~ TRUE,
        x_max_1 < x_max_2 & x_max_1 > x_min_2 ~ TRUE,
        x_min_2 < x_max_1 & x_min_2 > x_min_1 ~ TRUE,
        x_max_2 < x_max_1 & x_max_2 > x_min_1 ~ TRUE,
        TRUE ~ FALSE
      ),
      y_overlap = dplyr::case_when(
        y_min_1 < y_max_2 & y_min_1 > y_min_2 ~ TRUE,
        y_max_1 < y_max_2 & y_max_1 > y_min_2 ~ TRUE,
        y_min_2 < y_max_1 & y_min_2 > y_min_1 ~ TRUE,
        y_max_2 < y_max_1 & y_max_2 > y_min_1 ~ TRUE,
        TRUE ~ FALSE
      ),
      collision = x_overlap & y_overlap
    ) %>%
    dplyr::filter(collision)

  if (nrow(overlapping_labels) == 0) {
    return(NULL)
  }

  invalid_labels <- NULL
  while (nrow(overlapping_labels) != 0) {
    # remove the lower priority label of the first row since this will involve
    # the highest plotting priority label.
    # discard all overlaps involving the pruned label
    removed_label <- overlapping_labels$row_2[1]
    invalid_labels <- c(invalid_labels, removed_label)

    overlapping_labels <- overlapping_labels %>%
      dplyr::filter(
        row_1 != removed_label,
        row_2 != removed_label
      )
  }

  # obscured labels
  return(labels$ggtext_label[labels$row_id %in% invalid_labels])
}

select_score_overlay_palette <- function(score_palette, ...) {

  SCORE_PALETTE_NAMES <- c("indication scores", "log2 fold-change")
  if (is.null(score_palette)) {
    cli::cli_abort("Please provide a value for {.arg score_palette}. Valid palettes are {.field {SCORE_PALETTE_NAMES}}")
  }

  if (!(score_palette %in% SCORE_PALETTE_NAMES)) {
    cli::cli_abort("{.field {score_palette}} is not a valid value for {.arg score_palette}. Valid palettes are {.field {SCORE_PALETTE_NAMES}}")
  }

  if (score_palette == "indication scores") {
    score_palette_obj <- scale_color_gradientn(
      "Open Target's Indication Score",
      colors = c("gray90", "yellow", "orange", "orangered", "red"),
      limits = c(0, 1),
      breaks = c(0, 0.1, 0.4, 1),
      trans = "sqrt"
    )

  } else if (score_palette == "log2 fold-change") {
    score_palette_obj <- scale_color_gradient2(
      expression(log[2] ~ "fold-change"),
      low = "steelblue1",
      mid = "black",
      high = "yellow",
      midpoint = 0,
      ...
    )
  }

  return(score_palette_obj)
}
