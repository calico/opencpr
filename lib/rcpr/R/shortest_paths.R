#' Summarize Shortest Paths
#'
#' @param source_species_id species id to start with
#' @param dest_species_id species id to end at
#' @inheritParams create_neighborhood_table
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'
#'   # end_points <- sample(rownames(species_names), 2)
#'   # source_species_id <- end_points[1]
#'   # dest_species_id <- end_points[2]
#'
#'   # source_species_id <- "S00002615"
#'   # dest_species_id <- "S00004885"
#'
#'   source_species_id <- "S00002615"
#'   dest_species_id <- "S00000374"
#'
#'   summarize_shortest_paths(
#'     source_species_id,
#'     dest_species_id,
#'     consensus_model,
#'     consensus_graph,
#'     cpr = cpr
#'   )
#' }
#' @export
summarize_shortest_paths <- function(source_species_id, dest_species_id, consensus_model, consensus_graph, cpr) {
  checkmate::assertCharacter(source_species_id)
  checkmate::assertCharacter(dest_species_id)
  checkmate::assertClass(consensus_model, "cpr.sbml.SBML_dfs")
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertClass(cpr, "python.builtin.module")

  debugr::dwatch(msg = "Starting summarize_shortest_paths [cpr<utils.R>::summarize_shortest_paths]")
  print(glue::glue("Searching for path between {source_species_id} and {dest_species_id}"))

  target_species_paths <- cpr$network$compartmentalize_species_pairs(consensus_model, source_species_id, dest_species_id)

  if (nrow(target_species_paths) > 10) {
    print(glue::glue("{nrow(target_species_paths)} pairs of source and destination species; will only us the first 10"))
    target_species_paths <- target_species_paths %>%
      dplyr::slice(1:10)
  }

  shortest_paths_list <- try(
    cpr$network$find_all_shortest_reaction_paths(
      consensus_graph,
      consensus_model,
      target_species_paths,
      weight_var = "weights"
    ),
    silent = TRUE
  )

  debugr::dwatch(msg = "Shortest paths calculated [cpr<utils.R>::summarize_shortest_paths]")

  if (inherits(shortest_paths_list, "try-error")) {
    print(paste0("cpr$network$find_all_shortest_reaction_paths errored due to: ", attr(shortest_paths_list, "condition")))

    disconnected_plot <- ggplot(data.frame(x = 0, y = 0), aes(x = x, y = y)) +
      geom_text(label = glue::glue("No path exists between these species"), size = 10) +
      theme(text = element_blank(), line = element_blank())

    return(list(
      shortest_paths_plot = disconnected_plot,
      shortest_paths_table = tibble()
    ))
  }

  shortest_paths_grob <- plot_shortest_path_network(shortest_paths_list, consensus_graph, max_labeled_species = 10L)

  return(list(
    shortest_paths_plot = shortest_paths_grob,
    shortest_paths_table = shortest_paths_list[[1]] %>%
      dplyr::select(source, dest, path, step, node, node_type, label, weights, url) %>%
      dplyr::arrange(source, dest, path, step)
  ))
}

#' Plot Shortest Path Network
#'
#' @param shortest_paths_list results from cpr$network$find_all_shortest_reaction_paths
#' @inheritParams create_neighborhood_table
#' @inheritParams plot_one_neighborhood
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'
#'   end_points <- sample(rownames(species_names), 2)
#'   source_species_id <- end_points[1]
#'   dest_species_id <- end_points[2]
#'
#'   target_species_paths <- cpr$network$compartmentalize_species_pairs(
#'     consensus_model,
#'     source_species_id,
#'     dest_species_id
#'     )
#'
#'   shortest_paths_list <- try(cpr$network$find_all_shortest_reaction_paths(
#'     consensus_graph,
#'     consensus_model,
#'     target_species_paths,
#'     weight_var = "weights"
#'     ),silent = TRUE)
#'
#'   plot_shortest_path_network(
#'     shortest_paths_list,
#'     consensus_graph,
#'     max_labeled_species = 10L
#'     )
#' }
#' @export
plot_shortest_path_network <- function(shortest_paths_list, consensus_graph, max_labeled_species = 10L) {
  checkmate::assertList(shortest_paths_list, len = 4)
  checkmate::assertCount(max_labeled_species)

  # suppress some nodes labels to avoid overplotting
  nodes_to_label <-
    # preserve source and destination
    c(
      shortest_paths_list[[1]] %>%
        {
          .$node[.$node == .$source | .$node == .$dest]
        } %>%
        unique(),
      shortest_paths_list[[1]] %>%
        dplyr::filter(
          !(node == source | node == dest),
          node_type == "species"
        ) %>%
        dplyr::distinct(node, weights) %>%
        dplyr::arrange(weights) %>%
        dplyr::slice(1:max_labeled_species) %>%
        {
          .$node
        }
    )

  all_shortest_reaction_paths_df <- shortest_paths_list[[1]] %>%
    dplyr::relocate(node) %>%
    dplyr::group_by(node) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      special_node = factor(
        dplyr::case_when(
          node == source ~ "source",
          node == dest ~ "destination",
          TRUE ~ "mediator"
        ),
        levels = c("source", "destination", "mediator")
      ),
      node_name = label,
      label = dplyr::case_when(
        node %in% nodes_to_label ~ stringr::str_wrap(node_name, 40),
        TRUE ~ NA_character_
      )
    )

  all_shortest_reaction_path_edges_df <- shortest_paths_list[[2]] %>%
    dplyr::relocate(from, to)

  # summarize path
  path_labels <- all_shortest_reaction_path_edges_df %>%
    dplyr::group_by(source, dest, path) %>%
    dplyr::summarize(
      steps = n(),
      score = sum(weights),
      .groups = "drop"
    ) %>%
    dplyr::arrange(score) %>%
    dplyr::group_by(dest) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(source), !is.na(dest)) %>%
    dplyr::group_by(dest) %>%
    dplyr::mutate(
      path_qual_score = dplyr::case_when(
        score > 50 ~ "dubious",
        score > 30 ~ "possible",
        score > 10 ~ "likely",
        TRUE ~ "highly likely"
      ),
      path_label = glue::glue(
        "**{path_qual_score}** path<br>
        *steps* = {steps}<br>
        *score* = {score}"
      )
    )

  # use source info when available
  reaction_sources <- shortest_paths_list[[3]]
  if (!("NULL" %in% class(reaction_sources))) {
    all_shortest_reaction_path_edges_df <- all_shortest_reaction_path_edges_df %>%
      dplyr::bind_rows(
        reaction_sources %>%
          dplyr::select(from = "r_id", to = "pathway_id")
      )

    all_shortest_reaction_paths_df <- all_shortest_reaction_paths_df %>%
      dplyr::bind_rows(reaction_sources %>% dplyr::distinct(node = pathway_id))
  }

  shortest_paths_network <- igraph::graph_from_data_frame(all_shortest_reaction_path_edges_df,
    directed = consensus_graph$is_directed(),
    vertices = all_shortest_reaction_paths_df
  )

  debugr::dwatch(msg = "R igraph network created [cpr<shortest_paths.R>::summarize_shortest_paths]")

  # format title
  terminal_species <- all_shortest_reaction_paths_df %>%
    dplyr::filter(special_node != "mediator") %>%
    dplyr::group_by(special_node) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(species_name = stringr::str_split(node_name, " \\[([A-Za-z0-9 -]+)\\]$")[[1]][1])

  plot_title <- glue::glue(
    "Shortest path(s)<br>
     **from:** {terminal_species$species_name[terminal_species$special_node == 'source']}<br>
     **to:** {terminal_species$species_name[terminal_species$special_node == 'destination']}
     "
  )

  gg_network_layout <- layout_with_reaction_sources(shortest_paths_network)

  shortest_paths_grob <- ggraph::ggraph(graph = gg_network_layout)

  path_labels_layout <- gg_network_layout %>%
    dplyr::select(name, x, y) %>%
    dplyr::inner_join(path_labels, by = c("name" = "dest"))

  if (!("NULL" %in% class(reaction_sources))) {
    # add reaction sources if available
    pathway_coords <- layout_pathway_sources(gg_network_layout, reaction_sources)

    shortest_paths_grob <- shortest_paths_grob +
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

  shortest_paths_grob <- shortest_paths_grob +
    ggraph::geom_edge_link(color = "gray25", arrow = grid::arrow(type = "closed"), length = unit(0.15, "inches")) +
    ggraph::geom_node_point(aes(shape = factor(node_type), color = special_node), size = 3) +
    ggraph::geom_node_label(aes(label = label, color = special_node),
      repel = TRUE, alpha = 0.7, size = 4,
      label.padding = 0.1, fill = "gray90"
    ) +
    ggtext::geom_richtext(data = path_labels_layout, aes(x = x + 0.2, y = y, label = path_label), color = "RED", hjust = 0, size = 3) +
    scale_shape_manual("Type", values = c("reaction" = 15, "species" = 19)) +
    scale_color_manual(values = c("source" = "blue", "destination" = "red", "mediator" = "gray25")) +
    guides(
      shape = guide_legend(override.aes = list(size = 8, color = "gray25"), nrow = 2),
      fill = guide_legend(nrow = 2, byrow = TRUE),
      color = "none"
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

  return(shortest_paths_grob)
}
