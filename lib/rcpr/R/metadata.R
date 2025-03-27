layout_with_reaction_sources <- function(network_graph, network_layout = "fr") {
  checkmate::assertClass(network_graph, "igraph")

  # create a layout
  set.seed(1234)

  LAYOUT_OPTIONS <- c("fr", "kk", "drl")
  checkmate::assertChoice(network_layout, LAYOUT_OPTIONS)

  layout_fxn <- switch(
    network_layout,
    "fr" = igraph::layout_with_fr,
    "kk" = igraph::layout_with_kk,
    "drl" = igraph::layout_with_drl
  )

  network_layout_coords <- do.call(layout_fxn, list(graph = network_graph, dim = 2))
  network_layout_coords <- as.data.frame(network_layout_coords)
  colnames(network_layout_coords) <- c("x", "y")

  # remove pathways as layout nodes
  is_pathway_node <- is.na(igraph::vertex_attr(network_graph, "node_type"))
  network_layout_coords <- network_layout_coords[!is_pathway_node, ]
  network_graph <- igraph::delete_vertices(network_graph, igraph::vertex_attr(network_graph, "name")[is_pathway_node])

  # generate a ggraph object with node coordinates (for ggplot2-based plotting)
  gg_network_layout <- ggraph::create_layout(network_graph, layout = "manual", x = network_layout_coords$x, y = network_layout_coords$y)

  return(gg_network_layout)
}

layout_pathway_sources <- function(gg_network_layout, reaction_sources, max_pathways = 8L) {
  pathway_coords <- gg_network_layout %>%
    dplyr::select(name, x, y) %>%
    dplyr::inner_join(reaction_sources %>%
      dplyr::select(r_id, label = name),
    by = c("name" = "r_id")
    ) %>%
    dplyr::group_by(label) %>%
    dplyr::summarize(
      n_species = n(),
      x_min = min(x),
      x_max = max(x),
      y_min = min(y),
      y_max = max(y),
      w = x_max - x_min,
      h = y_max - y_min
    ) %>%
    dplyr::mutate(label_wrap = stringr::str_wrap(label, width = 50)) %>%
    #
    dplyr::arrange(dplyr::desc(n_species)) %>%
    dplyr::filter(n_species > 1) %>%
    dplyr::slice(1:max_pathways) %>%
    dplyr::mutate(label_wrap = factor(label_wrap, levels = label_wrap))

  pathway_coords
}
