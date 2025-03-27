#' Shiny Shortest Paths Test
#'
#' Demo for finding the shortest paths between pairs of entities.
#'
#' @inheritParams shortestPathsServer
#'
#' @returns a shiny app
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'   cpr <- reticulate::import("cpr")
#'
#'   shiny_shortestPaths_test(
#'     species_names,
#'     species_identifiers,
#'     consensus_model,
#'     consensus_graph,
#'     cpr = cpr
#'     )
#' }
#' @export
shiny_shortestPaths_test <- function(species_names,
                                     species_identifiers,
                                     consensus_model,
                                     consensus_graph,
                                     cpr) {
  checkmate::assertDataFrame(species_names)
  checkmate::assertDataFrame(species_identifiers)
  checkmate::assertClass(consensus_model, "cpr.sbml.SBML_dfs")
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertClass(cpr, "python.builtin.module")

  shiny::shinyApp(
    ui = shiny::fluidPage(
      shortestPathsInput("shortest_paths_app")
    ),
    server = function(input, output, session) {
      shortestPathsServer(
        "shortest_paths_app",
        species_names,
        species_identifiers,
        consensus_model,
        consensus_graph,
        cpr
      )
    }
  )
}

#' Shortest Paths Input
#'
#' UI components for the shortest paths shiny module
#'
#' @inheritParams shiny::moduleServer
#'
#' @returns Nothing; used for side-effects
#'
#' @export
shortestPathsInput <- function(id) {
  checkmate::assertCharacter(id, len = 1)

  ns <- shiny::NS(id)

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::tagList(
        shiny::h3("Choose source and target destination"),
        selectEntityInput(ns("source_entity"), "source reaction species"),
        selectEntityInput(ns("dest_entity"), "destination reaction species")
      )
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        type = "tabs",
        shiny::tabPanel(
          "Plot",
          romic::plotsaverInput(ns("ggsave_path"), "wide"),
          shiny::plotOutput(ns("shortest_paths_plot"))
        ),
        shiny::tabPanel(
          "Table",
          DT::dataTableOutput(ns("shortest_paths_table")),
          shiny::downloadButton(ns("download_paths_table"), "Download .tsv")
        )
      )
    )
  )
}

#' Shortest Paths Server
#'
#' Server-side components for the shortest paths shiny module
#'
#' @inheritParams shiny::moduleServer
#' @inheritParams selectEntityServer
#' @inheritParams neighborhoodServer
#'
#' @returns Nothing; used for side-effects
#'
#' @export
shortestPathsServer <- function(id,
                                species_names,
                                species_identifiers,
                                consensus_model,
                                consensus_graph,
                                cpr) {
  checkmate::assertCharacter(id, len = 1)
  checkmate::assertDataFrame(species_names)
  checkmate::assertDataFrame(species_identifiers)
  checkmate::assertClass(consensus_model, "cpr.sbml.SBML_dfs")
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertClass(cpr, "python.builtin.module")

  shiny::moduleServer(
    id,
    ## Below is the module function
    function(input, output, session) {
      ns <- session$ns

      selected_source_entity <- selectEntityServer("source_entity", species_names, species_identifiers)
      selected_dest_entity <- selectEntityServer("dest_entity", species_names, species_identifiers)

      shiny::observe({
        req(selected_source_entity(), selected_dest_entity())
        debugr::dwatch(msg = glue::glue("Finding path from {selected_source_entity()} to {selected_dest_entity()} [cpr<utils.R>::shortestPathsServer]"))

        # find shortest paths
        shortest_paths_summary <- summarize_shortest_paths(
          selected_source_entity(),
          selected_dest_entity(),
          consensus_model,
          consensus_graph,
          cpr
        )

        req(shortest_paths_summary)
        output$shortest_paths_plot <- renderPlot(
          shortest_paths_summary$shortest_paths_plot,
          height = 1000
        )
        output$shortest_paths_table <- DT::renderDataTable(
          add_url_links(shortest_paths_summary$shortest_paths_table),
          escape = FALSE
        )

        # download plot
        romic::plotsaverServer(
          "ggsave_path",
          shortest_paths_summary$shortest_paths_plot
        )

        # download table
        output$download_paths_table <- downloadHandler(
          filename = glue::glue("shortest_paths.tsv"),
          content = function(file) {
            vroom::vroom_write(shortest_paths_summary$shortest_paths_table, file)
          }
        )
      })
    }
  )
}

add_url_links <- function(df) {

  # convert url variable in df to html links

  checkmate::assertDataFrame(df)

  # find url variables

  url_variables <- df %>%
    dplyr::select_if(is.character) %>%
    purrr::map_lgl(~ any(stringr::str_detect(., "^http"), na.rm = TRUE)) %>%
    {
      names(.[.])
    }

  if (length(url_variables) == 0) {
    print("No urls found")
    return(df)
  }

  annotated_df <- df
  for (url_var in url_variables) {
    annotated_df <- annotated_df %>%
      # create a temporary variable containing the url of interest
      # (I couldn't get the tidy evaluation working with)
      dplyr::mutate(
        .tmp = !!rlang::sym(url_var),
        !!rlang::sym(url_var) := dplyr::case_when(
          .tmp == "" ~ glue::glue(""),
          is.na(.tmp) ~ glue::glue(""),
          TRUE ~ glue::glue('<a href="{.tmp}" target="_blank" class="btn btn-primary">Link</a>')
        )
      ) %>%
      dplyr::select(-.tmp)
  }

  return(annotated_df)
}
