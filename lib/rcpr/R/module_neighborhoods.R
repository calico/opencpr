#' Shiny Neighborhood Test
#'
#' Demo for interrogating the neighborhood of an entity.
#'
#' @inheritParams neighborhoodServer
#'
#' @returns A shiny app
#'
#' @examples
#'
#' if (interactive()) {
#'   interactive_initialization_wrapper()
#'   cpr <- reticulate::import("cpr")
#'
#'   shiny_neighborhood_test(
#'     species_names,
#'     species_identifiers,
#'     consensus_model,
#'     consensus_graph,
#'     cpr = cpr
#'   )
#' }
#' @export
shiny_neighborhood_test <- function(
  species_names,
  species_identifiers,
  consensus_model,
  consensus_graph,
  cpr
  ) {
  checkmate::assertDataFrame(species_names)
  checkmate::assertDataFrame(species_identifiers)
  checkmate::assertClass(consensus_model, "cpr.sbml.SBML_dfs")
  checkmate::assertClass(consensus_graph, "igraph.Graph")
  checkmate::assertClass(cpr, "python.builtin.module")

  shiny::shinyApp(
    ui = shiny::fluidPage(
      neighborhoodInput("neighborhood_app")
    ),
    server = function(input, output, session) {
      neighborhoodServer(
        "neighborhood_app",
        species_names,
        species_identifiers,
        consensus_model,
        consensus_graph,
        cpr
      )
    }
  )
}

#' Neighborhood Input
#'
#' UI components for the neighborhood shiny module
#'
#' @inheritParams shiny::moduleServer
#'
#' @returns Nothing; used for side-effects
#'
#' @export
neighborhoodInput <- function(id) {
  checkmate::assertCharacter(id, len = 1)

  ns <- shiny::NS(id)

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::h3("Choose a focal gene/molecule"),
      selectEntityInput(
        ns("neighborhood_central_entity"),
        "central reaction species"
      ),
      shiny::h3("Build a neighborhood"),
      shiny::sliderInput(
        ns("max_neighbors"),
        label = "Target Number of Gene/Metabolite Neighbors",
        min = 5L,
        max = 100L,
        value = 5L,
        step = 5L
      ),
      shiny::sliderInput(
        ns("max_steps_slider"),
        label = "Maximum Number of Steps to Find Neighbors",
        min = 1L,
        max = 25L,
        value = 3L,
        step = 2L
      ),
      shiny::radioButtons(
        ns("network_type"),
        label = "Neightborhood type",
        choices = c("hourglass", "downstream", "upstream"),
        inline = TRUE
      ),
      shiny::h3("Find neighbors' indications"),
      shiny::actionButton(
        ns("update_neighborhood"),
        "Query Open Targets (slow)"
      ),
      shiny::checkboxInput(
        ns("is_overlay_disease"),
        label = "Overlay a diseases?",
        FALSE
        ),
      shiny::uiOutput(ns("overlaid_disease_ui")),
      shiny::textOutput(ns('ot_overlay_error'))
    ),
    shiny::mainPanel(
      shiny::tabsetPanel(
        type = "tabs",
        shiny::tabPanel(
          "Plot",
          romic::plotsaverInput(ns("ggsave_neighborhood"), "wide"),
          shiny::plotOutput(ns("neighborhood_plot")),
        ),
        shiny::tabPanel(
          "Summary Table",
          DT::dataTableOutput(ns("neighborhood_table")),
          shiny::downloadButton(ns("download_neighborhood_table"), "Download .tsv")
        ),
        shiny::tabPanel(
          "Targets",
          DT::dataTableOutput(ns("targets")),
          shiny::downloadButton(ns("download_targets_table"), "Download .tsv")
        ),
        shiny::tabPanel(
          "Top Indications",
          DT::dataTableOutput(ns("indications")),
          shiny::downloadButton(ns("download_indications_table"), "Download .tsv")
        ),
        shiny::tabPanel(
          "Selected Indication",
          DT::dataTableOutput(ns("select_indication")),
          shiny::downloadButton(ns("download_select_indication"), "Download .tsv")
        )
      )
    )
  )
}

#' Neighborhood Server
#'
#' Server-side components of the neighborhood shiny module
#'
#' @inheritParams shiny::moduleServer
#' @inheritParams selectEntityServer
#' @param consensus_model A python pathway object with class cpr.sbml.SBML_dfs
#' @param consensus_graph A python igraph network of \code{consensus_model}
#' @param cpr reticulate connection to the cpr python package
#'
#' @returns Nothing; used for side-effects
#'
#' @export
neighborhoodServer <- function(id,
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

      shiny::observe({

        req(input$is_overlay_disease)

        output$overlaid_disease_ui <- shiny::renderUI({
          if (input$is_overlay_disease) {
            shiny::selectizeInput(
              ns("overlaid_disease"),
              "Disease:",
              selected = "loading",
              choices = NULL
            )
          } else {
            return()
          }
        })

        # add available diseases server-side for performance
        shiny::updateSelectizeInput(
          session,
          "overlaid_disease",
          selected = "diabetes mellitus",
          choices = efo_diseases$name,
          server = TRUE
        )
      })

      selected_neighborhood_central_entity <- selectEntityServer(
        "neighborhood_central_entity",
        species_names,
        species_identifiers
      )

      neighborhood_table <- reactive({
        shiny::req(selected_neighborhood_central_entity(), input$max_steps_slider)
        print(glue::glue("Creating neighborhood around {selected_neighborhood_central_entity()} within max steps {input$max_steps_slider}}"))

        # Find the neighbors of a species in each compartment and format as a table
        create_neighborhood_table(
          selected_neighborhood_central_entity(),
          consensus_model,
          consensus_graph,
          network_type = input$network_type,
          max_steps = input$max_steps_slider,
          max_neighbors = input$max_neighbors,
          cpr
        )
      })

      # create a summary table with one row per vertex
      neighborhood_summary_table <- shiny::reactive({
        shiny::req(neighborhood_table())
        create_neighborhood_summary_table(neighborhood_table())
      })

      # create a prefix for any saved files based on the focal nodes name
      # and neighborhood parameters

      filename_prefix <- shiny::reactive({
        shiny::req(neighborhood_summary_table())

        species_name <- neighborhood_summary_table() %>%
          dplyr::filter(node_orientation == "focal") %>%
          {
            .$sc_name[1]
          } %>%
          stringr::str_replace(" \\[[A-Za-z ]+\\]$", "") %>%
          stringr::str_replace_all(" ", "_")

        network_type_char <- stringr::str_extract(input$network_type, "^[a-z]")

        glue::glue("{species_name}_{network_type_char}{input$max_steps_slider}s{input$max_neighbors}n")
      })

      # query indication scores of neighbors for a disease
      indication_overlay <- shiny::reactive({
        if (!is.null(input$is_overlay_disease) && input$is_overlay_disease) {
          shiny::req(input$overlaid_disease)

          efo_id <- efo_diseases$id[efo_diseases$name == input$overlaid_disease]
          checkmate::assertString(efo_id)

          print(glue::glue(
            "Querying neighborhood for indications related to {efo_id}: {input$overlaid_disease}"
          ))

          indication_overlay <- summarize_indication(
            efo_id,
            consensus_model,
            neighborhood_summary_table(),
            species_identifiers
          )

          if ("ot_error" %in% names(attributes(indication_overlay))) {
            output$ot_overlay_error <- renderText({
              stop(safeError(attr(indication_overlay, "ot_error")))
            })
          }

          # update summary table

          output$select_indication <- DT::renderDataTable(
            indication_overlay %>%
              dplyr::rename(evidence = url) %>%
              add_url_links(),
            escape = FALSE
          )

          # download table
          disease_no_space <- stringr::str_replace_all(input$overlaid_disease, "[[:space:]]+", "_")

          output$download_select_indication <- shiny::downloadHandler(
            filename = glue::glue("{filename_prefix()}_{disease_no_space}.tsv"),
            content = function(file) {
              vroom::vroom_write(indication_overlay, file)
            }
          )

          # return indication overlay for broader consumption
          indication_overlay
        } else {
          NULL
        }
      })

      arranged_neighborhood_plots <- shiny::reactive({
        plot_neighborhoods(
          neighborhood_table(),
          consensus_graph,
          indication_overlay = indication_overlay(),
          indication_label = input$overlaid_disease
        )
      })

      # render plots
      shiny::observe({
        shiny::req(neighborhood_summary_table(), arranged_neighborhood_plots())

        # N neighborhoods
        n_neighborhoods <- length(unique(neighborhood_summary_table()$sc_name))
        if (n_neighborhoods > 4) {
          n_neighborhoods <- 4
        }

        output$neighborhood_plot <- shiny::renderPlot(
          graphics::plot(arranged_neighborhood_plots()),
          height = n_neighborhoods * 1000
        )
        output$neighborhood_table <- DT::renderDataTable(
          neighborhood_summary_table() %>%
            dplyr::select(-s_id) %>%
            dplyr::rename(reference = url, `gene info` = pharos) %>%
            add_url_links(),
          escape = FALSE
        )

        # save plot
        romic::plotsaverServer(
          "ggsave_neighborhood",
          arranged_neighborhood_plots(),
          glue::glue("{filename_prefix()}_plot.png")
        )

        # download table
        output$download_neighborhood_table <- shiny::downloadHandler(
          filename = glue::glue("{filename_prefix()}_neighborhood.tsv"),
          content = function(file) {
            vroom::vroom_write(neighborhood_summary_table(), file)
          }
        )
      })

      shiny::observeEvent(input$update_neighborhood, {
        shiny::req(neighborhood_summary_table())

        target_information <- summarize_diseases(
          consensus_model,
          neighborhood_summary_table(),
          species_identifiers
        )

        output$targets <- DT::renderDataTable(target_information$targets)
        output$indications <- DT::renderDataTable(target_information$indications)

        # download table
        output$download_targets_table <- shiny::downloadHandler(
          filename = glue::glue("{filename_prefix()}_targets.tsv"),
          content = function(file) {
            vroom::vroom_write(target_information$targets, file)
          }
        )

        output$download_indications_table <- shiny::downloadHandler(
          filename = glue::glue("{filename_prefix()}_indications.tsv"),
          content = function(file) {
            vroom::vroom_write(target_information$indications, file)
          }
        )
      })
    }
  )
}
