#' bulk_de UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_de_ui <- function(id, label) {
  ns <- NS(id)
  bslib::nav_panel(
    label,
    bslib::layout_sidebar(
      shiny::plotOutput(ns("volcano")),
      reactable::reactableOutput("de_results"),
      sidebar = bslib::sidebar(
        shiny::h3("Contrast"),
        shiny::selectizeInput(
          ns("contrast"),
          label = NULL,
          choices = NULL,
          multiple = FALSE
        )
      )
    ),
  )
}

#' de Server Functions
#'
#' @noRd
mod_de_server <- function(id, cached) {
  box::use(reactable[colDef, colFormat, reactable])

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cfg <- get_golem_config("de_viewer")
    # [2026-06-08 Mon] TODO: don't wanna have to read everything in...
    # is there a better way?
    de_results <- from_bfc(cached)
    all_contrasts <- c(names(de_results$ovr), names(de_results$pairwise)) |>
      sort()

    get_contrast <- function(contrast, indices = NULL) {
      if (contrast %in% names(de_results$ovr)) {
        df <- de_results$ovr[[contrast]]
      } else {
        df <- de_results$pairwise[[contrast]]
      }
      df <- df |>
        as.data.frame() |>
        tibble::rownames_to_column(var = "gene")
      if (!is.null(indices)) {
        df[indices, ]
      } else {
        df
      }
    }

    shiny::updateSelectizeInput(
      session,
      "contrast",
      choices = all_contrasts,
      server = TRUE
    )
    ## selection <- reactable::getReactableState(
    ##   "de_results",
    ##   name = "selected"
    ## )
    output$volcano <- shiny::renderPlot(
      volcano_plot(get_contrast(
        input$contrast,
        reactable::getReactableState(
          "de_results",
          name = "selected"
        )
      ))
    ) |>
      shiny::bindCache(input$contrast) |>
      shiny::bindEvent(input$contrast)

    output$de_results <- reactable::renderReactable(reactable(
      get_contrast(
        input$contrast
      ),
      columns = list()
    )) |>
      shiny::bindCache(input$contrast) |>
      shiny::bindEvent(input$contrast)
  })
}

## To be copied in the UI
# mod_de_ui("de_1")

## To be copied in the server
# mod_de_server("de_1")
