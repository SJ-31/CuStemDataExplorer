#' compare_expression UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_compare_expression_ui <- function(id) {
  ns <- NS(id)
  bslib::nav_panel(
    "Gene Expression",
    bslib::layout_sidebar(
      shiny::plotOutput(ns("expr_comparison")),
      sidebar = bslib::sidebar(
        shiny::h4("Gene selection"),
        shiny::selectizeInput(
          ns("gene_selection"),
          label = NULL,
          choices = NULL,
          multiple = TRUE,
          options = list(maxItems = 10)
        ),
        shiny::h4("Filters"),
        shiny::selectizeInput(
          ns("tumor_type"),
          "Tumor Type",
          choices = NULL,
          multiple = TRUE,
          options = list(items = NULL)
        ),
        shiny::selectizeInput(
          ns("cohort"),
          "Cohort",
          choices = NULL,
          multiple = TRUE,
          options = list(items = NULL)
        ),
      )
    )
  )
}

#' compare_expression Server Functions
#'
#' @noRd
mod_compare_expression_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    if (options()$golem.app.prod) {
      key <- "default"
    } else {
      key <- "dev"
    }
    cfg <- get_golem_config("expression_viewer", config = key)
    combined_expr <- from_bfc("bulk_expression")
    gene_ids <- purrr::discard(combined_expr$expr$gene_id, is.na) |>
      `names<-`(NULL)

    shiny::updateSelectizeInput(
      session,
      "gene_selection",
      choices = gene_ids,
      server = TRUE
    )
    shiny::updateSelectizeInput(
      session,
      "tumor_type",
      choices = unique(combined_expr$meta$tumor_type),
      server = TRUE
    )
    shiny::updateSelectizeInput(
      session,
      "cohort",
      choices = unique(combined_expr$meta$cohort),
      server = TRUE
    )

    output$expr_comparison <- shiny::renderPlot(
      do_heatmap(
        combined_expr,
        genes = input$gene_selection,
        cfg = cfg$palette %||% list(),
        tumor_types = input$tumor_type,
        cohorts = input$cohort
      ),
      res = 120
    ) |>
      shiny::bindCache(input$gene_selection, input$tumor_type, input$cohort) |>
      shiny::bindEvent(input$gene_selection, input$tumor_type, input$cohort)
  })
}

## To be copied in the UI
# mod_compare_expression_ui("compare_expression_1")

## To be copied in the server
# mod_compare_expression_server("compare_expression_1")
