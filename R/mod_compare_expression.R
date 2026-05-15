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
        shiny::h3("Filters"),
        shiny::selectInput(
          ns("gene_selection"),
          "Genes",
          choices = c(
            "ENSG00000142655",
            "ENSG00000171621",
            "ENSG00000173614",
            "ENSG00000171729",
            "ENSG00000157916"
          ),
          multiple = TRUE
        ),
        shiny::selectInput(
          ns("tumor_type"),
          "Tumor Type",
          choices = c("HCC", "PDAC", "CRC"),
          multiple = TRUE
        ),
        shiny::selectInput(
          ns("cohort"),
          "Cohort",
          choices = c("HCC", "CRC", "PHcase"),
          multiple = TRUE
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
    combined_expr <- read_all_expr(cfg)
    output$expr_comparison <- shiny::renderPlot(do_heatmap(
      combined_expr,
      genes = input$gene_selection,
      cfg = cfg$palette
    )) |>
      shiny::bindCache(input$gene_selection, input$tumor_type, input$cohort) |>
      shiny::bindEvent(input$gene_selection, input$tumor_type, input$cohort)
  })
}

## To be copied in the UI
# mod_compare_expression_ui("compare_expression_1")

## To be copied in the server
# mod_compare_expression_server("compare_expression_1")
