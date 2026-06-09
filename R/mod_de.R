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
  box::use(
    reactable[colDef, colFormat, reactable],
    glue[glue],
    DBI[dbGetQuery]
  )

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    cfg <- get_golem_config("de_viewer")
    # [2026-06-08 Mon] TODO: don't wanna have to read everything in...
    # is there a better way?
    keys <- from_bfc(cached)
    dbs <- list(
      ovr = from_bfc(keys$ovr),
      pairwise = from_bfc(keys$pairwise)
    )

    table_lookup <- "SELECT * FROM INFORMATION_SCHEMA.TABLES"

    all_contrasts <- c(
      dbGetQuery(dbs$ovr, table_lookup)$table_name,
      dbGetQuery(dbs$pairwise, table_lookup)$table_name
    ) |>
      sort()

    get_contrast <- function(contrast, indices = NULL) {
      if (stringr::str_ends(contrast, "vs. Rest")) {
        con <- dbs$ovr
      } else {
        con <- dbs$pairwise
      }
      df <- dbGetQuery(con, glue("SELECT * FROM '{contrast}'")) |>
        dplyr::arrange(log2FoldChange)
      if (!is.null(indices)) {
        df[indices, ]
      } else {
        df
      }
    }

    to_numeric <- c("stat", "pvalue", "lfcSE", "baseMean", "log2FoldChange")
    de_col_format <- lapply(
      rep(1, length(to_numeric)),
      \(x) {
        reactable::colDef(format = reactable::colFormat(digits = 3))
      }
    ) |>
      `names<-`(to_numeric)

    shiny::updateSelectizeInput(
      session,
      "contrast",
      choices = all_contrasts,
      selected = all_contrasts[1],
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
      columns = de_col_format,
      selection = "multiple",
      searchable = TRUE
    )) |>
      shiny::bindCache(input$contrast) |>
      shiny::bindEvent(input$contrast)
  })
}

## To be copied in the UI
# mod_de_ui("de_1")

## To be copied in the server
# mod_de_server("de_1")
