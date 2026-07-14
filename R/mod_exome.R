#' exome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_exome_ui <- function(id, label) {
  ns <- NS(id)
  bslib::nav_panel(
    label,
    bslib::layout_sidebar(
      reactable::reactableOutput(ns("table")),
      reactable::reactableOutput(ns("vtable")),
      sidebar = bslib::sidebar(
        shiny::h3("Tumor type"),
        shiny::selectizeInput(
          ns("ttype"),
          label = NULL,
          choices = NULL,
          multiple = FALSE
        ),
        variant_table_download_button("Download", ns("table"))
      )
    )
  )
}

#' exome Server Functions
#'
#' @noRd
mod_exome_server <- function(id, cached) {
  box::use(DBI[dbGetQuery])
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    con <- DBI::dbConnect(
      duckdb::duckdb(),
      dbdir = from_bfc(cached)$genes,
      read_only = TRUE
    )
    ttypes <- dbGetQuery(
      con,
      "SELECT * FROM INFORMATION_SCHEMA.TABLES"
    )$table_name

    shiny::updateSelectizeInput(
      inputId = "ttype",
      choices = ttypes,
      selected = ttypes[1],
      server = TRUE
    )

    make_table <- function(type) {
      data <- dbGetQuery(con, sprintf("SELECT * FROM %s", type))
      make_variant_table_at_level(data, level = "genes")
    }

    chosen_index <- reactive({
      reactable::getReactableState("table", name = "selected")
    })

    make_vtable <- function(type, index) {
      data <- dbGetQuery(con, sprintf("SELECT * FROM %s", type))
      make_variant_table_at_level(
        data,
        level = "variants",
        chosen = data[index, ]$Gene
      )
    }

    # Could also make a separate view of the variant-level details, rather than have it be a nested table
    # could make things faster

    output$table <- reactable::renderReactable({
      req(input$ttype)
      make_table(input$ttype)
    }) |>
      shiny::bindCache(input$ttype) |>
      shiny::bindEvent(input$ttype)

    output$vtable <- reactable::renderReactable({
      req(input$ttype)
      req(chosen_index())
      make_vtable(input$ttype, chosen_index())
    }) |>
      shiny::bindCache(input$ttype, chosen_index()) |>
      shiny::bindEvent(input$ttype, chosen_index())
  })
}

## To be copied in the UI
# mod_exome_ui("exome_1")

## To be copied in the server
# mod_exome_server("exome_1")
